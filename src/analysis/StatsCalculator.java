package analysis;

import model.SequenceRecord;
import model.StatsResult;

import java.util.*;
import java.util.concurrent.*;
import java.util.stream.Collectors;

public class StatsCalculator {
    private final int k = 5;
    private final int maxKmersToCount;
    private final int threads;

    public StatsCalculator() {
        this(Runtime.getRuntime().availableProcessors(), 100_000);
    }

    public StatsCalculator(int threads, int maxKmersToCount) {
        this.threads = Math.max(1, threads);
        this.maxKmersToCount = Math.max(1, maxKmersToCount);
    }

    public StatsResult compute(List<SequenceRecord> records) {
        StatsResult result = new StatsResult();
        if (records == null || records.isEmpty()) return result;

        List<List<SequenceRecord>> parts = partition(records, threads);
        ExecutorService ex = Executors.newFixedThreadPool(threads);
        try {
            List<Future<LocalStats>> futures = new ArrayList<>();
            for (List<SequenceRecord> p : parts) futures.add(ex.submit(() -> processPartition(p)));

            long totalLength = 0L;
            long totalUngapped = 0L;
            long ambiguousBases = 0L;
            long gapRuns = 0L;

            Map<Character, Long> baseCountGapped = initBaseCountMap();
            Map<Character, Long> baseCountUngapped = initBaseCountMap();

            List<Long> scaffoldLengths = new ArrayList<>();
            List<Long> contigLengths = new ArrayList<>();

            Map<String, Integer> seqCounts = new HashMap<>();

            double entropySum = 0.0;

            Map<Long, Long> kmerCounts = new HashMap<>();
            long kmerBudgetUsed = 0L;

            int seqTotal = 0;

            for (Future<LocalStats> f : futures) {
                LocalStats s = f.get();
                seqTotal += s.sequenceCount;
                totalLength += s.totalLength;
                totalUngapped += s.totalUngapped;
                ambiguousBases += s.ambiguousBases;
                gapRuns += s.gapRuns;
                entropySum += s.entropySum;

                for (var e : s.baseCountGapped.entrySet()) baseCountGapped.merge(e.getKey(), e.getValue(), Long::sum);
                for (var e : s.baseCountUngapped.entrySet()) baseCountUngapped.merge(e.getKey(), e.getValue(), Long::sum);

                scaffoldLengths.addAll(s.scaffoldLengths);
                contigLengths.addAll(s.contigLengths);

                for (var e : s.seqCounts.entrySet()) seqCounts.merge(e.getKey(), e.getValue(), Integer::sum);

                for (var e : s.kmerCounts.entrySet()) {
                    kmerCounts.merge(e.getKey(), e.getValue(), Long::sum);
                    kmerBudgetUsed += e.getValue();
                }
            }

            result.sequenceCount = seqTotal;
            result.totalLength = totalLength;
            result.totalUngappedLength = totalUngapped;
            result.ambiguousBaseCount = ambiguousBases;
            result.gapCount = gapRuns;

            Collections.sort(scaffoldLengths, Collections.reverseOrder());
            result.scaffoldN50 = computeNxx(scaffoldLengths, 0.5);
            result.scaffoldN90 = computeNxx(scaffoldLengths, 0.9);
            result.scaffoldL50 = computeNxx(scaffoldLengths, 0.5);
            result.scaffoldL90 = computeNxx(scaffoldLengths, 0.9);

            Collections.sort(contigLengths, Collections.reverseOrder());
            result.contigN50 = computeNxx(contigLengths, 0.5);
            result.contigN90 = computeNxx(contigLengths, 0.9);
            result.contigL50 = computeNxx(contigLengths, 0.5);
            result.contigL90 = computeNxx(contigLengths, 0.9);

            result.minLength = scaffoldLengths.isEmpty() ? 0 : (int) Math.min(scaffoldLengths.stream().mapToLong(Long::longValue).min().orElse(0L), Integer.MAX_VALUE);
            result.maxLength = scaffoldLengths.isEmpty() ? 0 : (int) Math.min(scaffoldLengths.stream().mapToLong(Long::longValue).max().orElse(0L), Integer.MAX_VALUE);
            result.avgLength = seqTotal == 0 ? 0.0 : (double) totalLength / seqTotal;

            long totalGapped = baseCountGapped.values().stream().mapToLong(Long::longValue).sum();
            result.baseFrequencies = new LinkedHashMap<>();
            for (char c : new char[]{'A', 'C', 'G', 'T', 'N'}) result.baseFrequencies.put(c, totalGapped == 0 ? 0.0 : baseCountGapped.getOrDefault(c, 0L) / (double) totalGapped);
            result.gcContent = (baseCountGapped.getOrDefault('G',0L) + baseCountGapped.getOrDefault('C',0L)) / (double) Math.max(1L, totalGapped);

            long totalUngappedBases = baseCountUngapped.values().stream().mapToLong(Long::longValue).sum();
            result.gcContentUngapped = (baseCountUngapped.getOrDefault('G',0L) + baseCountUngapped.getOrDefault('C',0L)) / (double) Math.max(1L, totalUngappedBases);

            int uniqueSeqs = seqCounts.size();
            long dupCount = seqCounts.values().stream().mapToLong(Integer::longValue).filter(v -> v > 1).sum();
            result.duplicateSequenceCount = (int) dupCount;
            result.duplicateFraction = result.sequenceCount == 0 ? 0.0 : (result.sequenceCount - uniqueSeqs) / (double) result.sequenceCount;

            result.avgShannonEntropy = result.sequenceCount == 0 ? 0.0 : entropySum / result.sequenceCount;

            result.kmerCounts = topKmersAsStrings(kmerCounts, 100, k);

            result.lengthHistogram = buildLengthHistogram(scaffoldLengths);

            return result;
        } catch (InterruptedException | ExecutionException e) {
            throw new RuntimeException(e);
        } finally {
            ex.shutdownNow();
        }
    }

    private LocalStats processPartition(List<SequenceRecord> records) {
        LocalStats s = new LocalStats();
        for (SequenceRecord rec : records) {
            String seq = rec.getSequence();
            if (seq == null) seq = "";
            String up = seq.toUpperCase(Locale.ROOT);
            int len = up.length();
            s.sequenceCount++;
            s.totalLength += len;
            s.scaffoldLengths.add((long) len);

            for (int i = 0; i < len; i++) {
                char c = up.charAt(i);
                if (c=='A' || c=='C' || c=='G' || c=='T' || c=='N') s.baseCountGapped.merge(c, 1L, Long::sum);
                else s.baseCountGapped.merge('N', 1L, Long::sum);
            }

            long curContig = 0L;
            int nRunCount = 0;
            long localUngapped = 0L;
            for (int i = 0; i < len; i++) {
                char c = up.charAt(i);
                if (c=='A' || c=='C' || c=='G' || c=='T') {
                    curContig++;
                    localUngapped++;
                    s.baseCountUngapped.merge(c, 1L, Long::sum);
                } else {
                    if (curContig > 0) { s.contigLengths.add(curContig); s.totalUngapped += curContig; curContig = 0L; }
                    if (c=='N') { s.ambiguousBases++; if (i==0 || up.charAt(i-1) != 'N') nRunCount++; }
                    else { s.ambiguousBases++; if (i==0 || up.charAt(i-1) != 'N') nRunCount++; }
                }
            }
            if (curContig > 0) { s.contigLengths.add(curContig); s.totalUngapped += curContig; }
            s.gapRuns += nRunCount;

            String canonical = canonicalize(up);
            s.seqCounts.merge(canonical, 1, Integer::sum);

            s.entropySum += shannonEntropy(up);

            if (s.kmerBudget < maxKmersToCount) {
                long mask = (1L << (2 * k)) - 1L;
                long code = 0L;
                int valid = 0;
                for (int i = 0; i < len; i++) {
                    int b = baseToBits(up.charAt(i));
                    if (b >= 0) {
                        code = ((code << 2) | b) & mask;
                        valid++;
                        if (valid >= k) {
                            s.kmerCounts.merge(code, 1L, Long::sum);
                            s.kmerBudget++;
                            if (s.kmerBudget >= maxKmersToCount) break;
                        }
                    } else {
                        valid = 0; code = 0L;
                    }
                }
            }
        }
        return s;
    }

    private static Map<Character, Long> initBaseCountMap() {
        Map<Character, Long> m = new HashMap<>();
        m.put('A', 0L); m.put('C', 0L); m.put('G', 0L); m.put('T', 0L); m.put('N', 0L);
        return m;
    }

    private static int baseToBits(char c) {
        return switch (c) {
            case 'A' -> 0; case 'C' -> 1; case 'G' -> 2; case 'T' -> 3; default -> -1;
        };
    }

    private static String bitsToKmer(long code, int k) {
        char[] buf = new char[k];
        for (int i = k - 1; i >= 0; i--) {
            int v = (int) (code & 0b11);
            buf[i] = switch (v) { case 0 -> 'A'; case 1 -> 'C'; case 2 -> 'G'; case 3 -> 'T'; default -> 'N'; };
            code >>= 2;
        }
        return new String(buf);
    }

    private static String reverseComplement(String s) {
        char[] out = new char[s.length()];
        for (int i = 0, j = s.length()-1; i < s.length(); i++, j--) out[i] = complement(s.charAt(j));
        return new String(out);
    }

    private static char complement(char c) {
        return switch (c) { case 'A' -> 'T'; case 'T' -> 'A'; case 'C' -> 'G'; case 'G' -> 'C'; default -> 'N'; };
    }

    private static String canonicalize(String up) {
        String rc = reverseComplement(up);
        return up.compareTo(rc) <= 0 ? up :rc;
    }

    private static String canonicalizeSequence(String up) { return canonicalize(up); }

    private static String canonicalize(SequenceRecord rec) { return canonicalize(rec.getSequence().toUpperCase(Locale.ROOT)); }

    private static String canonicalizeByContent(String seq) { return canonicalize(seq.toUpperCase(Locale.ROOT)); }

    private static String canonicalizeSimple(String up) { return up; }

    private static String canonicalizeForDup(SequenceRecord rec) { return canonicalize(rec.getSequence().toUpperCase(Locale.ROOT)); }

    private static double shannonEntropy(String seq) {
        long[] counts = new long[4]; int total = 0;
        for (int i = 0; i < seq.length(); i++) { int b = baseToBits(seq.charAt(i)); if (b>=0) { counts[b]++; total++; } }
        if (total == 0) return 0.0;
        double H = 0.0; for (long c : counts) if (c>0) { double p = c/(double)total; H -= p*Math.log(p)/Math.log(2); }
        return H;
    }

    private static int computeNxx(List<Long> lengths, double fraction) {
        if (lengths == null || lengths.isEmpty()) return 0;
        long total = lengths.stream().mapToLong(Long::longValue).sum();
        long threshold = (long)Math.ceil(total * fraction);
        long cum = 0L;
        for (long l : lengths) { cum += l; if (cum >= threshold) return (int)Math.min(l, Integer.MAX_VALUE); }
        return 0;
    }

    private static int computeLxx(List<Long> lengths, double fraction) {
        if (lengths == null || lengths.isEmpty()) return 0;
        long total = lengths.stream().mapToLong(Long::longValue).sum();
        long threshold = (long)Math.ceil(total * fraction);
        long cum = 0L; int cnt = 0;
        for (long l : lengths) { cum += l; cnt++; if (cum >= threshold) return cnt; }
        return 0;
    }

    private static Map<String, Long> topKmersAsStrings(Map<Long, Long> kmers, int topN, int k) {
        return kmers.entrySet().stream()
                .sorted(Map.Entry.<Long, Long>comparingByValue(Comparator.reverseOrder()))
                .limit(topN)
                .collect(Collectors.toMap(e -> bitsToKmer(e.getKey(), k), Map.Entry::getValue, (a,b)->a, LinkedHashMap::new));
    }

    private static Map<Integer, Integer> buildLengthHistogram(List<Long> lengths) {
        Map<Integer, Integer> hist = new TreeMap<>();
        for (long l : lengths) {
            int key = (l < 1_000) ? (int) l : (l < 10_000 ? (int)(l/1000)*1000 : (l < 100_000 ? (int)(l/10_000)*10_000 : (l < 1_000_000 ? (int)(l/100_000)*100_000 : (int)(Math.min(l, Integer.MAX_VALUE)/1_000_000)*1_000_000)));
            hist.merge(key, 1, Integer::sum);
        }
        return hist;
    }

    private static <T> List<List<T>> partition(List<T> list, int parts) {
        List<List<T>> out = new ArrayList<>(); int n = list.size(); int base = n/parts; int rem = n%parts; int idx=0;
        for (int i=0;i<parts;i++) { int sz = base + (i<rem?1:0); List<T> sub = new ArrayList<>(); for (int j=0;j<sz;j++) sub.add(list.get(idx++)); out.add(sub); }
        return out;
    }

    private static class LocalStats {
        int sequenceCount = 0;
        long totalLength = 0L;
        long totalUngapped = 0L;
        long ambiguousBases = 0L;
        long gapRuns = 0L;
        Map<Character, Long> baseCountGapped = initBaseCountMap();
        Map<Character, Long> baseCountUngapped = initBaseCountMap();
        List<Long> scaffoldLengths = new ArrayList<>();
        List<Long> contigLengths = new ArrayList<>();
        Map<String, Integer> seqCounts = new HashMap<>();
        double entropySum = 0.0;
        Map<Long, Long> kmerCounts = new HashMap<>();
        int kmerBudget = 0;
    }
}
