package analysis;

import model.SequenceRecord;
import model.StatsResult;
import parser.SequenceParser;

import java.util.*;
import java.util.stream.Collectors;

/**
 * High-performance streaming stats calculator.
 * Optimizations:
 * - Lookup tables instead of switch statements
 * - Single-pass processing per sequence
 * - Parallel sorting for large arrays
 * - No string case conversion (lookup handles both cases)
 */
public class StatsCalculator {
    private final int k = 5;
    private final int maxKmersToCount;

    // Lookup tables for O(1) character classification
    private static final int[] BASE_INDEX = new int[256]; // char -> 0-4 (A,C,G,T,N)
    private static final int[] BASE_BITS = new int[256]; // char -> 0-3 or -1
    private static final long[] HASH_MULT = new long[256]; // precomputed hash multipliers

    static {
        // Initialize lookup tables
        Arrays.fill(BASE_INDEX, 4); // Default = N
        Arrays.fill(BASE_BITS, -1); // Default = invalid

        // Map both upper and lower case
        BASE_INDEX['A'] = BASE_INDEX['a'] = 0;
        BASE_INDEX['C'] = BASE_INDEX['c'] = 1;
        BASE_INDEX['G'] = BASE_INDEX['g'] = 2;
        BASE_INDEX['T'] = BASE_INDEX['t'] = 3;

        BASE_BITS['A'] = BASE_BITS['a'] = 0;
        BASE_BITS['C'] = BASE_BITS['c'] = 1;
        BASE_BITS['G'] = BASE_BITS['g'] = 2;
        BASE_BITS['T'] = BASE_BITS['t'] = 3;

        // Precompute hash XOR values (FNV-1a inspired)
        for (int i = 0; i < 256; i++) {
            HASH_MULT[i] = i;
        }
    }

    public StatsCalculator() {
        this(100_000);
    }

    public StatsCalculator(int threads, int maxKmersToCount) {
        this.maxKmersToCount = Math.max(1, maxKmersToCount);
    }

    public StatsCalculator(int maxKmersToCount) {
        this.maxKmersToCount = Math.max(1, maxKmersToCount);
    }

    public StatsResult compute(SequenceParser parser) {
        StatsResult result = new StatsResult();

        // Counters
        int sequenceCount = 0;
        long totalLength = 0L;
        long totalUngapped = 0L;
        long ambiguousBases = 0L;
        long gapRuns = 0L;

        // Base counts - fixed size arrays
        long[] baseCounts = new long[5]; // A, C, G, T, N
        long[] baseCountsUngapped = new long[4]; // A, C, G, T

        // Lengths
        LongList scaffoldLengths = new LongList();
        LongList contigLengths = new LongList();

        // Hash-based duplicate detection
        Set<Long> seenHashes = new HashSet<>();
        int duplicateCount = 0;

        // Entropy - precompute log2 divisor
        final double LOG2 = Math.log(2);
        double entropySum = 0.0;

        // K-mer counting
        Map<Long, Long> kmerCounts = new HashMap<>();
        int kmerBudget = 0;
        final long kmerMask = (1L << (2 * k)) - 1L;

        // Min/max
        int minLength = Integer.MAX_VALUE;
        int maxLength = 0;

        while (parser.hasNext()) {
            SequenceRecord rec = parser.next();
            String seq = rec.getSequence();
            if (seq == null)
                seq = "";
            int len = seq.length();

            sequenceCount++;
            totalLength += len;
            scaffoldLengths.add(len);

            if (len < minLength)
                minLength = len;
            if (len > maxLength)
                maxLength = len;

            // === SINGLE-PASS PROCESSING ===
            long hash = 0xcbf29ce484222325L;
            long kmerCode = 0L;
            int kmerValid = 0;
            long curContig = 0L;
            char prevChar = 0;

            // Entropy accumulators
            long[] entropyBases = new long[4];
            int entropyTotal = 0;

            for (int i = 0; i < len; i++) {
                char c = seq.charAt(i);
                int idx = BASE_INDEX[c & 0xFF];
                int bits = BASE_BITS[c & 0xFF];

                // Base counting
                baseCounts[idx]++;

                // Hash computation (FNV-1a)
                hash ^= (c | 0x20); // lowercase
                hash *= 0x100000001b3L;

                if (bits >= 0) {
                    // Valid base (A, C, G, T)
                    curContig++;
                    baseCountsUngapped[bits]++;
                    entropyBases[bits]++;
                    entropyTotal++;

                    // K-mer counting
                    if (kmerBudget < maxKmersToCount) {
                        kmerCode = ((kmerCode << 2) | bits) & kmerMask;
                        kmerValid++;
                        if (kmerValid >= k) {
                            kmerCounts.merge(kmerCode, 1L, Long::sum);
                            kmerBudget++;
                        }
                    }
                } else {
                    // Invalid base (N or other)
                    if (curContig > 0) {
                        contigLengths.add(curContig);
                        totalUngapped += curContig;
                        curContig = 0L;
                    }
                    ambiguousBases++;

                    char prevUpper = (char) (prevChar & ~0x20);
                    if (prevUpper != 'N') {
                        gapRuns++;
                    }

                    // Reset k-mer state
                    kmerValid = 0;
                    kmerCode = 0L;
                }
                prevChar = c;
            }

            // Finalize contig
            if (curContig > 0) {
                contigLengths.add(curContig);
                totalUngapped += curContig;
            }

            // Duplicate detection
            if (!seenHashes.add(hash)) {
                duplicateCount++;
            }

            // Shannon entropy (computed from single-pass accumulators)
            if (entropyTotal > 0) {
                double H = 0.0;
                for (long cnt : entropyBases) {
                    if (cnt > 0) {
                        double p = cnt / (double) entropyTotal;
                        H -= p * Math.log(p) / LOG2;
                    }
                }
                entropySum += H;
            }
        }

        // Build result
        result.sequenceCount = sequenceCount;
        result.totalLength = totalLength;
        result.totalUngappedLength = totalUngapped;
        result.ambiguousBaseCount = ambiguousBases;
        result.gapCount = gapRuns;

        if (sequenceCount == 0) {
            result.minLength = 0;
            result.maxLength = 0;
        } else {
            result.minLength = minLength;
            result.maxLength = maxLength;
        }
        result.avgLength = sequenceCount == 0 ? 0.0 : (double) totalLength / sequenceCount;

        // Parallel sort for N50 calculation (faster on large arrays)
        long[] scaffoldArr = scaffoldLengths.toSortedDescendingParallel();
        long[] contigArr = contigLengths.toSortedDescendingParallel();

        result.scaffoldN50 = computeNxx(scaffoldArr, 0.5);
        result.scaffoldN90 = computeNxx(scaffoldArr, 0.9);
        result.scaffoldL50 = computeLxx(scaffoldArr, 0.5);
        result.scaffoldL90 = computeLxx(scaffoldArr, 0.9);

        result.contigN50 = computeNxx(contigArr, 0.5);
        result.contigN90 = computeNxx(contigArr, 0.9);
        result.contigL50 = computeLxx(contigArr, 0.5);
        result.contigL90 = computeLxx(contigArr, 0.9);

        // Base frequencies
        long totalBases = 0;
        for (long c : baseCounts)
            totalBases += c;
        result.baseFrequencies = new LinkedHashMap<>();
        result.baseFrequencies.put('A', totalBases == 0 ? 0.0 : baseCounts[0] / (double) totalBases);
        result.baseFrequencies.put('C', totalBases == 0 ? 0.0 : baseCounts[1] / (double) totalBases);
        result.baseFrequencies.put('G', totalBases == 0 ? 0.0 : baseCounts[2] / (double) totalBases);
        result.baseFrequencies.put('T', totalBases == 0 ? 0.0 : baseCounts[3] / (double) totalBases);
        result.baseFrequencies.put('N', totalBases == 0 ? 0.0 : baseCounts[4] / (double) totalBases);

        long gcGapped = baseCounts[1] + baseCounts[2];
        result.gcContent = totalBases == 0 ? 0.0 : gcGapped / (double) totalBases;

        long totalUngappedBases = 0;
        for (long c : baseCountsUngapped)
            totalUngappedBases += c;
        long gcUngapped = baseCountsUngapped[1] + baseCountsUngapped[2];
        result.gcContentUngapped = totalUngappedBases == 0 ? 0.0 : gcUngapped / (double) totalUngappedBases;

        result.duplicateSequenceCount = duplicateCount;
        result.duplicateFraction = sequenceCount == 0 ? 0.0 : duplicateCount / (double) sequenceCount;

        result.avgShannonEntropy = sequenceCount == 0 ? 0.0 : entropySum / sequenceCount;

        result.kmerCounts = topKmersAsStrings(kmerCounts, 100, k);
        result.lengthHistogram = buildLengthHistogram(scaffoldArr);

        return result;
    }

    private static String bitsToKmer(long code, int k) {
        char[] buf = new char[k];
        for (int i = k - 1; i >= 0; i--) {
            int v = (int) (code & 0b11);
            buf[i] = "ACGT".charAt(v);
            code >>= 2;
        }
        return new String(buf);
    }

    private static int computeNxx(long[] lengths, double fraction) {
        if (lengths == null || lengths.length == 0)
            return 0;
        long total = 0;
        for (long l : lengths)
            total += l;
        long threshold = (long) Math.ceil(total * fraction);
        long cum = 0L;
        for (long l : lengths) {
            cum += l;
            if (cum >= threshold)
                return (int) Math.min(l, Integer.MAX_VALUE);
        }
        return 0;
    }

    private static int computeLxx(long[] lengths, double fraction) {
        if (lengths == null || lengths.length == 0)
            return 0;
        long total = 0;
        for (long l : lengths)
            total += l;
        long threshold = (long) Math.ceil(total * fraction);
        long cum = 0L;
        int cnt = 0;
        for (long l : lengths) {
            cum += l;
            cnt++;
            if (cum >= threshold)
                return cnt;
        }
        return 0;
    }

    private static Map<String, Long> topKmersAsStrings(Map<Long, Long> kmers, int topN, int k) {
        return kmers.entrySet().stream()
                .sorted(Map.Entry.<Long, Long>comparingByValue(Comparator.reverseOrder()))
                .limit(topN)
                .collect(Collectors.toMap(
                        e -> bitsToKmer(e.getKey(), k),
                        Map.Entry::getValue,
                        (a, b) -> a,
                        LinkedHashMap::new));
    }

    private static Map<Integer, Integer> buildLengthHistogram(long[] lengths) {
        Map<Integer, Integer> hist = new TreeMap<>();
        for (long l : lengths) {
            int key;
            if (l < 1_000)
                key = (int) l;
            else if (l < 10_000)
                key = (int) (l / 1000) * 1000;
            else if (l < 100_000)
                key = (int) (l / 10_000) * 10_000;
            else if (l < 1_000_000)
                key = (int) (l / 100_000) * 100_000;
            else
                key = (int) (Math.min(l, Integer.MAX_VALUE) / 1_000_000) * 1_000_000;
            hist.merge(key, 1, Integer::sum);
        }
        return hist;
    }

    /**
     * Primitive long list with parallel sort support.
     */
    private static class LongList {
        private long[] data;
        private int size;

        LongList() {
            this.data = new long[1024];
            this.size = 0;
        }

        void add(long value) {
            if (size >= data.length) {
                data = Arrays.copyOf(data, data.length * 2);
            }
            data[size++] = value;
        }

        long[] toSortedDescendingParallel() {
            long[] result = Arrays.copyOf(data, size);
            // Use parallel sort for large arrays
            if (result.length > 10_000) {
                Arrays.parallelSort(result);
            } else {
                Arrays.sort(result);
            }
            // Reverse for descending
            for (int i = 0, j = result.length - 1; i < j; i++, j--) {
                long tmp = result[i];
                result[i] = result[j];
                result[j] = tmp;
            }
            return result;
        }
    }
}
