package analysis;

import model.SequenceRecord;
import model.StatsResult;
import parser.SequenceParser;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Streaming stats calculator - processes sequences one at a time.
 * Memory usage: O(unique_lengths + hash_set_size) instead of O(all_sequences).
 */
public class StatsCalculator {
    private final int k = 5;
    private final int maxKmersToCount;

    public StatsCalculator() {
        this(100_000);
    }

    public StatsCalculator(int threads, int maxKmersToCount) {
        // threads parameter kept for API compatibility, but streaming is single-pass
        this.maxKmersToCount = Math.max(1, maxKmersToCount);
    }

    public StatsCalculator(int maxKmersToCount) {
        this.maxKmersToCount = Math.max(1, maxKmersToCount);
    }

    public StatsResult compute(SequenceParser parser) {
        StatsResult result = new StatsResult();

        // Counters - O(1) memory
        int sequenceCount = 0;
        long totalLength = 0L;
        long totalUngapped = 0L;
        long ambiguousBases = 0L;
        long gapRuns = 0L;

        // Base counts - fixed size maps
        long[] baseCounts = new long[5]; // A, C, G, T, N
        long[] baseCountsUngapped = new long[4]; // A, C, G, T

        // Lengths - grows with sequence count but uses primitives
        LongList scaffoldLengths = new LongList();
        LongList contigLengths = new LongList();

        // Hash-based duplicate detection - O(unique_sequences) but only 8 bytes per
        // sequence
        Set<Long> seenHashes = new HashSet<>();
        int duplicateCount = 0;

        // Entropy accumulator
        double entropySum = 0.0;

        // K-mer counting with budget
        Map<Long, Long> kmerCounts = new HashMap<>();
        int kmerBudget = 0;
        long kmerMask = (1L << (2 * k)) - 1L;

        // Min/max tracking
        int minLength = Integer.MAX_VALUE;
        int maxLength = 0;

        while (parser.hasNext()) {
            SequenceRecord rec = parser.next();
            String seq = rec.getSequence();
            if (seq == null)
                seq = "";
            String up = seq.toUpperCase(Locale.ROOT);
            int len = up.length();

            sequenceCount++;
            totalLength += len;
            scaffoldLengths.add(len);

            if (len < minLength)
                minLength = len;
            if (len > maxLength)
                maxLength = len;

            // Count bases
            for (int i = 0; i < len; i++) {
                int idx = baseToIndex(up.charAt(i));
                baseCounts[idx]++;
            }

            // Process contigs and ungapped bases
            long curContig = 0L;
            int nRunCount = 0;
            for (int i = 0; i < len; i++) {
                char c = up.charAt(i);
                int b = baseToBits(c);
                if (b >= 0) {
                    curContig++;
                    baseCountsUngapped[b]++;
                } else {
                    if (curContig > 0) {
                        contigLengths.add(curContig);
                        totalUngapped += curContig;
                        curContig = 0L;
                    }
                    ambiguousBases++;
                    if (i == 0 || up.charAt(i - 1) != 'N')
                        nRunCount++;
                }
            }
            if (curContig > 0) {
                contigLengths.add(curContig);
                totalUngapped += curContig;
            }
            gapRuns += nRunCount;

            // Hash-based duplicate detection (64-bit hash)
            long hash = computeHash(up);
            if (!seenHashes.add(hash)) {
                duplicateCount++;
            }

            // Shannon entropy
            entropySum += shannonEntropy(up);

            // K-mer counting with budget limit
            if (kmerBudget < maxKmersToCount) {
                long code = 0L;
                int valid = 0;
                for (int i = 0; i < len && kmerBudget < maxKmersToCount; i++) {
                    int b = baseToBits(up.charAt(i));
                    if (b >= 0) {
                        code = ((code << 2) | b) & kmerMask;
                        valid++;
                        if (valid >= k) {
                            kmerCounts.merge(code, 1L, Long::sum);
                            kmerBudget++;
                        }
                    } else {
                        valid = 0;
                        code = 0L;
                    }
                }
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

        // Sort lengths for N50 calculation
        long[] scaffoldArr = scaffoldLengths.toSortedDescending();
        long[] contigArr = contigLengths.toSortedDescending();

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

        long gcGapped = baseCounts[1] + baseCounts[2]; // C + G
        result.gcContent = totalBases == 0 ? 0.0 : gcGapped / (double) totalBases;

        long totalUngappedBases = 0;
        for (long c : baseCountsUngapped)
            totalUngappedBases += c;
        long gcUngapped = baseCountsUngapped[1] + baseCountsUngapped[2];
        result.gcContentUngapped = totalUngappedBases == 0 ? 0.0 : gcUngapped / (double) totalUngappedBases;

        // Duplicates - using hash collision approximation
        result.duplicateSequenceCount = duplicateCount;
        result.duplicateFraction = sequenceCount == 0 ? 0.0 : duplicateCount / (double) sequenceCount;

        result.avgShannonEntropy = sequenceCount == 0 ? 0.0 : entropySum / sequenceCount;

        result.kmerCounts = topKmersAsStrings(kmerCounts, 100, k);
        result.lengthHistogram = buildLengthHistogram(scaffoldArr);

        return result;
    }

    // Hash function for sequence deduplication (FNV-1a inspired)
    private static long computeHash(String seq) {
        long hash = 0xcbf29ce484222325L;
        for (int i = 0; i < seq.length(); i++) {
            hash ^= seq.charAt(i);
            hash *= 0x100000001b3L;
        }
        return hash;
    }

    private static int baseToIndex(char c) {
        return switch (c) {
            case 'A' -> 0;
            case 'C' -> 1;
            case 'G' -> 2;
            case 'T' -> 3;
            default -> 4; // N and others
        };
    }

    private static int baseToBits(char c) {
        return switch (c) {
            case 'A' -> 0;
            case 'C' -> 1;
            case 'G' -> 2;
            case 'T' -> 3;
            default -> -1;
        };
    }

    private static String bitsToKmer(long code, int k) {
        char[] buf = new char[k];
        for (int i = k - 1; i >= 0; i--) {
            int v = (int) (code & 0b11);
            buf[i] = switch (v) {
                case 0 -> 'A';
                case 1 -> 'C';
                case 2 -> 'G';
                case 3 -> 'T';
                default -> 'N';
            };
            code >>= 2;
        }
        return new String(buf);
    }

    private static double shannonEntropy(String seq) {
        long[] counts = new long[4];
        int total = 0;
        for (int i = 0; i < seq.length(); i++) {
            int b = baseToBits(seq.charAt(i));
            if (b >= 0) {
                counts[b]++;
                total++;
            }
        }
        if (total == 0)
            return 0.0;
        double H = 0.0;
        for (long c : counts) {
            if (c > 0) {
                double p = c / (double) total;
                H -= p * Math.log(p) / Math.log(2);
            }
        }
        return H;
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
     * Simple primitive long list to avoid boxed Long objects.
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

        long[] toSortedDescending() {
            long[] result = Arrays.copyOf(data, size);
            Arrays.sort(result);
            // Reverse for descending order
            for (int i = 0, j = result.length - 1; i < j; i++, j--) {
                long tmp = result[i];
                result[i] = result[j];
                result[j] = tmp;
            }
            return result;
        }
    }
}
