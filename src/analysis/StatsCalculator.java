package analysis;

import model.SequenceRecord;
import model.StatsResult;

import java.util.*;

public class StatsCalculator {
    public StatsResult compute(List<SequenceRecord> records) {
        StatsResult result = new StatsResult();

        int totalLength = 0;
        int min = Integer.MAX_VALUE;
        int max = Integer.MIN_VALUE;

        Map<Character, Integer> baseCount = new HashMap<>();
        baseCount.put('A', 0);
        baseCount.put('C', 0);
        baseCount.put('G', 0);
        baseCount.put('T', 0);
        baseCount.put('N', 0);

        for (SequenceRecord rec : records) {
            int len = rec.length();
            totalLength += len;

            min = Math.min(min, len);
            max = Math.max(max, len);

            for (char c : rec.getSequence().toUpperCase().toCharArray()) {
                baseCount.put(c, baseCount.getOrDefault(c, 0) + 1);
            }
        }

        result.sequenceCount = records.size();
        result.minLength = min;
        result.maxLength = max;
        result.avgLength = (double) totalLength / records.size();

        int totalBases = baseCount.values().stream().mapToInt(Integer::intValue).sum();
        Map<Character, Double> freq = new HashMap<>();
        for (var e : baseCount.entrySet()) {
            freq.put(e.getKey(), e.getValue() / (double) totalBases);
        }
        result.baseFrequencies = freq;

        double gc = (baseCount.get('G') + baseCount.get('C')) / (double) totalBases;
        result.gcContent = gc;

        return result;
    }
}
