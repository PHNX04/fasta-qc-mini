package model;

import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;

public class StatsResult {
    public int sequenceCount;
    public int minLength;
    public int maxLength;
    public double avgLength;
    public Map<Character, Double> baseFrequencies;
    public double gcContent;

    @Override
    public String toString() {
        return "Sequences: " + sequenceCount + "\n" +
               "Min length: " + minLength + "\n" +
               "Max length: " + maxLength + "\n" +
               "Avg length: " + avgLength + "\n" +
               "GC content: " + gcContent + "\n" +
               "Base frequencies: " + baseFrequencies;
    }

    public String toJson() {
        StringBuilder sb = new StringBuilder();
        sb.append("{\n");
        sb.append(" \"sequenceCount\": ").append(sequenceCount).append(",\n");
        sb.append(" \"minLength\": ").append(minLength).append(",\n");
        sb.append(" \"maxLength\": ").append(maxLength).append(",\n");
        sb.append(String.format(Locale.ROOT, " \"avgLength\": %.8f,\n", avgLength));
        sb.append(String.format(Locale.ROOT, " \"gcContent\": %.8f,\n", gcContent));

        Map<Character, Double> ordered = new LinkedHashMap<>();
        ordered.put('A', baseFrequencies.getOrDefault('A', 0.0));
        ordered.put('C', baseFrequencies.getOrDefault('C', 0.0));
        ordered.put('G', baseFrequencies.getOrDefault('G', 0.0));
        ordered.put('T', baseFrequencies.getOrDefault('T', 0.0));
        ordered.put('N', baseFrequencies.getOrDefault('N', 0.0));

        sb.append(" \"baseFrequencies\": {");
        boolean first = true;
        for (Map.Entry<Character, Double> e : ordered.entrySet()) {
            if (!first) sb.append(", ");
            sb.append("\"").append(e.getKey()).append("\": ");
            sb.append(String.format(Locale.ROOT, "%.8f", e.getValue()));
            first = false;
        }
        sb.append("}\n");
        sb.append("}");
        return sb.toString();
    }
}
