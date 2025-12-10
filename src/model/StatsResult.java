package model;

import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;

public class StatsResult {
    public int sequenceCount;
    public long totalLength;
    public long totalUngappedLength;
    public long ambiguousBaseCount;
    public long gapCount;

    public int minLength;
    public int maxLength;
    public double avgLength;

    public Map<Character, Double> baseFrequencies;
    public double gcContent;
    public double gcContentUngapped;

    public int scaffoldN50, scaffoldN90, scaffoldL50, scaffoldL90;
    public int contigN50, contigN90, contigL50, contigL90;

    public int duplicateSequenceCount;
    public double duplicateFraction;

    public double avgShannonEntropy;

    public Map<Integer, Integer> lengthHistogram;

    public Map<String, Long> kmerCounts;

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Sequences: ").append(sequenceCount).append("\n");
        sb.append("Min length: ").append(minLength).append("\n");
        sb.append("Max length: ").append(maxLength).append("\n");
        sb.append(String.format(Locale.ROOT, "Avg length: %.8f", avgLength)).append("\n");
        sb.append(String.format(Locale.ROOT, "GC content: %.8f", gcContent)).append("\n");
        sb.append("Base frequencies: ").append(baseFrequencies).append("\n");
        sb.append("N50: ").append(scaffoldN50).append(" N90: ").append(scaffoldN90).append("\n");
        sb.append("Duplicate sequences: ").append(duplicateSequenceCount).append(" Fraction: ").append(String.format(Locale.ROOT, "%.8f", duplicateFraction)).append("\n");
        sb.append("Ambiguous bases (N): ").append(ambiguousBaseCount).append(" Fraction: ").append(String.format(Locale.ROOT, "%.8f", (totalLength==0?0.0:ambiguousBaseCount/(double)totalLength))).append("\n");
        sb.append("Length histogram (sample): ").append(lengthHistogram).append("\n");
        sb.append("Top k-mers (sample): ").append(kmerCounts).append("\n");
        return sb.toString();
    }

    public String toJson() {
        StringBuilder sb = new StringBuilder();
        sb.append("{\n");
        sb.append(" \"sequenceCount\": ").append(sequenceCount).append(",\n");
        sb.append(" \"minLength\": ").append(minLength).append(",\n");
        sb.append(" \"maxLength\": ").append(maxLength).append(",\n");
        sb.append(String.format(Locale.ROOT, " \"avgLength\": %.8f,\n", avgLength));
        sb.append(String.format(Locale.ROOT, " \"gcContent\": %.8f,\n", gcContent));
        sb.append(String.format(Locale.ROOT, " \"gcContentUngapped\": %.8f,\n", gcContentUngapped));

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
        sb.append("},\n");

        sb.append(" \"n50\": ").append(scaffoldN50).append(", \"n90\": ").append(scaffoldN90).append(",\n");
        sb.append(" \"duplicateFraction\": ").append(String.format(Locale.ROOT, "%.8f", duplicateFraction)).append(",\n");
        sb.append(" \"duplicateSequenceCount\": ").append(duplicateSequenceCount).append(",\n");
        sb.append(" \"ambiguousBaseCount\": ").append(ambiguousBaseCount).append(",\n");
        sb.append(String.format(Locale.ROOT, " \"ambiguousFraction\": %.8f, ", (totalLength==0?0.0:ambiguousBaseCount/(double)totalLength))).append("\n");

        sb.append(" \"lengthHistogram\": {");
        first = true;
        for (var e : lengthHistogram.entrySet()) {
            if (!first) sb.append(", ");
            sb.append("\"").append(e.getKey()).append("\": ").append(e.getValue());
            first = false;
        }
        sb.append("},\n");

        sb.append(" \"topKmers\": {");
        first = true;
        for (var e : kmerCounts.entrySet()) {
            if (!first) sb.append(", ");
            sb.append("\"").append(e.getKey()).append("\": ").append(e.getValue());
            first = false;
        }
        sb.append("}\n");

        sb.append("}");
        return sb.toString();
    }
}
