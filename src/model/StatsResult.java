package model;

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
}
