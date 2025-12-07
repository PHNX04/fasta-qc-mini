package model;

public class SequenceRecord {
    private final String id;
    private final String sequence;

    public SequenceRecord(String id, String sequence) {
        this.id = id;
        this.sequence = sequence;
    }

    public String getId() { return id; }

    public String getSequence() { return sequence; }

    public int length() { return sequence.length(); }
}
