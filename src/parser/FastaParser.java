package parser;

import model.SequenceRecord;
import java.io.*;
import java.util.NoSuchElementException;

/**
 * Streaming FASTA parser - reads one sequence at a time.
 */
public class FastaParser implements SequenceParser {
    private final BufferedReader reader;
    private String nextId;
    private StringBuilder nextSeq;
    private boolean finished;
    private SequenceRecord prefetched;

    public FastaParser(String filePath) throws IOException {
        this.reader = new BufferedReader(new FileReader(filePath), 1 << 16); // 64KB buffer
        this.nextSeq = new StringBuilder();
        this.finished = false;
        prefetch();
    }

    private void prefetch() {
        if (finished) return;
        try {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith(">")) {
                    if (nextId != null) {
                        prefetched = new SequenceRecord(nextId, nextSeq.toString());
                        nextId = line.substring(1).trim();
                        nextSeq.setLength(0);
                        return;
                    }
                    nextId = line.substring(1).trim();
                    nextSeq.setLength(0);
                } else {
                    nextSeq.append(line.trim());
                }
            }
            // End of file
            if (nextId != null) {
                prefetched = new SequenceRecord(nextId, nextSeq.toString());
                nextId = null;
            }
            finished = true;
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    @Override
    public boolean hasNext() {
        return prefetched != null;
    }

    @Override
    public SequenceRecord next() {
        if (prefetched == null) throw new NoSuchElementException();
        SequenceRecord result = prefetched;
        prefetched = null;
        prefetch();
        return result;
    }

    @Override
    public void close() {
        try {
            reader.close();
        } catch (IOException e) {
            // Ignore close errors
        }
    }
}
