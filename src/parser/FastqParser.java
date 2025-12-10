package parser;

import model.SequenceRecord;
import java.io.*;
import java.util.NoSuchElementException;

/**
 * Streaming FASTQ parser - reads one sequence at a time.
 */
public class FastqParser implements SequenceParser {
    private final BufferedReader reader;
    private SequenceRecord prefetched;
    private boolean finished;

    public FastqParser(String filePath) throws IOException {
        this.reader = new BufferedReader(new FileReader(filePath), 1 << 16); // 64KB buffer
        this.finished = false;
        prefetch();
    }

    private void prefetch() {
        if (finished) return;
        try {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("@")) {
                    String id = line.substring(1).trim();
                    String seq = reader.readLine();
                    reader.readLine(); // + line
                    reader.readLine(); // quality line
                    
                    if (seq != null) {
                        prefetched = new SequenceRecord(id, seq.trim());
                        return;
                    }
                }
            }
            finished = true;
            prefetched = null;
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
