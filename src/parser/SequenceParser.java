package parser;

import model.SequenceRecord;
import java.io.Closeable;
import java.util.Iterator;

/**
 * Streaming parser interface for FASTA/FASTQ files.
 * Processes sequences one-at-a-time to minimize memory usage.
 */
public interface SequenceParser extends Iterator<SequenceRecord>, Closeable {
    
    @Override
    boolean hasNext();
    
    @Override
    SequenceRecord next();
    
    @Override
    void close();
}
