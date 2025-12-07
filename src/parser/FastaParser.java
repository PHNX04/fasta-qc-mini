package parser;

import model.SequenceRecord;
import java.util.*;
import java.io.*;

public class FastaParser {

    public List<SequenceRecord> parse(String filePath) throws IOException {
        List<SequenceRecord> records = new ArrayList<>();
        BufferedReader br = new BufferedReader(new FileReader(filePath));

        String line;
        String currentId = null;
        StringBuilder currentSeq = new StringBuilder();

        while ((line = br.readLine()) != null) {
            if (line.startsWith(">")) {
                if (currentId != null) {
                    records.add(new SequenceRecord(currentId, currentSeq.toString()));
                }
                currentId = line.substring(1).trim();
                currentSeq = new StringBuilder();
            } else {
                currentSeq.append(line.trim());
            }
        }

        if (currentId != null) {
            records.add(new SequenceRecord(currentId, currentSeq.toString()));
        }

        return records;
    }
}
