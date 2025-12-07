package parser;

import model.SequenceRecord;
import java.io.*;
import java.util.*;

public class FastqParser {
    public List<SequenceRecord> parse(String filePath) throws IOException {
        List<SequenceRecord> records = new ArrayList<>();
        BufferedReader br = new BufferedReader(new FileReader(filePath));
        String line;

        while ((line = br.readLine()) != null) {
            if (!line.startsWith("@")) continue;

            String id = line.substring(1).trim();
            String seq = br.readLine().trim();
            br.readLine();
            br.readLine();

            records.add(new SequenceRecord(id, seq));
        }
        return records;
    }
}
