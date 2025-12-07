import analysis.StatsCalculator;
import model.SequenceRecord;
import model.StatsResult;
import parser.FastaParser;
import parser.FastqParser;

import java.util.*;

public class Main {
    public static void main(String[] args) throws Exception {
        if (args.length < 2 || !args[0].equals("--input")) {
            System.out.println("Usage: --input <file");
            return;
        }

        String file = args[1];

        List<SequenceRecord> records;
        if (file.endsWith(".fa") || file.endsWith(".fasta")) {
            records = new FastaParser().parse(file);
        } else if (file.endsWith(".fq") || file.endsWith(".fastq")) {
            records = new FastqParser().parse(file);
        } else {
            throw new IllegalArgumentException("Unsupported file format");
        }

        StatsResult result = new StatsCalculator().compute(records);
        System.out.println(result);
    }
}
