import analysis.StatsCalculator;
import model.SequenceRecord;
import model.StatsResult;
import parser.FastaParser;
import parser.FastqParser;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

public class Main {
    public static void main(String[] args) throws Exception {
        if (args.length < 2 || !args[0].equals("--input")) {
            System.out.println("Usage: --input <file> [--output <file>]");
            return;
        }

        String inputFile = null;
        String outputFile = null;

        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "--input":
                    if (i + 1 < args.length) {
                        inputFile = args[++i];
                    } else {
                        System.out.println("Missing value for --input");
                        return;
                    }
                    break;
                case "--output":
                    if (i + 1 < args.length) {
                        outputFile = args[++i];
                    } else {
                        System.out.println("Missing value for --output");
                        return;
                    }
                    break;
                default:
            }
        }

        if (inputFile == null) {
            System.out.println("Usage: --input <file> [--output <file>]");
            return;
        }

        List<SequenceRecord> records;
        if (inputFile.endsWith(".fa") || inputFile.endsWith(".fasta")) {
            records = new FastaParser().parse(inputFile);
        } else if (inputFile.endsWith(".fq") || inputFile.endsWith(".fastq")) {
            records = new FastqParser().parse(inputFile);
        } else {
            throw new IllegalArgumentException("Unsupported file format");
        }

        if (records.isEmpty()) {
            System.out.println("No sequence found in file.");
            return;
        }

        StatsResult result = new StatsCalculator().compute(records);
        System.out.println(result);

        if (outputFile != null) {
            String json = result.toJson();
            Files.writeString(Path.of(outputFile), json);
            System.out.println("Wrote JSON to " + outputFile);
        }
    }
}
