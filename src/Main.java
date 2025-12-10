import analysis.StatsCalculator;
import model.StatsResult;
import parser.FastaParser;
import parser.FastqParser;
import parser.SequenceParser;

import java.nio.file.Files;
import java.nio.file.Path;

public class Main {
    public static void main(String[] args) throws Exception {
        if (args.length < 2 || !args[0].equals("--input")) {
            System.out.println("Usage: --input <file> [--output <file>] [--kmers N]");
            return;
        }

        String inputFile = null;
        String outputFile = null;
        int maxKmers = 100_000;

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
                case "--kmers":
                    if (i + 1 < args.length) {
                        maxKmers = Integer.parseInt(args[++i]);
                    }
                    break;
                case "--threads":
                    // Kept for backwards compatibility but ignored (streaming is single-pass)
                    if (i + 1 < args.length)
                        i++;
                    break;
                default:
            }
        }

        if (inputFile == null) {
            System.out.println("Usage: --input <file> [--output <file>]");
            return;
        }

        StatsCalculator calc = new StatsCalculator(maxKmers);

        // Use try-with-resources for automatic cleanup
        try (SequenceParser parser = createParser(inputFile)) {
            StatsResult result = calc.compute(parser);

            if (result.sequenceCount == 0) {
                System.out.println("No sequence found in file.");
                return;
            }

            System.out.println(result);

            if (outputFile != null) {
                String json = result.toJson();
                Files.writeString(Path.of(outputFile), json);
                System.out.println("Wrote JSON to " + outputFile);
            }
        }
    }

    private static SequenceParser createParser(String inputFile) throws Exception {
        if (inputFile.endsWith(".fa") || inputFile.endsWith(".fasta") || inputFile.endsWith(".fna")) {
            return new FastaParser(inputFile);
        } else if (inputFile.endsWith(".fq") || inputFile.endsWith(".fastq")) {
            return new FastqParser(inputFile);
        } else {
            throw new IllegalArgumentException("Unsupported file format. Use .fa/.fasta/.fna or .fq/.fastq");
        }
    }
}
