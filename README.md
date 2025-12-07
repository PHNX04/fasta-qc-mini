# FASTA-QC-MINI

**FASTA-QC-MINI** is a lightweight command-line tool for **basic quality control of DNA/RNA sequences** in FASTA or FASTQ format.
It calculates key metrics such as sequence count, length statistics, base composition, and GC content.

This project is ideal for small bioinformatics pipelines, practicing sequence analysis, or as a professional portfolio project.

---

## Features 

- Supports FASTA and FASTQ files
- Computes:
  - Number of sequences / reads
  - Minimum, maximum, and average sequence length
  - Base frequencies (A, C, G, T, N)
  - Overall GC content
- Console output and optional JSON export (incoming)
- Modular, clean Java code structure
- Easily extendable (e.g., histograms, FASTQ quality score analysis)

---

## Installation & Compilation

### Requirements:
- Java 21 or higher

### Compile:
```bash
cd fasta-qc-mini/src
javac *.java
```

### Run:
- for FASTA:
```bash
cd out/production/fasta-qc-mini
java Main --input ../test.fasta
```
- for FASTQ
```bash
cd out/production/fasta-qc-mini
java Main --input ../test.fastq
```

## Example Output

```yaml
Sequences: 4
Min length: 32
Max length: 40
Avg length: 35.75
GC content: 0.34965034965034963
Base frequencies: {A=0.32167832167832167, C=0.17482517482517482, T=0.3076923076923077, G=0.17482517482517482, N=0.02097902097902098}
```

--- 
## Fyi
Feel free to use and extend the code.

It is a working progress and can take a while between updates.

I will try to improve the usability.