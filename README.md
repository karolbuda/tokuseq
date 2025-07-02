# TokuSeq pipeline documentation

The TokuSeq pipeline is inspired by levSeq from the Arnold lab (PUT LINK HERE), but exclusively uses python scripts and packages.

This documentation explains the usage of the TokuSeq pipeline for demultiplexing sequencing reads and analyzing amino acid mutations. The pipeline consists of two main scripts:

1. **`tokuseq.py`** - Demultiplexes nanopore reads and performs consensus sequence generation based on reference DNA
2. **`amino_acid_analysis.py`** - Analyzes amino acid mutations from consensus sequences based on a reference ORF

## Overview

The TokuSeq pipeline is designed to process long-read sequencing data (e.g., Oxford Nanopore) from 96-well plate experiments with dual barcoding. It extracts sequences between forward and reverse barcodes, aligns them to a reference, generates consensus sequences, and identifies mutations at both DNA and amino acid levels.

## Prerequisites

### Required Python Packages
- `biopython` - For sequence manipulation and file I/O
- `mappy` - For minimap2 alignment
- `pandas` - For data manipulation
- `tempfile`, `subprocess` - For running external tools

### External Tools
- **MAFFT** - Multiple sequence alignment tool (must be installed and available in PATH)
- **minimap2** - For sequence alignment (installed via mappy)

### Installation
```bash
pip install biopython mappy pandas
# Install MAFFT separately based on your system
```

## File Structure Requirements

Your project directory should contain:
```
project_directory/
├── barcodes.fasta          # Forward and reverse barcodes
├── reference.fasta         # Full reference sequence
├── reference_orf.fasta     # ORF region for amino acid analysis
└── input_reads.fastq       # Raw sequencing reads
```

## Script 1: tokuseq.py

### Purpose
Demultiplexes sequencing reads based on dual barcoding system and generates consensus sequences for each well.

### Usage
```bash
python tokuseq.py <reference_directory> <min_read_length> <max_read_length>
```

### Parameters
- `<reference_directory>`: Path to directory containing input files
- `<min_read_length>`: Minimum read length to include (e.g., 2000)
- `<max_read_length>`: Maximum read length to include (e.g., 3000)

### Example
```bash
python tokuseq.py /path/to/project/wt-p1 2000 3000
```

### Input Files

#### 1. barcodes.fasta
Contains forward and reverse barcodes for demultiplexing:
```fasta
>fwd_01
ATCGATCGATCG
>fwd_02
GCTAGCTAGCTA
...
>rev_01
TTGGAACCTTGG
>rev_02
AACCGGTTAACC
...
```

#### 2. reference.fasta
Reference sequence for alignment:
```fasta
>reference
ATCGATCGATCG...
```

#### 3. input_reads.fastq
Raw sequencing reads in FASTQ format.

### Output

#### Directory Structure
```
project_directory/
└── extracted_fastq_YYMMDD/
    ├── A1/
    │   ├── A1.fastq           # Demultiplexed reads for well A1
    │   └── A1_alignment.fasta # Multiple sequence alignment
    ├── A2/
    │   ├── A2.fastq
    │   └── A2_alignment.fasta
    └── ...
```

#### Summary File
`demultiplex_summary.csv` contains:
- **Well**: Well identifier (A1-H12)
- **Consensus**: Consensus sequence for the well
- **Mutations**: DNA-level mutations (e.g., "A123T;G456C" or "WT")
- **ReadCount**: Number of reads assigned to the well

### Algorithm Details

1. **Barcode Detection**: Searches for forward and reverse barcodes in reads
2. **Sequence Extraction**: Extracts sequence between barcodes
3. **Orientation Correction**: Handles forward and reverse orientations
4. **Reference Alignment**: Uses minimap2 for quality filtering
5. **Consensus Generation**: Uses MAFFT alignment and majority voting (50% threshold)
6. **Mutation Calling**: Compares consensus to reference sequence

### Key Parameters
- `MIN_ALIGNMENT_IDENTITY = 0.8`: Minimum alignment identity for read inclusion
- `threshold = 0.5`: Minimum fraction for consensus base calling

## Script 2: amino_acid_analysis.py

### Purpose
Translates consensus DNA sequences to proteins and identifies amino acid mutations.

### Usage
```bash
python amino_acid_analysis.py <reference_directory>
```

### Parameters
- `<reference_directory>`: Same directory used for tokuseq.py

### Example
```bash
python amino_acid_analysis.py /path/to/project/wt-p1
```

### Input Files

#### Required from tokuseq.py output:
- `demultiplex_summary.csv`

#### Additional input:
- `reference_orf.fasta`: Open reading frame sequence for translation

### Output

#### amino_acid_mutations_YYMMDD.csv
Contains detailed amino acid analysis:
- **Well**: Well identifier
- **ReadCount**: Number of reads
- **ORF_DNA**: Extracted ORF region from consensus
- **ConsensusProtein**: Translated protein sequence
- **AAMutations**: Amino acid mutations (e.g., "A123V;S193del" or "WT")
- **NumAAMutations**: Number of amino acid mutations

### Algorithm Details

1. **ORF Alignment**: Aligns consensus DNA to reference ORF using MAFFT
2. **ORF Extraction**: Extracts ORF region while preserving alignment gaps
3. **Translation**: Converts DNA to protein, handling gaps and ambiguities
4. **Protein Alignment**: Aligns translated proteins using MAFFT
5. **Mutation Detection**: Identifies substitutions, deletions, and insertions

### Mutation Notation
- **Substitutions**: `A123V` (Alanine at position 123 to Valine)
- **Deletions**: `S193del` (Serine at position 193 deleted)
- **Insertions**: `ins123V` (Valine inserted after position 123)
- **Wild-type**: `WT` (no mutations detected)

### Error Handling
- `No_consensus`: No consensus sequence available
- `Alignment_error`: Alignment failed
- `Translation_failed`: Translation unsuccessful
- `X`: Ambiguous amino acid due to gaps or ambiguous codons

## Complete Workflow Example

### 1. Prepare Input Files
```bash
# Ensure your directory has all required files
ls /path/to/project/wt-p1/
# Should show: barcodes.fasta, reference.fasta, reference_orf.fasta, input_reads.fastq
```

### 2. Run Demultiplexing
```bash
python tokuseq.py /path/to/project/wt-p1 2000 3000
```

### 3. Run Amino Acid Analysis
```bash
python amino_acid_analysis.py /path/to/project/wt-p1
```

### 4. Results
Check the output files:
- `demultiplex_summary.csv` - DNA-level analysis
- `amino_acid_mutations_YYMMDD.csv` - Protein-level analysis
- `extracted_fastq_YYMMDD/` - Demultiplexed reads and alignments

## Quality Control and Troubleshooting

### Common Issues

1. **No reads assigned to wells**
   - Check barcode sequences in `barcodes.fasta`
   - Verify read length cutoffs are appropriate
   - Ensure barcodes are present in sequencing reads

2. **Poor consensus quality**
   - Increase minimum read count threshold
   - Check alignment quality in `*_alignment.fasta` files
   - Verify reference sequence accuracy

3. **Translation errors**
   - Ensure `reference_orf.fasta` contains complete ORF
   - Check for frameshifts in consensus sequences
   - Verify ORF alignment quality

### Quality Metrics
- **ReadCount**: Higher counts improve consensus reliability
- **Alignment Identity**: Should be ≥0.8 for included reads
- **Consensus Quality**: Check alignment files for consistency

## Output Interpretation

### DNA-level Results (demultiplex_summary.csv)
- Wells with high read counts (>10) provide more reliable consensus
- "WT" indicates no mutations detected
- Multiple mutations separated by semicolons

### Protein-level Results (amino_acid_mutations_YYMMDD.csv)
- `NumAAMutations = 0`: Wild-type protein
- `NumAAMutations > 0`: Mutant protein
- Check `ConsensusProtein` for translation quality

### Summary Statistics
The amino acid analysis script provides summary statistics:
- Total wells analyzed
- Wells with sequencing reads
- Wells with translatable proteins
- Wells with amino acid mutations

## Best Practices

1. **Read Length Filtering**: Use appropriate length cutoffs for your amplicon
2. **Quality Control**: Review alignment files for consensus quality
3. **Reference Preparation**: Ensure accurate reference and ORF sequences
4. **Barcode Design**: Use distinct barcodes to minimize cross-contamination
5. **Data Validation**: Cross-check DNA and protein-level results

## Version Information

This documentation is current as of July 2, 2025. The pipeline uses:
- BioPython for sequence manipulation
- MAFFT for multiple sequence alignment
- minimap2 for reference alignment
- Consensus calling with majority voting (≥50%)
