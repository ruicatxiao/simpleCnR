![Image](https://github.com/user-attachments/assets/68139213-5a1d-4948-8d40-ec32ae434382)


# Simple CnR QC

A Simple CUT&RUN Processing Pipeline specifically tailored for Cryptosporidium parvum

## Overview

This is a pipeline for processing CUT&RUN sequencing data for C.parvum. It automates the workflow from raw FASTQ files to reports and figures, including read trimming, alignment, deduplication, filter for primary alignment, coverage analysis, BAM signal calculation, BigWig intergration and data visualization.

## Pipeline Author

Rui Xiao, PhD  

## Pipeline Contributor

 - Katelyn Walzer, PhD
 - Christopher Noetzel, PhD

## Prerequisites

The pipeline requires the following tools to be installed and available in your PATH:

| Software      | Tested Version | Purpose                        |
|---------------|----------------|--------------------------------|
| trim_galore   | 0.6.7          | Read quality trimming          |
| bowtie2       | 2.4.4          | Sequence alignment             |
| samtools      | 1.13           | SAM/BAM manipulation           |
| GATK          | 4.6            | Duplicate marking and metrics  |
| deeptools     | 3.5.5          | Coverage analysis and plotting |
| multiqc       | 1.27.1         | Quality control reporting      |

Tested on Ubuntu-server 22.04 LTS.

## Installation

```bash
# Clone the repository
git clone https://github.com/username/CNR2Results.git
cd CNR2Results

# No additional installation required - the script is ready to use
```

## Usage

```bash
python cnr2results.py \
    --fqinput <PATH_TO_RAW_FASTQ> \
    --genomeindex cpbgfbt2idx \
    --output <PATH_TO_OUTPUT> \
    --annotation CpBGF_genome_V1.bed \
    --threads <NUM_OF_CPU>
```

### Parameters

- `--fqinput`: Directory containing raw FASTQ files
- `--genomeindex`: Bowtie2 index base path (provided in this repository)
- `--output`: Directory for output files
- `--annotation`: BED file with genomic annotations (provided in this repository)
- `--threads`: Number of CPU cores to use (default: 8)

## Input Structure

The pipeline expects Illumina paired-end sequencing data with standard naming convention:

```
<PATH_TO_RAW_FASTQ>/
├── Sample1_S1_L001_R1_001.fastq.gz
├── Sample1_S1_L001_R2_001.fastq.gz
├── Sample2_S2_L001_R1_001.fastq.gz
├── Sample2_S2_L001_R2_001.fastq.gz
└── ...
```

## Output Structure

The pipeline generates the following directory structure:

```
<PATH_TO_OUTPUT>/
├── trimmed/                   # Quality-trimmed FASTQ files
├── bam/                       # Initial BAM alignments
├── dedup_bams/                # Deduplicated BAM files, USED FOR MACS3!!!
├── primary_bams/              # BAMs with only primary alignments
│   └── *.bam and *.bam.bai    # Indexed BAM files
├── coverage_reports/          # Coverage statistics
├── histograms/                # Insert size distributions
├── bigwigs/                   # Normalized coverage tracks
├── deeptools/                 # Matrix files for heatmaps/profiles
├── plots/                     # Publication-ready figures
│   ├── heatmap.png            # Gene-centered heatmap
│   └── profile.png            # Aggregated signal profile
├── multiqc/                   # Comprehensive QC report
│   └── cnr_report.html        # MultiQC report
└── logs/                      # Detailed logs for each step
    ├── Trimming_Reads.log
    ├── Alignment.log
    ├── BAM_Processing.log
    └── ...
```

## Pipeline Steps

1. **Trimming Reads**: Adapter and quality trimming with Trim Galore
2. **Alignment**: Mapping to reference genome with Bowtie2
3. **BAM Processing**:
   - Marking duplicates
   - Filtering for primary alignments
   - Indexing BAM files
   - Generating coverage reports
   - Collecting insert size metrics
   - Creating normalized bigWig files
4. **Matrix Generation**: Computing gene-centered matrices
5. **Plotting**: Creating heatmaps and profile plots
6. **Final Report**: Generating MultiQC report

## Example Visualization Outputs

Upon successful completion, the pipeline generates:

- Heatmaps showing protein binding distribution around genes
- Aggregated signal profiles
- Insert size distributions
- Comprehensive MultiQC report with quality metrics

## Troubleshooting

Check the log files in the `logs/` directory for detailed information about each step. Common issues include:

- Missing paired-end files
- Insufficient disk space
- Memory limitations for large datasets

## License

[MIT License](LICENSE)

## Citation

If you use this pipeline in your research, please cite:
```
<PLACE HOLDER FOR MANUSCRIPT>
```