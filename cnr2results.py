#!/usr/bin/env python
"""
CNR2Results - Integrated CUT&RUN Processing Pipeline
Rui Xiao, PhD
Bioinformatician C
School of Veterinary Medicine
University of Pennsylvania
ruicatx@upenn.edu

Required tools in PATH:
- trim_galore
- bowtie2
- samtools
- gatk
- bamCoverage (deeptools)
- computeMatrix/computeMatrixOperations (deeptools)
- plotHeatmap/plotProfile (deeptools)
- multiqc
"""

import argparse
import subprocess
import sys
from pathlib import Path
from tqdm import tqdm

def validate_tools():
    required_tools = [
        'trim_galore', 'bowtie2', 'samtools',
        'gatk', 'bamCoverage', 'computeMatrix',
        'computeMatrixOperations', 'plotHeatmap',
        'plotProfile', 'multiqc'
    ]
    missing = []
    for tool in required_tools:
        if not subprocess.run(['which', tool], capture_output=True).returncode == 0:
            missing.append(tool)
    if missing:
        sys.exit(f"Error: Missing required tools: {', '.join(missing)}")

def run_stage(command, stage_name, progress, dirs):
    progress.set_description(f"{stage_name: <25}")
    
    # Create unique log filename
    log_file = dirs['logs'] / f"{stage_name.replace(' ', '_')}.log"
    
    with open(log_file, 'w') as f:
        # Write command header
        f.write(f"COMMAND: {command}\n{'='*50}\n")
        
        # Run with real-time logging
        result = subprocess.run(
            ['bash', '-c', command],
            stdout=f,
            stderr=subprocess.STDOUT,  # Merge stderr into stdout
            text=True
        )
    
    if result.returncode != 0:
        sys.exit(
            f"ERROR in {stage_name}\n"
            f"Full log: {log_file}\n"
            f"Last 5 lines:\n{open(log_file).readlines()[-5:]}"
        )

def main():
    parser = argparse.ArgumentParser(
        description='CUT&RUN Processing Pipeline (Revised)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--fqinput', required=True, help='Path to raw FASTQ directory')
    parser.add_argument('--genomeindex', required=True, help='Bowtie2 index base path')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--threads', type=int, default=8, help='Number of threads')
    parser.add_argument('--annotation', required=True, help='Path to annotation BED file')
    args = parser.parse_args()

    validate_tools()
    
    # Setup directory structure
    base = Path(args.output)
    dirs = {
        'trimmed': base/'trimmed',
        'bam': base/'bam',
        'dedup': base/'dedup_bams',
        'primary': base/'primary_bams',
        'coverage': base/'coverage_reports',
        'hist': base/'histograms',
        'bw': base/'bigwigs',
        'matrix': base/'deeptools',
        'plots': base/'plots',
        'multiqc': base/'multiqc',
        'logs': base/'logs'
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    # Processing pipeline (revised commands)
    stages = [
        ('Trimming Reads', f"""
        find "{args.fqinput}" -name '*_R1_001.fastq.gz' | while read read1; do
            read2="${{read1/_R1_001.fastq.gz/_R2_001.fastq.gz}}"
            if [[ -f "$read2" ]]; then
                trim_galore --cores {args.threads} --paired --gzip \\
                    -o "{dirs['trimmed']}" "$read1" "$read2"
            else
                echo "Missing paired-end file for $read1"
                exit 1
            fi
        done
        """),
        
        ('Alignment', f"""
        find "{dirs['trimmed']}" -name '*_R1_001_val_1.fq.gz' | while read read1; do
            read2="${{read1/_R1_001_val_1.fq.gz/_R2_001_val_2.fq.gz}}"
            samplename=$(basename "$read1" _R1_001_val_1.fq.gz)
            bamname="${{samplename%_S*}}"
            bowtie2 --very-sensitive -k 10 -x "{args.genomeindex}" \\
                -1 "$read1" -2 "$read2" \\
                --rg-id "$samplename" --rg "SM:$samplename" --rg "PL:ILLUMINA" \\
                -p {args.threads} | samtools sort -@ {args.threads} \\
                -O BAM -o "{dirs['bam']}/$bamname.bam"
        done
        """),
        
        ('BAM Processing', f"""
        for bam in "{dirs['bam']}"/*.bam; do
            [[ ! -f "$bam" ]] && continue
            sample=$(basename "$bam" .bam)
            
            # Fix 1: Add proper curly braces for variable expansion
            # Mark duplicates
            gatk MarkDuplicates \\
                -I "$bam" \\
                -O "{dirs['dedup']}/${{sample}}_dedup.bam" \\
                -M "{dirs['dedup']}/${{sample}}_dup_report.txt" \\
                --VALIDATION_STRINGENCY SILENT
            
            # Primary alignments
            samtools view -@ {args.threads} -F 256 -b \\
                "{dirs['dedup']}/${{sample}}_dedup.bam" | \\
                samtools sort -@ {args.threads} -O BAM \\
                -o "{dirs['primary']}/${{sample}}_dedup_primary.bam"
            
            samtools index -@ {args.threads} \\
                "{dirs['primary']}/${{sample}}_dedup_primary.bam"
            
            # Coverage report
            samtools coverage -d 0 \\
                -o "{dirs['coverage']}/${{sample}}_cov_report.txt" \\
                "{dirs['primary']}/${{sample}}_dedup_primary.bam"
            
            # Fragment metrics
            gatk CollectInsertSizeMetrics \\
                -I "{dirs['primary']}/${{sample}}_dedup_primary.bam" \\
                -O "{dirs['hist']}/${{sample}}_histogram.txt" \\
                -H "{dirs['hist']}/${{sample}}_histogram.pdf"
            
            # BigWig generation
            bamCoverage -b "{dirs['primary']}/${{sample}}_dedup_primary.bam" \\
                -o "{dirs['bw']}/${{sample}}.bw" -p {args.threads} -bs 1 \\
                --normalizeUsing RPGC --effectiveGenomeSize 7068952
        done
        """),
        
        ('Matrix Generation', f"""
        gz_files=()
        for bw in "{dirs['bw']}"/*.bw; do
            [[ ! -f "$bw" ]] && continue
            sample=$(basename "$bw" .bw)
            computeMatrix scale-regions --metagene \\
                -S "$bw" -R "{args.annotation}" \\
                -bs 10 -m 1200 -a 300 -b 300 \\
                -o "{dirs['matrix']}/$sample.gz" \\
                --outFileNameMatrix "{dirs['matrix']}/${{sample}}_matrix.tab" \\
                --outFileSortedRegions "{dirs['matrix']}/${{sample}}_regions.tab" \\
                --missingDataAsZero -p {args.threads}
            gz_files+=("{dirs['matrix']}/$sample.gz")
        done
        [[ ${{#gz_files[@]}} -gt 0 ]] && computeMatrixOperations cbind \\
            -m "${{gz_files[@]}}" \\
            -o "{dirs['matrix']}/combined_matrix.gz"
        """),
        
        ('Plotting', f"""
        [[ -f "{dirs['matrix']}/combined_matrix.gz" ]] && {{
            plotHeatmap -m "{dirs['matrix']}/combined_matrix.gz" \\
                -out "{dirs['plots']}/heatmap.png" \\
                --colorMap viridis --zMax 5 \\
                --heatmapHeight 30 --heatmapWidth 20 \\
                --perGroup --whatToShow 'heatmap and colorbar' \\
                --outFileSortedRegions "{dirs['plots']}/heatmap_regions.tab"
            
            plotProfile -m "{dirs['matrix']}/combined_matrix.gz" \\
                -out "{dirs['plots']}/profile.png" \\
                --plotType se --perGroup \\
                --plotWidth 20 --plotHeight 15
        }}
        """),
        
        ('Final Report', f"""
        multiqc --outdir "{dirs['multiqc']}" \\
            --filename cnr_report \\
            --clean-up "{base}"
        """)
    ]

    with tqdm(total=len(stages), desc="Overall Progress") as main_progress:
        for stage_name, command in stages:
            run_stage(command, stage_name, main_progress, dirs)
            main_progress.update(1)

if __name__ == "__main__":
    main()