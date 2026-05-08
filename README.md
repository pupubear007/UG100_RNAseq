# UG100_RNAseq

A bioinformatics pipeline for differential expression analysis of **Sclerotinia sclerotiorum** (white mold) infecting **soybean** (*Glycine max*) and **sunflower** (*Helianthus annuus*).

## Experimental Design

Dual RNA-seq experiment profiling host and pathogen transcriptomes across two host species infected with six fungal isolates of varying aggressiveness, sampled at four time points.

| Factor | Levels |
|---|---|
| Host species | Soybean (Gm), Sunflower (Ha) |
| Fungal isolates | WISS47, JS659 (low generalist), Xtra7, MNSS6 (high generalist), SsPotter (soybean-specific), BN172 (sunflower-specific), NC (negative control) |
| Time points | 0h/NC, 24h, 48h, 96h post-infection |
| Replicates | 3–6 per condition |
| Total samples | 120 |

## Pipeline

All steps are implemented as standalone SLURM array jobs designed for HPC (Minnesota Supercomputing Institute).

### 1. Raw Data Preparation (`Raw_cram_rename_n_relocate/`)
Copy 120 CRAM files from shared storage and rename to a clean format (`{time}_{rep}_{host}_{isolate}.cram`). Generates `sample_log.csv` metadata.

### 2. CRAM to FASTQ + QC (`Clean_fastq_preparation/`)
Convert CRAM to paired-end FASTQ via `samtools fastq`, then trim with Trimmomatic (ILLUMINACLIP, LEADING:3, TRAILING:3, SLIDINGWINDOW:4:15, MINLEN:36) and Cutadapt (Q20, custom adapters). FASTQ validation and repair with BBMap `reformat.sh`/`repair.sh`. Quality reports via FastQC and MultiQC.

### 3. FASTQ Repair (`fix_fastq/`)
Optional validation and repair of corrupted FASTQ files that fail STAR alignment.

### 4. Reference Preparation (`ref_and_indice_preparation/`)
Build 15 indices (5 per aligner: host, pathogen, host+pathogen combined) for **STAR** (custom build with increased `seqReadLengthMax`), **HISAT2**, and **Salmon**. References: Phytozome V13 (soybean Wm82.a6.v1), Phytozome V12 (sunflower r1.2), JGI MycoCosm (S. sclerotiorum). Combined references add organism prefixes to FASTA headers.

### 5. Alignment (`alignment/`)
- **STAR**: Spliced alignment to combined host+pathogen indices with `--quantMode GeneCounts`
- **Salmon**: Mapping-based transcript quantification with `--validateMappings` and tx2gene mapping

### 6. Count Collection (`collect_DEG_counts/`)
- **STAR counts**: `featureCounts` on BAM files, separately quantifying host and pathogen genes. Merged into unified count matrices.
- **Salmon counts**: Collect gene-level counts from `quant.genes.sf`, separate by organism prefix.

### 7. Statistical Analysis (`count_analysis/`)
R Markdown notebooks using **DESeq2** for differential expression, **clusterProfiler** for GO enrichment, **vegan** for PERMANOVA, and various visualization packages (ComplexHeatmap, EnhancedVolcano, UpSetR).

Analyses include: per-host isolate comparisons, cross-species comparisons, ANOVA/clustering, pairwise isolate contrasts, and GO enrichment (host and pathogen separately). Parallel analysis sets for STAR and Salmon counts.

## Requirements

- **Conda environment**: `pu_RNAseq` (STAR, Salmon, HISAT2, Trimmomatic, Cutadapt, FastQC, MultiQC, subread, BBMap, samtools, pigz, gffread, GNU parallel)
- **R packages**: DESeq2, clusterProfiler, vegan, ComplexHeatmap, EnhancedVolcano, UpSetR, tidyverse, etc.
- **HPC**: SLURM workload manager

## Usage

Each step is submitted as a SLURM array job. Modify input/output paths and account parameters before running.

```bash
# Step 1: Copy and rename CRAM files
cd Raw_cram_rename_n_relocate
sbatch --array=0-119 copy_cram_array.sh
bash cram_rename_csv.sh

# Step 2: CRAM to FASTQ
cd ../Clean_fastq_preparation
bash cram_to_rawfastq_R1R2.sh

# Step 3: QC and trimming
sbatch raw_to_1st_trimmed_n_qc.sh
sbatch trimmed_to_clean_fastq.sh

# Step 4: Build reference indices
cd ../ref_and_indice_preparation
sbatch build_all_jgi_indices_FIXED.sh

# Step 5: Alignment
cd ../alignment
sbatch alignment_STAR_20260326.sh
sbatch alignment_Salmon_20260326.sh

# Step 6: Count collection
cd ../collect_DEG_counts/STAR_counts
sbatch --array=0-7 dual_count_v5.sh
bash create_count_matrix_FIXED.sh
cd ../Salmon_count
bash collect_salmon_counts_ONLY.sh

# Step 7: Statistical analysis (in RStudio)
# Open and knit the appropriate .Rmd files in count_analysis/
```

## Output

- Unified count matrices (host + pathogen genes)
- Separated host-only and pathogen-only count files
- Pathogen load statistics per sample
- DEG lists (CSV) for each comparison
- GO enrichment results
- MultiQC quality reports (HTML)

## License

Apache License 2.0
