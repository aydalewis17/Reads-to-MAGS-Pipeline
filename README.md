# Reads-to-MAGS-Pipeline_within_sample

Metagenomic pipeline that takes raw paired-end reads from multiple samples and produces dereplicated, taxonomy-annotated MAGs with **within-sample** relative abundance estimates.

Written by AL and Claude Sonnet 4.6, incorporating steps from the Borton Lab metagenomics workshop.

---

## Files

| File | Description |
|------|-------------|
| `reads_to_MAGs_pipeline_within_sample.py` | Main pipeline script |
| `RUN_reads_to_MAGs.sh` | SLURM submission script |

---

## Before You Run

**1. Copy and prepare raw reads**

```bash
# Copy raw reads to the project directory
cp -r /ORG-Data-phoenix/Manure_Lagoons/seq_CSU_Dec2025 /home/projects-phoenix/Manure_Lagoons/MetaG

# Rename to raw_reads
mv /home/projects-phoenix/Manure_Lagoons/MetaG/seq_CSU_Dec2025 /home/projects-phoenix/Manure_Lagoons/MetaG/raw_reads

# Unzip all FASTQ files
gunzip /home/projects-phoenix/Manure_Lagoons/MetaG/raw_reads/*.fastq.gz
```

Raw reads must be named: `{sample}_R1_001.fastq` and `{sample}_R2_001.fastq`

**2. Update USER CONFIG in `reads_to_MAGs_pipeline_within_sample.py`**

At minimum, check `PROJECT_DIR` and `RAW_DIR` point to the correct locations for your project. See the full list of configurable parameters below.

**3. Update the sample list in `RUN_reads_to_MAGs.sh`**

Edit the `-s` argument to include your sample names, comma-separated with no spaces:

```bash
python3 /path/to/reads_to_MAGs_pipeline_UPDATED.py \
-s "SAMPLE1,SAMPLE2,SAMPLE3"
```

---

## Running the Pipeline

Submit to SLURM:

```bash
sbatch RUN_reads_to_MAGs.sh
```

The job requests 1 node, 40 tasks, 500 GB memory, and a 4-day walltime on the `borton-hi` or `borton-low` partition. Email notifications are sent on job start, end, and failure.

To run interactively (not recommended for full datasets):

```bash
python3 reads_to_MAGs_pipeline_UPDATED.py -s "SAMPLE1,SAMPLE2"
```

---

## Configurable Parameters

All parameters are set at the top of `reads_to_MAGs_pipeline_UPDATED.py` in the `USER CONFIG` block. No changes should be needed anywhere else in the script.

### Paths & Resources

| Parameter | Default | Description |
|-----------|---------|-------------|
| `PROJECT_DIR` | `/home/projects-phoenix/Manure_Lagoons/MetaG` | Root project directory. Sample subdirectories are created here. |
| `RAW_DIR` | `{PROJECT_DIR}/raw_reads` | Directory containing raw FASTQ files. |
| `threads` | `50` | Global thread count used by most tools. |
| `memory` | `400G` | Memory cap passed to BBMap and reformat.sh. |

### Assembly

| Parameter | Default | Description |
|-----------|---------|-------------|
| `KMER_MIN` | `31` | Smallest k-mer size used by MEGAHIT. Smaller values help assemble low-coverage regions. |
| `KMER_MAX` | `121` | Largest k-mer size used by MEGAHIT. Larger values improve assembly of high-coverage, repetitive regions. |
| `KMER_STEP` | `10` | Step size between k-mer iterations. Smaller steps = more thorough but slower assembly. |
| `MIN_CONTIG_LEN` | `2500` | Minimum contig length (bp) kept after assembly. Shorter contigs are discarded before binning. |
| `POLYG_TRIM` | `50` | Number of terminal G bases bbduk trims from the right end of reads. Relevant for NovaSeq/NextSeq poly-G artifact. |

### Mapping Identity Filters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `BBMAP_MIN_ID` | `0.99` | Minimum read identity (fraction) for BBMap when mapping reads back to contigs for binning. |
| `BOWTIE2_MIN_ID` | `0.99` | Minimum read identity (fraction) for reads mapped to the final MAG database (Bowtie2 post-filter). |
| `COVERM_MIN_ID` | `0.97` | Minimum per-read-pair identity for CoverM to count a read as mapping. Used in all three CoverM calls. |
| `COVERM_MIN_BREADTH` | `0.75` | Fraction of a genome that must be covered before CoverM reports abundance (breadth-of-coverage cutoff). |

### MAG Quality Thresholds

| Parameter | Default | Description |
|-----------|---------|-------------|
| `CHECKM2_MIN_COMP` | `50` | Minimum completeness (%) for a bin to be kept as a MAG. >50% = medium quality, >90% = high quality (MIMAG). |
| `CHECKM2_MAX_CONT` | `10` | Maximum contamination (%) allowed in a kept MAG. Higher values indicate sequences from multiple organisms. |
| `DREP_MIN_COMP` | `50` | dRep completeness filter applied before clustering. Should match `CHECKM2_MIN_COMP` unless re-filtering is desired. |
| `DREP_MAX_CONT` | `10` | dRep contamination filter applied before clustering. Should match `CHECKM2_MAX_CONT` unless re-filtering is desired. |

### Per-Tool Thread Overrides

These are set independently from the global `threads` variable and will **not** auto-scale with it.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `CHECKM2_THREADS` | `10` | Threads for CheckM2 quality assessment. |
| `GTDBTK_CPUS` | `20` | CPUs for GTDB-Tk taxonomy classification. |
| `BOWTIE2_THREADS` | `15` | Threads for Bowtie2 index building and read mapping. |

---

## Pipeline Steps

### Step 1 — Quality trimming (sickle)
Trims low-quality bases from raw paired-end reads using sickle in paired-end mode with Sanger quality scores. Discarded/unpaired reads are written to a separate file.

### Step 1.2 — Adapter & poly-G removal (bbduk)
Removes adapter sequences and poly-G tails from the right end of reads. Important for NovaSeq/NextSeq data where poly-G is a common sequencing artifact.

### Step 2 — Assembly (MEGAHIT)
Assembles cleaned paired-end reads into contigs using a k-mer range of 31–121 (step 10), iteratively building from short to long k-mers to maximize assembly quality across coverage depths.

### Step 3 — Contig length filtering (pullseq)
Filters out contigs shorter than `MIN_CONTIG_LEN` (default 2500 bp). Short contigs are too small to bin reliably and add noise.

### Step 4 — Read mapping back to contigs (BBMap)
Maps trimmed reads back to the filtered contigs at ≥90% identity, producing a SAM file with per-contig coverage depth needed for binning.

### Step 5 — SAM to BAM conversion (samtools view)
Converts SAM to binary BAM format to reduce file size and enable sorting.

### Step 6 — Sort BAM (samtools sort)
Sorts the BAM file by genomic coordinate, required for coverage-based tools.

### Step 6.1 — High-quality mapping filter (reformat.sh)
Filters the sorted BAM to keep only read pairs mapping at ≥`BBMAP_MIN_ID` identity, requiring both reads in a pair to map (pairedonly) and only primary alignments (primaryonly). Produces a cleaner depth signal for binning.

### Step 7 — Binning (MetaBAT)
Clusters contigs into MAGs based on coverage depth and tetranucleotide frequency using the high-quality filtered BAM.

### Step 8 — Assembly statistics (contig_stats.pl)
Calculates summary statistics (N50, total length, etc.) for the raw MEGAHIT assembly.

### Step 9 — MAG quality assessment (CheckM2)
Evaluates completeness and contamination of each bin using CheckM2 v1.1.0. Results are written to `quality_report.tsv`.

### Step 10 — Filter and copy MQ/HQ MAGs
Reads CheckM2 results and copies only bins passing quality thresholds (completeness >`CHECKM2_MIN_COMP`%, contamination <`CHECKM2_MAX_CONT`%) into the `MAGs/` folder, prepending the sample name to each filename.

### Step 11 — Taxonomy assignment (GTDB-Tk)
Assigns taxonomy to all passing MAGs using GTDB-Tk v2.7.1 against the r232 database. Uses the pre-sketched skani database included with the package — no manual sketch building required.

### Step 12 — Cleanup of contig-mapping files
Deletes intermediate SAM, sorted BAM, and filtered BAM files from Steps 4–6 to free disk space. The unsorted BAM is gzip-compressed and kept as an archive.

### Step 13 — Within-sample dereplication (dRep)
Dereplicates MAGs within the sample using dRep, selecting representative genomes at the default ANI threshold. Only MAGs passing the quality thresholds are considered.

### Step 14 — Rename contigs and build MAG database FASTA
Renames contig headers in dereplicated MAGs using `rename_bins_like_dram.py` (required for DRAM compatibility), then concatenates all renamed MAGs into a single FASTA file (`{sample}_derep_MAGs.fa`) used as the mapping database.

### Step 15 — Build Bowtie2 index
Builds a Bowtie2 index from the concatenated MAG database FASTA for fast read alignment.

### Step 16 — Map reads to MAG database (Bowtie2)
Maps trimmed reads to the MAG database, converts to BAM, sorts, and filters for ≥`BOWTIE2_MIN_ID` identity paired primary alignments. Counts mapped read pairs and calculates the percent of total trimmed reads mapped, appending results to `sample_MAG_database_mapping_summary.tsv` in `PROJECT_DIR`.

### Step 17 — Relative abundance (CoverM)
Runs CoverM three times on the same sorted BAM:
- **reads_per_base** — raw read depth normalized by genome length, no breadth filter
- **coverm_min75** — same but only reports genomes where ≥`COVERM_MIN_BREADTH` of the genome is covered (filters spurious hits)
- **trimmed_mean** — trimmed mean coverage depth, robust to coverage outliers

---

## Output Structure

For each sample, the following directories are created under `{PROJECT_DIR}/{sample}/`:

```
{sample}/
├── trimmed_reads/               # Sickle and bbduk trimmed reads
├── megahit_out/                 # MEGAHIT assembly and all intermediate mapping files
│   ├── checkm2_v1.1.0_fa/       # CheckM2 quality results
│   └── {sample}_2500.fa         # Filtered contigs (>2500 bp)
├── MAGs/                        # MQ/HQ MAGs passing CheckM2 thresholds
├── gtdb_v2.7.0_r232/            # GTDB-Tk taxonomy output
├── dRep_v3.0.0_MAGs/            # dRep dereplication output
│   └── dereplicated_genomes/
│       └── genome_renamed/      # Contigs renamed for DRAM compatibility
├── bowtie_DB/                   # MAG database FASTA and Bowtie2 index
├── MAG_db_mapping/              # BAM files from mapping reads to MAG database
└── coverm_output/               # CoverM relative abundance tables
```

At the project level:
```
{PROJECT_DIR}/
└── sample_MAG_database_mapping_summary.tsv   # Mapping rate summary across all samples
```

---

## Tools Required

| Tool | Version | Conda env |
|------|---------|-----------|
| sickle | any | base |
| bbduk / bbmap / reformat | BBTools | base |
| MEGAHIT | any | base |
| pullseq | any | base |
| samtools | any | base |
| MetaBAT | any | base |
| contig_stats.pl | — | base |
| CheckM2 | v1.1.0 | `checkm2_v1.1.0` |
| GTDB-Tk | v2.7.1 | `gtdbtk_v2.7.0` |
| dRep | v3.0.0 | base |
| rename_bins_like_dram.py | — | `scripts` (Miniconda2) |
| Bowtie2 | any | base |
| CoverM | any | base |

---

## Notes

- Errors for individual samples are printed to stderr but do not stop the pipeline — remaining samples continue processing.
- The dRep step dereplicates **within** each sample, not across samples.
- `CHECKM2_MIN_COMP`, `DREP_MIN_COMP`, and `DREP_MAX_CONT` should be kept in sync with `CHECKM2_MAX_CONT` unless you intentionally want dRep to apply different thresholds than CheckM2.
- The three per-tool thread overrides (`CHECKM2_THREADS`, `GTDBTK_CPUS`, `BOWTIE2_THREADS`) do **not** inherit from the global `threads` variable — update them separately if scaling up or down.