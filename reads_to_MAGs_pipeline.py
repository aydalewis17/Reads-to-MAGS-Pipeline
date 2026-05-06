#!/usr/bin/env python3

# ---------------------------------------------------------------------------------------------------------------------------
## AL Notes: this script was edited by AL and Claude Sonnet 4.6 to include steps from the workshop, so all QC steps are included.
### Steps 1-14 run per-sample (assembly, binning, QC, taxonomy, within-sample dRep).
### Steps 15-17 run ONCE after all samples complete, using a CROSS-SAMPLE dereplicated MAG database
### for BETWEEN-SAMPLE relative abundance (Bowtie2 + CoverM).
### Things this script has that the original does not:
# 1. Assembly statistics (contig_stats.pl)
# 2. MAG quality control with CheckM2
# 3. Filtering for HQ/MQ MAGs
# 4. Copying HQ/MQ MAGs into a MAGs/ folder
# 5. GTDB-Tk taxonomy assignment
# 6. Within-sample dRep of MAGs (Step 13) - removes redundancy within each sample
# 7. Cross-sample dRep (Step 15) - builds a single shared MAG database across all samples
# 8. rename_bins_like_dram.py to fix contig headers
# 9. Building a single shared Bowtie2 index (Step 16)
# 10. Mapping ALL samples' reads to the shared database (Step 16)
# 11. Calculating percent reads mapped per sample
# 12. Between-sample CoverM abundance table (Step 17) - all samples as columns, MAGs as rows
# 13. Creating additional folders:
    # MAGs
    # dRep_v3.0.0_MAGs (per-sample within-sample dRep)
    # cross_sample_dRep/  (project-level between-sample dRep database)
    # MAG_db_mapping (per-sample BAMs mapped to cross-sample DB)
    # gtdb_v2.7.0_r232
    # coverm_output (project-level, between-sample abundance tables)
    # Cleanup of intermediate mapping files
# ----------------------------------------------------------------------------------------------------------------------------

import subprocess
import os
import sys
import shutil
import gzip
import csv
import argparse

parser = argparse.ArgumentParser(description="Metagenomic Pipeline")
parser.add_argument("-s", "--sample", required=True, help="Samples separated with commas if more than 1")
args = parser.parse_args()

# =====================
# USER CONFIG
# =====================
PROJECT_DIR = "/home/projects-phoenix/Manure_Lagoons/MetaG" # added to make sample directories OUTSIDE of the raw_reads directory
RAW_DIR = os.path.join(PROJECT_DIR, "raw_reads")
threads = "50"
memory = "400G"

# NOTE: contig_stats.pl and rename_bins_like_dram.py are Borton Lab shared utility scripts stored in ORG-Data, not conda-installed tools. If this path breaks, check that ORG-Data is mounted and the script hasn't moved.
CONTIG_STATS = "/ORG-Data/scripts/quicklooks/contig_stats.pl"
RENAME_BINS_SCRIPT = "rename_bins_like_dram.py"

# --- Assembly ---
# Smallest k-mer size used by MEGAHIT. Smaller values help assemble low-coverage regions.
KMER_MIN = "31"
# Largest k-mer size used by MEGAHIT. Larger values improve assembly of high-coverage, repetitive regions.
KMER_MAX = "121"
# Step size between k-mer iterations. Smaller steps = more thorough but slower assembly.
KMER_STEP = "10"
# Minimum contig length (bp) kept after assembly. Shorter contigs are discarded before binning.
MIN_CONTIG_LEN = "2500"
# Number of terminal G bases bbduk trims from the right end of reads. Relevant for NovaSeq/NextSeq poly-G artifact.
POLYG_TRIM = "50"

# --- Mapping identity filters ---
# Minimum read identity (fraction) for BBMap when mapping reads back to contigs for binning.
BBMAP_MIN_ID = "0.99"
# Minimum read identity (fraction) for reads mapped to the final MAG database (Bowtie2 post-filter).
BOWTIE2_MIN_ID = "0.99"
# Minimum per-read-pair identity for CoverM to count a read as mapping. Used in all three CoverM calls.
COVERM_MIN_ID = "0.97"
# Fraction of a genome that must be covered before CoverM reports abundance (breadth-of-coverage cutoff).
COVERM_MIN_BREADTH = "0.75"

# --- MAG quality thresholds ---
# Minimum completeness (%) for a bin to be kept as a MAG. >50% = medium quality, >90% = high quality (MIMAG).
CHECKM2_MIN_COMP = 50
# Maximum contamination (%) allowed in a kept MAG. Higher values indicate sequences from multiple organisms.
CHECKM2_MAX_CONT = 10
# dRep completeness filter applied before clustering. Should match CHECKM2_MIN_COMP unless re-filtering is desired.
DREP_MIN_COMP = "50"
# dRep contamination filter applied before clustering. Should match CHECKM2_MAX_CONT unless re-filtering is desired.
DREP_MAX_CONT = "10"

# --- Per-tool thread overrides ---
# NOTE: these are hardcoded separately from the global `threads` variable and will NOT auto-scale with it.
# Threads for CheckM2 quality assessment.
CHECKM2_THREADS = "10"
# CPUs for GTDB-Tk taxonomy classification.
GTDBTK_CPUS = "20"
# Threads for Bowtie2 index building and read mapping.
BOWTIE2_THREADS = "15"

samples = args.sample.split(',')

# =====================
# HELPERS
# =====================
def join_raw(filename: str) -> str:
    if RAW_DIR:
        return os.path.abspath(os.path.join(RAW_DIR, filename))
    return os.path.abspath(filename)

def count_fastq_reads(fastq_path: str) -> int:
    line_count = 0
    with open(fastq_path, "r") as f:
        for _ in f:
            line_count += 1
    return line_count // 4

def parse_checkm2_results(results_file: str):
    good_bins = []
    with open(results_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            genome = row["Name"]
            completeness = float(row["Completeness"])
            contamination = float(row["Contamination"])
            if completeness >= CHECKM2_MIN_COMP and contamination <= CHECKM2_MAX_CONT:
                bin_name = genome.replace(".fa", "")
                good_bins.append((bin_name, completeness, contamination))
    return good_bins

# =====================
# MAIN
# =====================
start_cwd = os.getcwd()
mapping_summary = os.path.join(PROJECT_DIR, "sample_MAG_database_mapping_summary.tsv")

with open(mapping_summary, "w") as out_summary:
    out_summary.write("sample\ttotal_trimmed_pairs\tmapped_pairs\tpercent_pairs_mapped\n")

for sample in samples:
    print(f"\n=== Processing {sample} ===", flush=True)

    r1 = join_raw(f"{sample}_R1_001.fastq")
    r2 = join_raw(f"{sample}_R2_001.fastq")

    if not os.path.exists(r1) or not os.path.exists(r2):
        print(f"ERROR: Missing reads for {sample}. Expected: {r1} and {r2}", file=sys.stderr)
        continue

    sample_dir = os.path.join(PROJECT_DIR, sample) #changed to make sample directories OUTSIDE of raw reads
    trimmed_dir = os.path.join(sample_dir, "trimmed_reads")
    megahit_dir = os.path.join(sample_dir, "megahit_out")
    mags_dir = os.path.join(sample_dir, "MAGs")

    os.makedirs(trimmed_dir, exist_ok=True)
    os.makedirs(mags_dir, exist_ok=True)

    sickle_r1 = os.path.join(trimmed_dir, f"{sample}_R1_sickle_trimmed.fastq")
    sickle_r2 = os.path.join(trimmed_dir, f"{sample}_R2_sickle_trimmed.fastq")
    trimmed_r1 = os.path.join(trimmed_dir, f"{sample}_R1_trimmed_noplyG.fastq")
    trimmed_r2 = os.path.join(trimmed_dir, f"{sample}_R2_trimmed_noplyG.fastq")
    discarded = os.path.join(trimmed_dir, f"{sample}_discarded.fastq")

    # ========== Step 1: Trimming with sickle ==========
    try:
        subprocess.run([
            "sickle", "pe",
            "-f", r1,
            "-r", r2,
            "-t", "sanger",
            "-o", sickle_r1,
            "-p", sickle_r2,
            "-s", discarded
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: sickle failed for {sample}: {e}", file=sys.stderr)
        continue

    # ========== Step 1.2: Remove poly G tails & adapters with bbduk ==========
    try:
        subprocess.run([
            "bbduk.sh",
            "threads=10",
            "overwrite=t",
            f"in1={sickle_r1}",
            f"in2={sickle_r2}",
            "ref=/opt/bbtools/bbmap/resources/adapters.fa",
            "tpe",
            "tbo",
            f"trimpolygright={POLYG_TRIM}",
            f"out1={trimmed_r1}",
            f"out2={trimmed_r2}"
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: bbduk failed for {sample}: {e}", file=sys.stderr)
        continue

    # ========== Step 2: MEGAHIT assembly ==========
    try:
        subprocess.run([
            "megahit",
            "-1", trimmed_r1,
            "-2", trimmed_r2,
            "--k-min", KMER_MIN,
            "--k-max", KMER_MAX,
            "--k-step", KMER_STEP,
            "-m", "0.4",
            "-t", threads,
            "-o", megahit_dir
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: MEGAHIT failed for {sample}: {e}", file=sys.stderr)
        continue

    os.chdir(megahit_dir)

    # ========== Step 3: pullseq (>2500 bp) ==========
    filtered_scaffolds = os.path.abspath(os.path.join(megahit_dir, f"{sample}_{MIN_CONTIG_LEN}.fa"))
    try:
        with open(filtered_scaffolds, "w") as out_fa:
            subprocess.run([
                "pullseq",
                "-i", "final.contigs.fa",
                "-m", MIN_CONTIG_LEN
            ], stdout=out_fa, check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: pullseq failed for {sample}: {e}", file=sys.stderr)
        os.chdir(start_cwd)
        continue

    # ========== Step 4: BBMap mapping to contigs ==========
    try:
        subprocess.run([
            "bbmap.sh",
            f"-Xmx{memory}",
            f"threads={threads}",
            "minid=90",
            "overwrite=t",
            f"ref={filtered_scaffolds}",
            f"in1={trimmed_r1}",
            f"in2={trimmed_r2}",
            f"out={sample}_B_mapped.sam"
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: BBMap failed for {sample}: {e}", file=sys.stderr)
        os.chdir(start_cwd)
        continue

    # ========== Step 5: SAM -> BAM ==========
    try:
        with open(f"{sample}_B_mapped.bam", "wb") as out_bam:
            subprocess.run([
                "samtools", "view",
                "-@", threads,
                "-bS", f"{sample}_B_mapped.sam"
            ], stdout=out_bam, check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: samtools view failed for {sample}: {e}", file=sys.stderr)
        os.chdir(start_cwd)
        continue

    # ========== Step 6: Sort BAM ==========
    try:
        subprocess.run([
            "samtools", "sort",
            "-T", f"{sample}.sorted",
            "-o", f"{sample}_B_mapped.sorted.bam",
            f"{sample}_B_mapped.bam",
            "-@", threads
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: samtools sort failed for {sample}: {e}", file=sys.stderr)
        os.chdir(start_cwd)
        continue

    # ========== Step 6.1: Filter for high quality mapping ==========
    try:
        subprocess.run([
            "reformat.sh",
            f"-Xmx{memory}",
            f"minidfilter={BBMAP_MIN_ID}",
            f"in={sample}_B_mapped.sorted.bam",
            f"out={sample}_B_mapped99per.sorted.bam",
            "pairedonly=t",
            "primaryonly=t"
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: reformat.sh failed for {sample}: {e}", file=sys.stderr)
        os.chdir(start_cwd)
        continue

    # ========== Step 7: MetaBAT ==========
    try:
        subprocess.run([
            "runMetaBat.sh",
            filtered_scaffolds,
            f"{sample}_B_mapped99per.sorted.bam"
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: MetaBAT failed for {sample}: {e}", file=sys.stderr)
        os.chdir(start_cwd)
        continue

    # ========== Step 8: Assembly stats ==========
    try:
        subprocess.run([
            CONTIG_STATS,
            "-i", "final.contigs.fa",
            "-o", f"{sample}_final.contigs_STATS"
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: contig stats failed for {sample}: {e}", file=sys.stderr)
        os.chdir(start_cwd)
        continue

    # ========== Step 9: CheckM2 QC ==========
    bins_dir = f"{filtered_scaffolds}.metabat-bins"
    checkm2_out = os.path.join(megahit_dir, "checkm2_v1.1.0_fa")
    checkm_results = os.path.join(checkm2_out, "quality_report.tsv")

    try:
        subprocess.run([
            "bash", "-c",
            f"""
            source /home/opt/Miniconda3/miniconda3/bin/activate checkm2_v1.1.0
            checkm2 predict \
                -x fa \
                --input {bins_dir} \
                --output-directory {checkm2_out} \
                --threads {CHECKM2_THREADS}
            """
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: CheckM2 failed for {sample}: {e}", file=sys.stderr)
        os.chdir(start_cwd)
        continue

    # ========== Step 10: Copy/rename MQ + HQ MAGs ==========
    try:
        good_bins = parse_checkm2_results(checkm_results)

        for bin_name, completeness, contamination in good_bins:
            src_bin = os.path.join(bins_dir, f"{bin_name}.fa")
            dst_bin = os.path.join(mags_dir, f"{sample}_{bin_name}.fa")
            if os.path.exists(src_bin):
                shutil.copy2(src_bin, dst_bin)

        if len(good_bins) == 0:
            print(f"WARNING: No MQ/HQ MAGs found for {sample}", file=sys.stderr)

    except Exception as e:
        print(f"ERROR: MAG selection/renaming failed for {sample}: {e}", file=sys.stderr)
        os.chdir(start_cwd)
        continue

    # ========== Step 11: GTDB-Tk taxonomy ==========
    try:
        mag_files = [x for x in os.listdir(mags_dir) if x.endswith(".fa")]
        if len(mag_files) > 0:

            gtdb_out = os.path.join(sample_dir, "gtdb_v2.7.0_r232")

            os.makedirs(gtdb_out, exist_ok=True)

            subprocess.run([
                "bash", "-c",
                f"""
                source /home/opt/Miniconda3/miniconda3/bin/activate gtdbtk_v2.7.0
                gtdbtk classify_wf \
                    -x fa \
                    --genome_dir {mags_dir} \
                    --out_dir {gtdb_out} \
                    --cpus {GTDBTK_CPUS}
                """
            ], check=True)

    except subprocess.CalledProcessError as e:
        print(f"ERROR: GTDB-Tk failed for {sample}: {e}", file=sys.stderr)
        os.chdir(start_cwd)
        continue

    # ========== Step 12: Cleanup contig-mapping files ==========
    try:
        sam_file = os.path.join(megahit_dir, f"{sample}_B_mapped.sam")
        sorted_bam = os.path.join(megahit_dir, f"{sample}_B_mapped.sorted.bam")
        filt_bam = os.path.join(megahit_dir, f"{sample}_B_mapped99per.sorted.bam")
        bam_file = os.path.join(megahit_dir, f"{sample}_B_mapped.bam")
        bam_gz = bam_file + ".gz"

        if os.path.exists(sam_file):
            os.remove(sam_file)
        if os.path.exists(sorted_bam):
            os.remove(sorted_bam)
        if os.path.exists(filt_bam):
            os.remove(filt_bam)

        if os.path.exists(bam_file) and not os.path.exists(bam_gz):
            with open(bam_file, "rb") as f_in, gzip.open(bam_gz, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(bam_file)

    except Exception as e:
        print(f"ERROR: cleanup failed for {sample}: {e}", file=sys.stderr)
        os.chdir(start_cwd)
        continue

    # ========== Step 13: dRep within sample ==========
    drep_outdir = os.path.join(sample_dir, "dRep_v3.0.0_MAGs")
    try:
        mag_files = [x for x in os.listdir(mags_dir) if x.endswith(".fa")]
        if len(mag_files) > 0:
            subprocess.run([
                "bash", "-c",
                f"dRep dereplicate {drep_outdir} -p {threads} -comp {DREP_MIN_COMP} -con {DREP_MAX_CONT} -g {mags_dir}/*.fa"
            ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: dRep failed for {sample}: {e}", file=sys.stderr)
        os.chdir(start_cwd)
        continue

    # ========== Step 14: Build sample MAG database from within-sample dereplicated HQ/MQ MAGs ==========
    derep_genomes_dir = os.path.join(drep_outdir, "dereplicated_genomes")
    renamed_db_dir = os.path.join(derep_genomes_dir, "genome_renamed")
    bowtie_db_dir = os.path.join(sample_dir, "bowtie_DB")                    # new
    mag_db_fa = os.path.join(bowtie_db_dir, f"{sample}_derep_MAGs.fa")       # moved here from renamed_db_dir

    try:
        if os.path.exists(derep_genomes_dir):

            # Remove renamed_db_dir if it exists from a previous run
            if os.path.exists(renamed_db_dir):
                shutil.rmtree(renamed_db_dir)

            os.makedirs(bowtie_db_dir, exist_ok=True)                        # new

            # Single call with quoted glob - lets rename_bins_like_dram.py handle expansion internally
            subprocess.run([
                "bash", "-c",
                f"source /opt/Miniconda2/miniconda2/bin/activate scripts && "
                f"rename_bins_like_dram.py -i '{derep_genomes_dir}/*.fa' -o {renamed_db_dir}"
            ], check=True)

            # Build file list before opening output to avoid self-concatenation
            db_basename = f"{sample}_derep_MAGs.fa"
            renamed_fas = [
                os.path.join(renamed_db_dir, x) for x in os.listdir(renamed_db_dir)
                if x.endswith(".fa") and x != db_basename
            ]

            with open(mag_db_fa, "wb") as out_f:                             # now writes to bowtie_DB/
                for fa_file in renamed_fas:
                    with open(fa_file, "rb") as in_f:
                        shutil.copyfileobj(in_f, out_f)

        else:
            print(f"ERROR: dereplicated_genomes dir not found for {sample}", file=sys.stderr)
            continue

    except subprocess.CalledProcessError as e:
        print(f"ERROR: MAG database build failed for {sample}: {e}", file=sys.stderr)
        os.chdir(start_cwd)
        continue

    os.chdir(start_cwd)
    print(f"=== Finished {sample} ===", flush=True)

print("\nAll per-sample steps complete. Beginning cross-sample steps.", flush=True)

# =============================================================================
# CROSS-SAMPLE STEPS (run once after all samples complete)
# Steps 15-17 build a single shared MAG database from all samples and compute
# between-sample relative abundance — every sample maps to the same reference
# so abundance values are directly comparable across samples.
# =============================================================================

# ========== Step 15: Cross-sample dRep and build shared MAG database ==========
# Collect all per-sample MAGs into one pool, dereplicate across samples,
# rename contig headers, and concatenate into one shared FASTA database.

cross_drep_dir = os.path.join(PROJECT_DIR, "cross_sample_dRep")
cross_derep_genomes_dir = os.path.join(cross_drep_dir, "dereplicated_genomes")
cross_renamed_dir = os.path.join(cross_derep_genomes_dir, "genome_renamed")
cross_bowtie_dir = os.path.join(PROJECT_DIR, "cross_sample_bowtie_DB")
cross_mag_db_fa = os.path.join(cross_bowtie_dir, "all_samples_derep_MAGs.fa")
cross_bowtie_prefix = os.path.join(cross_bowtie_dir, "all_samples_derep_MAG_DB")

os.makedirs(cross_bowtie_dir, exist_ok=True)

# Collect all per-sample MAG .fa files from every sample's MAGs/ directory
all_mag_files = []
for sample in samples:
    mags_dir = os.path.join(PROJECT_DIR, sample, "MAGs")
    if os.path.isdir(mags_dir):
        for f in os.listdir(mags_dir):
            if f.endswith(".fa"):
                all_mag_files.append(os.path.join(mags_dir, f))

if len(all_mag_files) == 0:
    print("ERROR: No MAG .fa files found across any sample MAGs/ directories. Cannot build cross-sample database.", file=sys.stderr)
else:
    print(f"Found {len(all_mag_files)} total MAGs across all samples for cross-sample dRep.", flush=True)

    # Build a space-separated glob string for dRep -g
    mag_glob_str = " ".join(all_mag_files)

    try:
        subprocess.run([
            "bash", "-c",
            f"dRep dereplicate {cross_drep_dir} -p {threads} -comp {DREP_MIN_COMP} -con {DREP_MAX_CONT} -g {mag_glob_str}"
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Cross-sample dRep failed: {e}", file=sys.stderr)
        cross_derep_genomes_dir = None

    if cross_derep_genomes_dir and os.path.exists(cross_derep_genomes_dir):
        # Rename contig headers for DRAM compatibility
        if os.path.exists(cross_renamed_dir):
            shutil.rmtree(cross_renamed_dir)
        try:
            subprocess.run([
                "bash", "-c",
                f"source /opt/Miniconda2/miniconda2/bin/activate scripts && "
                f"rename_bins_like_dram.py -i '{cross_derep_genomes_dir}/*.fa' -o {cross_renamed_dir}"
            ], check=True)
        except subprocess.CalledProcessError as e:
            print(f"ERROR: rename_bins_like_dram.py failed for cross-sample DB: {e}", file=sys.stderr)
            cross_renamed_dir = None

        if cross_renamed_dir and os.path.exists(cross_renamed_dir):
            # Concatenate all renamed MAGs into one shared database FASTA
            db_basename = "all_samples_derep_MAGs.fa"
            renamed_fas = [
                os.path.join(cross_renamed_dir, x) for x in os.listdir(cross_renamed_dir)
                if x.endswith(".fa") and x != db_basename
            ]
            with open(cross_mag_db_fa, "wb") as out_f:
                for fa_file in renamed_fas:
                    with open(fa_file, "rb") as in_f:
                        shutil.copyfileobj(in_f, out_f)
            print(f"Cross-sample MAG database written to: {cross_mag_db_fa}", flush=True)

# ========== Step 16: Build one shared Bowtie2 index; map ALL samples to it ==========
# Each sample's trimmed reads are mapped to the shared cross-sample database.
# This produces one BAM per sample, all using the same reference, so coverage
# values are directly comparable across samples.

if os.path.exists(cross_mag_db_fa):
    # Build the shared Bowtie2 index once
    try:
        subprocess.run([
            "bowtie2-build",
            cross_mag_db_fa,
            cross_bowtie_prefix,
            "--threads", BOWTIE2_THREADS
        ], check=True)
        print("Cross-sample Bowtie2 index built.", flush=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Bowtie2 index build failed for cross-sample DB: {e}", file=sys.stderr)

    # Map each sample's trimmed reads to the shared index
    all_filt_bams = []  # collected for CoverM in Step 17

    for sample in samples:
        trimmed_dir = os.path.join(PROJECT_DIR, sample, "trimmed_reads")
        trimmed_r1 = os.path.join(trimmed_dir, f"{sample}_R1_trimmed_noplyG.fastq")
        trimmed_r2 = os.path.join(trimmed_dir, f"{sample}_R2_trimmed_noplyG.fastq")

        if not os.path.exists(trimmed_r1) or not os.path.exists(trimmed_r2):
            print(f"WARNING: Trimmed reads not found for {sample}, skipping mapping.", file=sys.stderr)
            continue

        sample_map_dir = os.path.join(PROJECT_DIR, sample, "MAG_db_mapping")
        os.makedirs(sample_map_dir, exist_ok=True)

        sam_file    = os.path.join(sample_map_dir, f"{sample}_mapped_crossDB.sam")
        bam_file    = os.path.join(sample_map_dir, f"{sample}_mapped_crossDB.bam")
        filt_bam    = os.path.join(sample_map_dir, f"{sample}_mapped_crossDB_99id.bam")
        possort_bam = os.path.join(sample_map_dir, f"{sample}_mapped_crossDB_99id_POSSORT.bam")

        try:
            # Map reads to the shared cross-sample database
            subprocess.run([
                "bowtie2",
                "-D", "10",
                "-R", "2",
                "-N", "0",
                "-L", "22",
                "-i", "S,0,2.50",
                "-p", BOWTIE2_THREADS,
                "-x", cross_bowtie_prefix,
                "-S", sam_file,
                "-1", trimmed_r1,
                "-2", trimmed_r2
            ], check=True)

            # SAM -> BAM
            with open(bam_file, "wb") as out_bam:
                subprocess.run([
                    "samtools", "view",
                    "-bS", sam_file
                ], stdout=out_bam, check=True)

            # Filter for high-quality paired primary alignments (filter before sort)
            subprocess.run([
                "reformat.sh",
                f"-Xmx{memory}",
                f"minidfilter={BOWTIE2_MIN_ID}",
                f"in={bam_file}",
                f"out={filt_bam}",
                "pairedonly=t",
                "primaryonly=t"
            ], check=True)

            # Position sort the filtered BAM (POSSORT) — this is what goes into CoverM
            subprocess.run([
                "samtools", "sort",
                "-@", threads,
                "-o", possort_bam,
                filt_bam
            ], check=True)

            # Delete SAM, unsorted BAM, and unfiltered BAM; keep only POSSORT for CoverM
            if os.path.exists(sam_file):
                os.remove(sam_file)
            if os.path.exists(bam_file):
                os.remove(bam_file)
            if os.path.exists(filt_bam):
                os.remove(filt_bam)

            # Count mapped read pairs and write to project-level summary
            total_pairs = count_fastq_reads(trimmed_r1)
            mapped_reads = int(subprocess.check_output([
                "samtools", "view", "-c", "-F", "4", possort_bam
            ]).decode().strip())
            mapped_pairs = mapped_reads // 2
            percent_pairs_mapped = (mapped_pairs / total_pairs * 100) if total_pairs > 0 else 0.0
            with open(mapping_summary, "a") as out_summary:
                out_summary.write(
                    f"{sample}\t{total_pairs}\t{mapped_pairs}\t{percent_pairs_mapped:.2f}\n"
                )

            all_filt_bams.append(possort_bam)
            print(f"  {sample}: {percent_pairs_mapped:.2f}% pairs mapped to cross-sample DB", flush=True)

        except subprocess.CalledProcessError as e:
            print(f"ERROR: Cross-sample mapping failed for {sample}: {e}", file=sys.stderr)

# ========== Step 17: CoverM between-sample relative abundance ==========
# CoverM receives ALL samples' BAMs in a single call against the shared
# cross-sample genome directory. Output tables have one column per sample
# and one row per MAG — enabling direct between-sample comparisons.

if all_filt_bams:
    coverm_out_dir = os.path.join(PROJECT_DIR, "coverm_output")
    os.makedirs(coverm_out_dir, exist_ok=True)

    try:
        # reads per base — all samples, no breadth filter
        subprocess.run([
            "coverm", "genome",
            "--proper-pairs-only",
            "--genome-fasta-extension", "fa",
            "--genome-fasta-directory", cross_renamed_dir,
            "--bam-files", *all_filt_bams,
            "--threads", BOWTIE2_THREADS,
            "--min-read-percent-identity-pair", COVERM_MIN_ID,
            "--min-covered-fraction", "0",
            "-m", "reads_per_base",
            "--output-file", os.path.join(coverm_out_dir, "coverm_reads_per_base.txt")
        ], check=True,
           stderr=open(os.path.join(coverm_out_dir, "reads_per_base_stats.txt"), "w"))

        # min covered fraction — only report MAGs with >= COVERM_MIN_BREADTH coverage
        subprocess.run([
            "coverm", "genome",
            "--proper-pairs-only",
            "--genome-fasta-extension", "fa",
            "--genome-fasta-directory", cross_renamed_dir,
            "--bam-files", *all_filt_bams,
            "--threads", BOWTIE2_THREADS,
            "--min-read-percent-identity-pair", COVERM_MIN_ID,
            "--min-covered-fraction", COVERM_MIN_BREADTH,
            "--output-file", os.path.join(coverm_out_dir, "coverm_min75.txt")
        ], check=True,
           stderr=open(os.path.join(coverm_out_dir, "min75_stats.txt"), "w"))

        # trimmed mean — robust to coverage outliers
        subprocess.run([
            "coverm", "genome",
            "--proper-pairs-only",
            "--genome-fasta-extension", "fa",
            "--genome-fasta-directory", cross_renamed_dir,
            "--bam-files", *all_filt_bams,
            "--threads", BOWTIE2_THREADS,
            "--min-read-percent-identity-pair", COVERM_MIN_ID,
            "-m", "trimmed_mean",
            "--output-file", os.path.join(coverm_out_dir, "coverm_trimmed_mean.txt")
        ], check=True,
           stderr=open(os.path.join(coverm_out_dir, "trimmed_mean_stats.txt"), "w"))

        print(f"\nBetween-sample CoverM tables written to: {coverm_out_dir}", flush=True)

    except subprocess.CalledProcessError as e:
        print(f"ERROR: CoverM between-sample abundance failed: {e}", file=sys.stderr)

else:
    print("WARNING: No BAM files collected — CoverM between-sample step skipped.", file=sys.stderr)

print("\nAll samples processed (errors per-sample will be reported above).")

