#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --time=4-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=500gb
#SBATCH --mail-user=ayda.lewis@colostate.edu
#SBATCH --partition=borton-hi,borton-low

# Copy raw reads directory
#  cp -r /ORG-Data-phoenix/Manure_Lagoons/seq_CSU_Dec2025 /home/projects-phoenix/Manure_Lagoons/MetaG

# Rename directory to raw_reads
#  mv /home/projects-phoenix/Manure_Lagoons/MetaG/seq_CSU_Dec2025 /home/projects-phoenix/Manure_Lagoons/MetaG/raw_reads

# Unzip all FASTQ files
#  gunzip /home/projects-phoenix/Manure_Lagoons/MetaG/raw_reads/*.fastq.gz

cd /home/projects-phoenix/Manure_Lagoons/MetaG/raw_reads

echo "Beginning pipeline at $(date)"

python3 /home/projects-phoenix/Manure_Lagoons/MetaG/scripts/reads_to_MAGs_pipeline_UPDATED.py \
-s "M10_S50,M11_S5,M12_S45,M13_S30,M14_S46,M15_S6,M16_S31,M3_S38,M4_S26,M5_S27,M6_S3,M7_S28,M8_S4,M9_S29"

echo "Pipeline completed at $(date)"