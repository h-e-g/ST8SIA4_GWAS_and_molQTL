#!/bin/bash
################################################################################
################################################################################
# File name: 02_run_ggsashimi_per_sample.sh
# Author: M. Rotival
################################################################################
################################################################################

#SBATCH --job-name=run_ggsashimi_per_sample
#SBATCH --output=logfile_run_ggsashimi_per_sample_%A_%a.log
#SBATCH --array=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1         # Number of CPUs per task
#SBATCH --mem=100G                # Memory per task
#SBATCH --time=01:00:00           # Time limit
#SBATCH --qos=fast


SCRIPT_DIR="./scripts"
OUT_DIR="./data/04_coverage"
ST8SIA4_b38_coords="chr5:10796933-100913282"
GTF_annot="./data/07_ST8SIA4_gene_structure/genes_hg38.gtf"

module load R/4.1.0
module load Python/3.8.3
chmod 770 $SCRIPT_DIR/ggsashimi.py
python3 $SCRIPT_DIR/ggsashimi.py --help
export GGSASHIMI_DEBUG=yes


BAMfile=$(head -n $SLURM_ARRAY_TASK_ID ./data/03_bam/all_970_bams.tsv | tail -n 1 | cut -f 1)
IID=$(echo ${BAMfile##*/})
IID=${IID%%_aligned.sortedByCoord.out*}

python3 $SCRIPT_DIR/ggsashimi.py --bam ${BAMfile} \
  --coordinates ${ST8SIA4_b38_coords} --gtf $GTF_annot \
  --out-prefix ${OUT_DIR}/sashimi_ST8SIA4_${IID}_ \
  --min-coverage=5 \
  --out-format pdf --out-resolution 400 --base-size 14 \
  --width 10 --height 2 --ann-height 1.5
