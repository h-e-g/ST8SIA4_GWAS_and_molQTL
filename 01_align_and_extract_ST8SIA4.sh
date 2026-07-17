#!/bin/bash
################################################################################
################################################################################
# File name: 01_align_and_extract_ST8SIA4.sh
# Author: M. Rotival
################################################################################
################################################################################

#SBATCH --job-name=align_and_extract
#SBATCH --output=logfile_alignST8SIA4_%A_%a.log
#SBATCH --array=1-970
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4         # Number of CPUs per task
#SBATCH --mem=40G                # Memory per task
#SBATCH --time=02:00:00           # Time limit
#SBATCH --qos=fast
#SBATCH --mail-type=END

cd ./data || exit

ST8SIA4_b38=chr5:99806933-101903282
### directory containg the hg38 reference geneome for STAR
REFDIR=HG38__2024/star/ 

FASTQ_DIR=/data/02_fastq
mapfile -t FASTQ_LIST <<<$(ls $FASTQ_DIR | grep ".fq.gz$")

BAM_DIR=./data/03_bam
module load STAR/2.7.3a

NTHREADS=${SLURM_CPUS_PER_TASK:-1}
TMP_DIR=$(mktemp -d -p .//tmp)

FASTQ=${FASTQ_LIST[$SLURM_ARRAY_TASK_ID]}
echo "$FASTQ"
gunzip -c ${FASTQ_DIR}/"${FASTQ}" > "$TMP_DIR"/"${FASTQ%.gz}"
STAR --runThreadN "$NTHREADS" --genomeDir ${REFDIR} --readFilesIn "$TMP_DIR"/"${FASTQ%.gz}" --outFileNamePrefix "${TMP_DIR}"/"${FASTQ%.fq.gz}"_ --outSAMtype BAM SortedByCoordinate
rm "$TMP_DIR"/"${FASTQ%.gz}"

module load samtools/1.21
BAMFILE=${TMP_DIR}/${FASTQ%.fq.gz}.sortedByCoord.out.bam
samtools index "${BAMFILE}"
samtools view -b "${BAMFILE}" ${ST8SIA4_b38}  > ${BAM_DIR}/"${FASTQ%.fq.gz}".sortedByCoord.out.ST8SIA4_b38.bam
rm "$BAMFILE"
rm "${BAMFILE}".bai
samtools index ${BAM_DIR}/"${FASTQ%.fq.gz}".sortedByCoord.out.ST8SIA4_b38.bam

echo 'All done\n'
exit 0