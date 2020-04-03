#!/bin/bash

#SBATCH --time=0-04:00:00
#SBATCH --mail-user=grossfel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=/home/dtg-00/Groups/Hjelm_lab/Grossfel/slurm_logs/%x\_%j.out
#SBATCH --job-name=picard_mark_duplicates
##SBATCH --nodes=1
#SBATCH --ntasks=8
##SBATCH --exclusive
##SBATCH --partition=main
#SBATCH --account=lc_dwc
#SBATCH --mem=30G
##SBATCH --exclude=hpc3969,hpc3970,hpc3971,hpc3972,hpc3973,hpc3974,hpc3975,hpc3976,hpc3977,hpc3978



BAMFILE=$1
THREADS=8


PICARDPATH=/home/dtg-00/Packages/picard-tools/picard-2.18.22
REF=/home/dtg-00/Resources/dtgref/pecan/bwa_index/hs37d5_plusRibo_plusOncoViruses_plusERCC.fa

WES_BAMS=/home/dtg-scratch-01/Projects/wesRepo/

cd $WES_BAMS

OUTDIR=/staging/dwc/grossfel/wes_bams_working/

METRICS="METRICS"
mkdir -p $OUTDIR/$METRICS

echo $BAMFILE

java -Xmx22g -XX:ParallelGCThreads=$THREADS -jar ${PICARDPATH}/picard.jar MarkDuplicates \
ASSUME_SORTED=true REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT \
TMP_DIR=/tmp INPUT=${BAMFILE} OUTPUT=$OUTDIR/${BAMFILE} METRICS_FILE=$OUTDIR/$METRICS/${BAMFILE}.picStats.MarkDupMetrics \
MAX_RECORDS_IN_RAM=18000000 CREATE_INDEX=true