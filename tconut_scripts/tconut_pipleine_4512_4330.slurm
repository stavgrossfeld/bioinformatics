#!/bin/bash

#SBATCH --time=3-0:00:00
#SBATCH --mail-user=grossfel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=/home/dtg-00/Groups/Hjelm_lab/Grossfel/slurm_logs/%x\_%j.out
#SBATCH --job-name=tc_no_hc_4512
##SBATCH --nodes=1
#SBATCH --ntasks=3
##SBATCH --exclusive
#SBATCH --partition=dtg
#SBATCH --account=lc_dwc
#SBATCH --mem=12G
#SBATCH --exclude=hpc3969,hpc3970,hpc3971,hpc3972,hpc3973,hpc3974,hpc3975,hpc3976,hpc3977,hpc3978

PROJECT_BACKUP_DIR="/home/dtg-00/Groups/Hjelm_lab/Grossfel/"


SAMPLE_BAM=46__4330.bam
CTRL_BAM=65__4512.bam
GENDER_FLAG=M

BAMFILES_DIR="/staging/dwc/grossfel/wes_bams_working/"
cd $BAMFILES_DIR


SAMPLE_NAME=`echo $SAMPLE_BAM | cut -d"_" -f3 | cut -d"." -f1`


CTRL_NAME_NO_FLAG=`echo $CTRL_BAM | cut -d"_" -f3 | cut -d"." -f1`
CTRL_NAME=$CTRL_NAME_NO_FLAG$GENDER_FLAG

# HAPLOTYPE CALLING

DBSNP=/home/dtg-00/Resources/dtgref/pecan/annotations/public_data_bases/DBSNP_V147/dbsnp_147.b37_norm.vcf
REF=/home/dtg-00/Resources/dtgref/pecan/bwa_index/hs37d5_plusRibo_plusOncoViruses_plusERCC.fa
GATKPATH=/home/dtg-00/Packages/gatk/4.0.10.1
tCoNuTdir=/home/dtg-00/Packages/tCoNuT/1.0
SAMTOOLS=/home/dtg-00/Packages/samtools/1.2/bin/samtools
MCRPATH=/home/dtg-00/Packages/MCR/9.0/v90





TRK=$SAMPLE_NAME"_"$CTRL_NAME
echo $TRK

OUTDIR_CTRL="$BAMFILES_DIR/TCONUT/CTRL_"$CTRL_NAME
OUTDIR_VCF="$BAMFILES_DIR/TCONUT/"$CTRL_NAME
OUTDIR_SAMPLE="$BAMFILES_DIR/TCONUT/CTRL_"$CTRL_NAME/$SAMPLE_NAME

BACKUP_DIR_CTRL="$PROJECT_BACKUP_DIR/TCONUT/CTRL_"$CTRL_NAME
BACKUP_DIR_VCF="$PROJECT_BACKUP_DIR/TCONUT/"$CTRL_NAME
BACKUP_OUTDIR_SAMPLE="$PROJECT_BACKUP_DIR/TCONUT/CTRL_"$CTRL_NAME/$SAMPLE_NAME


mkdir -p $OUTDIR_CTRL
mkdir -p $OUTDIR_SAMPLE
mkdir -p $OUTDIR_VCF

mkdir -p $BACKUP_DIR_CTRL
mkdir -p $BACKUP_DIR_VCF
mkdir -p $BACKUP_OUTDIR_SAMPLE

# HAPLOTYPE CALLING 

${GATKPATH}/gatk --java-options "-Xmx10g" HaplotypeCaller \
    -I ${CTRL_BAM} -I ${SAMPLE_BAM}\
    --reference ${REF} \
    --dbsnp ${DBSNP} \
    --min-base-quality-score 10 \
    --output $OUTDIR_VCF/${TRK}.HC.vcf


# ANNOTATION
cd $OUTDIR_VCF

SNPEFFPATH=/home/dtg-00/Resources/dtgref/pipeline_v0.4/bin/snpEff/snpEff_3.5h
VCF=$TRK.HC.vcf
DB=GRCh37.74
## snpEff
java -Xmx8g -jar ${SNPEFFPATH}/snpEff.jar \
	GRCh37.74 \
	$VCF \
   -verbose \
    -noLog \
    -stats ${summaryOut} \
    -noDownload \
    -config ${SNPEFFPATH}/snpEff.config \
    > $OUTDIR_VCF/$TRK.snpeff.txt
## snpsift
java -Xmx8g -jar $SNPEFFPATH/SnpSift.jar annotate $DBSNP $OUTDIR_VCF/$TRK.snpeff.txt > $OUTDIR_VCF/$TRK.HC_ALL.vcf

rm -f ./-noDownload* $TRK.snpeff.txt #$TRK.HC.vcf

# STEP 2: tconut clone covariance for both control.bam and sample.bam - produces .dat files
cd $BAMFILES_DIR

SAMTOOLS=/home/dtg-00/Packages/samtools/1.2/bin/samtools
OUTFILE=$OUTDIR_SAMPLE/$CTRL_NAME
perl ${tCoNuTdir}/tgen_CloneCov.pl I=$BAMFILES_DIR/${CTRL_BAM} O=${OUTFILE} M=RG: S=${SAMTOOLS}
OUTFILE=$OUTDIR_SAMPLE/$SAMPLE_NAME
perl ${tCoNuTdir}/tgen_CloneCov.pl I=$BAMFILES_DIR/${SAMPLE_BAM} O=${OUTFILE} M=RG: S=${SAMTOOLS}


# STEP 3: parse merge vcf

cd $OUTDIR_SAMPLE

HCVCF=$TRK.HC_ALL.vcf
CTRL_NAME_VCF=$CTRL_NAME_NO_FLAG"_R1_001"
SAMPLE_NAME_VCF=$SAMPLE_NAME"_R1_001"
SNPDEPTH=50 	#	<<<< This can be modified for specific cases depending on target coverage. Another option would be to use ${TUMORX} from hsMetrics above

perl ${tCoNuTdir}/parseMergeVCF.pl $OUTDIR_VCF/${HCVCF} ${CTRL_NAME_VCF} ${SAMPLE_NAME_VCF} ${SNPDEPTH}
mv $CTRL_NAME_VCF"-"$SAMPLE_NAME_VCF".baf.txt" $TRK.baf.txt
mv merged.vcf.txt $TRK.merged.vcf



# constants for next steps
OFILE=$TRK
BAF=$TRK.baf.txt



# STEP 4:

assayID="Exome"

##
## Parameters for tCoNuT.  Currently set for EXOME data. Please ngs_cna2015_WG.pbs for suggested whole genome parameters.
##
NORMALDAT=$CTRL_NAME.cln.dat
TUMORDAT=$SAMPLE_NAME.cln.dat
HETFILE=$TRK.merged.vcf   #   from parseMergeVCF.pl
TARGETSFILE=/home/dtg-00/Resources/dtgref/pipeline_v0.4/cna/exome.overlap.bed # Copy number BED file of Agilent Clinic Research exome targets

smWin=6                 #   <<<< THIS CAN BE ADJUSTED - smoothing window size >>>>
fcThresh=0.75           #   <<<< THIS CAN BE ADJUSTED - fold-change threshold - plot >>>>
res=2                   #   <<<< THIS CAN BE ADJUSTED - min resolution >>>>
maxGap=1000             #   <<<< THIS CAN BE ADJUSTED - max distance between consecutive postions >>>>

##
## Average Coverage calculated Picard hsMetrics can be used for minimum depth (${NORMX} and ${TUMORX})
##
hetDepthN=50 #${NORMX}      #   <<<< THIS CAN BE ADJUSTED - min depth of diploid het position >>>>
hetDepthT=50 #${TUMORX}     #   <<<< THIS CAN BE ADJUSTED - min depth of diploid het position >>>>
hetDev=0.025            #   <<<< THIS CAN BE ADJUSTED - allowable deviation from ref allele frequency of 0.5+/-0.025

readDepth=$( echo "$hetDepthN * 3" | bc )  # Set max read depth to 3 x control sample's average target coverage (from Picard HS metrics)

${tCoNuTdir}/tCoNuT/run_tCoNuT.sh ${MCRPATH} ${NORMALDAT} ${TUMORDAT} ${OFILE} ${HETFILE} ${smWin} ${fcThresh} ${assayID} ${res} ${readDepth} ${maxGap} ${hetDepthN} ${hetDepthT} ${hetDev} ${TARGETSFILE}




# STEP 5:

### Copy Number ###
Rscript --vanilla ${tCoNuTdir}/segmentation/runDNAcopyExomeV2.R ${OFILE}.cna.tsv ${OFILE}.seg

### BAF ###
 #4512_R1_001-4464_R1_001.baf.txt
Rscript --vanilla ${tCoNuTdir}/segmentation/runDNAcopyBAF.R ${BAF} ${OFILE}.baf


# STEP 6:

###COPY NUMBER####
##Annotate and convert SEG file to gVCF


ZTABLE=ztable.txt
CCDSLIST=/home/dtg-00/Groups/Hjelm_lab/references/Homo_sapiens.GRCh37.74.gtf.hs37d5.EGFRvIII_CCDS.gtf
DUPTHRESH=0.58     #   <<<< THIS CAN BE ADJUSTED - Amplification Threshold - log2 fold-change >>>>
DELTHRESH=-0.99    #   <<<< THIS CAN BE ADJUSTED - Deletion Threshold - log2 fold-change >>>>

perl ${tCoNuTdir}/annotSeg.pl ${CCDSLIST} ${OFILE}.cna.seg ${DUPTHRESH} ${DELTHRESH}

perl ${tCoNuTdir}/validateCNAVariantsVCF.pl ${OFILE}.cna.seg.vcf ${BAF} ${ZTABLE}




# step 7:

## Plotting of copy number



Rscript --vanilla ${tCoNuTdir}/plotting/plotCGH_EXOME.R ${OFILE}.cna.tsv ${OFILE}.amp.tsv ${OFILE}.del.tsv ${OFILE}

## Plotting of copy number with hets superimposed
if [ -f ${OFILE}.hets.tsv ];then
        Rscript --vanilla ${tCoNuTdir}/plotting/plotCGHwithHets.R ${OFILE}.cna.tsv ${OFILE}.amp.tsv ${OFILE}.del.tsv ${OFILE}.hets.tsv ${OFILE}_withhets
fi

## Plotting of BAF
Rscript --vanilla ${tCoNuTdir}/plotting/plotBAF.R $BAF ${OFILE}.baf

## Linear Genome Plotting of CNA and BAF (both scripts are MATLAB code)
${tCoNuTdir}/plotting/plotLinearCNA/run_plotLinearCNAandBAF.sh ${MCRPATH} ${OFILE}.cna.tsv $BAF ${OFILE}.cnaBAF.png
${tCoNuTdir}/plotting/plotLinearCNAabsBAF/run_plotLinearCNAandAbsBAF.sh ${MCRPATH} ${OFILE}.cna.tsv $BAF ${OFILE}.cnaAbsBAF.png


# copy the vcf that takes 14 hours to produce - haplotype caller

cp -r $OUTDIR_VCF/$TRK.HC_ALL.vcf $BACKUP_DIR_VCF/.

# copy the tconut produced files for the sample
cp -r $OUTDIR_SAMPLE/* $BACKUP_OUTDIR_SAMPLE/.
