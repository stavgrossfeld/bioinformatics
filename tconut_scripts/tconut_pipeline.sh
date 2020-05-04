CTRL_BAM=$1
GENDER_FLAG=$2
BAMFILES_DIR=$PWD
#for BAMFILE in `ls *.bam`
for BAMFILE in `cat ../samples_to_run.txt`
	do
		if [ "$BAMFILE" != "$CTRL_BAM" ];
	then sbatch ~/bioinformatics/tconut_scripts/tconut_pipeline.slurm $BAMFILE $CTRL_BAM $BAMFILES_DIR $GENDER_FLAG;
		fi
done