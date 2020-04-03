CTRL_BAM=$1
GENDER_FLAG=$2
BAMFILES_DIR=$PWD
for BAMFILE in `ls *.bam`
	do
		if [ "$BAMFILE" != "$CTRL_BAM" ];
	then sbatch /home/rcf-40/grossfel/tconut_scripts/tconut_pipeline.slurm $BAMFILE $CTRL_BAM $BAMFILES_DIR $GENDER_FLAG;
		fi
done