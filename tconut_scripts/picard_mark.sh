for BAMFILE in `ls *.bam`
	do sbatch ~/tconut_scripts/picard_mark_dup.slurm $BAMFILE
done