import os
for decimal in [.01, .1, .5, .9, .99]:
    print(decimal)
    file_var = "%d_sample.bam" % decimal

    # extract sample
    sample_cmd = "samtools view ../SIGAF1_144520_mRNA.bam -s $arg_var -b -@12 > %s" % file_var


for decimal in [.01, .1, .5, .9, .99]:
    print(decimal)
    file_var = "%d_sample.bam" % decimal
    sample_cmd = "samtools view ../SIGAF1_144520_mRNA.bam -s $arg_var -b -@12 > %s" % file_var

    cmd_1 = "printf '. % s \t' >> log.txt" % decimal
    cmd_2 = "samtools flagstat %s - @12 | awk 'NR==4{print $1}' | cat | xargs > file1.tmp" % file_var
    cmd_2 = "samtools flagstat %s - @12 | awk 'NR==5{print $1}' | cat | xargs > file2.tmp" % file_var
    cmd_3 = "paste file1.tmp file2.tmp >> log.txt"
