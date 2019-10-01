import os


def create_sample_files(bam_file_name):
    """create samples using samtools view command"""

    print("running sample retreival: \n\n")
#    os.mkdir("trial")

    for decimal in [.01, .1, .5, .9, .99]:
        print("%s: complete") % decimal
        file_var = "%d_sample.bam" % decimal

        # extract sample
        sample_cmd = "samtools view ../%s -s $arg_var -b -@12 > %s" % bam_file_name, file_var

        os.system(sample_cmd)


def create_flag_stats():
    """run flag stats command and create log.txt"""

    print("running flag stats: \n\n")
    for decimal in [.01, .1, .5, .9, .99]:
        print("%s: complete") % decimal
        file_var = "%d_sample.bam" % decimal

        cmd_1 = "printf '. % s \t' >> log.txt" % decimal
        cmd_2 = "samtools flagstat %s - @12 | awk 'NR==4{print $1}' | cat | xargs > file1.tmp" % file_var
        cmd_3 = "samtools flagstat %s - @12 | awk 'NR==5{print $1}' | cat | xargs > file2.tmp" % file_var
        cmd_4 = "paste file1.tmp file2.tmp >> log.txt"

        os.system(cmd_1)
        os.system(cmd_2)
        os.system(cmd_3)
        os.system(cmd_4)


create_sample_files("SIGAF1_144520_mRNA.bam")
create_flag_stats()
