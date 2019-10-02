import os


def index_bam_file(bam_file):
    index_cmd = "samtools index ./%s" % bam_file
    os.system(index_cmd)


def create_sample_files(bam_file_name):
    """create samples using samtools view command"""

    print("running sample retreival: \n\n")

    sample_folder = "./%s_sample_stats/" % bam_file_name.replace(
        "_mRNA.bam", "")
    if not os.path.exists(sample_folder):
        os.mkdir(sample_folder)
        os.chdir(sample_folder)
    else:
        os.chdir(sample_folder)
    for decimal in [.01, .1, .5, .9, .99]:
        print("%s: complete" % decimal)
        file_var = "%s_sample.bam" % str(decimal).replace("0.", "")

        # extract sample
        print(bam_file_name, file_var)
        sample_cmd = "samtools view ../%s -s $arg_var -b -@12 > %s" % (
            bam_file_name, file_var)
        print(sample_cmd)
        os.system(sample_cmd)


def create_flag_stats():
    """run flag stats command and create log.txt"""

    print("running flag stats: \n\n")

    create_log_cmd = """printf "decimal\tduped\tall\n" > log.txt"""
    os.system(create_log_cmd)
    for decimal in [.01, .1, .5, .9, .99]:
        print("%s: complete" % decimal)
        file_var = "%s_sample.bam" % str(decimal).replace("0.", "")

        cmd_1 = "printf '. % s \t' >> log.txt" % decimal
        cmd_2 = "samtools flagstat %s - @12 | awk 'NR==4{print $1}' | cat | xargs > file1.tmp" % file_var
        cmd_3 = "samtools flagstat %s - @12 | awk 'NR==5{print $1}' | cat | xargs > file2.tmp" % file_var
        cmd_4 = """paste file1.tmp file2.tmp >> log.txt"""

        os.system(cmd_1)
        os.system(cmd_2)
        os.system(cmd_3)
        os.system(cmd_4)

        # remove all tmp files
        os.system("rm *.tmp")


if __name__ == "__main__":
    # include standard modules
    import argparse
    # initiate the parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--filename", "-f", help="bam_file_name")
    args = parser.parse_args()
    # read arguments from the command line
    bam_file = args.filename
    create_sample_files(bam_file)
    create_flag_stats()
