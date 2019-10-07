import glob
import os
import sys
from multiprocessing import Pool


def main(folder, multi):

    os.chdir(folder)
    bam_files = glob.glob('./*.bam')
    os.system("touch ./log.txt")

    # rint(multi)
    if multi.lower() == "true":
        print("multi processing on these files: ", bam_files)

        p = Pool(5)
        p.map(call_cmd, bam_files)

    else:
        for filename in bam_files:
            call_cmd(filename)


def call_cmd(filename):

    cmd = "samtools view %s | python ~/bioinformatics_scripts/umi_dict.py %s " % (
        filename, filename)

    os.system(cmd)


if __name__ == "__main__":
    folder = sys.argv[1]
    multi = sys.argv[2]
    main(folder, multi)
