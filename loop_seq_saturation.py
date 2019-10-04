import glob
import os
import sys
from multiprocessing import Pool


def main(folder):

    os.chdir(folder)
    bam_files = glob.glob('./*.bam')

    number_of_workers = 4
    with Pool(number_of_workers) as p:
        # Do something with pool here
        for filename in bam_files:

            os.system("touch ./log.txt")
            cmd = "samtools view %s | python ~/bioinformatics_scripts/umi_dict.py %s " % (
                filename, filename)

            os.system(cmd)


if __name__ == "__main__":
    folder = sys.argv[1]
    main(folder)
