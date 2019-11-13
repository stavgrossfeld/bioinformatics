import glob
import os
import sys
import subprocess
import re
import pandas as pd

regex_dict = {'mapped': "(\d+) \+ (\d+) mapped \((.+):(.+)\)"}


def main(folder):
    #file_path = os.path.abspath(folder)
    file_path = '/Users/stav/Desktop/bam_split'
    os.chdir(file_path)

    bam_files = glob.glob('./*.bam')

    mapped_dict = {}
    for sample in bam_files:
        print(sample)

        cmd = "samtools flagstat %s" % (sample)
        flag_stats = subprocess.check_output(cmd, shell=True).decode("utf-8")

        mapped_reads = re.search(regex_dict["mapped"], flag_stats)

        mapped = mapped_reads[0].split(" ")[0]

        mapped_dict[sample] = mapped

    mapped_df = pd.DataFrame(pd.Series(mapped_dict))

    print(mapped_df)


if __name__ == "__main__":
    #folder = sys.argv[1]
    folder = "~/Desktop/bam_split/"
    #multi = sys.argv[2]
    main(folder)
