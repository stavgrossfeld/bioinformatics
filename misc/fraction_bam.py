""" split bam files into fractions in a sample directory

# command for fractioning bams
python ~/bioinformatics_scripts/misc/fraction_bam.py -f SIGAB2_144100_mRNA.bam -s 

 """

import os

import pandas as pd
import numpy as np


# def index_bam_file(bam_file):
#     index_cmd = "samtools index ./%s" % bam_file
#     os.system(index_cmd)


def create_sample_files(bam_file_name, fractions):
    """create samples using samtools view command"""

    print("running sample retreival: \n\n")

    sample_folder = "./%s_sample_stats/" % bam_file_name.replace(
        "_mRNA.bam", "")

    if not os.path.exists(sample_folder):
        os.mkdir(sample_folder)
        os.chdir(sample_folder)
    else:
        os.chdir(sample_folder)
    print("fractions: ", fractions)
    for decimal in fractions:
        print("%s: started" % decimal)
        file_var = "%s_sample.bam" % str(decimal).replace("0.", "")

        # extract sample
        print(bam_file_name, file_var)
        sample_cmd = "samtools view ../%s -s %s -b -@12 > %s" % (
            bam_file_name, '{:f}'.format(decimal), file_var)
        print(sample_cmd)
        os.system(sample_cmd)


if __name__ == "__main__":
    # include standard modules
    import argparse
    # initiate the parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--filename", "-f", help="bam_file_name")

    parser.add_argument("--fractions", "-s", type=float, default=[.01, .1, .5, .9, .99], nargs="+",
                        help="create random samples - input: --fractions .001 .0001")

    args = parser.parse_args()
    # read arguments from the command line
    bam_file = args.filename.strip()
    fractions = args.fractions

    create_sample_files(bam_file, fractions)
