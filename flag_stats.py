import os

import pandas as pd
import numpy as np
from bam_reader_simple import main as bam_reader_simple


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
        sample_cmd = "samtools view ../%s -s %s -@12 > %s" % (
            bam_file_name, decimal, file_var)
        print(sample_cmd)
        os.system(sample_cmd)


def loop_over_samples():
    for i in os.listdir():
        if "sample" in i:
            ct_umis(i)


def ct_umis(bam_file_name):
    # | python ~/bioinformatics_scripts/bam_reader_simple.py"
    ct_umis_cmd = "samtools view 01_sample.bam &"
    os.system(ct_umis_cmd)
    bam_reader_simple()


if __name__ == "__main__":
    # include standard modules
    import argparse
    # initiate the parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--filename", "-f", help="bam_file_name")
    parser.add_argument("--sample", "-s", default=False, type=bool,
                        help="create random samples")
    parser.add_argument("--run-stats", "-r", default=False, type=bool,
                        help="run flagstats")

    args = parser.parse_args()
    # read arguments from the command line
    bam_file = args.filename.strip()

    os.chdir("/Users/stav/Desktop/bam_split")
    if args.sample == True:
        create_sample_files(bam_file)
    if args.run_stats == True:
        ct_umis(bam_file)
