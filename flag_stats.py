import os

import pandas as pd
import numpy as np


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


def create_barcodes_cut(file_name):
    print("creating barcodes cut file: ")
    create_urz_cmd = "samtools view %s | grep UR:Z: | cut -f19-27 >> urz_%s.txt" % (
        file_name, file_name.replace(".bam", ""))
    os.system(create_urz_cmd)
    urz_file_name = "urz_%s.txt" % (file_name.replace(".bam", ""))

    return urz_file_name


def ct_umis(urz_file_name):
    umi_cts = {}

    print("creating seq saturation numbers cut file: ")

    with open(urz_file_name, 'r') as f:
        for line in f:
            read = line.replace("\n", "")
            if "UR:Z:" in read:

                if read not in umi_cts:
                    umi_cts[read] = 1
                else:
                    umi_cts[read] += 1

    umi_df = pd.DataFrame(list(umi_cts.items()), columns=['barcode', 'ct'])
    umi_df.to_csv("./umi_df_sample.csv", index=False)

    sequence_saturation = 1 - (umi_df.shape[0] / np.sum(umi_df.ct))
    print(sequence_saturation)


def get_saturation(file_name):
    urz_file_name = create_barcodes_cut(file_name)

    ct_umis(urz_file_name)


if __name__ == "__main__":
    # include standard modules
    import argparse
    # initiate the parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--filename", "-f", help="bam_file_name")
    parser.add_argument("--sample", "-s", default=False,
                        help="create random samples")
    parser.add_argument("--run-stats", "-r", default=False,
                        help="run flagstats")

    args = parser.parse_args()
    # read arguments from the command line
    bam_file = args.filename
    if args.sample == True:
        create_sample_files(bam_file)
    if args.run_stats == True:
        get_saturation(bam_file)
