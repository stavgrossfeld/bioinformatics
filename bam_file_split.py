
import re
import os
import pandas as pd
import pysam
import numpy as np
umi_cts = {}


def create_barcodes_cut(file_name):

    create_urz_cmd = "samtools view %s | grep UR:Z: | cut -f19-27 >> urz_%s.txt" % (
        file_name, file_name.replace(".bam", ""))
    os.system(create_urz_cmd)
    urz_file_name = "urz_%s.txt" % (file_name.replace(".bam", ""))


def ct_umis():

    with open('/Users/stav/Desktop/bam_split/urz_01_sample.txt', 'r') as f:
        for line in f:
            read = line.replace("\n", "")
            if "UR:Z:" in read:

                if read not in umi_cts:
                    umi_cts[read] = 1
                else:
                    umi_cts[read] += 1

    umi_df = pd.DataFrame(list(umi_cts.items()), columns=['barcode', 'ct'])
    umi_df.to_csv("umi_df_sample.csv", index=false)

    sequence_saturation = 1 - (umi_df.shape[0] / np.sum(umi_df.ct))

    print(sequence_saturation)


if __name__ == "__main__":
    import argparse
    # initiate the parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--filename", "-f", help="bam_file_name")
    parser.add_argument("--sample", "-s", default=False,
                        help="create random samples")
    parser.add_argument("--run-stats", "-r", default=False,
                        help="run flagstats")

    create_barcodes_cut(file_name)
    ct_umis()
    # unique_barcodes = len(urz_df.barcode.unique())
    # mapped_reads = float(urz_df.shape[0])
    # print(unique_barcodes, mapped_reads)

    # file_name = '/Users/stav/Desktop/bam_split/01_sample.bam'

    # samfile = pysam.AlignmentFile(file_name, "rb")

    # # os.chdir("/Users/stav/Desktop/bam_split")

    # pure_bam = pybam.bgunzip('/Users/stav/Desktop/bam_split/01_sample.bam')
    # parser = pybam.compile_parser(
    #     ['pos', 'flag', 'rname', 'mapq', 'rnext', 'pnext', 'tlen'])
    # for read in parser(pure_bam):
    #     read[0]
    #     read[1]
    #     read[2]
    #     read[3]
    #     read[4]
    #     read[5]
    #     read[6]
