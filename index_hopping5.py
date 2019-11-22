"""Script to create count index hopped / multi mapped in bams
cb_umi but different read name = index hopping
"""

import sys
import os
import re
from tqdm import tqdm
import subprocess

import numpy as np
import mmap
import pandas as pd
pd.set_option('display.width', 1000)
pd.options.display.max_colwidth = 100


def main(number_of_lines):
    """ two dictionaries:
    cb_umi : each key is a cb_umi and a list of its read names
    multi_map_index_hop_dict: each key is a cb_umi and contains its index_hopped_reads and multi_mapped_reads

    index_hop: different readnames for cb_umi
    multi_map: same readname for cb_umi
    """
    number_of_lines = float(number_of_lines)

    # original dict
    cb_umi_dict = {}

    multi_map_index_hop_dict = {}

    for line in tqdm(sys.stdin, total=number_of_lines):
        bam_line = line.split()

        chr_num = bam_line[2]
        chr_loc = bam_line[3]

        chr_loc = str(int(round(float(chr_loc), -6)))
        read_name = bam_line[0]

        chr_location_mega = "".join([chr_num, "_", chr_loc])

        # chr_location = "".join([bam_line[2], "_", bam_line[3]])

        cb = re.findall(r"CB:Z:\w*", line)
        umi = re.findall(r"UR:Z:\w*", line)

        if len(cb) == 0:
            continue

        umi = umi[0].strip()
        cb = cb[0].strip()

        cb_umi = "".join([cb, "_", umi])  # , "_", chr_location_mega])

        # umi_chr = "".join([umi, ":", chr_location])
        # cb_chr = "".join([cb, ":", chr_location])
        # attributes = [umi, cb, chr_location]

        if cb_umi not in cb_umi_dict:
            # if read name not in cb_umi dict, create  a list
            cb_umi_dict[cb_umi] = [read_name]

        else:

            if len(cb_umi_dict[cb_umi]) == 1:
                # check at least one read in cb_umi, if it is create a key,value pair for it in multi_map_index_hop_dict
                multi_map_index_hop_dict[cb_umi] = {"multi_mapped_reads": {}, "index_hop_reads": {},
                                                    }
            if len(cb_umi_dict[cb_umi]) >= 1:
                # if the length of reads for a cb_umi is > 1
                if read_name in cb_umi_dict[cb_umi]:
                     # multimap
                    if read_name not in multi_map_index_hop_dict[cb_umi]["multi_mapped_reads"]:
                        multi_map_index_hop_dict[cb_umi]["multi_mapped_reads"][read_name] = 1
                    else:
                        multi_map_index_hop_dict[cb_umi]["multi_mapped_reads"][read_name] += 1

                    multi_map_index_hop_dict[cb_umi]["unique_multi_map_ct"] = len(
                        multi_map_index_hop_dict[cb_umi]["multi_mapped_reads"])
                else:
                    # index hop
                    if read_name not in multi_map_index_hop_dict[cb_umi]["index_hop_reads"]:
                        multi_map_index_hop_dict[cb_umi]["index_hop_reads"][read_name] = 1
                    else:
                        multi_map_index_hop_dict[cb_umi]["index_hop_reads"][read_name] += 1

                    multi_map_index_hop_dict[cb_umi]["unique_index_hop_ct"] = len(
                        multi_map_index_hop_dict[cb_umi]["index_hop_reads"])

            # once the read has gone through the flow append it to the original cb_umi_dict
            cb_umi_dict[cb_umi].append(read_name)

    print("total cb_umi combos in bam", len(cb_umi_dict))

    # create a dataframe
    jumped_index = pd.DataFrame(multi_map_index_hop_dict).T

    print("number of cb_umi seen more than once \n", jumped_index.shape)
    # create distribution

    create_distriubtion_dataframes(jumped_index, number_of_lines)

    jumped_index.to_csv("./jumped_index.csv", index=True)


def create_distriubtion_dataframes(jumped_index, number_of_lines):
    multi_mapped_dist = pd.DataFrame(
        jumped_index.unique_index_hop_ct.value_counts()).reset_index()
    multi_mapped_dist.columns = ["read_number", "ct"]
    multi_mapped_dist["total_reads"] = multi_mapped_dist.read_number * \
        multi_mapped_dist.ct
    multi_mapped_dist["pct_reads"] = multi_mapped_dist.total_reads / \
        number_of_lines * 100

    index_hop_dist = pd.DataFrame(
        jumped_index.unique_multi_map_ct.value_counts()).reset_index()
    index_hop_dist.columns = ["read_number", "ct"]
    index_hop_dist["total_reads"] = index_hop_dist.read_number * \
        index_hop_dist.ct
    index_hop_dist["pct_reads"] = index_hop_dist.total_reads / \
        number_of_lines * 100

    print("distributions: \n\n")
    print("index_hop_dist: \n", index_hop_dist.head(10))
    print("multi_map_dist: \n", multi_mapped_dist.head(10))

    index_hop_dist.to_csv("./index_hop_dist.csv", index=True)
    multi_mapped_dist.to_csv("./multi_mapped_dist.csv", index=True)


if __name__ == "__main__":
    # filename = sys.argv[1]
    number_of_lines = sys.argv[1]
    # filename = "~/Desktop/bam_split/01_sample.bam"
    main(number_of_lines)
