"""Script to filter index hopped or seen once in bams
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
    cb_umi_read_dict = {}
    cb_umi_line_dict = {}

    index_hop_bam = open("index_hop.bam", "a")
    seen_once_bam = open("seen_once.bam", "a")
    multi_map_bam = open("mutli_map.bam", "a")
    pcr_replicate_bam = open("pcr_replicate.bam", "a")

    for line in tqdm(sys.stdin, total=number_of_lines):

        bam_line = line.split()

        read_name = bam_line[0]
        is_pcr_replicate = int(
            bam_line[1]) / 1024 or int(bam_line[1]) / 1040 == 1  # samtools pcr flag

        cb = re.findall(r"CB:Z:\w*", line)
        umi = re.findall(r"UR:Z:\w*", line)

        if len(cb) == 0:
            continue

        umi = umi[0].strip()
        cb = cb[0].strip()

        cb_umi = "".join([cb, "_", umi])  # , "_", chr_location_mega])

        if cb_umi not in cb_umi_read_dict:
            # if read name not in cb_umi dict, create  a list
            cb_umi_read_dict[cb_umi] = [read_name]
            cb_umi_line_dict[cb_umi] = [line]

        else:
            if len(cb_umi_read_dict[cb_umi]) >= 1:
                # if the length of reads for a cb_umi is > 1

                if is_pcr_replicate:
                    pcr_replicate_bam.write(line)

                else:  # multimap
                    if read_name in cb_umi_read_dict[cb_umi]:
                        multi_map_bam.write(line)

                    else:   # index hop
                        index_hop_bam.write(line)

                        # once the read has gone through the flow append it to the original cb_umi_dict
        cb_umi_read_dict[cb_umi].append(read_name)
        cb_umi_line_dict[cb_umi].append(line)

    print("\n creating seen once bam: ")
    for cb_umi in tqdm(cb_umi_read_dict, total=len(cb_umi_read_dict)):
        if len(cb_umi_read_dict[cb_umi]) == 1:
            seen_once_bam.write(cb_umi_line_dict[cb_umi])
        else:
            print(cb_umi_line_dict[cb_umi])

    seen_once_bam.close()
    index_hop_bam.close()
    multi_map_bam.close()
    pcr_replicate_bam.close()

    print("total cb_umi combos in bam", len(cb_umi_read_dict))


if __name__ == "__main__":
    number_of_lines = sys.argv[1]
    main(number_of_lines)
