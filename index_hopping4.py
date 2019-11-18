"""Script to create filter out index hopped bams
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
# cell barcodes dictionary of umi dictionary - deduepd reads
CELL_BARCODE_UMI_CT_DICT = {}
# cell barcode read_ct
CELL_BARCODE_RD_CT_DICT = {}

SEQ_SATURATION_LIST = []
DEDUPED_READS_LIST = []
TOTAL_READS_LIST = []


def main(number_of_lines):

    # number_of_lines = get_num_lines(filename)

    # print("\n counting umis in: %s " % filename)

    # parse_bam_cmd = "samtools view %s" % filename
    # bam_output = subprocess.run(parse_bam_cmd)
    number_of_lines = float(number_of_lines)
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

            cb_umi_dict[cb_umi] = [read_name]

        else:
            if len(cb_umi_dict[cb_umi]) == 1:
                multi_map_index_hop_dict[cb_umi] = {"multi_mapped_reads": [], "index_hop_reads": [],
                                                    }
            if len(cb_umi_dict[cb_umi]) >= 1:
                if read_name in cb_umi_dict[cb_umi]:
                     # multimap
                    multi_map_index_hop_dict[cb_umi]["multi_mapped_reads"].append(
                        read_name)
                    multi_map_index_hop_dict[cb_umi]["multi_map_ct"] = len(
                        multi_map_index_hop_dict[cb_umi]["multi_mapped_reads"])
                else:
                    multi_map_index_hop_dict[cb_umi]["index_hop_reads"].append(
                        read_name)
                    multi_map_index_hop_dict[cb_umi]["index_hop_ct"] = len(
                        multi_map_index_hop_dict[cb_umi]["index_hop_reads"])

            cb_umi_dict[cb_umi].append(read_name)

    # print(multi_map_index_hop_dict)
    # exit()
    # if chr_num not in cb_umi_dict[cb_umi]["chr_num"]:
    #     cb_umi_dict[cb_umi]["chr_num"].append(chr_num)
    #     cb_umi_dict[cb_umi]["chr_ct"] = len(
    #         cb_umi_dict[cb_umi]["chr_num"])

    jumped_index = pd.DataFrame(multi_map_index_hop_dict).T

    print(jumped_index.shape)

    multi_mapped_dist = pd.DataFrame(
        jumped_index.multi_map_ct.value_counts() / number_of_lines * 100)

    index_hop_dist = pd.DataFrame(
        jumped_index.index_hop_ct.value_counts() / number_of_lines * 100)

    print("index_hop_dist: \n\n", index_hop_dist)
    print("multi_map_dist: \n\n", multi_mapped_dist)

    index_hop_dist.to_csv("./index_hop_dist.csv", index=True)
    multi_mapped_dist.to_csv("./multi_mapped_dist.csv", index=True)

    jumped_index.to_csv("./jumped_index.csv", index=True)


if __name__ == "__main__":
    # filename = sys.argv[1]
    number_of_lines = sys.argv[1]
    # filename = "~/Desktop/bam_split/01_sample.bam"
    main(number_of_lines)
