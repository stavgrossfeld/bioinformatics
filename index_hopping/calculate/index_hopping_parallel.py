import multiprocessing
from multiprocessing import Pool
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


def worker(line):

    print(line)
    return(line)


def check_read_type(line):
    bam_line = line.split()
    chr_num = bam_line[2]
    chr_loc = bam_line[3]
    chr_loc = str(int(round(float(chr_loc), -6)))
    read_name = bam_line[0]
    is_pcr_replicate = int(
        bam_line[1]) / 1024 or int(bam_line[1]) / 1040 == 1  # samtools pcr flag
    # chr_location_mega = "".join([chr_num, "_", chr_loc])
    # chr_location = "".join([bam_line[2], "_", bam_line[3]])
    cb = re.findall(r"CB:Z:\w*", line)
    umi = re.findall(r"UR:Z:\w*", line)

    if len(cb) == 0:
        return

    umi = umi[0].strip()
    cb = cb[0].strip()
    cb_umi = "".join([cb, "_", umi])  # , "_", chr_location_mega])
    if cb_umi not in cb_umi_dict:
        # if read name not in cb_umi dict, create  a list
        cb_umi_dict[cb_umi] = [read_name]
    else:
        if len(cb_umi_dict[cb_umi]) == 1:
            # check at least one read in cb_umi, if it is create a key,value pair for it in multi_map_index_hop_dict
            multi_map_index_hop_dict[cb_umi] = {
                "multi_mapped_reads": {}, "index_hop_reads": {}, "pcr_replicates": {}}
        if len(cb_umi_dict[cb_umi]) >= 1:
            # if the length of reads for a cb_umi is > 1
            if is_pcr_replicate:
                if read_name not in multi_map_index_hop_dict[cb_umi]["pcr_replicates"]:
                    multi_map_index_hop_dict[cb_umi]["pcr_replicates"][read_name] = 1
                else:
                    multi_map_index_hop_dict[cb_umi]["pcr_replicates"][read_name] += 1
                multi_map_index_hop_dict[cb_umi]["pcr_replicate_ct"] = len(
                    multi_map_index_hop_dict[cb_umi]["pcr_replicates"])
            else:  # multimap
                if read_name in cb_umi_dict[cb_umi]:
                    if read_name not in multi_map_index_hop_dict[cb_umi]["multi_mapped_reads"]:
                        multi_map_index_hop_dict[cb_umi]["multi_mapped_reads"][read_name] = 1
                    else:
                        multi_map_index_hop_dict[cb_umi]["multi_mapped_reads"][read_name] += 1
                    multi_map_index_hop_dict[cb_umi]["unique_multi_map_ct"] = len(
                        multi_map_index_hop_dict[cb_umi]["multi_mapped_reads"])
                else:   # index hop
                    if read_name not in multi_map_index_hop_dict[cb_umi]["index_hop_reads"]:
                        multi_map_index_hop_dict[cb_umi]["index_hop_reads"][read_name] = 1
                    else:
                        multi_map_index_hop_dict[cb_umi]["index_hop_reads"][read_name] += 1
                    multi_map_index_hop_dict[cb_umi]["unique_index_hop_ct"] = len(
                        multi_map_index_hop_dict[cb_umi]["index_hop_reads"])
        # once the read has gone through the flow append it to the original cb_umi_dict
            cb_umi_dict[cb_umi].append(read_name)

    # if len(cb_umi_dict) % 10000 == 0:
    #    print(len(cb_umi_dict))
    return cb_umi_dict, multi_map_index_hop_dict


def main(number_of_lines):
    """ two dictionaries:
    cb_umi : each key is a cb_umi and a list of its read names
    multi_map_index_hop_dict: each key is a cb_umi and contains its index_hopped_reads and multi_mapped_reads

    index_hop: different readnames for cb_umi
    multi_map: same readname for cb_umi
    """
    number_of_lines = float(number_of_lines)

    # original dict
    global cb_umi_dict, multi_map_index_hop_dict
    manager = multiprocessing.Manager()

    cb_umi_dict = manager.dict()
    multi_map_index_hop_dict = manager.dict()

    pool = multiprocessing.Pool(processes=4)

    # [line for line in enumerate(
    #   pool.map(check_read_type, args=, iterable=sys.stdin))]
    for ix, line in enumerate(tqdm(pool.imap_unordered(func=check_read_type, iterable=sys.stdin), total=number_of_lines)):
        continue

    #    cb_umi_dict, multi_map_index_hop_dict = check_read_type
    #    print(len(cb_umi_dict))

    # cb_umi_dict, multi_map_index_hop_dict = check_read_type(line)

    # for line in tqdm(pool.imap(worker(), sys.stdin)):
    # check_read_type()
    # print(line)

    print("total cb_umi combos in bam:", len(cb_umi_dict))
    # create a dataframe
    jumped_index = pd.DataFrame(multi_map_index_hop_dict).T

    print("number of cb_umi seen more than once \n", jumped_index.shape)
    # create distribution

    create_distriubtion_dataframes(jumped_index, number_of_lines)

    jumped_index.to_csv("./jumped_index.csv", index=True)


def create_distriubtion_dataframes(jumped_index, number_of_lines):
    multi_mapped_dist = pd.DataFrame(
        jumped_index.unique_index_hop_ct.value_counts()).reset_index()
    multi_mapped_dist.columns = [
        "cb_umis_with_number_of_reads", "total_read_ct"]
    multi_mapped_dist["total_reads"] = multi_mapped_dist.cb_umis_with_number_of_reads * \
        multi_mapped_dist.total_read_ct
    multi_mapped_dist["pct_reads"] = multi_mapped_dist.total_reads / \
        number_of_lines * 100

    index_hop_dist = pd.DataFrame(
        jumped_index.unique_multi_map_ct.value_counts()).reset_index()
    index_hop_dist.columns = ["cb_umis_with_number_of_reads", "total_read_ct"]
    index_hop_dist["total_reads"] = index_hop_dist.cb_umis_with_number_of_reads * \
        index_hop_dist.total_read_ct
    index_hop_dist["pct_reads"] = index_hop_dist.total_reads / \
        number_of_lines * 100

    pcr_dist = pd.DataFrame(
        jumped_index.pcr_replicate_ct.value_counts()).reset_index()
    pcr_dist.columns = ["cb_umis_with_number_of_reads", "total_read_ct"]
    pcr_dist["total_reads"] = pcr_dist.cb_umis_with_number_of_reads * \
        pcr_dist.total_read_ct
    pcr_dist["pct_reads"] = pcr_dist.total_reads / number_of_lines * 100

    print("distributions: \n\n")
    print("index_hop_dist: \n", index_hop_dist.head(10))
    print("multi_map_dist: \n", multi_mapped_dist.head(10))
    print("pcr_dist: \n", pcr_dist.head(10))

    index_hop_dist.to_csv("./index_hop_dist.csv", index=True)
    multi_mapped_dist.to_csv("./multi_mapped_dist.csv", index=True)
    pcr_dist.to_csv("./pcr_dist.csv", index=True)


if __name__ == "__main__":
    number_of_lines = sys.argv[1]
    main(number_of_lines)
