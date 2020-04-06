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
from pymongo import MongoClient
import time
pd.set_option('display.width', 1000)
pd.options.display.max_colwidth = 100


def main_mongo():

    myclient = MongoClient("mongodb://localhost:27017/")

    mydb = myclient["bam_db"]
    mycol = mydb["bam_file"]
    mycol.drop()
    return mycol
# create tables


def insert_reads_mongo(mycol, data):
    # mongo_query = """INSERT INTO tmp_bam_reads VALUES("blahblahblah","yallahyallahyallah")"""

    x = mycol.insert_many(data)
    # print(x.inserted_ids)


def query_read_mongo(mycol, read_names):

    # bam_line = mycol.find({"read_name": read_name}, {"line": 1})
    bam_lines = mycol.find({"read_name": {"$in": read_names}})
    return bam_lines


def main(number_of_lines, mycol):
    """ two dictionaries:
    cb_umi: each key is a cb_umi and a list of its read names
    multi_map_index_hop_dict: each key is a cb_umi and contains its index_hopped_reads and multi_mapped_reads

    index_hop: different readnames for cb_umi
    multi_map: same readname for cb_umi
    """
    number_of_lines = float(number_of_lines)

    # original dict
    cb_umi_read_dict = {}
    cb_umi_line_dict = []

    # os.chdir("/Users/stav/Desktop/bam_split/")

    index_hop_bam = open("index_hop.bam", "a")
    seen_once_bam = open("seen_once.bam", "a")
    multi_map_bam = open("multi_map.bam", "a")
    pcr_replicate_bam = open("pcr_replicate.bam", "a")

    for ix, line in enumerate(tqdm(sys.stdin, total=number_of_lines)):
        ix = ix + 1
        if ix % 10000 == 0 and len(cb_umi_line_dict) > 0:

            insert_reads_mongo(mycol, cb_umi_line_dict)

            cb_umi_line_dict = []
            continue

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
            cb_umi_read_dict[cb_umi] = {read_name: 1}
            cb_umi_line_dict.append({"read_name": read_name, "line": line})

        else:

            # this adds another .3 % to to the total
            #cb_umi_line_dict.append({"read_name": read_name, "line": line})

            if len(cb_umi_read_dict[cb_umi]) >= 1:
                # if the length of reads for a cb_umi is > 1

                if is_pcr_replicate:
                    pcr_replicate_bam.write(line)
                    pass
                else:
                    # multimap
                    if read_name in cb_umi_read_dict[cb_umi]:
                        cb_umi_read_dict[cb_umi][read_name] += 1
                        multi_map_bam.write(line)
                        pass
                    else:
                        # index hop
                        cb_umi_read_dict[cb_umi][read_name] = 1
                        index_hop_bam.write(line)
                        pass
            else:
                print(ix, bam_line)

            # once the read has gone through the flow append it to the original cb_umi_dict

    # insert last one
    try:
        insert_reads_mongo(mycol, cb_umi_line_dict)
    except:
        pass
    print("index mongodb collection on read name:")

    start_time = time.time()
    mycol.create_index("read_name")

    print("--- %s seconds to index ---" % (time.time() - start_time))

    print("\n creating seen once bam query mongo: ")
    # print(cb_umi_line_dict)
    reads_to_query = []
    for ix, cb_umi in enumerate(tqdm(cb_umi_read_dict, total=len(cb_umi_read_dict))):
        if len(cb_umi_read_dict[cb_umi]) <= 1:
            #  print(cb_umi_line_dict[cb_umi][line])
            read_name = list(cb_umi_read_dict[cb_umi].items())[0][0]
            reads_to_query.append(read_name)

        if ix % 1000 == 0:
            # read mongo and write to bam
            reads_requested = query_read_mongo(mycol, reads_to_query)

            for read in reads_requested:
                seen_once_bam.write(read["line"])
            reads_to_query = []

    seen_once_bam.close()
    index_hop_bam.close()
    multi_map_bam.close()
    pcr_replicate_bam.close()

    print("total cb_umi combos in bam", len(cb_umi_read_dict))


if __name__ == "__main__":
    number_of_lines = sys.argv[1]
    mycol = main_mongo()
    main(number_of_lines, mycol)
