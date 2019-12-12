import sqlite3
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


def main_sql():
    db_file = "./bam.db"

    create_table_sql = """ DROP TABLE IF EXISTS tmp_bam_reads; CREATE TABLE IF NOT EXISTS tmp_bam_reads (
                                        read text text,
                                        bam_line text
                                    ); """

    # create a database connection
    conn = sqlite3.connect(db_file)
   # create tables
    if conn is not None:
        # create projects table

        c = conn.cursor()
        c.executescript(create_table_sql)
        conn.commit()
        print("created sql reads table")
        return conn
    else:
        print("Error! cannot create the database connection.")


def insert_reads_sql(conn, data):
    # sql_query = """INSERT INTO tmp_bam_reads VALUES("blahblahblah","yallahyallahyallah")"""

    c = conn.cursor()
    c.executemany("""
    INSERT INTO
        tmp_bam_reads
        (read, bam_line)
    VALUES
        (:read, :bam_line)""", data)

    conn.commit()


def query_read_sql(conn, read_name):

    c = conn.cursor()
    qry = """SELECT bam_line FROM tmp_bam_reads WHERE read = "%s"; """ % read_name
    c.execute(qry)
    bam_line = c.fetchone()[0]

    return bam_line


def main(number_of_lines, conn):
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

        # for ix, line in enumerate(read_bam_sample.readlines()):
        # print(line)
       # print(ix)
        if ix % 100000 == 0:
            print(ix)
            insert_reads_sql(conn, cb_umi_line_dict)
            cb_umi_line_dict = []

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
            cb_umi_line_dict.append((read_name, line))

        else:

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

            # once the read has gone through the flow append it to the original cb_umi_dict

    print("\n creating seen once bam: ")
    # print(cb_umi_line_dict)
    for cb_umi in tqdm(cb_umi_read_dict, total=len(cb_umi_read_dict)):
        if len(cb_umi_read_dict[cb_umi]) == 1:
            #  print(cb_umi_line_dict[cb_umi][line])
            read_name = list(cb_umi_read_dict[cb_umi].items())[0][0]
            returned_bam_line = query_read_sql(conn, read_name)
            seen_once_bam.write(returned_bam_line)

    conn.close()
    seen_once_bam.close()
    index_hop_bam.close()
    multi_map_bam.close()
    pcr_replicate_bam.close()

    print("total cb_umi combos in bam", len(cb_umi_read_dict))


if __name__ == "__main__":
    number_of_lines = sys.argv[1]
    #number_of_lines = 1000
    conn = main_sql()
    main(number_of_lines, conn)
