"""Script to create sequence saturation log for a file
according to 10x: sequence saturaiton is the mean of the following equation:

1-(deduped_reads/total_reads)*100
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


def get_num_lines(filename):
    count_lines_cmd = "samtools view %s | wc -l" % filename

    direct_output = subprocess.check_output(count_lines_cmd, shell=True)
    number_of_lines = int(direct_output.decode(sys.stdout.encoding))
    return number_of_lines


def main(filename):

    number_of_lines = get_num_lines(filename)

    print("\n counting umis in: %s " % filename)

    # parse_bam_cmd = "samtools view %s" % filename
    # bam_output = subprocess.run(parse_bam_cmd)

    umi_dict = {}
    umi_dict_ct = {}
    for line in tqdm(sys.stdin, total=number_of_lines):
        bam_line = line.split()

        chr_num = bam_line[2]
        chr_loc = bam_line[3]

        #chr_location = "".join([bam_line[2], "_", bam_line[3]])

        cb = re.findall(r"CB:Z:\w*", line)
        umi = re.findall(r"UR:Z:\w*", line)

        if len(cb) == 0:
            continue

        umi = umi[0].strip()
        cb = cb[0].strip()

        #umi_chr = "".join([umi, ":", chr_location])
        #cb_chr = "".join([cb, ":", chr_location])
        #attributes = [umi, cb, chr_location]

        if umi not in umi_dict:
            # umi_dict[umi] = [cb_chr]

            umi_dict[umi] = {"chr_num": [chr_num], "chr_ct": 1}

            #umi_dict[umi] = {"cb_code": [cb_chr], "ct": 1}
            # umi_dict_ct[umi] = 1
        else:
            if chr_num not in umi_dict[umi]["chr_num"]:
                umi_dict[umi]["chr_num"].append(chr_num)
                umi_dict[umi]["chr_ct"] += 1
            else:
                pass
    # print(umi_dict)
    # umi_dict_ct[umi] += 1

    umi_df = pd.DataFrame(umi_dict).T
    # umi_df.reset_index()
    # umi_df.columns = ["umi", "cb_chr", "ct"]
    #umi_df["cb_code"] = umi_df.cb_code.astype(str)

    # print(umi_df)
    jumped_index = umi_df[umi_df.chr_ct > 1]
    print(jumped_index)


if __name__ == "__main__":
    # filename = sys.argv[1]
    filename = "~/Desktop/bam_split/01_sample.bam"
    main(filename)
