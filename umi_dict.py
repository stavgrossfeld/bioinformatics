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
    for line in tqdm(sys.stdin, total=number_of_lines):

        cb = re.findall(r"CB:Z:\w*", line)
        umi = re.findall(r"UR:Z:\w*", line)

        if len(cb) == 0:
            continue

        umi = umi[0].strip()
        cb = cb[0].strip()

        if cb not in CELL_BARCODE_UMI_CT_DICT:
            CELL_BARCODE_UMI_CT_DICT[cb] = {}
            CELL_BARCODE_RD_CT_DICT[cb] = 1
        if umi not in CELL_BARCODE_UMI_CT_DICT[cb]:
            CELL_BARCODE_UMI_CT_DICT[cb][umi] = 1
            CELL_BARCODE_RD_CT_DICT[cb] += 1
        else:
            CELL_BARCODE_UMI_CT_DICT[cb][umi] += 1
            # add umi to total reads
            CELL_BARCODE_RD_CT_DICT[cb] += 1

    for cell in CELL_BARCODE_UMI_CT_DICT:
        seq_saturation_for_cell = 1 - \
            (float(len(CELL_BARCODE_UMI_CT_DICT[cell])
                   ) / float(CELL_BARCODE_RD_CT_DICT[cell]))

        TOTAL_READS_LIST.append(CELL_BARCODE_RD_CT_DICT[cb])
        DEDUPED_READS_LIST.append(len(CELL_BARCODE_UMI_CT_DICT[cell]))

        SEQ_SATURATION_LIST.append(seq_saturation_for_cell)

    mean_seq_saturation = round(np.mean(SEQ_SATURATION_LIST)*100, 6)
    mean_deduped_reads = round(np.mean(DEDUPED_READS_LIST), 6)
    mean_reads_per_spot = round(np.mean(TOTAL_READS_PER_CELL), 6)

    metrics_dict = {"file": filename.replace("_sample.bam", ""),
                    "mean_deduped_reads": mean_deduped_reads,
                    "mean_seq_saturation": mean_seq_saturation,
                    "mean_reads_per_spot": mean_reads_per_spot}

    print(metrics_dict)

    with open('log.txt', 'a') as f:
        f.write(str(metrics_dict)+"\n")

    # return umi_dict


if __name__ == "__main__":

    filename = sys.argv[1]
    main(filename)
