import sys
import os
import re
from tqdm import tqdm
import subprocess

import numpy as np
import mmap


def get_num_lines(filename):
    count_lines_cmd = "samtools view %s | wc -l" % filename
    # print(count_lines_cmd)
    # could be anything here.
    direct_output = subprocess.check_output(count_lines_cmd, shell=True)
    number_of_lines = int(direct_output.decode(sys.stdout.encoding))
    return number_of_lines


def main(filename):

    # cell barcodes dictionary of umi dictionary - deduepd reads
    cb_umi_dict = {}
    # cell barcode read_ct
    cb_umi_ct_dict = {}

    number_of_lines = get_num_lines(filename)

    print("\n counting umis in: %s " % filename)
    for line in tqdm(sys.stdin, total=number_of_lines):

        cb = re.findall(r"CB:Z:\w*", line)
        umi = re.findall(r"UR:Z:\w*", line)

        if len(cb) == 0:
            continue

        umi = umi[0].strip()
        cb = cb[0].strip()

        if cb not in cb_umi_dict:
            cb_umi_dict[cb] = {}
            cb_umi_ct_dict[cb] = 1
        if umi not in cb_umi_dict[cb]:
            cb_umi_dict[cb][umi] = 1
            cb_umi_ct_dict[cb] += 1
        else:
            cb_umi_dict[cb][umi] += 1
            # add umi to total reads
            cb_umi_ct_dict[cb] += 1

    seq_saturation_list = []
    deduped_reads_list = []
    reads_per_spot_list = []
    for cell in cb_umi_dict:
        seq_saturation_for_cell = 1 - \
            (len(cb_umi_dict[cell]) / cb_umi_ct_dict[cell])

        reads_per_spot_list.append(cb_umi_ct_dict[cb])
        deduped_reads_list.append(len(cb_umi_dict[cell]))

        seq_saturation_list.append(seq_saturation_for_cell)

    mean_seq_saturation = round(np.mean(seq_saturation_list)*100, 6)
    mean_deduped_reads = round(np.mean(deduped_reads_list), 6)
    mean_reads_per_spot = round(np.mean(reads_per_spot_list), 6)

    metrics_dict = {"file": filename.replace("_sample.bam"),
                    "deduped_reads": mean_deduped_reads,
                    "mean_seq_saturation": mean_seq_saturation,
                    "mean_reads_per_spot": mean_reads_per_spot}

    print(metrics_dict)
    # sequence_saturation = 1 -
    #     float(umi_dict["deduped_reads"] / float(umi_dict["rd_ct"]))
    # print("\n metrics: ", umi_dict)
    # print("\n sequence_saturation: ", sequence_saturation)

    # umi_dict["sequence_saturation"] = sequence_saturation

    with open('log.txt', 'a') as f:
        f.write(str(metrics_dict)+"\n")

    # return umi_dict


if __name__ == "__main__":

    filename = sys.argv[1]
    main(filename)
