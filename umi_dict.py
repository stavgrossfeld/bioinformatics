import sys
import os
import re
from tqdm import tqdm
import subprocess

import mmap


def get_num_lines(filename):
    count_lines_cmd = "samtools view %s | wc -l" % filename
    # print(count_lines_cmd)
    # could be anything here.
    direct_output = subprocess.check_output(count_lines_cmd, shell=True)
    number_of_lines = int(direct_output.decode(sys.stdout.encoding))
    return number_of_lines


def main(filename):

    umi_ct_dict = {}
    rd_ct = 0

    number_of_lines = get_num_lines(filename)

    # to maybe read in python later

    #cmd_samtools_view = "samtools view ~/Desktop/bam_split/01_sample.bam | less"
    #samtools_output = subprocess.Popen(cmd_samtools_view, shell=False)
    # for line in tqdm(samtools_output, number_of_lines):
    #    print(line)

    print("\n counting umis in: %s " % filename)
    for line in tqdm(sys.stdin, total=number_of_lines):
        umi = re.findall(r"UR:Z:\w*", line)

        if len(umi) > 1:
            print(umi, len(umi))
            break
        if umi[0] not in umi_ct_dict:
            umi_ct_dict[umi[0]] = 1
            rd_ct += 1
        else:
            umi_ct_dict[umi[0]] += 1
            rd_ct += 1

    deduped_reads = len(umi_ct_dict)

    umi_dict = {"file": filename,
                "deduped_reads": deduped_reads,
                "rd_ct": rd_ct}

    sequence_saturation = 1 - \
        float(umi_dict["deduped_reads"] / float(umi_dict["rd_ct"]))
    print("\n metrics: ", umi_dict)
    print("\n sequence_saturation: ", sequence_saturation)

    umi_dict["sequence_saturation"] = sequence_saturation

    with open('log.txt', 'a') as f:
        f.write(str(umi_dict)+"\n")

    return umi_dict


if __name__ == "__main__":

    filename = sys.argv[1]
    main(filename)
