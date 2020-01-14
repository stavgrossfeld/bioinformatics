"""script to filter bam files
ex run: python ~/bioinformatics_scripts/loop_seq_saturation.py ./ True
"""

import glob
import os
import sys
from multiprocessing import Pool


def main(bam_file):

    # rint(multi)
    cwd = os.getcwd()

    read_types_dir = "%s/%s_read_types" % (cwd, bam_file)

    if not os.path.isdir(read_types_dir):
        os.mkdir(read_types_dir)
    os.chdir(read_types_dir)

    file_types = ["index_hop.bam", "seen_once.bam",
                  "multi_map.bam", "pcr_replicate.bam"]

    # p = Pool(4)
    # p.map(create_files, file_type, bam_file)

    for file_name in file_types:
        create_files(file_name, bam_file)

    call_cmd(bam_file)


def create_files(file_type, bam_file):
    print(bam_file)
    file_type = file_type

    cmd = "samtools view -H ../%s > ./%s" % (bam_file, file_type)
   # print(cmd)
    os.system(cmd)


def call_cmd(bam_file):

    number_of_lines = os.system("samtools view -c ../%s" % bam_file)

    cmd = "samtools view ../%s | head -10000 | python3 ~/bioinformatics_scripts/index_hopping/filter_bams/filter_bam_mongo.py %s " % (
        bam_file, number_of_lines)

    os.system(cmd)


if __name__ == "__main__":

    bam_file = sys.argv[1]
    main(bam_file)
