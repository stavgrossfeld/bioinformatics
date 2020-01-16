
import sys
from tqdm import tqdm


def main(number_of_lines, not_seen_once_reads):
    print("creating seen once bam: ")
    not_seen_once_reads = not_seen_once_reads.split(" ")
    not_seen_once_reads = [int(i) for i in not_seen_once_reads]

    seen_once_bam = open("seen_once.bam", "a")

    for bam_ix, line in enumerate(tqdm(sys.stdin, total=number_of_lines)):
        if bam_ix not in not_seen_once_reads:
            # print(line)
            seen_once_bam.write(line)
    seen_once_bam.close()


if __name__ == "__main__":

    number_of_lines = float(sys.argv[1])
    not_seen_once_reads = sys.argv[2]

    main(number_of_lines, not_seen_once_reads)
