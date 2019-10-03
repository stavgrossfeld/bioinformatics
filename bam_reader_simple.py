import sys
import re
from tqdm import tqdm


def main():

    umi_ct_dict = {}
    rd_ct = 0
    for line in tqdm(sys.stdin):
        # print(line)
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

        # if rd_ct % 100000:
        #    print("read_ct:  ", rd_ct)

    # print(rd_ct)
    # print(len(umi_cts))

    deduped_reads = len(umi_ct_dict)

    umi_dict = {"deduped_reads": deduped_reads, "rd_ct": rd_ct}

    sequence_saturation = 1 - (deduped_reads / rd_ct)
    print("\n metrics: ", umi_dict)
    print("\n sequence_saturation: ", sequence_saturation)


if __name__ == "__main__":
    main()
