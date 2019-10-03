
import re
import os
import pandas as pd

file_name = "99_sample.bam"
# os.chdir("/Users/stav/Desktop/bam_split")

create_urz_cmd = "samtools view %s | cut -f19-25 >> urz_%s.txt" % (
    file_name, file_name.replace(".bam", ""))

urz_file_name = "urz_%s.txt" % (file_name.replace(".bam", ""))

# os.system(create_urz_cmd)

df = pd.read_csv(urz_file_name, delimiter="\t")
df.columns = ["barcode"]
urz_df = df[df.barcode.str.contains("UR:Z:")]

unique_barcodes = len(urz_df.barcode.unique())
mapped_reads = float(urz_df.shape[0])

print(unique_barcodes, mapped_reads)
