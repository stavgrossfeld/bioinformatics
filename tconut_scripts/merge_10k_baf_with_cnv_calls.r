setwd("~/Desktop/hjelm_lab/cnv/")

library(GenomicRanges)
library(dplyr)
baf_df <- read.table("./consolidated/baf/consolidate_cnvs_bp_10000.seg", header = T)
cnv_df <- read.table("./consolidated/cnv/consolidate_cnv_bp_1e6.seg", header = T)

cnv_df <- cnv_df[cnv_df$amp_del == -1,]


baf_gr <- GenomicRanges::makeGRangesFromDataFrame(baf_df,keep.extra.columns=TRUE)
cnv_gr <- GenomicRanges::makeGRangesFromDataFrame(cnv_df, keep.extra.columns=TRUE)


#hits <- findOverlaps(baf_gr, cnv_gr)


baf_overlap <- findOverlaps(query = cnv_gr, subject = baf_gr, type = "within")
#baf_overlap <- findOverlaps(query = baf_gr, subject = cnv_gr, type = "within")

baf_overlap <- data.frame(baf_gr[subjectHits(baf_overlap),]) %>% select(name, seqnames, start, end, width, freq, BAF)
colnames(baf_overlap) <- c("name", "Chr", "start", "end", "width","freq","BAF")

#write.table(baf_overlap, "./consolidated/filtered_baf_cnv_1e5_sz.seg", row.names = F, quote=FALSE, sep ="\t")
write.table(baf_overlap, "./consolidated/filtered_10000_baf_cnv_1e6_bp.seg", row.names = F, quote=FALSE, sep ="\t")

