require(BiocManager)
#BiocManager::install("CNVRanger")
library(dplyr)
library(ggplot2)
library(data.table)
library(GenomicRanges)
library(rlist)


setwd("~/Desktop/hjelm_lab/cnv/")
# investigate

df <- read.table("./consolidated/consolidate_cnv_sz_5e5.seg", sep = "\t", header = T)
df_ctrl <- read.table("./consolidated/consolidate_cnv_ctrl_5e5.seg", sep = "\t", header = T)

df_all <- rbind(df, df_ctrl)
df_all <- df_all %>% filter(start!=-100)

colnames(df_all) <- c("name","Chr", "start","end","width", "freq","state")


df_all$state <- df_all$state + 2

df_all <- df_all %>% select(-"freq")
grl_amp <- GenomicRanges::makeGRangesListFromDataFrame(df_all,
                                                       split.field="name", keep.extra.columns=TRUE)
##grl_amp_overlaps <- CNVRanger::populationRanges(grl_amp, mode="RO", ro.thresh=.9, verbose, est.recur=TRUE)
cnvrs <- CNVRanger::populationRanges(grl_amp, density=0.1, ro.thresh=.9, est.recur=TRUE)

CNVRanger::plotRecurrentRegions(cnvrs, genome="hg19", chr="Chr9")#, from =0, to = 1e6)

cnvrs[cnvrs$pvalue < .05,]

plotCoverage(cnvrs)
library(Gviz)

Gviz::plotTracks(cnvrs[cnvrs$pvalue < .05,])
