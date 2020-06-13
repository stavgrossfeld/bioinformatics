setwd("~/Desktop/hjelm_lab/cnv/")
#https://bioconductor.org/packages/devel/bioc/vignettes/CNVRanger/inst/doc/CNVRanger.html
#install.packages('DescTools')
#require(DescTools)
#install.packages("CNVRanger")
require(BiocManager)
#BiocManager::install("CNVRanger")
library(dplyr)
library(ggplot2)
library(data.table)
library(GenomicRanges)
library(rlist)

parse_name <- function(name) {
  
  stub_list <- unlist(strsplit(name[1], split="_"))
  sample_name <- paste0(stub_list[1], "_",stub_list[2])
  return(sample_name)
}



pre_process <- function(df_all) {
  df_all$sample_name <- as.character(lapply(as.character(df_all$name), FUN=function(i) parse_name(i)))
  # df_all <- df_all[df_all$Position_start > 0,]
  df_all$BAF[is.na(df_all$BAF)] <- 0
  df_all["loss_het"] <- ifelse(df_all$BAF >= 0 & df_all$BAF <= .1, -1, 0)
  df_all["loss_het"] <- ifelse(df_all$BAF >= .9 & df_all$BAF <= 1, -1, df_all$loss_het)
  df_all["loss_het"] <- as.integer(df_all$loss_het)
  return(df_all)
}

find_overlapping_cnv <- function(df_all) {
  cnv_list_df <- list()
  for (i in unique(df_all$sample_name)) {
    df <- df_all[df_all$sample_name == i,]
    df$name <- as.factor(as.character(df$name))
    df <- df[,c("name", "Chr", "Position_start", "Position", "loss_het")]
    colnames(df) <- c("name","Chr", "start","end","state")
    
    # df <- df[!df$Chr %in% c(23,24),]
    
    df_amp <- df[df$state == 1,]
    df_del <- df[df$state == -1,]
    
    
    grl_amp_overlaps <- c()
    grl_del_overlaps <- c()
    if (nrow(df_amp) > 0){
      grl_amp <- GenomicRanges::makeGRangesListFromDataFrame(df_amp,
                                                             split.field="name", keep.extra.columns=TRUE)
      grl_amp_overlaps <- CNVRanger::populationRanges(grl_amp, mode="RO", ro.thresh=.9, verbose)
      
    }
    
    if (nrow(df_del) > 0) {
      grl_del <- GenomicRanges::makeGRangesListFromDataFrame(df_del,
                                                             split.field="name", keep.extra.columns=TRUE)
      grl_del_overlaps <- CNVRanger::populationRanges(grl_del, mode="RO", ro.thresh=.9, verbose)
      
    }
    
    cnv_list_df[[i]] <- list("amp"=grl_amp_overlaps, "del"=grl_del_overlaps)
    
  }
  return(cnv_list_df)
}


post_process_cnv <- function(cnv_list_df) {
  new_cnv_list <- list.flatten(cnv_list_df, use.names = TRUE, classes = "ANY")
  
  new_cnv_list_df <- lapply(new_cnv_list, data.frame)
  new_df <- do.call("rbind", new_cnv_list_df)
  new_df$BAF <- ifelse(grepl("BAF",rownames(new_df)), 1, -1)
  
  new_df$sample_name <- as.character(rownames(new_df))
  new_df$sample_name <- as.character(lapply((new_df$sample_name), FUN=function(i) unlist(strsplit(i,"[.]"))[1]))
  
  rownames(new_df) <- 1:nrow(new_df)
  new_df <- new_df %>% select(sample_name, seqnames, start, end, width, freq, BAF)
  new_df <- new_df[new_df$freq >= 5,]
  colnames(new_df) <- c("name", "Chr", "start", "end", "width","freq","BAF")
  
  return(new_df)
  
  
}



run_overlap_ct <- function(read_path, write_path) {
  df_all <- read.table(read_path, header = T)
  df_all <- pre_process(df_all)
  cnv_list_df <- find_overlapping_cnv(df_all)
  new_df <- post_process_cnv(cnv_list_df)
  write.table(new_df, write_path, col.names = TRUE, row.names = F, quote = FALSE, sep = "\t")
  
}

#1e6
run_overlap_ct(read_path = "./found_bafs/SZ_found_cnvs_1e+06.seg", write_path ="./consolidated/baf/consolidate_cnvs_sz_1e6.seg")
run_overlap_ct(read_path = "./found_bafs/BP_found_cnvs_1e+06.seg", write_path ="./consolidated/baf/consolidate_cnvs_bp_1e6.seg")
run_overlap_ct(read_path = "./found_bafs/CTRL_found_cnvs_1e+06.seg", write_path ="./consolidated/baf/consolidate_cnvs_ctrl_1e6.seg")


#1e5

run_overlap_ct(read_path = "./found_bafs/SZ_found_cnvs_1e+05.seg", write_path ="./consolidated/baf/consolidate_cnvs_sz_1e5.seg")
run_overlap_ct(read_path = "./found_bafs/BP_found_cnvs_1e+05.seg", write_path ="./consolidated/baf/consolidate_cnvs_bp_1e5.seg")
run_overlap_ct(read_path = "./found_bafs/CTRL_found_cnvs_1e+05.seg", write_path ="./consolidated/baf/consolidate_cnvs_ctrl_1e5.seg")

# 5e5

run_overlap_ct(read_path = "./found_bafs/SZ_found_cnvs_5e+05.seg", write_path ="./consolidated/baf/consolidate_cnvs_sz_5e5.seg")
run_overlap_ct(read_path = "./found_bafs/BP_found_cnvs_5e+05.seg", write_path ="./consolidated/baf/consolidate_cnvs_bp_5e5.seg")
run_overlap_ct(read_path = "./found_bafs/CTRL_found_cnvs_5e+05.seg", write_path ="./consolidated/baf/consolidate_cnvs_ctrl_5e5.seg")


# 10000, 1000
run_overlap_ct(read_path = "./found_bafs/CTRL_found_cnvs_10000.seg", write_path ="./consolidated/baf/consolidate_cnvs_ctrl_10000.seg")
run_overlap_ct(read_path = "./found_bafs/CTRL_found_cnvs_1000.seg", write_path ="./consolidated/baf/consolidate_cnvs_ctrl_1000.seg")


run_overlap_ct(read_path = "./found_bafs/SZ_found_cnvs_10000.seg", write_path ="./consolidated/baf/consolidate_cnvs_sz_10000.seg")
run_overlap_ct(read_path = "./found_bafs/SZ_found_cnvs_1000.seg", write_path ="./consolidated/baf/consolidate_cnvs_sz_1000.seg")

run_overlap_ct(read_path = "./found_bafs/BP_found_cnvs_10000.seg", write_path ="./consolidated/baf/consolidate_cnvs_bp_10000.seg")
run_overlap_ct(read_path = "./found_bafs/BP_found_cnvs_1000.seg", write_path ="./consolidated/baf/consolidate_cnvs_bp_1000.seg")

