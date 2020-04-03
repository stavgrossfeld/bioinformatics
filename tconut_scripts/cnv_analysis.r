#!/usr/bin/env Rscript

#args = commandArgs(trailingOnly=TRUE)
#ctrl_to_analyze <- args[1]

ctrl_to_analyze <- "CTRL_4512"
workdir <- paste0(paste0("/home/dtg-00/Groups/Hjelm_lab/Grossfel/wes_bams_dup/TCONUT/",ctrl_to_analyze,""),"/","")
setwd(workdir)

files <- list.files(path=".", pattern="*.cna.tsv", full.names=TRUE, recursive=TRUE)

df_list = list()
for (file_name in files) {
  df <- read.table(file_name, header = T)


  df_list <- append(list(df), df_list)
 }

names(df_list) <- lapply(strsplit(files, "/"), function(x) x[2])




create_amp_del <- function(df) {
  df["amp_del"] <- ifelse(df$Fold.Change >= .75, 1, 0)
  df["amp_del"] <- ifelse(df$Fold.Change <= -.75, -1, df$amp_del)
  df["amp_del"] <- as.integer(df$amp_del)


  return (df)
}


find_cnvs <- function(df, min_cnv_length) {

  position_keeper <- data.frame()

  df["cnv_length"] <- 0
  df["cnv_start_end"] <- 0
  for (ix in seq_along(df$amp_del[-1])) {

    #print(ix)

    prev_amp_del <- df[ix-1, "amp_del"]
    next_amp_del <- df[ix+1,"amp_del"]
    amp_del <- df[ix, "amp_del"]
    position <- df[ix, "Position"]



    cnv_check <- (((amp_del == 1) && (prev_amp_del == 1)) || (amp_del == -1) && (prev_amp_del == -1)) && df[ix, "Chr"] == df[ix+1, "Chr"]
    if (is.na(cnv_check) == FALSE && (cnv_check == TRUE))

    {
      position_keeper <- rbind(position_keeper, df[ix,])


      if ((nrow(position_keeper) > 1) && (next_amp_del != amp_del)) {
        cnv_length <- max(position_keeper$Position) - min(position_keeper$Position)

        if ((cnv_length > min_cnv_length) && (df[ix,"Chr"] == df[ix+1,"Chr"])) {

          df[ix, "cnv_length"] <- cnv_length
          df[ix-(nrow(position_keeper)),"cnv_start_end"] <- 1
          df[ix,"cnv_start_end"] <- 2


        }

      }

    }


    else  {position_keeper <- data.frame() }
  }
  return(df)
}



library(parallel)

# Calculate the number of cores
no_cores <- 12

# Initiate cluster
cl <- makeCluster(no_cores)



process_parallel <- function(i, df_list) {
  df <- df_list[[i]]

  df['name'] <- i
  df <- create_amp_del(df)
  df <- find_cnvs(df, min_cnv_length = 1e6)
  print(paste0("completed:", i))
  return(df)
}


df_list <- mclapply(X=names(df_list), FUN=process_parallel, df_list)

stopCluster(cl)




#rbind(df_list)

df_final <- do.call("rbind", df_list)
found_cnvs <- df_final[df_final$cnv_start_end==2,]

found_cnvs$Position_start <- found_cnvs$Position - found_cnvs$cnv_length

found_bed <- found_cnvs[,c("Chr", "Position_start", "Position","name", "amp_del")]

filename_bed <- paste("../found_cnvs_",ctrl_to_analyze,"_1e6.bed",sep="")

write.table(found_bed, file = filename_bed, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

found_seg <- found_cnvs[,c("name","Chr", "Position_start", "Position","cnv_length","Fold.Change")]

filename_seg <- paste("../found_cnvs_",ctrl_to_analyze,"_1e6.seg",sep="")

write.table(found_seg, file = filename_seg, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)