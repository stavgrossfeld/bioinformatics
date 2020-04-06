setwd("/auto/dtg-00/Groups/Hjelm_lab/Grossfel")

files <-list.files(path="./tconut_cna_files", pattern="*.tsv", full.names=TRUE, recursive=FALSE)


#SZ
files_to_process <- c("3993", "4213", "4301", "4385", "4404", "4413", "4463", "4469", "4812", "4114", "4656", "3985", "4142", "4284", "4357", "4508", "4984", "4100", "4564", "4619", "4646", "4661", "4735", "4938", "4429", "4296", "4801")

#grepl(files_to_process)
#grepl(files_to_process[[1]],files)


df_list = list()
file_name_list = list()

for (single_file in files_to_process) {
  
  files_like <-  files[grepl(single_file, files)]
    
  for (file_name in files_like) {
      df <- read.table(file_name, header = T)
      file_name<-gsub(".cna.tsv", x=gsub("./tconut_cna_files",x=file_name,""),"")
      file_name_list <- append(file_name, file_name_list)
      df_list <- append(list(df), df_list)
    }
}

names(df_list) <- file_name_list


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



library("parallel")
no_cores <- 12

cl <- makeCluster(no_cores)

create_cnv_calls <- function(i, df_list) {
  df <- df_list[[i]]
  df$name <- i
  df <- create_amp_del(df)
  df <- find_cnvs(df, min_cnv_length = 1e6)
  return(df)
}
df_new_list <- lapply(X=names(df_list), FUN=create_cnv_calls, df_list)


stopCluster(cl)


#rbind(df_list)

df_final <- do.call("rbind", df_new_list)
found_cnvs <- df_final[df_final$cnv_start_end==2,]

found_cnvs$Position_start <- found_cnvs$Position - found_cnvs$cnv_length

found_bed <- found_cnvs[,c("Chr", "Position_start", "Position","name", "amp_del")]

write.table(found_bed, file = "./found_cnvs/multi_found_cnvs.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

found_seg <- found_cnvs[,c("name","Chr", "Position_start", "Position","cnv_length","Fold.Change")]

write.table(found_seg, file = "./found_cnvs/multi_found_cnvs.seg", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


     