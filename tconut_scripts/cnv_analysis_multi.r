setwd("/auto/dtg-00/Groups/Hjelm_lab/Grossfel")



#SZ
sz_files_to_process <- c("3993", "4213", "4301", "4385", "4404", "4413", "4463", "4469", "4812", "4114", "4656", "3985", "4142", "4284", "4357", "4508", "4984", "4100", "4564", "4619", "4646", "4661", "4735", "4938", "4429", "4296", "4801")

# BP
bp_files_to_process <- c("3927","4330","4345","4794","3772","3967","4087","4242","4262","4359","4366","4405","4639","4306","4383","4419","4857","3618","3711","4400","4707","4714","4741","4848","4944","4978")

# CTRL
ctrl_files_to_process <- c("3965","3969","4048","4063","4236","4286","4350","4520","3933","3686","3774","3776","3952","4111","3947","4015","4191","4207","4512","4638","4652","4706","4839","4844","4905","3878","3896","4069","4088","4250","4327","4744","4754","4796","3850","4235","4387","4537","4635","4623","4464","4314","4302")


# controls to process
ref_ctrls <- c("4069F","4191M","4314F","4387M","4512M","4744M","4754F")
#files_to_process <- append(files_to_process,ref_ctrls)
faulty_ctrl <- "3774F"
# 3774 reference sample to remove



diagnosis_list <- list("SZ"=sz_files_to_process, "CTRL"=ctrl_files_to_process)

df_list = list()
file_name_list = list()


files <-list.files(path="./tconut_cna_files", pattern="*.tsv", full.names=TRUE, recursive=FALSE)

# files vs ctrls
  
for (diagnosis_tag in names(diagnosis_list)) {
  files_to_process <- diagnosis_list[[diagnosis_tag]]
  for (single_file in files_to_process) {
    files_like <-  files[grepl(single_file, files)]
    for (file_name in files_like) {
      new_file_name<-gsub(".cna.tsv", x=gsub("./tconut_cna_files/",x=file_name,""),"")
      new_file_name <- paste0(diagnosis_tag, "_", new_file_name)
      print(new_file_name)
      if (grepl(faulty_ctrl, file_name) == FALSE) {
          df <- read.table(file_name, header = T)
          df_list <- append(list(df), df_list)
          file_name_list <- append(new_file_name, file_name_list)
      }
    }
  }
}


# ctrls vs ctrls

for (single_ctrl in ref_ctrls) {
  files_like <- files[grepl(single_ctrl, files)]
  for (file_name in files_like) {
    new_file_name<-gsub(".cna.tsv", x=gsub("./tconut_cna_files/",x=file_name,""),"")
    split_file <- strsplit(new_file_name,"_")[[1]]
    if (any(grepl(split_file[1], ref_ctrls)) && any((grepl(split_file[2], ref_ctrls))))
     {
        stub_1 <- ref_ctrls[grepl(split_file[1],ref_ctrls)] 
        stub_2 <- ref_ctrls[grepl(split_file[2],ref_ctrls)] 
        new_file_name <- paste0("ref_vs_ref_",stub_1, "_", stub_2)
        df <- read.table(file_name, header = T)
        df_list <- append(list(df), df_list)
        file_name_list <- append(new_file_name, file_name_list)
     }

  }
}
print("finished reading all files")


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
  # use amp_del[-1] because last row
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

create_cnv_calls <- function(i, df_list, min_cnv_length) {
  df <- df_list[[i]]
  df$name <- i
  df <- create_amp_del(df)
  df <- find_cnvs(df, min_cnv_length)
  found_cnvs <- df[df$cnv_start_end==2,]

  found_cnvs$Position_start <- found_cnvs$Position - found_cnvs$cnv_length
  if (nrow(found_cnvs) == 0) {
    found_cnvs <- df[1,]
    found_cnvs$Chr <- 1
    found_cnvs$Positition=-100
    found_cnvs$Position_start=-100
    found_cnvs$cnv_legnth=-100

    return(found_cnvs)
  }
  #print(found_cnvs)
}

run_cnv_and_write_file <- function(df_list, min_cnv_length, tag) {
  df_found_cnv_list <- lapply(X=names(df_list), FUN=create_cnv_calls, df_list, min_cnv_length)
  #return(df_found_cnv_list)
  found_cnvs <- do.call("rbind", df_found_cnv_list)


  found_bed <- found_cnvs[,c("Chr", "Position_start", "Position","name", "amp_del")]
  write.table(found_bed, file = paste0("./found_cnvs/",tag,"_found_cnvs_",min_cnv_length,".bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  found_seg <- found_cnvs[,c("name","Chr", "Position_start", "Position","cnv_length","Fold.Change")]
  write.table(found_seg, file = paste0("./found_cnvs/",tag,"_found_cnvs_",min_cnv_length,".seg"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

}
run_cnv_and_write_file(df_list, 1e5, "sz") #100k
print("finished 100k")
run_cnv_and_write_file(df_list, 5e5, "sz") #500k
print("finished 500k")
run_cnv_and_write_file(df_list, 1e6, "sz") # 1mil
print("finished 1M")






stopCluster(cl)


#rbind(df_list)




     