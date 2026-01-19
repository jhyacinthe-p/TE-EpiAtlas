#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#input directory, output file, starting directory
# test if there is at least one argument: if not, return an error
defaultW <- getOption("warn") 
options(warn = -1) 
if (length(args)==0) {
  #stop("At least one argument must be supplied (input file).n", call.=FALSE)
  args[1] = "/analyze_peaks_results/output"
  args[2] = paste0("core_","mergedtables_",sample(1000:9999, 1),"_objects.RData")
  args[3] = "./files"
  args[4] = 1000
  args[5] = 0.005
} else if (length(args)==1) {
  # default output file
  args[2] = paste0("core_","mergedtables_",sample(1000:9999, 1),"_objects.RData")
  args[3] = "./files"
  args[4] = 1000
  args[5] = 0.005
}else if (length(args)==2) 
{
  args[3] = "./files"
  args[4] = 1000
  args[5] = 0.005
}else if (length(args)==3) 
{
  args[4] = 1000
  args[5] = 0.005
}else if (length(args)==4) 
{
  args[5] = 0.005
}else if (length(args)==5) 
{
}else
{
  stop("Only 5 arguments are needed at most.", call.=FALSE)
}
print("TE Table merge Start")
print("Arguments: ")
print(args)
#Set up
library(readr)
library(readxl)
#library(ggplot2)
#library(ggrepel)
library(tibble)
library(reshape2)
#library(RColorBrewer)
library(dplyr)

# Initiating
## Arguments
start_path <- args[3]#"~/scratch"

group_folders = args[1]#list(paste0(start_path,args[1]))#"/code/analyze_peaks_results/epiatlas2022_100"#list(paste0(start_path,"/My Drive/code/analyze_peaks_results/ihec_set1"))
group_names = list("ihec")
savefilename=args[2] #paste0("core_","mergedtables_",sample(1000:9999, 1),"_objects.RData")
# Load needed data
source(paste0(start_path,"/lib/repeat_analysis_helper.R"))

exp_path = paste0(start_path,"/data/Joint_experiment_xmlparseout.csv")
samp_path = paste0(start_path,"/data/Joint_sample_xmlparseout.csv")
megadata_path = paste0(start_path,"/data/megadata.csv")
ihec_path = paste0(start_path,"/new_data/all_ihec/samples-17495.csv")
epidata_path = paste0(start_path,"/data/ihec_metadata_17052022.csv")#"/data/epiatlas_metadata_23022023.csv")#ihec_metadata.csv#epiatlas_metadata_07012020.csv
epidata2022_path = paste0(start_path,"/data/IHEC_metadata_harmonization.v0.8.csv")#"/data/IHEC_metadata_harmonization.v1.2.csv")#
ihecjson_path = paste0(start_path,"/new_data/all_ihec/metadata_17496.json")
ihec19json_path = paste0(start_path,"/data/hg19_portal_json.json")
cell_rename_table_path = paste0(start_path,"/data/cell_renaming_table.csv")
run_name=paste0(group_names,collapse="_")
#- extra info
exp_table <- read_csv(exp_path)
sample_table <- read_csv(samp_path)
mega_table <- read_csv(megadata_path)
ihec_table <- read_csv(ihec_path)
epiatlas_table <- read_csv(epidata_path)
epiatlas_metadata <- read_csv(epidata2022_path)
ihecjson_table <- jsonlite::fromJSON( ihecjson_path )
ihec19json_table <- jsonlite::fromJSON( ihec19json_path )
cell_type_to_cat <- make_type_to_cat()
cell_rename_table <- read_csv(cell_rename_table_path)
renaming = FALSE
if(identical(sort(rownames(cell_type_to_cat)),sort(cell_rename_table$name)))
{
  renaming = TRUE
  cell_type_to_cat <- column_to_rownames(select(cell_rename_table,"New name","New category")[!duplicated(cell_rename_table$`New name`),],"New name")
  #print out the duplicated entries
  #cell_rename_table[ which(cell_rename_table$"New name" %in% cell_rename_table[duplicated(cell_rename_table$`New name`),]$"New name"),]
  #setNames(cell_rename_table$`New name`,cell_rename_table$`New category`)
} else {
  print("Cell Types involved not the same as those accounted for by the curation. Using raw cell types")
}
#epirr_table <- make_epirr()

#- params
save_tables=FALSE
genome_size= 3234830000
trials = as.numeric(args[4])#1000
threshold = as.numeric(args[5])#0.005
stringent=TRUE#FALSE
calculate_gc = FALSE
calculate_dist = FALSE

#Build Table and Data Objects

ptm <- proc.time()


global_metrics_list = list()
global_GC_list = list()
global_all_samp_list = list()
global_samp_info_table = list()
samp_name_list = list()
metrics_list = list()
global_wind_table = list()
global_dist_table = list()
wind_count = NA
k=1
i= 1
for (i in 1:length(group_folders))
{
  group_folder <- group_folders[[i]]
  group_name <- group_names[[i]]
  GC_file_dir <- paste0(group_folder,"/GC")
  
  
  ##========================== GC
  GC_list = list()
  GC_file_names <- dir(GC_file_dir, pattern ="_GC.bed")
  batch_name= group_name
  joint_gc_file <- paste(GC_file_dir,"joint_gc.csv",sep="/")
  if (file.exists(joint_gc_file))
  {
  gc_table_i <- read.csv(joint_gc_file,header = F, sep = "", col.names = c("name","num","pct_at","pct_gc","num_A","num_C","num_G","num_T","num_oth","seq_len"))
  gc_table <- mutate(select(gc_table_i, c("name","num","pct_gc")),name = as.character(name), mean = pct_gc, batch = batch_name)
  gc_table <- mutate(gc_table, verdict = 0)
  
  gc_table <-cbind(gc_table, bind_rows(mapply(ihec_info_fun, gc_table$name ,SIMPLIFY = F)))
  
  gc_table$sample_name <- factor(gc_table$sample_name, levels = unique(gc_table$sample_name[order(as.factor(gc_table$batch),gc_table$assay)] ))
  }else
  {
    gc_table <- c()
  }
  global_GC_list[[i]] = gc_table
  
  #- analysis [IN DEV]
  file.names <- dir(group_folder, pattern = paste0("t.csv"))#"_1000t.csv"
  file.names <- file.names[!grepl("*_distribution*", file.names)]
  
  all_samp_list = list()
  all_samp_info_list = list()
  dist_table_list = list()
  wind_table_list = list()
  skipped_files = c()
  skippedwind_files = c()
  for(j in 1:length(file.names))
  {
    name1 <- file.names[j]
    full_sample_name = name1
    if ((j<100 & j%%10==0) | j%%100==0)
    {
      print(j)
      print(name1)
    }
    path_raw_repeats <- paste(group_folder,file.names[j],sep="/")
    print(path_raw_repeats)
    if (file.info(path_raw_repeats)$size==0)
    {
      print(j)
      print(name1)
      print("Skipped")
      skipped_files <- c(skipped_files, path_raw_repeats)
      next
    }
    suppressWarnings(
    raw_table <- read_csv(path_raw_repeats, show_col_types = FALSE)
    )
    print("trials")
    print(trials)
    print("tresh")
    print(threshold)
    the_samp <- analyse_samp_all(raw_table, trials, threshold)
    the_samp <- arrange(the_samp,name)
    the_samp <- mutate(the_samp, samp_name=full_sample_name)
    all_samp <- the_samp #%>% rowwise() %>% mutate(sample_name = make_sample_name(samp_name,batch))
    
    
    #sample_info <- data.frame(t(mapply(info_from_name, all_samp$sample_name[[1]], batch_name)))
    all_samp <- cbind(all_samp, ihec_info_fun(all_samp$samp_name))
    ##colnames(sample_info) <- c("assay","condition","donor")
    sample_info <- all_samp
    assay = as.character(sample_info$assay)[[1]]
    condition = as.character(sample_info$condition)[[1]]
    donor = as.character(sample_info$donor)[[1]]
    mybatch = as.character(sample_info$batch)[[1]]
    samp_name = all_samp$samp_name[[1]]
    short_name = as.character(all_samp$sample_name[[1]])
    cell_full = as.character(all_samp$cell_full[[1]])
    samp_name_list[k] <- paste(group_name,short_name,sep=".")
    verdict = 0 ##metrics_table[match(short_name, metrics_table$sample_name),]$"verdict"
    
    
    #Rearrange the data and save a summary copy(info)
    all_samp <- select(
      mutate(all_samp, 
             sample=samp_name, 
             sample_short=sample_name, 
             assay, 
             condition, 
             donor, 
             batch,
             cell_full,
             tissue,
             verdict = verdict,##metrics_table[match(short_name, metrics_table$sample_name),]$"verdict",
             med_name=samp_name_list[[k]])
      ,sample,everything())
    
    #samp_info_cols <- c("sample","sample_short","assay","condition","donor","batch","cell_full","verdict","med_name")
    #a_samp_info <- setNames(data.frame(list( samp_name, short_name, assay, condition, donor, mybatch, cell_full, verdict, samp_name_list[k])), samp_info_cols)
    
    a_samp_info <- all_samp %>% group_by(sample) %>% top_n(1,name) %>% select(sample, sample_short, sample_name, epirr_id, assay, condition, donor, batch,cell_full, tissue, biomaterial, read_count, read_bases, repeat_read_count, med_name, life_stage, age, age_unit, sex, disease) %>% ungroup()
    all_samp <- select(all_samp, -c( life_stage, age, age_unit, sex, disease) )
    #Store within full list to merge at end
    all_samp_list[[j]] = all_samp
    all_samp_info_list[[j]] = a_samp_info
    
    #=================== Distribution
    
    dist_name <- gsub(".csv","_distribution.csv",file.names[j])
    filedir_dist <- paste(group_folder,dist_name,sep="/")
    
    sampledist_table <- read.csv(filedir_dist)
    onerow <- c(short_name,samp_name, short_name, as.character(assay), as.character(condition), as.character(donor), mybatch, samp_name_list[k], as.numeric(sampledist_table$ratio))
    dist_cols <- c("short_name","sample", "name", "assay", "condition", "donor", "batch","mid_name", as.character(unlist(sampledist_table[,1])))
    dist_table <- setNames(data.frame(onerow), dist_cols)#onerow#data.frame(onerow)
    dist_table_list[[j]] = dist_table
    
    #=================== Windows
    wind_name <- gsub("_bed.*","_10kwin_overlap.bed",file.names[j])#gsub("_bed.*","_10kwin_overlap.bed",file.names[j])
    filedir_wind <- paste(paste0(group_folder,"/windows"),wind_name,sep="/")
    if (file.exists(filedir_wind))
    {
      samplewind_table <- read.table(filedir_wind, sep ="\t")
      wind_table <- setNames(data.frame(samplewind_table[[4]]), short_name)
      if (is.na(wind_count) | wind_count==nrow(wind_table))
      {
      wind_table_list[[j]] = wind_table
      wind_count <- nrow(wind_table)
      }else{
        print(paste0("skipping windows: ",filedir_wind))
        skippedwind_files <- c(skippedwind_files, filedir_wind)
      }
      
    }else
    {
      print(paste0("missing windows: ",filedir_wind))
      skippedwind_files <- c(skippedwind_files, filedir_wind)
    }
  }
  print("almost done")
  print("all sample")
  all_samp_table <- bind_rows(all_samp_list)
  global_all_samp_list[[k]] = all_samp_table
  print("all sample info")
  all_samp_info_table <- bind_rows(all_samp_info_list)
  global_samp_info_table[[k]] <- all_samp_info_table
  print("dist table")
  dist_table <- bind_rows(dist_table_list)
  global_dist_table[[k]] <- dist_table
  print("wind table")
  wind_table <- bind_cols(wind_table_list)
  global_wind_table[[k]] <- wind_table
  print("done")
  k=k+1
}
gc_table <- bind_rows(global_GC_list)
#new
all_samp_table <- bind_rows(global_all_samp_list)
all_samp_info_table <- bind_rows(global_samp_info_table)
dist_table <- bind_rows(global_dist_table)
wind_table <- bind_cols(global_wind_table)
#Long Version
dist_table.long <- melt(select(dist_table,-c(all)),id.vars = c("name","batch","assay") ,measure.vars = c(9,10,11,12,13,14), variable.name = "zone")
dist_table.long$assay <- factor(dist_table.long$assay)
dist_table.long$name <- factor(dist_table.long$name, levels = unique(dist_table.long$name[order(dist_table.long$assay)] ))
dist_table.long <- arrange(dist_table.long, assay)

options(warn = defaultW)
print(proc.time()-ptm)

remove(global_all_samp_list, global_wind_table)

# Output Data
Sys.time()
print(getwd())
print(savefilename)
save(all_samp_table, all_samp_info_table, dist_table, wind_table, file = savefilename)
print("Skipped files")
print(skipped_files)
print("Skipped windfiles")
print(skippedwind_files)

