# Library ####
{
  library(tidyverse)
  library(openxlsx)
  library(stringi)
  library(matrixStats)
  library(rstudioapi)
  library(VennDiagram)
  
  code_dir = dirname(rstudioapi::getSourceEditorContext()$path)
  main_dir = dirname(code_dir)
  main_data_dir = paste(main_dir, "data", sep = .Platform$file.sep)
}

# Parameter setup ####
{
  # Run peak picking using MZmine first 
  filename = './raw_data/mzmine.csv'
  
  # put the raw peak table file directory here
  analyzed_folder = "QTOF_6600_demo/neg"
  setwd(paste(main_data_dir, analyzed_folder, sep = .Platform$file.sep))
  analyzed_dir = getwd()
  
  # Define duplicated peaks
  mz_abs_tol = 0.002
  mz_rel_tol = 20e-6
  RT_abs_tol = 0.1 
  RT_rel_tol = 0.01 
  
  # filter out peaks that are 2x less than blanks
  high_blank_ratio = 2 
  
  # filter out a peak that its mean inten < threshold 
  min_average_inten = 500
  
  # filter out a peak that its QC peak has variation > threshold
  QC_variation = 0.5
  
  # Sample info
  sample_info = read_csv("sample_info.csv") %>%
    # left_join(read_csv('sample_order.csv')) %>%
    arrange(Cohort)
  
  cohort = sample_info$Cohort
  names(cohort) = sample_info$Samples
  
  sample_names = names(cohort)[!cohort %in% c("quality_control","blank")]
  blank_names = names(cohort)[cohort == "blank"]
  QC_names = names(cohort)[cohort == "quality_control"]
  
  select_samples = c(sample_names) # which samples are used as representatives
  
  # EVA
  eva_original_dir = 'C:/Users/Li/AppData/Local/Packages/11888HuanLab.EVAanalysistool_8tmgtde19meca/LocalState'
  
  # new folders
  dir.create('eva_output',showWarnings = F) 
  dir.create('processing_file',showWarnings = F)
} 

# Function ####
{
  
  fuzzy_group = function(data, abs_tol, rel_tol){
    if(length(data) == 1){
      return(0)
    }
    data_order = order(data)
    data_sort = data[data_order]
    count = 1
    data_group = rep(1,(length(data_sort)))
    for(i in 2:length(data_sort)){
      temp_tol = max(abs_tol, data_sort[i] * rel_tol)
      if(data_sort[i]-data_sort[i-1]>temp_tol){
        count = count+1
      }
      data_group[data_order[i]]=count
    }
    data_group
  }
  
  mz_RT_group = function(mzs, rts,
                         mz_abs_tol = 0.001, mz_rel_tol = 10e-6,
                         RT_abs_tol = 0.1, RT_rel_tol = 0.01)
  {
    group1 = fuzzy_group(mzs,abs_tol = mz_abs_tol, rel_tol = mz_rel_tol) 
    
    group2_split = split(rts, group1)
    group2_ls = sapply(group2_split, fuzzy_group, RT_abs_tol, RT_rel_tol) %>% unlist()
    seq_split_by_group1 = split(1:length(mzs), group1) %>% unlist()
    group2 = numeric(length(mzs))
    group2[seq_split_by_group1] = group2_ls 
    
    group12_split = split(mzs, paste(group1, group2))
    group12_ls = sapply(group12_split, fuzzy_group, mz_abs_tol, mz_rel_tol) %>% unlist()
    seq_split_by_group12 = split(1:length(mzs), paste(group1, group2)) %>% unlist()
    group12 = numeric(length(mzs))
    group12[seq_split_by_group12] = group12_ls 
    
    clu_group = split(1:length(mzs), paste(group1, group2, group12))
    
    group_mapping = rep(c(1:length(clu_group)), times=sapply(clu_group, length))
    names(group_mapping) = unlist(clu_group)
    
    return(group_mapping)
  }
  
  plot_filter_peak_venn = function(s, folder = './'){
    filter_venn = venn.diagram(
      x = list(
        EVA_remove = s %>% filter(EVA_remove) %>% pull(id),
        duplicated_peaks = s %>% filter(is_duplicate) %>% pull(id),
        high_blank = s %>% filter(high_blank) %>% pull(id),
        high_QC_var = s %>% filter(high_QC_var) %>% pull(id),
        low_inten = s %>% filter(low_inten) %>% pull(id)
      ),
      filename = paste0(folder,"processing_file/VENN.png"),
      imagetype = "png",
      fill = RColorBrewer::brewer.pal(5, "Dark2"),
      cat.just = list(l1=c(0.5,1.5),l2=c(0,-5),l3=c(0,0),l4=c(0,0),l5=c(1,-5)),
      disable.logging = F) %>% suppressMessages()
  }
  
}
   

# Read files #### 
{
  if(grepl('mzmine', basename(filename))){
    summary = read_csv(filename) %>%
      select(c('id','mz','rt','rt_range:min','rt_range:max'),contains(':height')) %>%
      rename_all(str_remove_all,'datafile:|.mzXML:height') %>%
      rename_all(str_remove_all,'[:blank:]') %>%
      rename_at(vars(4:5),~c('rtmin', 'rtmax')) %>%
      select(1:5, gtools::mixedsort(colnames(.)[6:ncol(.)])) %>%
      filter(T)
   
  } else {
    stop('Please set filename to either elmaven or mzmine')
  }
  
  # EVA format 
  eva_input_format = summary %>%
    select(-id) %>%
    mutate(rt = rt * 60,
           rtmin = rtmin * 60,
           rtmax = rtmax * 60) %>%
    select(mz, rt, rtmin, rtmax, all_of(gtools::mixedsort(select_samples))) %>%
    # position 5 is the first sample
    # give a small signal to avoid all NA in selected samples
    mutate_at(5, function(x)ifelse(is.na(x),1,x)) %>% 
    filter(T)
  
  write_csv(eva_input_format, './eva_output/eva_input_format.csv')
}

# Move EVA output to current folder ####
{
  output_files = list.files(eva_original_dir, full.names = T,recursive = T)
  output_files = output_files[grepl('outcome|Rplots|pie',output_files)]
  if(length(output_files) == 0){
    stop('Run EVA first.')
  } 
  if(!any(grepl('outcome', output_files))){
    stop('folder does not contain outcome.csv')
  }
  
  file.rename(
    from = output_files,
    to = paste('./eva_output',basename(output_files),sep='/'))
  
}

# Data cleaning ####
{
  # Flag eva results
  # Flag low mean inten 
  # Flag high blank 
  # Flag duplicated entries
  # Flag QC variation 
  # Remove above
  
  s = summary
  s[is.na(s)] = 0
  
  # Flag EVA reults
  eva_outcome = read_csv('./eva_output/outcome.csv') %>%
    mutate(id = as.numeric(str_extract(image, '[:digit:]+'))) %>%
    arrange(id) %>%
    select(image, prediction)
  s = s %>%
    arrange(mz) %>%
    bind_cols(eva_outcome) %>%
    mutate(EVA_remove = !prediction)
  
  # Flag duplicated 
  peak_groups = mz_RT_group(s$mz, s$rt, 
                          mz_abs_tol, mz_rel_tol,
                          RT_abs_tol, RT_rel_tol)
  
  
  s_duplicated = s %>%
    mutate(peak_group = peak_groups[as.character(1:nrow(.))]) %>%
    group_by(peak_group) %>%
    mutate(n=n()) %>%
    filter(n>1) 
  
  if(nrow(s_duplicated)!=0){
    s_duplicated = s_duplicated %>%
      mutate_at(all_of(sample_info$Samples), function(x){
        max(x, na.rm=T)
      }) %>% 
      arrange(EVA_remove, id) %>%
      mutate(is_duplicate = c(F,rep(T,times=n()-1))) %>%
      ungroup()
    
    s = s_duplicated %>%
      bind_rows(s) %>%
      distinct(id, .keep_all=T) %>%
      mutate(is_duplicate = replace_na(is_duplicate, F))
  } else {
    s = s %>%
      mutate(is_duplicate = F)
  }
  
  # Flag min_intensity 
  s = s %>%
    mutate(mean_inten = rowMeans(.[,select_samples], na.rm = T),
           low_inten =  mean_inten < min_average_inten)
  
  
  # Flag high blank
  if(sum(cohort == "blank")>0){
    s = s %>%
      mutate(high_blank = rowMeans(.[,select_samples, na.rm=T]) < 
               high_blank_ratio*rowMeans(.[,blank_names], na.rm=T))
  }else {
    s = s %>%
      mutate(high_blank = F)
  }
  
  
  # Flag QC variation or QC no signal
  if(sum(cohort == "quality_control")>0){
    s = s %>%
      mutate(rel_QC_variation = sqrt(rowVars(.[,QC_names] %>% as.matrix(), na.rm = FALSE)) / 
               rowMeans(.[,QC_names],na.rm=T)) %>%
      mutate(high_QC_var = rel_QC_variation > QC_variation | rowMeans(.[,QC_names]) == 0)
  } else {
    s = s %>%
      mutate(high_QC_var = F)
  }
  
  summary_cleanup = s %>%
    filter(!EVA_remove) %>%
    filter(!low_inten) %>%
    filter(!high_blank) %>%
    filter(!is_duplicate) %>%
    filter(!high_QC_var) %>%
    select(colnames(summary), mean_inten) %>%
    # arrange(-mean_inten) %>%
    filter(T)

  
}


# Output ####
{
  
  plot_filter_peak_venn(s)
  
  write_csv(summary_cleanup, "peak_table.csv", na="")
  write_csv(s, "./processing_file/peak_table_raw.csv", na="")

}





