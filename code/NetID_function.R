# NetID Fucntions ####
library(tidyverse)
library(openxlsx)
library(stringi)
library(janitor)

library(lc8)
library(enviPat)
library(slam)
library(pracma)
library(igraph)
library(lpsymphony)
library(cplexAPI)

# Read_Related_Files read related files ####
Read_Related_Files = function(
  ion_mode = -1,
  neg_MS2_library_file = "",
  pos_MS2_library_file = "",
  HMDB_library_file = "",
  known_library_file = "",
  LC_method = "",
  manual_library_file = "",
  empirical_rule_file = "",
  propagation_rule_file = "",
  MS_type = "", 
  instrument_parameter_type = "", 
  instrument_parameter_custom = NULL
){
  Related_files = list()
  Related_files[["HMDB_library"]] = readRDS(HMDB_library_file) %>% filter(grepl("HMDB", accession))
  
  # known_library contains RT information of documented metabolites
  {
    Related_files[["known_library"]] = read.xlsx(known_library_file) %>%
      left_join(Related_files$HMDB_library %>%
                  dplyr::select(accession, mass) %>%
                  dplyr::rename(mass_hmdb=mass), by = "accession") %>%
      mutate(mass = ifelse(is.na(mass_hmdb), mass,mass_hmdb)) %>%
      dplyr::select(-mass_hmdb)
     
    if(LC_method %in% colnames(Related_files[["known_library"]])){
      Related_files[["known_library"]] = Related_files[["known_library"]] %>%
        filter(!is.na(.[,eval(LC_method)]))
    } else {
      Related_files[["known_library"]] = Related_files[["known_library"]] %>%
        filter(FALSE)
    }
  }
  # Read in manual_library, which is data dependent, unlike known_library is global applicable.
  # manual_library provides a more flexible way to temporarily add annotations
  {
    if(!file.exists(manual_library_file)){
      # warning("No manual library file found in data folder.")
      manual_library = NULL
    } else {
      data("isotopes")
      manual_library = read.csv(manual_library_file, stringsAsFactors = FALSE)
      if(ion_mode!=0){
        manual_library = manual_library %>% filter(mode %in% c(ion_mode, 0))
      }
      check_formula = check_chemform(isotopes, manual_library$formula) 
      if(any(check_formula$warning)){
        stop(paste("Check manual library for formula error:", 
                   paste(check_formula$new_formula[check_formula$warning], collapse = ", ")))
      }
      manual_library = manual_library %>%
        mutate(formula = check_formula$new_formula) %>%
        mutate(formula = my_calculate_formula(formula, "C1"),
               formula = my_calculate_formula(formula, "C-1")) %>%
        mutate(formula = as.character(formula)) %>%
        mutate(mass = formula_mz(formula),
               rdbe = formula_rdbe(formula)) %>%
        mutate(note = as.character(note))
      
    }
    
    
    Related_files[["manual_library"]] = manual_library
  }
  
  # Read empirical_rules
  {
    data("isotopes")
    Connect_rules = read.csv(empirical_rule_file,stringsAsFactors = FALSE)
    if(nrow(Connect_rules) == 0){return(Connect_rules)}
    for(i in 1: nrow(Connect_rules)){
      Connect_rules$formula[i] = check_chemform(isotopes,Connect_rules$formula[i])$new_formula
      Connect_rules$formula[i] = my_calculate_formula(Connect_rules$formula[i], "C1")
      Connect_rules$formula[i] = my_calculate_formula(Connect_rules$formula[i], "C1", -1)
      Connect_rules$mass[i] = formula_mz(Connect_rules$formula[i])
    }
    Related_files[["empirical_rules"]] = Connect_rules
  }
  
  # Read global_parameter
  {
    propagation_rule = read.csv(propagation_rule_file, row.names = 1)
    propagation_rule_ls = list()
    for(i in colnames(propagation_rule)){
      propagation_rule_ls[[i]] = colnames(propagation_rule)[propagation_rule[,i]]
    }
    
    if(grepl('Custom',instrument_parameter_type,ignore.case = T)){
      instrument_parameter = instrument_parameter_custom
    } else if(MS_type == "Orbitrap"){
      instrument_parameter = data.frame(ppm = 10e-6,
                                        ms1_ms2_match_ppm = 20e-6,
                                        propagation_ppm = 10e-6,
                                        record_ppm = 5e-6,
                                        edge_expand_inten_cutoff = 5e4,
                                        score_mz_ppm_error_rate = -0.5)
    } else if(MS_type == "QTOF"){
      instrument_parameter = data.frame(ppm = 10e-6,
                                        ms1_ms2_match_ppm = 20e-6,
                                        propagation_ppm = 10e-6,
                                        record_ppm = 10e-6,
                                        edge_expand_inten_cutoff = 1e3,
                                        score_mz_ppm_error_rate = -0.125)
    }
    
    Related_files[["global_parameter"]] = list(mode = ion_mode,
                                      LC_method = LC_method,
                                      MS_type = MS_type,
                                      instrument_parameter = instrument_parameter,
                                      propagation_rule = propagation_rule_ls)
  }
  
  # Read MS2_library
  {
    if(ion_mode == -1){
      Related_files[["neg_MS2_library"]] = readRDS(neg_MS2_library_file)
    } else if(ion_mode == 1){
      Related_files[["pos_MS2_library"]] = readRDS(pos_MS2_library_file)
    } else {
      stop(paste("Allowed ion_mode = -1 (negative ionization) or 1 (positive ionization)"))
    }
  }
  dir.create(paste(MASTER_ADDRESS, "/log", sep = ""))
  return(Related_files)
  
}

## Fill in missing values ####
detail_fill = function(detail){
  
  all_names=unique(detail$sample)
  
  if(length(grep("blank|blk", all_names, ignore.case = TRUE))!=0){
    sample_names=all_names[-grep("blank|blk", all_names, ignore.case = TRUE)]
  } else {
    sample_names=all_names
  } 
  sample_len = length(sample_names)
  detail =  detail %>%
    group_by(groupId) %>% mutate(group_label = cur_group_id()) %>% ungroup()
  detail = apply(detail, 1, function(peak){
    groupId = as.numeric(peak[["groupId"]])
    sample = peak[["sample"]]
    if(is.na(peak[["peakMz"]])){
      group_label = as.numeric(peak[["group_label"]])
      group_index = ((group_label-1)*sample_len+1) : (group_label*sample_len)
      peakMz = median(detail[group_index,]$peakMz, na.rm=TRUE)
      rt = median(detail[group_index,]$rt, na.rm=TRUE)
      rtmax = median(detail[group_index,]$rtmax, na.rm=TRUE)
      rtmin = median(detail[group_index,]$rtmin, na.rm=TRUE)
      peakIntensity = median(detail[group_index,]$peakIntensity, na.rm=TRUE)
      mzmax = median(detail[group_index,]$mzmax, na.rm=TRUE)
      mzmin = median(detail[group_index,]$mzmin, na.rm=TRUE)
      peakArea = median(detail[group_index,]$peakArea, na.rm=TRUE)
      peak = data.frame(groupId = groupId, sample = sample, peakMz = peakMz, mzmax = mzmax,
                        mzmin = mzmin, peakArea = peakArea, peakIntensity = peakIntensity,
                        rt = rt, rtmax = rtmax, rtmin = rtmin)
      
      # unmissing = which(!is.na(detail$peakMz[group_index]))[1]
      # id = group_index[unmissing]
      # peak = detail[id,,drop=FALSE] %>% dplyr::select(-group_label)
    }else{
      peak = data.frame(groupId = groupId,
                        sample = sample,
                        peakMz = as.numeric(peak[["peakMz"]]),
                        mzmax = as.numeric(peak[["mzmax"]]),
                        mzmin = as.numeric(peak[["mzmin"]]),
                        peakArea = as.numeric(peak[["peakArea"]]),
                        peakIntensity = as.numeric(peak[["peakIntensity"]]),
                        rt = as.numeric(peak[["rt"]]),
                        rtmax = as.numeric(peak[["rtmax"]]),
                        rtmin = as.numeric(peak[["rtmin"]]))
    }
    peak
    
  }) %>% bind_rows() 
}

## For netid input ####
Write_NetID_Input <- function(info, em_info, mode){
  inds <- unique(em_info$id) %in% info@id
  sample_names <- unique(em_info$sample)
  medMz <- sapply(info@Mz, mean)
  medRt <- (sapply(info@Rt, mean))/60
  ints <-t(matrix(em_info$peakIntensity, nrow = length(sample_names)))[inds, , drop = FALSE]
  
  base::colnames(ints) = sample_names
  # maxQuality <- matrix(em_info$quality, nrow = length(sample_names)) %>% apply(., 2, max)
  a <- rep(NA, length(info@id))
  output <-data.frame(label = a, metaGroupId = a,
                      id = info@id, goodPeakCount = a,
                      medMz = medMz, medRt = medRt,
                      maxQuality = a, isotopeLabel = a,
                      compound = a, compoundId = a,
                      formula = a, expectedRtDiff = a,
                      ppmDiff = a, parent = a) %>% cbind(., ints)
  
}

# Read_Files read pos and neg files ####
Read_Files = function(
  Mset,
  Related_files,
  peak_table_files,
  untarget_raw_files,
  sample_filter
){
  H_mass = 1.00782503224
  e_mass = 0.00054857990943
  ion_mode = Related_files$global_parameter$mode
  
  # Read detail peak table
  table_label = case_when(
    grepl("pos", peak_table_files) ~ "pos",
    grepl("neg", peak_table_files) ~ "neg")
    
  Mset[[paste(table_label, "peak_table", sep = "_")]] = read_csv(peak_table_files) %>%
    filter(grepl(sample_filter, sample)) %>% detail_fill()
  
  peak_table_colnames = unique(Mset[[paste(table_label, "peak_table", sep = "_")]]$sample)
  mzXML_file_names = names(untarget_raw_files)
  intersect_names = intersect(peak_table_colnames, mzXML_file_names)
  mzXML_file_path = untarget_raw_files[intersect_names]
  
  
  Mset[[paste(table_label, "sample_names", sep='_')]]=intersect_names
  
  if(length(intersect_names) == 0){
    stop(paste('None of the mzXML files exist in peak table. Check files and filter settings.'))
  }
  
  if(!all(mzXML_file_names %in% intersect_names) | 
     !all(peak_table_colnames %in% intersect_names)){
    warning(paste('mzXML files and samples in peak table do not fully match.',
                  'Analysis is performed using the common files.'))
  }
    
  # Get EICs
  cat(crayon::green("Getting EICs...\n"))
  Mset = Get_EICs(Mset, Related_files,
                  mzXML_file_path = mzXML_file_path)
  
  # delete peaks with bad chromatogram
 if(ion_mode==1){
    Mset$pos_peak_table = Mset$pos_peak_table %>% 
      filter(groupId %in% Mset$pos_msdata@id) %>%
      mutate(peakMz = peakMz - (H_mass-e_mass)*1) %>% 
      dplyr::rename(id = groupId)
    Mset[["pos_raw_data"]] = Write_NetID_Input(Mset$pos_msdata, Mset$pos_peak_table, 1) %>% 
      mutate(medMz = medMz - (H_mass-e_mass)*1)
  } else{
    Mset$neg_peak_table = Mset$neg_peak_table %>%
      filter(groupId %in% Mset$neg_msdata@id) %>%
      mutate(peakMz = peakMz - (H_mass-e_mass)*(-1)) %>% 
      dplyr::rename(id = groupId)
    Mset[["neg_raw_data"]] = Write_NetID_Input(Mset$neg_msdata, Mset$neg_peak_table, -1) %>% 
      mutate(medMz = medMz - (H_mass-e_mass)*(-1))
  }
  return(Mset)
}

### fuzzy_match ####
### More memory-efficient matching
fuzzy_match = function(data, ref, abs_tol = 0.002, rel_tol = 10e-6, error_type = c("abs", "relative")){
  ref = c(ref, -Inf, Inf)
  ref_sort = sort(ref)
  data_sort = sort(data)
  ref_order = order(ref)
  data_order = order(data)
  i_min=1
  j=1
  match_position = list()
  
  while(j<=length(data_sort)){
    if(error_type == "abs"){temp_data_tol = abs_tol}
    else{temp_data_tol = data_sort[j]*rel_tol}
    
    while(ref_sort[i_min+1] < (data_sort[j]-temp_data_tol)){
      i_min=i_min+1
    }
    k=1
    while(ref_sort[i_min+k]<=(data_sort[j]+temp_data_tol)){
      k=k+1
    }
    if(k!=1){
      match_position[[length(match_position)+1]] = data.frame(data_position = data_order[j],
                                                              ref_position = ref_order[(i_min+1):(i_min+k-1)])
    }
    j=j+1
  }
  return(match_position)
}


## MS1_MS2_match ####
MS1_MS2_Match = function(ms1_data, ms2_data, ms1_spec, mz_ppm){
  ms1_mz = ms1_data$peakMz
  # ms1_rt = ms1_data$rt
  ms1_groupId = ms1_data$id
  ms2_mz = unlist(lapply(ms2_data, function(x){
    return(x[[1]][1])
  }))
  ms2_rt = unlist(lapply(ms2_data, function(x){
    return(x[[1]][2])
  }))
  
  mz_match = fuzzy_match(ms1_mz, ms2_mz, abs_tol = 0.0002, rel_tol = mz_ppm, error_type = "relative") %>%
    bind_rows() %>%
    apply(., 1, function(x){
      x[1] = ms1_data[x[1], ]$id
      return(x)
    }) %>% t() %>% as.data.frame() %>% distinct(.,.keep_all = TRUE)
  group = apply(mz_match, 1, function(x){
    ms1_groupId = x[1]
    ms2_index = x[2]
    rt_ms2 = ms2_rt[ms2_index]
    ms1_spec_index = which(ms1_spec@id == ms1_groupId)
    ms1_rt_max = max(unlist(ms1_spec@rtime_list[[ms1_spec_index]]))/60
    ms1_rt_min = min(unlist(ms1_spec@rtime_list[[ms1_spec_index]]))/60
    if(rt_ms2 <= ms1_rt_max & rt_ms2 >= ms1_rt_min){
      return(x)
    } else{
      x = NULL
      return(x)
    }
  })
  group = group[!sapply(group, is.null)] %>% bind_rows()

  temp_ms1_idx = unlist(unique(group[,1]))

  match_result <- lapply(temp_ms1_idx, function(idx) {
    idx2 <- unlist(group[which(group[,1] == idx), 2])
    if (length(idx2) == 1) {
      return(c(idx, idx2))
    } else{
      temp_ms2_info <- ms2_data[idx2]
      return(c(idx, idx2[which.max(unlist(lapply(temp_ms2_info, function(y) {
        y = y[[2]]
        y <- y[order(y[, 2], decreasing = TRUE), , drop = FALSE]
        if (nrow(y) > 5){
          y <- y[1:5, ]
        }
        return(sum(y[, 2]))
      })))]))
    }
  }) %>% lapply(.,t) %>% lapply(., function(x){
    x = as.data.frame(x)
    colnames(x) = c("idx1","idx2")
    return(x)
  }) %>% bind_rows()
  
  return(match_result)
}
## read_mgf ####
read_mgf = function(file) {
  ListMGF = function(file) {
    mgf.data <- readLines(file)
    nl.rec.new <- 1
    idx.rec <- 1
    rec.list <- list()
    for (nl in seq_along(mgf.data))
    {
      if (mgf.data[nl] == "END IONS")
      {
        rec.list[idx.rec] <- list(Compound = mgf.data[nl.rec.new:nl])
        nl.rec.new <- nl + 1
        idx.rec <- idx.rec + 1
      }
    }
    rec.list
  }
  # main
  pbapply::pboptions(style = 1)
  cat(crayon::green("Reading mgf data...\n"))
  # mgf.data.list <- pbapply::pblapply(file, ListMGF)
  ms2 <- purrr::map(
    .x = file,
    .f = function(mgf.data) {
      mgf.data <- ListMGF(mgf.data)
      # nl.spec <- grep('^\\d', mgf.data)
      nl.spec <-
        lapply(mgf.data, function(x)
          grep('^\\d', x))
      info.mz <-
        lapply(mgf.data, function(x)
          grep('^PEPMASS', x, value = T))
      info.rt <-
        lapply(mgf.data, function(x)
          grep('^RTINSECONDS', x, value = T))
      
      info.mz <- unlist(info.mz)
      #for orbitrap data, the intensity of precursor ion should be removed
      info.mz <-
        unlist(lapply(strsplit(x = info.mz, split = " "), function(x)
          x[1]))
      info.mz <-
        as.numeric(gsub(pattern = "\\w+=", "", info.mz))
      info.rt <- unlist(info.rt)
      info.rt <-
        as.numeric(gsub(pattern = "\\w+=", "", info.rt))
      
      if (length(mgf.data) == 1) {
        spec <- mapply(function(x, y) {
          temp <- do.call(rbind, strsplit(x[y], split = " "))
          list(temp)
        },
        x = mgf.data,
        y = nl.spec)
      } else{
        spec <- mapply(function(x, y) {
          do.call(rbind, strsplit(x[y], split = " "))
        },
        x = mgf.data,
        y = nl.spec)
      }
      
      spec <- lapply(spec, function(x) {
        temp <- cbind(as.numeric(x[, 1]), as.numeric(x[, 2]))
        temp <- matrix(temp, ncol = 2)
        # if(nrow(temp) > 0) temp <- temp[temp[,2] >= max(temp[,2])*0.01,]
        temp <- matrix(temp, ncol = 2)
        colnames(temp) <- c("mz", "intensity")
        temp
      })
      
      ms2 <- mapply(function(x, y, z) {
        info <- c(y, z)
        names(info) <- c("mz", "rt")
        spectrum <- as.matrix(x)
        temp <- list(info, spectrum)
        names(temp) <- c("info", "spec")
        list(temp)
      },
      x = spec,
      y = info.mz,
      z = info.rt)
      
      ms2
      
    }
  )
  
  spec.info <- ms2[[1]]
  if (length(ms2) > 1) {
    for (i in 2:length(ms2)) {
      spec.info <- c(spec.info, ms2[[i]])
    }
  }
  
  remove.idx <-
    which(unlist(lapply(spec.info, function(x)
      nrow(x[[2]]))) == 0)
  if (length(remove.idx) != 0)
    spec.info <- spec.info[-remove.idx]
  
  spec.info <- spec.info
  spec.info <- lapply(spec.info, function(x){
    x[[1]][2] = x[[1]][2]/60
    return(x)
  })
  return(spec.info)
}
# Read_MS2data read_MS2data ####
Read_MS2data = function(Mset, Related_files, MS2_folder){
 
  # ms1_ms2_flow workflow
  ms1_ms2_flow = function(ms1_data, ms, ms2_filenames, mode, mz_ppm){
    H_mass = 1.00782503224
    e_mass = 0.00054857990943
    ms2_data = read_mgf(ms2_filenames)
    if(mode == "pos"){
      ms1_data$peakMz = ms1_data$peakMz + (H_mass-e_mass)*1
    } else{
      ms1_data$peakMz = ms1_data$peakMz + (H_mass-e_mass)*(-1)
    }
    ms1_ms2_matches = MS1_MS2_Match(ms1_data, ms2_data, ms, mz_ppm)
    ms2_data = ms2_data[as.numeric(ms1_ms2_matches[,2])]
    names(ms2_data) = as.character(ms1_ms2_matches[,1])
    ms2_data = ms2_data %>% lapply(.,function(x){
      x[[1]] = mode
      x[[2]] = as.data.frame(x[[2]])
      return(x)
    })
    names(ms2_data) = as.character(ms1_ms2_matches[,1])
    return(ms2_data)
  }

  # main
  input_mode = Related_files$global_parameter$mode
  neg_ms2_filenames = dir(path = MS2_folder$neg, pattern = "mgf", full.names = TRUE)
  pos_ms2_filenames = dir(path = MS2_folder$pos, pattern = "mgf", full.names = TRUE)
  
  if(input_mode >= 0){
    if(length(pos_ms2_filenames)==0){
      stop('No pos MS2 files (.mgf) are found. Please check files and data path.')
    }
    pos_ms1_data = Mset$pos_peak_table
    pos_ms = Mset$pos_ms
    Mset[["MS2_ls"]] = ms1_ms2_flow(pos_ms1_data, pos_ms, pos_ms2_filenames, mode = "pos", 
                                    mz_ppm = Related_files$global_parameter$instrument_parameter$ms1_ms2_match_ppm)
  } 
  if(input_mode <= 0){
    if(length(neg_ms2_filenames)==0){
      stop('No neg MS2 files (.mgf) are found. Please check files and data path.')
    }
    neg_ms1_data = Mset$neg_peak_table
    neg_ms = Mset$neg_ms
    Mset[["MS2_ls"]] = ms1_ms2_flow(neg_ms1_data, neg_ms, neg_ms2_filenames, mode = "neg", 
                                    mz_ppm = Related_files$global_parameter$instrument_parameter$ms1_ms2_match_ppm)
  } 
  return(Mset)
}





# Initiate NodeSet ####
## NodeSet = c(single_pos, both, single_neg)
## each node has two mode information(if exists): pos_info, neg_info
## each mode information has three kinds of peak info: summary data(raw_data), detail data(peak), chromatographic peak data(ms)
Initiate_Nodeset = function(Mset){
 
  if(!is.null(Mset$neg_raw_data) & is.null(Mset$pos_raw_data)){
    neg = Mset$neg_msdata
    neg_raw_data = Mset$neg_raw_data %>% 
      dplyr::select(-c(label, metaGroupId, goodPeakCount, isotopeLabel, compound, compoundId, formula, expectedRtDiff, ppmDiff, parent))
    neg_inten = log10(rowMeans(neg_raw_data[, Mset$neg_sample_names, drop=FALSE], na.rm = TRUE))
    node_neg =list()
    for(i in 1:nrow(Mset$neg_raw_data)){
      nid = Mset$neg_raw_data$id[i]
      node = list(
          pos_info = NULL,
          neg_info = list(
            neg_raw = as.list(neg_raw_data[i,] %>%
                                mutate(inten = neg_inten[i],
                                      sample_inten = .[Mset$neg_sample_names]
                                ) %>%
                                dplyr::select(-c(Mset$neg_sample_names))),
            neg_peak = as.list(Mset$neg_peak_table[Mset$neg_peak_table$id == nid,]),
            neg_ms = list(
              id = nid,
              intensity_list = neg@intensity_list[[i]],
              rtime_list = neg@rtime_list[[i]]
            )
                   
          ),
          MS2 = NULL
                 
        )
        node_neg[[i]] = node
        
    }
    node_neg = node_neg[!sapply(node_neg, is.null)]
    NodeSet = node_neg
    names(NodeSet) = 1:length(NodeSet)
    {
      # map MS2 in nodeset
      MS2_ls = Mset$MS2_ls 
      # MS2_ls = lapply(MS2_ls, function(x){as.data.frame(x[[2]])})
      if(!is.null(MS2_ls)){
        MS2_ls_names = names(MS2_ls)
        id_map = 1:nrow(Mset$neg_raw_data)
        names(id_map) = as.character(Mset$neg_raw_data$id)
        for(i in 1:length(MS2_ls_names)){
          node_id = id_map[MS2_ls_names[i]]
          if(is.na(node_id)){next}
          MS2 = list(pos_MS2 = NULL,
                     neg_MS2 = MS2_ls[[i]])
          NodeSet[[node_id]][["MS2"]] = MS2
          
        }
        
      }
    }
    print(paste(length(NodeSet), "negative nodes"))
    return(NodeSet)
  } else{
    pos = Mset$pos_msdata
    pos_raw_data = Mset$pos_raw_data %>% 
      dplyr::select(-c(label, metaGroupId, goodPeakCount, isotopeLabel, compound, compoundId, formula, expectedRtDiff, ppmDiff, parent))
    pos_inten = log10(rowMeans(pos_raw_data[, Mset$pos_sample_names, drop=FALSE], na.rm = TRUE))
    node_pos =list()
    for(i in 1:nrow(Mset$pos_raw_data)){
      pid = Mset$pos_raw_data$id[i]
      
      node = list(
        pos_info = list(
          pos_raw = as.list(pos_raw_data[i,] %>%
                              mutate(inten = pos_inten[i],
                                     sample_inten = .[Mset$pos_sample_names]
                              ) %>%
                              dplyr::select(-c(Mset$pos_sample_names))),
          pos_peak = as.list(Mset$pos_peak_table[Mset$pos_peak_table$id == pid,]),
          pos_ms = list(
            id = pid,
            intensity_list = pos@intensity_list[[i]],
            rtime_list = pos@rtime_list[[i]]
          )
        ),
        neg_info = NULL,
        MS2 = NULL
      )
      node_pos[[i]] = node
     
    }
    node_pos = node_pos[!sapply(node_pos, is.null)]
    NodeSet = node_pos
    names(NodeSet) = 1:length(NodeSet)
    {
      # map MS2 in nodeset
      MS2_ls = Mset$MS2_ls 
      # MS2_ls = lapply(MS2_ls, function(x){as.data.frame(x[[2]])})
      if(!is.null(MS2_ls)){
        MS2_ls_names = names(MS2_ls)
        id_map = 1:nrow(Mset$pos_raw_data)
        names(id_map) = as.character(Mset$pos_raw_data$id)
        for(i in 1:length(MS2_ls_names)){
          node_id = id_map[MS2_ls_names[i]]
          if(is.na(node_id)){next}
          MS2 = list(pos_MS2 = MS2_ls[[i]],
                     neg_MS2 = NULL)
          NodeSet[[node_id]][["MS2"]] = MS2
         
          
        }
        
      }
    }
    print(paste(length(NodeSet), "positive nodes"))
    return(NodeSet)
  }
  
}

## match_edge ####
match_edge = function(temp_fg,temp_mz_list,mz_tol_ppm,mz_tol_abs,temp_RT_list,temp_deltaRT){
  if(temp_fg>0){
    temp_mz_list = temp_mz_list %>% sort()
  } else{temp_mz_list = temp_mz_list %>% sort(decreasing = TRUE)}
  temp_id_list = names(temp_mz_list)
  temp_edge_list = list()
  # Find the i,j combination gives potential transformation 
  ## Memeory efficient, but may be slower
  ## matrix calculation could be a faster & simpler approach 
  i=j=j_pos=1
  merge_nrow = length(temp_mz_list)
  while(i<=merge_nrow){
    
    while(TRUE){
      j=j+1
      if(j>merge_nrow){break}
      temp_ms = temp_mz_list[j]-temp_mz_list[i]
      
      mass_tol = max(max(temp_mz_list[i],temp_mz_list[j])*mz_tol_ppm,mz_tol_abs)
      if(temp_fg>0){
        if(temp_ms < (temp_fg - mass_tol)){
          j_pos = j # locate the last j that has smaller ms
          next
        }
      } else{
        if(temp_ms > (temp_fg + mass_tol)){
          j_pos = j # locate the last j that has smaller ms
          next
        }
      }
      # Criteria to entry
      if(abs(temp_ms-temp_fg)<mass_tol){
        delta_RT = temp_RT_list[names(temp_mz_list[j])] - temp_RT_list[names(temp_mz_list[i])]
        if(abs(delta_RT) < temp_deltaRT){
          temp_edge_list[[length(temp_edge_list)+1]]= list(node1=temp_id_list[i],
                                                           node2=temp_id_list[j], 
                                                           mass_dif=(temp_ms-temp_fg)/temp_mz_list[j]*1E6)
        }
      }
      if(temp_fg>0){
        if(temp_ms > (temp_fg + mass_tol)){break}
      } else{
        if(temp_ms < (temp_fg - mass_tol)){break}
      }
    }
    i = i + 1
    j = j_pos - 1
  }
  return(temp_edge_list)
}
# Initiate_edgeset ####
Initiate_edgeset = function(Mset, Related_files, NodeSet, mz_tol_abs = 0, mz_tol_ppm = 10, 
                            rt_tol_bio = Inf, rt_tol_nonbio = 0.2){
  mz_tol_ppm = mz_tol_ppm
  temp_mz_list = NodeSet %>% sapply(get_mz)
  temp_RT_list = NodeSet %>% sapply(get_rt) 
  
  temp_rules = Related_files$empirical_rules %>% arrange(mass)
  
  {
    edge_ls = list()
    for (k in 1:nrow(Related_files$empirical_rules)){
      temp_fg=temp_rules$mass[k]
      temp_deltaRT = ifelse(temp_rules$category[k] == "Biotransform", rt_tol_bio, rt_tol_nonbio)
      temp_edge_list = match_edge(temp_fg,temp_mz_list,mz_tol_ppm,mz_tol_abs,temp_RT_list,temp_deltaRT)
      
      edge_ls[[k]]= bind_rows(temp_edge_list) %>%
        mutate(category = temp_rules$category[k],
               linktype = temp_rules$formula[k],
               direction = temp_rules$direction[k],
               rdbe = temp_rules$rdbe[k])
      # print(paste(temp_rules$category[k], temp_rules$formula[k], nrow(edge_ls[[k]]),"found."))
    }
  }
  
  
  edge_list = bind_rows(edge_ls) %>%
    filter(!linktype == "") %>% # remove data-data isomer connection
    mutate(node1 = as.numeric(node1),
           node2 = as.numeric(node2)) %>%
    mutate(edge_id = 1:nrow(.)) %>%
    filter(TRUE)
  
  print(table(edge_list$category))
  
  
  return(edge_list)
}

# Initiate_libraryset ####
Initiate_Libraryset = function(Mset, Related_files){
  Metabolites_HMDB = Related_files$HMDB_library %>%
    dplyr::rename(note = accession) %>%
    mutate(origin = "HMDB_library")
  Metabolites_HMDB[which(is.na(Metabolites_HMDB$PubMed_Count)),]$PubMed_Count = 0
  Metabolites_known = Related_files$known_library %>%
    mutate(category = "Metabolite") %>%
    dplyr::rename(note = accession) %>%
    mutate(origin = "known_library")
  if(Related_files$global_parameter$LC_method %in% colnames(Related_files[["known_library"]])){
    Metabolites_known = Metabolites_known %>%
      mutate(rt = .[,eval(Related_files$global_parameter$LC_method)]) %>%
      mutate(rt = as.numeric(unlist(rt)))
  }
  
  Adducts = Related_files$empirical_rules %>%
    filter(category == "Adduct") %>%
    # mutate(category = "Artifact") %>%
    dplyr::select(-direction) %>%
    mutate(origin = "empirical_rules")
  
  if(!is.null(Related_files$manual_library)){
    Manual = Related_files$manual_library %>%
      mutate(origin = "manual_library")
  } else {
    Manual = NULL
  }
  
  # Remove entries in HMDB that are adducts 
  Metabolites_HMDB = Metabolites_HMDB %>%
    filter(!formula %in% Adducts$formula)
 
  if(!is.null(Mset$pos_raw_data) & is.null(Mset$neg_raw_data)){
    count = nrow(Mset$pos_raw_data)
  } else if(is.null(Mset$pos_raw_data) & !is.null(Mset$neg_raw_data)){
    count = nrow(Mset$neg_raw_data)
  }
  Metabolites_known = Metabolites_known %>% left_join(Metabolites_HMDB %>%
                                                        dplyr::select(note, PubMed_Count, FirstBlock), by="note")
  Metabolites_known[which(is.na(Metabolites_known$PubMed_Count)),]$PubMed_Count = 0
  Metabolites = bind_rows(Metabolites_known, Metabolites_HMDB) %>%
    arrange(note)
  if(!is.null(Metabolites$rt)){
    Metabolites_HMDB = Metabolites_HMDB %>% left_join(Metabolites_known %>% dplyr::select(note, rt), by = "note")
    Metabolites = bind_rows(Metabolites_known, Metabolites_HMDB)
    Metabolites_1 = Metabolites %>% filter(!is.na(rt))
    Metabolites_2 = Metabolites %>% filter(is.na(rt)) %>% arrange(note)
    Metabolites_1 = Metabolites_1 %>% left_join(Metabolites_1 %>%
                                            group_by(FirstBlock) %>%
                                            summarise_at(vars(rt), list(rt = mean), na.rm = TRUE) %>%
                                            mutate(rt = ifelse(is.nan(rt), NA, rt)), by = "FirstBlock") %>%
      dplyr::select(-rt.x) %>%
      dplyr::rename(rt = rt.y) %>%
      arrange(note)
    Metabolites = bind_rows(Metabolites_1, Metabolites_2) %>% arrange(note)
  }
   
  LibrarySet = bind_rows(Metabolites, Manual, Adducts) %>%
    distinct(SMILES, formula, .keep_all = TRUE) %>%
    mutate(library_id = (1+count) : (count+nrow(.))) %>%
    # group_by(SMILES) %>%
    # filter(n()>1) %>%
    filter(TRUE)
  
  return(LibrarySet)
}

## expand_library ####
expand_library = function(lib, rule, direction, category){
  if(nrow(rule) == 0){
    return(NULL)
  }
  # initial_lib_adduct_1 = expand_library(lib_adduct, rule_1, direction = 1, category = "Artifact")
  mz_lib = lib$mass
  mz_rule = rule$mass
  mass_exp = outer(mz_lib, mz_rule * direction, FUN = "+")
  
  rdbe_lib = lib$rdbe
  rdbe_rule = rule$rdbe
  rdbe_exp = outer(rdbe_lib, rdbe_rule * direction, FUN = "+")
  
  formula_lib = lib$formula
  formula_rule = rule$formula
  formula_exp = my_calculate_formula(formula_lib, formula_rule, sign = direction)
  
  parent_lib = lib$formula
  transformation_rule = rule$formula
  
  expansion = data.frame(parent_formula = rep(parent_lib, length(transformation_rule)),
                         transform = rep(transformation_rule, each = length(parent_lib))) %>%
    mutate(parent_formula = as.character(parent_formula),
           transform = as.character(transform),
           mass = as.vector(mass_exp),
           rdbe = as.vector(rdbe_exp),
           formula = as.vector(formula_exp),
           direction = direction,
           parent_id = rep(lib$node_id, nrow(rule)),
           category = category) %>%
    dplyr::select(formula, mass, rdbe, category, parent_id, parent_formula, transform, direction)
  return(expansion)
}
# Expand_libraryset ####
Expand_Libraryset = function(LibrarySet){
  ## initialize
  seed_library = LibrarySet %>% 
    mutate(node_id = library_id,
           formula = formula,
           mass = mass,
           parent_id = library_id,
           parent_formula = formula,
           transform = "",
           direction = 1,
           rdbe = rdbe,
           category = category,
           status = status,
           PubMed_Count = PubMed_Count,
           accession = note,
           FirstBlock = FirstBlock
    ) %>%
    dplyr::select(node_id, formula, mass, rdbe, parent_id, parent_formula, transform, direction, category, status, PubMed_Count, accession, FirstBlock) %>%
    filter(TRUE)
  
  initial_rule = Related_files$empirical_rules %>%
    filter(category %in% c("Biotransform", "Adduct"))
  
  # Expanding known metabolites 
  {
    lib_known_id = LibrarySet %>%
      filter(origin == "known_library") %>%
      pull(library_id)
    
    lib_known = seed_library %>%
      filter(node_id %in% lib_known_id) %>%
      filter(category == "Metabolite")
    
    if(nrow(lib_known) > 0){
      rule_1 = initial_rule %>% filter(category == "Biotransform") %>% filter(direction %in% c(0,1))
      rule_2 = initial_rule %>% filter(category == "Biotransform") %>% filter(direction %in% c(0,-1))
      
      lib_known_1 = expand_library(lib_known, rule_1, direction = 1, category = "Metabolite")
      lib_known_2 = expand_library(lib_known, rule_2, direction = -1, category = "Metabolite")
      lib_met = bind_rows(lib_known_1, lib_known_2) %>%
        filter(!grepl("-|NA", formula)) 
    } else {
      lib_met = lib_known
    }
    
    }
  
  # Expanding known adducts
  {
    lib_adduct = seed_library %>%
      filter(category == "Adduct") %>%
      filter(TRUE)
    
    rule_1 = initial_rule %>% filter(category == "Adduct") %>% filter(direction %in% c(0,1))
    rule_2 = initial_rule %>% filter(category == "Adduct") %>% filter(direction %in% c(0,-1))
    
    lib_adduct_1 = expand_library(lib_adduct, rule_1, direction = 1, category = "Adduct")
    lib_adduct_2 = expand_library(lib_adduct, rule_2, direction = -1, category = "Adduct")
    
    lib_adduct = bind_rows(lib_adduct_1, lib_adduct_2) %>%
      filter(!str_detect(formula, "((?<!H)-)|(-(?!1))"))
    # filter((!grepl("-|NA", formula)) | grepl("H-(1|2)", formula)) 
  }
  
  expanded_lib = bind_rows(lib_met, lib_adduct) %>%
    arrange(mass) %>%
    mutate(node_id = 1:nrow(.)+nrow(LibrarySet))
  
  return(expanded_lib)
}
 
## get_mz get mz from NodeSet ####
get_mz = function(x){
 if(! is.null(x$pos_info) & is.null(x$neg_info)){
    mass = x$pos_info$pos_raw$medMz
  }else if(! is.null(x$neg_info) & is.null(x$pos_info)){
    mass = x$neg_info$neg_raw$medMz
  }
  return(mass)
}

# Initilize_empty_structureset ####
Initilize_Empty_Structureset = function(NodeSet){
  
  node_mass = sapply(NodeSet, get_mz)
  sf = lapply(1:length(NodeSet), function(x) {
    data.frame(
      node_id = x,
      formula = "Unknown",
      mass = node_mass[x], 
      rdbe = 0, 
      category = "Unknown", 
      parent_id = x, 
      parent_formula = "Unknown",
      transform = "",
      direction = 0,
      steps = 0,
      status = "",
      PubMed_Count = 0,
      accession = "",
      FirstBlock = "",
      stringsAsFactors = FALSE
    )
  })
  return(sf)
}

### get_rt get rt from NodeSet ####
get_rt = function(x){
 if(! is.null(x$pos_info) & is.null(x$neg_info)){
    rt = x$pos_info$pos_raw$medRt
  }else if(! is.null(x$neg_info) & is.null(x$pos_info)){
    rt = x$neg_info$neg_raw$medRt
  }
  return(rt)
  
}

## match_library ####
match_library = function(lib, sf, record_ppm_tol, record_RT_tol, current_step, NodeSet){
  lib = lib %>%
    arrange(mass)
  lib_mass = lib$mass
  length_lib = length(lib_mass)
  node_mass = sapply(NodeSet, get_mz) %>% sort()
  node_RT = sapply(NodeSet, get_rt)[names(node_mass)]
  temp_id = as.numeric(names(node_mass)) # numeric
  
  i=i_min=i_max=1
  while(i <= length(node_mass)){
    mass_tol = node_mass[i]*record_ppm_tol
    # Move i_min to a position larger than lower threshold
    while(lib_mass[i_min] < node_mass[i] - mass_tol & i_min < length_lib){
      i_min = i_min + 1
    }
    # if i_min's position larger than upper threhold, then move up i
    if(lib_mass[i_min] > node_mass[i] + mass_tol){
      i = i+1
      next
    }
    # Arriving here, means i_min reaches maximum or/and i_min's position is larger than lower threhold
    i_max = i_min
    
    # Move i_max to a position that is below upper threhold
    while(lib_mass[i_max] - node_mass[i] < mass_tol & i_max < length_lib){
      i_max = i_max + 1
    }
    i_max = i_max - 1
    
    # if there is no overlap of between i_min above lower threhold and i_max below upper threhold
    # it means i_min is maximum and i_max is maximum - 1, then break
    if(i_min > i_max){break}
    
    # Otherwise, record
    candidate = lib[i_min:i_max,]
    if(record_RT_tol < 999){
      candidate_filter = abs(node_RT[i] - node_RT[as.character(candidate$parent_id)]) < record_RT_tol
      candidate = candidate[candidate_filter, ]
      
      if(nrow(candidate) == 0){
        i = i + 1
        next
      }
      # We may implement another filter to break a cycle a->b->a, but it is not necessary
      # Actually, keeping both propagation direction makes it easier to trace parents
    }
    
    candidate["node_id"] = temp_id[i]
    candidate["steps"] = current_step
    
    sf[[temp_id[i]]] = bind_rows(sf[[temp_id[i]]], candidate)
    
    i = i+1
  }
  return(sf)
}
# Match_library_structureset ####
Match_Library_Structureset = function(LibrarySet, 
                                      expanded_lib,
                                      StructureSet, 
                                      NodeSet, 
                                      ppm_tol = 5e-6){
  ## initialize
  seed_library = LibrarySet %>% 
    mutate(node_id = library_id,
           formula = formula,
           mass = mass,
           parent_id = library_id,
           parent_formula = formula,
           transform = "",
           direction = 1,
           rdbe = rdbe,
           category = category,
           status = status,
           PubMed_Count = PubMed_Count,
           accession = note,
           FirstBlock = FirstBlock
    ) %>%
    dplyr::select(node_id, formula, mass, rdbe, parent_id, parent_formula, transform, direction, category, status, PubMed_Count, accession, FirstBlock) %>%
    filter(TRUE)
  
  if(!is.null(expanded_lib)){
    seed_library = bind_rows(seed_library, expanded_lib)
  }
  sf = StructureSet
  sf = match_library(lib = seed_library, sf, 
                     record_ppm_tol = ppm_tol, 
                     record_RT_tol = Inf, 
                     current_step = 0, 
                     NodeSet)
  sf = lapply(sf, function(x){
    unknown = x[1, ,drop = FALSE]
    if(nrow(x)>1){
      known_match = x[2:nrow(x), ,drop = FALSE]
      known_match  = known_match %>% arrange(accession)
      x = bind_rows(unknown, known_match)
    }
    return(x)
  })
 
  return(sf)
  
}



## get_inten get intensity from NodeSet ####
get_inten = function(x){
  if(! is.null(x$pos_info) & is.null(x$neg_info)){
    inten = x$pos_info$pos_raw$inten
  }else if(! is.null(x$neg_info) & is.null(x$pos_info)){
    inten = x$neg_info$neg_raw$inten
  }
  return(inten)
}

## Peak_grouping ####
Peak_grouping = function(NodeSet, RT_cutoff = 0.2, inten_cutoff = 1e4)
{
  node_mass = sapply(NodeSet, get_mz)
  node_RT = sapply(NodeSet, get_rt) 
  node_inten = sapply(NodeSet, get_inten)
  temp_id = as.numeric(names(node_RT)) # numeric
  
  peak_group_ls = list()
  
  for(i in 1:length(node_RT)){
    if(node_inten[i] < log10(inten_cutoff)){next}
    RT_min = node_RT[i] - RT_cutoff
    RT_max = node_RT[i] + RT_cutoff
    partner_id = temp_id[which(node_RT <= RT_max & node_RT >= RT_min)]
    peak_group_ls[[length(peak_group_ls)+1]] = list(node1 = rep(temp_id[i], length(partner_id)),
                                                    node2 = partner_id)
    
  }
  
  peak_group = bind_rows(peak_group_ls) %>%
    mutate(mass1 = node_mass[node1],
           mass2 = node_mass[node2]) %>%
    mutate(RT1 = node_RT[node1],
           RT2 = node_RT[node2]) %>%
    mutate(inten1 = node_inten[node1],
           inten2 = node_inten[node2]) %>%
    filter(TRUE)
  
  return(peak_group)
}

## Ring_artifact_connection ####
Ring_artifact_connection = function(peak_group,
                                    ppm_range_lb = 50, ppm_range_ub = 1000, ring_fold = 50, inten_threshold = 1e6){
  
  ring_artifact = peak_group %>%
    filter(inten1 > log10(inten_threshold)) %>%
    mutate(mz_dif = mass2 - mass1) %>%
    mutate(ppm_mz_dif = mz_dif / mass2 * 1e6) %>%
    filter(abs(ppm_mz_dif) < ppm_range_ub & abs(ppm_mz_dif) > ppm_range_lb) %>%
    mutate(inten_ratio = inten2 - inten1) %>%
    filter(inten_ratio < log10(1/ring_fold))
  
  
  EdgeSet_ring_artifact = apply(ring_artifact, 1, function(x){
    list(
      node1 = as.numeric(x["node1"]),
      node2 = as.numeric(x["node2"]),
      category = "Ring_artifact",
      linktype = "Ring_artifact",
      direction = 1
    )
  })
  
  return(EdgeSet_ring_artifact)
}

## Oligomer_multicharge_connection ####
Oligomer_multicharge_connection = function(peak_group, ppm_tol = 10){
  oligomer = peak_group %>% 
    mutate(mz_ratio12 = mass1/mass2,
           mz_ratio21 = mass2/mass1) %>% 
    filter(mz_ratio12 > 1.5 | mz_ratio21 > 1.5) %>%
    mutate(mz_ratio12_dif = mz_ratio12 - round(mz_ratio12),
           mz_ratio21_dif = mz_ratio21 - round(mz_ratio21),
           ratio = round(pmax(mz_ratio12, mz_ratio21))) %>%
    filter(abs(mz_ratio12_dif) < (ratio * ppm_tol / 1e6) | abs(mz_ratio21_dif) < (ratio * ppm_tol / 1e6)) %>%
    mutate(direction = ifelse(mz_ratio12 > mz_ratio21, -1, 1))
  
  oligomer1 = oligomer %>%
    filter(direction == 1)
  oligomer2 = oligomer %>%
    filter(direction == -1) %>%
    dplyr::rename(node1 = node2, node2 = node1)
  oligomer = bind_rows(oligomer1, oligomer2) %>%
    distinct(node1, node2, .keep_all = TRUE)
  
  EdgeSet_oligomer = apply(oligomer, 1, function(x){
    list(
      node1 = as.numeric(x["node1"]),
      node2 = as.numeric(x["node2"]),
      category = "Oligomer",
      linktype = as.character(x["ratio"]),
      direction = 1
    )
  })
  
  EdgeSet_multicharge = apply(oligomer, 1, function(x){
    list(
      node1 = as.numeric(x["node1"]),
      node2 = as.numeric(x["node2"]),
      category = "Multicharge",
      linktype = as.character(x["ratio"]),
      direction = -1
    )
  })
  
  EdgeSet_oligomer_multicharge = bind_rows(EdgeSet_oligomer, EdgeSet_multicharge)
  
  return(EdgeSet_oligomer_multicharge)
}

## Multicharge_isotope_connection ####
Multicharge_isotope_connection = function(EdgeSet, EdgeSet_oligomer_multicharge){
  
  EdgeSet_natural_abundance = bind_rows(EdgeSet) %>%
    filter(category == "Natural_abundance")
  if(is.null(EdgeSet_oligomer_multicharge)|nrow(EdgeSet_oligomer_multicharge)==0){
    return (EdgeSet_natural_abundance[0,])
  }
  
  
  EdgeSet_multicharge = EdgeSet_oligomer_multicharge %>%
    filter(category == "Multicharge") %>%
    mutate(linktype = as.numeric(linktype)) %>%
    dplyr::select(node1, node2, linktype)
  
  edge_multicharge_isotope = EdgeSet_natural_abundance %>%
    filter(node1 %in% EdgeSet_multicharge$node2, node2 %in% EdgeSet_multicharge$node2) %>% 
    merge(EdgeSet_multicharge, by.x = "node1", by.y="node2", suffixes = c("",".parent")) %>%
    merge(EdgeSet_multicharge, by.x = "node2", by.y="node2", suffixes = c("",".isotope")) %>%
    filter(linktype.parent == linktype.isotope) %>%
    mutate(linktype = mapply(my_calculate_formula, linktype, linktype, 
                             sign = -(linktype.parent-1)/linktype.parent)) %>%
    mutate(linktype = as.character(linktype)) %>%
    mutate(node1 = node1.parent,
           node2 = node1.isotope,
           category = "Multicharge_isotope") %>%
    dplyr::select(colnames(EdgeSet_natural_abundance))
  
  return(edge_multicharge_isotope)
}
### Chrom_Hetero_Sim
Chrom_Hetero_Sim = function(chrom1, chrom2, chrom_link){
  rts1 = t(as.data.frame(chrom1$rtime_list))
  intens1 = t(as.data.frame(chrom1$intensity_list))
  rts_link = t(as.data.frame(chrom_link$rtime_list))
  intens_link = t(as.data.frame(chrom_link$intensity_list))
  rts2 = t(as.data.frame(chrom2$rtime_list))
  intens2 = t(as.data.frame(chrom2$intensity_list))
  res = c()
  for(i in 1:nrow(rts1)){
    inten1 = Get_Ints(rts2[i,], 301, rts1[i,] , intens1[i,])
    inten_link = Get_Ints(rts2[i,], 301, rts_link[i,] , intens_link[i,])
    inten_merge = inten1*inten_link
    cos = Cos_Calculate(inten_merge, intens2[i,])
    cos = ifelse(is.na(cos), 0, cos)
    res = c(res, cos)
  }
  return(median(res, na.rm = TRUE))

}
## Heterodimer_connection ####
Heterodimer_connection = function(Related_files, NodeSet, ppm_tol = 10, inten_threshold = 1e5){
  peak_group = Peak_grouping(NodeSet, RT_cutoff = 0.4, inten_cutoff = Related_files$global_parameter$instrument_parameter$edge_expand_inten_cutoff)
  peak_group_ls = peak_group %>%
    split(.$node2)
  
  node_inten = sapply(NodeSet, get_inten)
  
  hetero_dimer_ls = list()
 
  hetero_dimer_ls = Heterodimer_connection_core(peak_group_ls, ppm_tol)
  if(length(hetero_dimer_ls) == 0){
    return(NULL)
  }
  
  hetero_dimer_df1 = bind_rows(hetero_dimer_ls) %>% 
    mutate(inten1 = node_inten[node1], 
           inten2 = node_inten[node2]) %>%
    filter(inten1 > log10(inten_threshold)) %>% # retain only high intensity as node1
    filter(inten1 > inten2) %>% # Heterodimer parent has higher intensity
    filter(node1 != linktype) %>% # remove homo-dimer
    arrange(node1)
  hetero_dimer_df2 = hetero_dimer_df1 %>%
    mutate(temp = linktype,
           linktype = node1,
           node1 = temp) %>%
    dplyr::select(-temp)
  

  hetero_dimer_df = bind_rows(hetero_dimer_df1, hetero_dimer_df2) %>%
    distinct(node1, linktype, .keep_all =TRUE)
  mode = ifelse(Related_files$global_parameter$mode==1, "pos_info","neg_info")
  ms_mode = ifelse(Related_files$global_parameter$mode==1, "pos_ms","neg_ms")
  EdgeSet_heterodimer = apply(hetero_dimer_df, 1, function(x){
    list(
      node1 = as.numeric(x["node1"]),
      node2 = as.numeric(x["node2"]),
      category = "Heterodimer",
      linktype = as.character(x["linktype"]),
      direction = 1
    )
  }) %>% lapply(.,function(x){
    node1 = x$node1
    node2 = x$node2
    linktype = as.numeric(x$linktype)
    chrom1 = NodeSet[[node1]][[mode]][[ms_mode]]
    chrom2 = NodeSet[[node2]][[mode]][[ms_mode]]
    chrom_link = NodeSet[[linktype]][[mode]][[ms_mode]]
    score_chrom = Chrom_Hetero_Sim(chrom1, chrom2, chrom_link)
    if(score_chrom>=0.75){
      return(x)
    } else{
      return(NULL)
    }
  })
  
  EdgeSet_heterodimer = EdgeSet_heterodimer[!sapply(EdgeSet_heterodimer, is.null)]
  return(EdgeSet_heterodimer)
}

## Experiment_MS2_fragment_connection ####
Experiment_MS2_fragment_connection = function(Related_files, peak_group, NodeSet, ppm_tol = 10, inten_threshold = 1e5){
  # It requires the parent node has msr MS2, and build connection between nodes
  # See comparison to Library_MS2_fragment_connection
  node_MS2_valid = sapply(NodeSet, function(x){!is.null(x$MS2)})
  node_MS2 = (1:length(NodeSet))[node_MS2_valid]
  mz = sapply(NodeSet, get_mz)
  if(length(node_MS2)==0){
    return(NULL)
  }
  
  H_mass = 1.00782503224
  e_mass = 0.00054857990943
  ion_mode = Related_files$global_parameter$mode
 if(ion_mode==1){
    pos_count = 1:nrow(Mset$pos_raw_data)
    neg_count = NULL
    both_count =  NULL
  }else{
    neg_count = 1:nrow(Mset$neg_raw_data)
    pos_count = NULL
    both_count =  NULL
  }
  
  peak_group_MS2 = peak_group %>%
    filter(node2 %in% node_MS2 & inten2 > log10(inten_threshold),
           node1 != node2,
           mass1 < mass2) %>%
    split(.$node2)
  
  experiment_MS2_fragment_ls = list()
  for(i in 1:length(peak_group_MS2)){
    
    node2 = as.numeric(names(peak_group_MS2)[i])
    mz_parent = mz[node2] 
    
    MS2 = NodeSet[[node2]]$MS2
    if(is.null(MS2$pos_MS2)){
      pos_MS2_mz = NULL
      neg_MS2_mz = MS2$neg_MS2$spec[,1] - (H_mass-e_mass)*(-1)
    }
    if(is.null(MS2$neg_MS2)){
      pos_MS2_mz = MS2$pos_MS2$spec[,1] - (H_mass-e_mass)*1
      neg_MS2_mz = NULL
    }
    if(!is.null(MS2$pos_MS2) & !is.null(MS2$neg_MS2)){
      pos_MS2_mz = MS2$pos_MS2$spec[,1] - (H_mass-e_mass)*1
      neg_MS2_mz = MS2$neg_MS2$spec[,1] - (H_mass-e_mass)*(-1)
    }
    MS2_mz = c(pos_MS2_mz, neg_MS2_mz)
    
    temp_e = peak_group_MS2[[i]]
    node1_mz = temp_e$mass1 
    node1 = temp_e$node1 %>% lapply(.,function(x){
      if(x %in% pos_count){
        return(1)
      }else{
          return(-1)
        }
      
    }) %>% unlist()
    MS2_mode = c(rep(1, length(pos_MS2_mz)), rep(-1, length(neg_MS2_mz)))
    temp_matrix = abs(outer(node1_mz, MS2_mz, FUN = "-")) < (ppm_tol * mz_parent / 1e6)
    temp_mode = matrix(outer(node1, MS2_mode, FUN="*") %in% c(0,1), ncol=length(MS2_mode)) 
    # temp_index1 = which(abs(temp_matrix) < ppm_tol * mz_parent / 1e6, arr.ind = TRUE)
    # temp_index2 = which(temp_mode==1, arr.ind = TRUE)
    temp_index = which(temp_matrix&temp_mode, arr.ind = TRUE)
    temp_matrix_value = abs(outer(node1_mz, MS2_mz, FUN = "-")) 
    # print(i)
    # print(temp_index)
    
    if(dim(temp_index)[1]>0){
      temp_dif = temp_matrix_value[temp_index]
      temp_node_1 = temp_e$node1[temp_index[,1]]
      temp_node_2 = node2
      
      # temp_linktype = mz_parent - MS2_mz[temp_index[,2]]
      
      temp_df = data.frame(node1 = temp_node_1, linktype = "Experiment_MS2_fragment", 
                           node2 = temp_node_2, mass_dif = temp_dif)
      experiment_MS2_fragment_ls[[length(experiment_MS2_fragment_ls)+1]] = temp_df
    }
    
  }
  
  experiment_MS2_fragment_df = bind_rows(experiment_MS2_fragment_ls)
  
  EdgeSet_experiment_MS2_fragment = apply(experiment_MS2_fragment_df, 1, function(x){
    list(
      node1 = as.numeric(x["node1"]),
      node2 = as.numeric(x["node2"]),
      category = "Experiment_MS2_fragment",
      linktype = as.character(x["linktype"]),
      mass_dif = as.numeric(x["mass_dif"]),
      direction = -1
    )
  })
  
  return(EdgeSet_experiment_MS2_fragment)
}

## Library_MS2_fragment_connection ####
Library_MS2_fragment_connection = function(Mset, NodeSet, LibrarySet, peak_group, StructureSet,Related_files,
                                           inten_threshold = 1e5,
                                           ppm_tol = 10, abs_tol = 1e-4){
  # It only requires the parent node has potential database MS2
  
  StructureSet_df = bind_rows(StructureSet)
  pos_MS2_library = Related_files$pos_MS2_library
  neg_MS2_library = Related_files$neg_MS2_library
  MS2_library = c(pos_MS2_library, neg_MS2_library)
  pos_MS2_library_external_id = sapply(pos_MS2_library, "[[", 'external_id')
  neg_MS2_library_external_id = sapply(neg_MS2_library, "[[", 'external_id')
  MS2_library_external_id = sapply(MS2_library, "[[", 'external_id')
  mz = sapply(NodeSet, get_mz)
  H_mass = 1.00782503224
  e_mass = 0.00054857990943
  ion_mode = Related_files$global_parameter$mode
 if(ion_mode==1){
    pos_count = 1:nrow(Mset$pos_raw_data)
    neg_count = NULL
    both_count =  NULL
  }else{
    neg_count = 1:nrow(Mset$neg_raw_data)
    pos_count = NULL
    both_count =  NULL
  }
  
  StructureSet_df_MS2 = StructureSet_df %>%
    # filter(steps == 0) %>%
    filter(steps == 0, transform == "") %>%
    merge(LibrarySet %>%
            dplyr::select(library_id, note, mass), 
          by.x = "parent_id",
          by.y = "library_id") %>%
    mutate(contain_lib_MS2 = ((note %in% pos_MS2_library_external_id) & (node_id %in% c(pos_count, both_count)))|
             ((note %in% neg_MS2_library_external_id) & (node_id %in% c(both_count,neg_count)))) %>%
    filter(contain_lib_MS2)
  
  peak_group_MS2 = peak_group %>%
    filter(node2 %in% StructureSet_df_MS2$node_id & inten2 > log10(inten_threshold), #Only look at high-inten peaks
           inten2 > inten1, # Fragment should have lower inten comparaed to parent
           mass2 > mass1) %>%
    split(.$node2)
  
  Library_MS2_fragment = list()
  for(i in 1:length(peak_group_MS2)){
    id_node1 = peak_group_MS2[[i]]$node1
    mz_node1 = mz[id_node1]
    
    id_node2 = peak_group_MS2[[i]]$node2[1]
    StructureSet_df_MS2_filter = StructureSet_df_MS2 %>%
      filter(node_id == id_node2)
    
    node2 = as.numeric(names(peak_group_MS2)[i])
    mz_parent = mz[node2] 
    
    for(j in 1:nrow(StructureSet_df_MS2_filter)){
      match_MS2_IDs = which(StructureSet_df_MS2_filter$note[j] == MS2_library_external_id)
      if(length(match_MS2_IDs) == 0){next}
      match_MS2 = MS2_library[match_MS2_IDs] %>%
        lapply(., function(x){
          if(x[["polarity"]] == "positive"){
            x[["spectrum"]][,1] = x[["spectrum"]][,1] - (H_mass-e_mass)*1
          }else if(x[["polarity"]] == "negative"){
            x[["spectrum"]][,1] = x[["spectrum"]][,1] - (H_mass-e_mass)*(-1)
          }
          return(x)
        })
      node1 = id_node1 %>% lapply(.,function(x){
        if(x %in% pos_count){
          return(1)
        }else if(x %in% both_count){return(0)}
        else{return(-1)}
        
      }) %>% unlist()
      MS2_mode = lapply(match_MS2, function(x){
        if(x[["polarity"]] == "positive"){
          return(rep(1,nrow(x[["spectrum"]])))
        }else{
          return(rep(-1,nrow(x[["spectrum"]])))
        }
      }) %>% unlist()
      match_MS2 = bind_rows(lapply(match_MS2, "[[","spectrum"))
      if(nrow(match_MS2)==0){next}
      mz_tol = max(mz_node1 * ppm_tol * 1e-6, abs_tol)
      temp_matrix = abs(outer(mz_node1, match_MS2$mz, FUN = "-")) < mz_tol
      temp_mode = matrix(outer(node1, MS2_mode, FUN = "*") %in% c(0,1), ncol = length(MS2_mode))
      
      temp_index = which(temp_matrix & temp_mode, arr.ind = TRUE)
      
      formula_node1 = match_MS2$formula[temp_index[,2]]
      # Some MS2 library contains NULL formulas
      # Clean MS2 library to make sure each fragment has formula
      if(any(is.null(formula_node1))){next}
      
      length_formula_node1 = length(formula_node1)
      
      formula_node2 = StructureSet_df_MS2_filter$formula[j]
      
      temp_node1 = id_node1[temp_index[,1]]
      
      
      Library_MS2_fragment[[length(Library_MS2_fragment)+1]] = list(node1 = temp_node1,
                                                                    node2 = rep(id_node2,length_formula_node1),
                                                                    formula1 = formula_node1,
                                                                    formula2 = rep(formula_node2,length_formula_node1))
      
    }
    
  }
  
  
  fragment_df = bind_rows(Library_MS2_fragment) %>%
    # group_by(node1, node2) %>%
    # filter(n()>1) %>%
    distinct() %>% 
    # filter(!is.na(formula1)) %>%
    # mutate(rdbe = formula_rdbe(formula1)) %>%
    mutate(linktype = "Library_MS2_fragment",
           category = "Library_MS2_fragment",
           direction = -1)
}

# Expand_edgeset ####
Expand_edgeset = function(Mset,
                          Related_files,
                          NodeSet,
                          LibrarySet,
                          StructureSet,
                          EdgeSet,
                          RT_cutoff = 0.2, inten_cutoff = 1e4,
                          heterodimer_option = F,
                          types = c("ring_artifact",
                                    "oligomer_multicharge",
                                    "heterodimer",
                                    "experiment_MS2_fragment",
                                    "library_MS2_fragment"
                          ))
{
  ppm_tol = Related_files$global_parameter$instrument_parameter$ppm * 1e6
  peak_group = Peak_grouping(NodeSet, RT_cutoff = RT_cutoff, inten_cutoff = inten_cutoff)
  EdgeSet_ring_artifact = EdgeSet_oligomer_multicharge = EdgeSet_multicharge_isotope = 
    EdgeSet_heterodimer = EdgeSet_experiment_MS2_fragment = EdgeSet_library_MS2_fragment = NULL
  if("ring_artifact" %in% types & Related_files$global_parameter$MS_type == "Orbitrap"){
    EdgeSet_ring_artifact = Ring_artifact_connection(peak_group, 
                                                     ppm_range_lb = 0, ppm_range_ub = 1000, 
                                                     ring_fold = 50, inten_threshold = 1e6)
    print(paste("ring_artifact", length(EdgeSet_ring_artifact)))
  }
  if("oligomer_multicharge" %in% types){
    EdgeSet_oligomer_multicharge = Oligomer_multicharge_connection(peak_group, ppm_tol = ppm_tol)
    EdgeSet_multicharge_isotope = Multicharge_isotope_connection(EdgeSet, EdgeSet_oligomer_multicharge)
    
    print(paste("oligomer_multicharge", 
                nrow(EdgeSet_oligomer_multicharge)+nrow(EdgeSet_multicharge_isotope)))
  }
  if(heterodimer_option){
    EdgeSet_heterodimer = Heterodimer_connection(Related_files, NodeSet, ppm_tol = ppm_tol, inten_threshold = 100*inten_cutoff)
    print(paste("heterodimer", length(EdgeSet_heterodimer)))
  }
  if("experiment_MS2_fragment" %in% types){
    EdgeSet_experiment_MS2_fragment = Experiment_MS2_fragment_connection(Related_files, peak_group, NodeSet, ppm_tol = ppm_tol, 
                                                                         inten_threshold = Related_files$global_parameter$instrument_parameter$edge_expand_inten_cutoff)
    print(paste("experiment_MS2_fragment", length(EdgeSet_experiment_MS2_fragment)))
  }
  if("library_MS2_fragment" %in% types){
    EdgeSet_library_MS2_fragment = Library_MS2_fragment_connection(Mset, NodeSet, LibrarySet, peak_group, StructureSet, Related_files,
                                                                   inten_threshold = Related_files$global_parameter$instrument_parameter$edge_expand_inten_cutoff,
                                                                   ppm_tol = ppm_tol, abs_tol = 1e-4)
    print(paste("library_MS2_fragment", length(EdgeSet_library_MS2_fragment)))
  }
  
  EdgeSet_expand = bind_rows(EdgeSet_ring_artifact, 
                             EdgeSet_oligomer_multicharge,
                             EdgeSet_multicharge_isotope,
                             EdgeSet_heterodimer,
                             EdgeSet_experiment_MS2_fragment,
                             EdgeSet_library_MS2_fragment
  )
  # print(table(EdgeSet_expand$category))
  return(EdgeSet_expand)
}
# merge_edgeset ####
merge_edgeset = function(EdgeSet, ...){
  
  EdgeSet_all = bind_rows(EdgeSet, ...) %>%
    mutate(edge_id = 1:nrow(.))
  return(EdgeSet_all)
}

## propagate_ring_artifact ####
propagate_ring_artifact = function(new_nodes_df, sf, EdgeSet_ring_artifact, NodeSet, current_step){
  
  if(nrow(EdgeSet_ring_artifact) == 0){return (NULL)}
  node_mass = sapply(NodeSet, get_mz)
  
  node1 = EdgeSet_ring_artifact$node1
  node2 = EdgeSet_ring_artifact$node2
  node1_node2_mapping = bind_rows(EdgeSet_ring_artifact) %>%
    dplyr::select(-c("category", "linktype", "direction"))
  
  new_nodes_df_ring_artifact = new_nodes_df %>% 
    filter(node_id %in% node1) %>% 
    merge(node1_node2_mapping, by.x="node_id", by.y="node1") %>%
    mutate(parent_id = node_id,
           parent_formula = formula,
           category = "Ring_artifact",
           transform = "Ring_artifact", 
           direction = 1,
           steps = current_step,
           formula = paste0("Ring_artifact_", formula)) %>%
    mutate(node_id = node2, 
           mass = node_mass[node_id]) %>%
    dplyr::select(-node2)
  
  # for(i in unique(new_nodes_df_ring_artifact$node_id)){
  #   sf[[i]] = bind_rows(sf[[i]], new_nodes_df_ring_artifact[new_nodes_df_ring_artifact$node_id == i, ])
  # }
  
  return(new_nodes_df_ring_artifact)
}
## propagate_oligomer ####
propagate_oligomer = function(new_nodes_df, sf, EdgeSet_oligomer, NodeSet, current_step){
  
  node_mass = sapply(NodeSet, get_mz)
  
  node1_node2_mapping = bind_rows(EdgeSet_oligomer) %>%
    dplyr::select(node1, node2, linktype, edge_id) %>%
    # dplyr::select(-c("category","direction")) %>%
    mutate(linktype = as.numeric(linktype))
  
  new_nodes_df1 = new_nodes_df %>% 
    filter(node_id %in% node1_node2_mapping$node1) %>% 
    dplyr::select(-c("category")) %>%
    merge(node1_node2_mapping, by.x="node_id", by.y="node1") %>%
    mutate(parent_id = node_id,
           parent_formula = formula,
           category = "Oligomer",
           direction = 1,
           steps = current_step,
           transform = as.character(linktype)) %>%
    mutate(node_id = node2, 
           mass = mass * linktype,
           formula = as.character(mapply(my_calculate_formula, formula, formula, sign = linktype - 1)),
           rdbe = rdbe * linktype) %>%
    dplyr::select(-c("node2","linktype"))
  
  new_nodes_df_oligomer = new_nodes_df1
  
  return(new_nodes_df_oligomer)
}
## propagate_multicharge ####
propagate_multicharge = function(new_nodes_df, sf, EdgeSet_oligomer, NodeSet, current_step){
  
  node_mass = sapply(NodeSet, get_mz)
  
  node1_node2_mapping = bind_rows(EdgeSet_oligomer) %>%
    dplyr::select(node1, node2, linktype, edge_id) %>%
    # dplyr::select(-c("category","direction")) %>%
    mutate(linktype = as.numeric(linktype))
  
  new_nodes_df2 = new_nodes_df %>% 
    filter(node_id %in% node1_node2_mapping$node2) %>% 
    dplyr::select(-c("category")) %>%
    merge(node1_node2_mapping, by.x="node_id", by.y="node2") %>%
    mutate(parent_id = node_id,
           parent_formula = formula,
           category = "Multicharge",
           direction = -1,
           steps = current_step,
           transform = as.character(linktype)) %>%
    mutate(node_id = node1, 
           mass = mass / linktype,
           formula = as.character(mapply(my_calculate_formula, formula, formula, sign = -(linktype-1)/linktype)),
           rdbe = rdbe / linktype) %>%
    dplyr::select(-c("node1","linktype"))
  
  new_nodes_df_oligomer = new_nodes_df2
  
  return(new_nodes_df_oligomer)
}

## propagate_heterodimer ####
propagate_heterodimer = function(new_nodes_df, sf, EdgeSet_heterodimer, NodeSet, 
                                 current_step, propagation_ppm_threshold, temp_propagation_category){
  # heterodimer = propagate_heterodimer(temp_new_nodes_df, sf, 
  #                                     temp_edgeset, NodeSet, current_step, 
  #                                     propagation_ppm_threshold, temp_propagation_category)
  node_mass = sapply(NodeSet, get_mz)
  node1_node2_mapping = bind_rows(EdgeSet_heterodimer) %>%
    dplyr::select(-c("category", "direction")) %>%
    dplyr::select(node1, node2, linktype, edge_id) %>%
    filter(TRUE)
  
  new_nodes_df_heterodimer = new_nodes_df %>% 
    filter(node_id %in% node1_node2_mapping$node1) %>% 
    merge(node1_node2_mapping, by.x="node_id", by.y="node1") %>%
    mutate(parent_id = node_id,
           parent_formula = formula,
           category = "Heterodimer",
           direction = 1,
           steps = current_step,
           transform = linktype) %>%
    dplyr::select(-c("linktype")) 
  
  if(nrow(new_nodes_df_heterodimer) == 0){
    return(NULL)
  }
  
  heterodimer_ls = list()

  # for(i in 1:nrow(new_nodes_df_heterodimer)){
  #   temp = new_nodes_df_heterodimer[i,]
  #   transform = sf[[as.numeric(temp$transform)]] %>%
  #     arrange(steps) %>%
  #     distinct(formula, .keep_all = TRUE) %>%
  #     # limits the choice of heterodimer partner
  #     filter(steps <= 0.01 | (steps >= 1  & steps <= 1.01)) %>%
  #     filter(steps < current_step) %>%
  #     filter(category %in% temp_propagation_category) %>%
  #     filter(abs(node_mass[as.character(node_id)] - mass) < propagation_ppm_threshold * mass) # propagate from accurate formulas
  # 
  #   if(nrow(transform) == 0){next}
  # 
  #   heterodimer_ls[[length(heterodimer_ls)+1]] = list(parent_id = rep(temp$parent_id, length(transform$formula)),
  #                                                     transform = rep(temp$transform, length(transform$formula)),
  #                                                     node2 = rep(temp$node2, length(transform$formula)),
  #                                                     parent_formula = rep(temp$parent_formula, length(transform$formula)),
  #                                                     formula = as.character(my_calculate_formula(temp$formula, transform$formula)),
  #                                                     rdbe = temp$rdbe + transform$rdbe,
  #                                                     mass = temp$mass + transform$mass)
  # 
  # }
  
  new_nodes_df_heterodimer[,8] = as.integer(new_nodes_df_heterodimer[,8])
  heterodimer_ls = propagate_heterodimer_core(new_nodes_df_heterodimer, sf, propagation_category = temp_propagation_category,
                                              node_mass = node_mass, ppm_threshold = propagation_ppm_threshold)
  heterodimer_ls = lapply(heterodimer_ls, function(x){
    x[["formula"]] = as.character(x[["formula"]])
    return(x)
  })
  
  heterodimer = bind_rows(heterodimer_ls) %>%
    merge(new_nodes_df_heterodimer %>% dplyr::select(-c("formula", "mass", "rdbe")), all.x = TRUE) %>%
    mutate(node_id = node2) %>%
    dplyr::select(-node2) %>% 
    distinct(formula, node_id, transform, parent_formula, parent_id, .keep_all = TRUE)
 
  return(heterodimer)
}
## propagate_experiment_MS2_fragment ####
propagate_experiment_MS2_fragment = function(new_nodes_df, sf, 
                                             EdgeSet_experiment_MS2_fragment, 
                                             NodeSet, current_step){
  
  
  if(nrow(EdgeSet_experiment_MS2_fragment)==0){
    return(NULL)
  }
  node_mass = sapply(NodeSet, get_mz)
  
  node1 = EdgeSet_experiment_MS2_fragment$node1
  node2 = EdgeSet_experiment_MS2_fragment$node2
  node1_node2_mapping = bind_rows(EdgeSet_experiment_MS2_fragment) %>%
    dplyr::select(-c("category", "linktype", "direction"))
  
  new_nodes_df_MS2_fragment = new_nodes_df %>% 
    filter(node_id %in% node2) %>% 
    merge(node1_node2_mapping, by.x="node_id", by.y="node2") %>%
    mutate(parent_id = node_id,
           parent_formula = formula,
           category = "Experiment_MS2_fragment",
           transform = "Experiment_MS2_fragment", 
           direction = -1,
           steps = current_step,
           formula = paste0("MS2_fragment_", formula)) %>%
    mutate(node_id = node1, 
           mass = node_mass[node_id] - mass_dif) %>%
    dplyr::select(-node1, -mass_dif)
  
  return(new_nodes_df_MS2_fragment)
}

## propagate_library_MS2_fragment ####
propagate_library_MS2_fragment = function(new_nodes_df, sf, EdgeSet_library_MS2_fragment, 
                                          NodeSet, current_step){
  
  # (temp_new_nodes_df, sf, temp_edgeset, NodeSet, current_step)
  # new_nodes_df = temp_new_nodes_df
  # EdgeSet_library_MS2_fragment = temp_edgeset
  if(nrow(EdgeSet_library_MS2_fragment) == 0){return (NULL)}
  
  node_mass = sapply(NodeSet, get_mz)
  
  node1 = EdgeSet_library_MS2_fragment$node1
  node2 = EdgeSet_library_MS2_fragment$node2
  node1_node2_mapping = EdgeSet_library_MS2_fragment %>%
    dplyr::select(node1, node2, formula1, formula2, edge_id) %>%
    # dplyr::select(-c("category", "linktype", "direction", "rdbe")) %>%
    dplyr::rename(formula = formula2, node_id = node2)
  
  new_nodes_df_library_MS2_fragment = new_nodes_df %>% 
    filter(node_id %in% node2) %>%
    merge(node1_node2_mapping)
  if(nrow(new_nodes_df_library_MS2_fragment) == 0){
    return(NULL)
  }
 
  
  new_nodes_df_library_MS2_fragment = new_nodes_df_library_MS2_fragment%>% 
    # merge(node1_node2_mapping) %>%
    mutate(parent_id = node_id,
           parent_formula = formula,
           category = "Library_MS2_fragment",
           transform = "Library_MS2_fragment", 
           direction = -1,
           steps = current_step,
           formula = formula1) %>%
    mutate(node_id = node1, 
           mass = formula_mz(formula)) %>%
    dplyr::select(-node1, -formula1,)
  
  return(new_nodes_df_library_MS2_fragment)
}
# Propagate_structureset ####
Propagate_structureset = function(Mset, 
                                  NodeSet,
                                  StructureSet,
                                  EdgeSet_all,
                                  biotransform_step = 2,
                                  artifact_step = 3,
                                  propagation_ppm_threshold = 5e-6,
                                  propagation_abs_threshold = 2e-4,
                                  record_RT_tol = 0.1,
                                  record_ppm_tol = 5e-6)
{
  
  
  sf = StructureSet
  empirical_rules = Related_files$empirical_rules
  propagation_rule = Related_files$global_parameter$propagation_rule
  node_mass = sapply(NodeSet, get_mz) %>% sort()
  node_RT = sapply(NodeSet, get_rt)[names(node_mass)]
  temp_id = as.numeric(names(node_mass)) # numeric
  
  
  
  ## Expansion 
  timer = Sys.time()
  
  # Handle biotransform
  step_count = 0
  while(step_count < biotransform_step){
    
    
    
    {
      all_nodes_df = bind_rows(sf)
      new_nodes_df = all_nodes_df %>%
        filter(category == "Metabolite") %>% # only metabolites go to biotransformation, also garantee it is not filtered.
        distinct(node_id, formula, .keep_all=TRUE) %>% # This garantee only new formulas will go to next propagation
        filter(steps == step_count) %>% # only formulas generated from the step go to next propagation
        filter(rdbe > -1) %>% # filter out formula has ring and double bind equal or less than -1
        filter(abs(node_mass[as.character(node_id)] - mass) < pmax(propagation_ppm_threshold * mass, 
                                                                   propagation_abs_threshold)) # propagate from accurate formulas
      
      step_count = step_count+1
      current_step = step_count
      
      if(nrow(new_nodes_df)==0){break}
      
      rule_1 = empirical_rules %>% filter(category == "Biotransform") %>% filter(direction %in% c(0,1)) #direction whether + or - atom diff
      rule_2 = empirical_rules %>% filter(category == "Biotransform") %>% filter(direction %in% c(0,-1))
      
      new_nodes_df = new_nodes_df %>%
        mutate(parent_id = node_id)
      lib_1 = expand_library(new_nodes_df, rule_1, direction = 1, category = "Metabolite")
      lib_2 = expand_library(new_nodes_df, rule_2, direction = -1, category = "Metabolite")
      
      lib_met = bind_rows(lib_1, lib_2) %>%
        filter(!grepl("-|NA", formula)) %>% 
        arrange(mass)
      
      sf = match_library(lib_met,
                         sf,
                         record_ppm_tol,
                         record_RT_tol=Inf,
                         current_step,
                         NodeSet)
      
      print(paste("Step",current_step,"elapsed="))
      print(Sys.time()-timer)
    }
    
  }
  
  # Handle artifacts
  step_count = 0
  for(step_count in 0:(biotransform_step-1)){
    sub_step = 0
    while(sub_step < 0.01 * artifact_step ){
      
      node_step = step_count + sub_step
      sub_step = sub_step+0.01
      current_step = step_count + sub_step
      
      
      all_nodes_df = bind_rows(sf)
      new_nodes_df = all_nodes_df %>%
        distinct(node_id, formula, category, .keep_all = TRUE) %>%
        filter(category != "Unknown") %>%
        filter(abs(node_mass[as.character(node_id)] - mass) <
                 pmax(propagation_ppm_threshold * mass, propagation_abs_threshold)) # propagate from accurate formulas
      
      if(nrow(new_nodes_df)==0){break}
      
      lib_artifact_ls = list()
      # to Adduct ##
      {
        temp_new_nodes_df = new_nodes_df %>%
          filter(category %in% propagation_rule[["Adduct"]]) %>%
          arrange(steps) %>%
          distinct(node_id, formula, .keep_all = TRUE) %>%
          filter(steps==node_step)
        
        if(nrow(temp_new_nodes_df) != 0){
          temp_rule = empirical_rules %>% filter(category == "Adduct")
          temp_propagation = expand_library(temp_new_nodes_df, temp_rule, direction = 1, category = "Adduct")
          lib_artifact_ls[[length(lib_artifact_ls)+1]] = temp_propagation
        }
      }
      
      # to Fragment ##
      {
        temp_new_nodes_df = new_nodes_df %>% 
          filter(category %in% propagation_rule[["Fragment"]]) %>%
          arrange(steps) %>%
          distinct(node_id, formula, .keep_all = TRUE) %>%
          filter(steps==node_step) 
        
        if(nrow(temp_new_nodes_df) != 0){
          temp_rule = empirical_rules %>% filter(category == "Fragment")
          temp_propagation = expand_library(temp_new_nodes_df, temp_rule, direction = -1, category = "Fragment")
          lib_artifact_ls[[length(lib_artifact_ls)+1]] = temp_propagation
        }
        
      }
      
      # to Natural_abundance ##
      {
        temp_new_nodes_df = new_nodes_df %>%
          filter(category %in% propagation_rule[["Natural_abundance"]]) %>%
          arrange(steps) %>%
          distinct(node_id, formula, .keep_all = TRUE) %>%
          filter(steps==node_step) 
        
        if(nrow(temp_new_nodes_df) != 0){
          
          temp_rule = empirical_rules %>% filter(category == "Natural_abundance") %>% filter(direction == 1)
          temp_propagation = expand_library(temp_new_nodes_df, temp_rule, direction = 1, category = "Natural_abundance")
          lib_artifact_ls[[length(lib_artifact_ls)+1]] = temp_propagation
          
          # # [10]B has a direction of -1
          # temp_rule = empirical_rules %>% filter(category == "Natural_abundance") %>% filter(direction == -1)
          # temp_propagation = expand_library(temp_new_nodes_df, temp_rule, direction = -1, category = "Natural_abundance")
          # lib_artifact_ls[[length(lib_artifact_ls)+1]] = temp_propagation
        }
      }
      
      # to Radical ##
      {
        temp_new_nodes_df = new_nodes_df %>%
          filter(category %in% propagation_rule[["Radical"]]) %>%
          arrange(steps) %>%
          distinct(node_id, formula, .keep_all = TRUE) %>%
          filter(steps==node_step)
        
        if(nrow(temp_new_nodes_df) != 0){
          temp_rule = empirical_rules %>% filter(category == "Radical")
          
          temp_propagation = expand_library(temp_new_nodes_df, temp_rule, direction = -1, category = "Radical")
          
          lib_artifact_ls[[length(lib_artifact_ls)+1]] = temp_propagation
        }
      }
      
      # to Heterodimer ##
      {
        temp_new_nodes_df = new_nodes_df %>%
          filter(category %in% propagation_rule[["Heterodimer"]]) %>%
          arrange(steps) %>%
          distinct(node_id, formula, .keep_all = TRUE) %>%
          filter(steps==node_step)
        
        temp_edgeset = EdgeSet_all %>%
          filter(category == "Heterodimer")
        
        temp_propagation_category = propagation_rule[["Heterodimer"]]
        heterodimer = propagate_heterodimer(temp_new_nodes_df, sf,
                                            temp_edgeset, NodeSet, current_step,
                                            propagation_ppm_threshold, temp_propagation_category)
        if(!is.null(heterodimer)){
          heterodimer$transform = as.character(heterodimer$transform)
        }
      }
      
      # to Oligomer ##
      {
        temp_new_nodes_df = new_nodes_df %>%
          filter(category %in% propagation_rule[["Oligomer"]]) %>%
          arrange(steps) %>%
          distinct(node_id, formula, .keep_all = TRUE) %>%
          filter(steps==node_step) 
        
        temp_edgeset = EdgeSet_all %>%
          filter(category == "Oligomer")
        
        oligomer = propagate_oligomer(temp_new_nodes_df, sf, temp_edgeset, NodeSet, current_step)
      }
      
      # to Multicharge ##
      {
        temp_new_nodes_df = new_nodes_df %>%
          filter(category %in% propagation_rule[["Multicharge"]]) %>%
          arrange(steps) %>%
          distinct(node_id, formula, .keep_all = TRUE) %>%
          filter(steps==node_step)
        
        temp_edgeset = EdgeSet_all %>%
          filter(category == "Multicharge")
        
        multicharge = propagate_multicharge(temp_new_nodes_df, sf, temp_edgeset, NodeSet, current_step)
      }
      
      # to Library_MS2_fragment ##
      {
        temp_new_nodes_df = new_nodes_df %>%
          filter(category %in% propagation_rule[["Library_MS2_fragment"]]) %>%
          arrange(steps) %>%
          distinct(node_id, formula, .keep_all = TRUE) %>%
          filter(steps==node_step) 
        
        temp_edgeset = EdgeSet_all %>%
          filter(category == "Library_MS2_fragment")
        
        library_MS2_fragment = propagate_library_MS2_fragment(temp_new_nodes_df, sf, temp_edgeset, NodeSet, current_step)
      }
      
      # to Experiment_MS2_fragment ##
      {
        temp_new_nodes_df = new_nodes_df %>%
          filter(category %in% propagation_rule[["Experiment_MS2_fragment"]]) %>%
          arrange(steps) %>%
          distinct(node_id, formula, .keep_all = TRUE) %>%
          filter(steps==node_step) 
        
        temp_edgeset = EdgeSet_all %>%
          filter(category == "Experiment_MS2_fragment")
        
        experiment_MS2_fragment = propagate_experiment_MS2_fragment(temp_new_nodes_df, sf, temp_edgeset, NodeSet, current_step)
      }
      
      # to Ring_artifact ##
      {
        temp_new_nodes_df = new_nodes_df %>%
          filter(category %in% propagation_rule[["Ring_artifact"]]) %>%
          arrange(steps) %>%
          distinct(node_id, formula, .keep_all = TRUE) %>%
          filter(steps==node_step) 
        
        temp_edgeset = EdgeSet_all %>%
          filter(category == "Ring_artifact")
        
        ring_artifact = propagate_ring_artifact(temp_new_nodes_df, sf, temp_edgeset, NodeSet, current_step)
      }
      
      
      # Handle common edge list, adding lib_artifact_ls to library ##
      {
        lib_artifact = bind_rows(lib_artifact_ls) %>%
          filter(!str_detect(formula, "((?<!H)-)|(-(?!1))")) %>% # any "-" not preceded by H, or not followed by 1 is removed
          # mutate(RT = node_RT[as.character(parent_id)]) %>%
          arrange(mass)
        sf = match_library(lib_artifact,
                           sf,
                           record_ppm_tol,
                           record_RT_tol,
                           current_step,
                           NodeSet)
        
      }
      
      # Handle expanded list #
      {
        sf_add_mass_filter = bind_rows(oligomer, heterodimer, multicharge, library_MS2_fragment) %>%
          mutate(msr_mass = node_mass[as.character(node_id)],
                 mass_dif = abs(msr_mass-mass)) %>%
          mutate(mass_retain = mass_dif / mass <= record_ppm_tol) %>%
          filter(mass_retain) %>%
          dplyr::select(colnames(sf[[1]]))
        
        sf_add = bind_rows(sf_add_mass_filter, ring_artifact, experiment_MS2_fragment) %>%
          dplyr::select(colnames(sf[[1]]))
        
        for(i in unique(sf_add$node_id)){
          sf[[i]] = bind_rows(sf[[i]], sf_add[sf_add$node_id == i, ])
        }
      }
      print(paste("Step", current_step, "elapsed="))
      print((Sys.time()-timer))
    }
    
  }
  
  return(sf)
}


## Score_structureset_database_origin ####
Score_structureset_database_origin = function(StructureSet_df, LibrarySet,
                                              manual_match = 1){
  # Score HMDB and known adduct match
  
  manual_library_id = LibrarySet %>%
    filter(origin == "manual_library") %>%
    pull(library_id)
  
  structureset_database_origin = StructureSet_df %>%
    mutate(score_database_origin = case_when(
      category == "Unknown" ~ 0,
      parent_id %in% manual_library_id & transform == "" ~ manual_match,
      steps == 0 & transform == ""  & status=="quantified" ~ 0.4 + log10(PubMed_Count+1) * 0.05,
      steps == 0 & transform == ""  & status=="detected" ~ 0.3 + log10(PubMed_Count + 1) * 0.05,
      steps == 0 & transform == ""  & !(status %in% c("quantified","detected")) ~ log10(PubMed_Count + 1) * 0.05,
      TRUE ~ 0 # Everything else
    )) %>%
    dplyr::select(struct_set_id, score_database_origin)
  
}
## Score_structureset_mz ####
Score_structureset_mz = function(StructureSet_df, NodeSet, 
                                 score_ppm_error_rate = -0.5){
  
  node_mass = sapply(NodeSet, get_mz)
  
  structureset_mz = StructureSet_df %>%
    mutate(msr_mass = node_mass[node_id],
           ppm_error = (mass - msr_mass) / mass * 1e6) %>%
    mutate(score_mz = abs(ppm_error) * score_ppm_error_rate) %>%
    # mutate(score_mz = dgamma(abs(ppm_error), 1, scale = mass_dist_gamma_rate) * mass_dist_gamma_rate) %>%
    # mutate(score_mz = log10(score_mz)) %>%
    dplyr::select(struct_set_id, score_mz)
  
}

## Score_structureset_RT ####
Score_structureset_RT = function(StructureSet_df, NodeSet, LibrarySet, 
                                 rt_match = 1, known_rt_tol = 0.5){
  if(is.null(LibrarySet$rt)){
    return(NULL)
  }
  node_RT = sapply(NodeSet, get_rt)
  library_RT = LibrarySet$rt
  names(library_RT) = LibrarySet$library_id
  
  structureset_RT = StructureSet_df %>%
    mutate(node_rt = node_RT[node_id], 
           library_rt = library_RT[as.character(parent_id)]) %>%
    mutate(score_known_rt = ifelse(abs(node_rt - library_rt) < known_rt_tol & 
                                     steps == 0 & transform == "", 
                                   rt_match-abs(node_rt - library_rt), 0)) %>%
    dplyr::select(struct_set_id, score_known_rt)
}

### mergeMzIntensity ####
mergeMzIntensity = function(spec1_df, spec2_df, ppmTol = 10E-6, absTol = 1e-3){
  
  
  # absTol = 1e-3
  # ppmTol = 10E-6
  # true_merge = merge(as.data.frame(spec1_df),spec2_df, by="mz", all=TRUE)
  # Assuming mz came in sorted
  mz1 = spec1_df[,1]
  mz2 = spec2_df[,1]
  mzs = c(mz1, mz2)
  mzs = sort(mzs)
  inten1 = spec1_df[,2]
  inten2 = spec2_df[,2]
  keeps = c(1)
  count = 1
  MS_group = rep(1,(length(mzs)))
  for(i in 2:length(mzs)){
    # Determine when to treat two mz as same, and when they are different
    # Currently, when two masses are at least ppmTol different, and absTol different, they are considered different.
    if(mzs[i]-mzs[i-1]>mzs[i-1]*ppmTol & mzs[i]-mzs[i-1]>absTol){
      
      count = count+1
      keeps=c(keeps,i)
    }
    MS_group[i]=count
  }
  
  mz_result = mzs[keeps]
  count_mzs = 2
  count_keep = 1
  i = j = 1
  intens_mat = matrix(0, ncol=3, nrow = length(mz_result))
  while(count_mzs <= length(mzs)){
    if(MS_group[count_mzs] == MS_group[count_mzs-1]){
      intens_mat[count_keep,2] = inten1[i]
      intens_mat[count_keep,3] = inten2[j]
      i = i+1
      j = j+1
      count_mzs = count_mzs+2
    } else if(is.na(mz1[i])){
      intens_mat[count_keep,3] = inten2[j]
      count_mzs = count_mzs+1
      j=j+1
    } else if(mz1[i] == mzs[count_mzs-1]){
      intens_mat[count_keep,2] = inten1[i]
      i=i+1
      count_mzs = count_mzs+1
    } else {
      intens_mat[count_keep,3] = inten2[j]
      j=j+1
      count_mzs = count_mzs+1
    }
    count_keep = count_keep+1
  }
  # Handle boundary
  if(!is.na(mz1[i])){intens_mat[count_keep,2] = inten1[i]}
  if(!is.na(mz2[j])){intens_mat[count_keep,3] = inten2[j]}
  
  intens_mat[,1] = mz_result
  return(intens_mat)
}
mergeMzIntensity_backup = function(spec1_df, spec2_df, ppmTol = 10E-6, absTol = 1e-3){
  # if(identical(spec1_df, spec2_df))
  
  colnames(spec1_df)[2] = "inten1"
  colnames(spec2_df)[2] = "inten2"
  spec_df = merge(spec1_df, spec2_df, all=TRUE, by = "mz")
  
  MS_group = rep(1,(nrow(spec_df)))
  mzs = spec_df$mz
  if(any(is.na(mzs))){
    return(NULL)
  }
  
  keeps = c(1)
  count = 1
  for(i in 2:length(mzs)){
    if(mzs[i]-mzs[i-1]>mzs[i-1]*ppmTol & mzs[i]-mzs[i-1]>absTol){
      count = count+1
      keeps=c(keeps,i)
    }
    MS_group[i]=count
  }
  spec_df[is.na(spec_df)]=0
  
  k_max=k_min=1
  while (k_max <= length(MS_group)){
    k_min = k_max
    while (MS_group[k_min] == MS_group[k_max]){
      k_max = k_max+1
      if(k_max > length(MS_group)){break}
    }
    # spec_df$mz[k_min]=max(spec_df$mz[k_min:(k_max-1)], na.rm = TRUE)
    spec_df$inten1[k_min]=max(spec_df$inten1[k_min:(k_max-1)], na.rm = TRUE)
    spec_df$inten2[k_min]=max(spec_df$inten2[k_min:(k_max-1)], na.rm = TRUE)
  }
  spec_df = spec_df[keeps,]
  return(spec_df)
}
### Score_merge_MS2 ####
Score_merge_MS2 = function(spec_merge_df, mz_parent = 0){
  
  # A more realistic evaluation of two MS2 spectra simialrity is needed.
  # idea1 discount the weight of parent peak matched
  # idea2 add weights to the number of shared fragment
  # idea3 make sqrt of the intensity - 
  # which adds more weight to low-signal shared peaks
  
  # Implementing idea 3
  # a = a/max(a)
  # b = b/max(b)
  # a = sqrt(a)
  # b = sqrt(b)
  # Needs to use experimental data to validate the score cutoff for match and similarity
  # Implementing idea 1
  if(is.null(spec_merge_df)){
    return(0)
  }
  
  a = spec_merge_df[,2]
  b = spec_merge_df[,3]
  
  # 5e-4 and divide by 5 is arbitrary, needs experimental data to validate.
  # If shared parent, then discount the parent intensity during scoring
  mz_parent_position = which(abs(spec_merge_df[,1] - mz_parent) < 5e-4)
  if(length(mz_parent_position) != 0){
    a[mz_parent_position] = a[mz_parent_position]/5
    b[mz_parent_position] = b[mz_parent_position]/5
  }
  
  spec_score = a%*%b / (sqrt(a%*%a) * sqrt(b%*%b))
  
  return(as.numeric(spec_score))
}

#### cal_spec_entropy ####
cal_spec_entropy = function(inten, weight = 0.25){
  # smaller weight means mz match is more important
  # default weight comes from paper 
  inten = inten/sum(inten)
  S_raw = -sum(ifelse(inten > 0, inten * log(inten), 0))
  
  if(S_raw < 3 & weight <= 0.25){
    w = weight + S_raw*weight
    inten = inten^w
    inten = inten/sum(inten)
    S = -sum(ifelse(inten > 0, inten * log(inten), 0))
  } else {
    S = S_raw
  }
  return(list(S = S, inten = inten))
}

### Score_merge_MS2_entropy ####
Score_merge_MS2_entropy = function(spec_merge_df, weight = 0.25){
  # mz1 = mz; inten1=inten
  # mz2 = mz_sim; inten2 = inten_sim
  # weight = 1
  spec_merge_df[is.na(spec_merge_df)] = 0
  inten1 = spec_merge_df[,2]
  inten2 = spec_merge_df[,3]
  inten1 = inten1/max(inten1) * 100
  inten2 = inten2/max(inten2) * 100
  inten1 = inten1/sum(inten1)
  inten2 = inten2/sum(inten2)
  
  S1 = cal_spec_entropy(inten1, weight)
  S2 = cal_spec_entropy(inten2, weight)
  
  inten_m = S1$inten + S2$inten
  inten_m = inten_m/2
  Sm = -sum(ifelse(inten_m > 0, inten_m * log(inten_m), 0))
  
  spec_similarity = 1-(Sm*2-S2$S-S1$S)/log(4)
  spec_similarity = ifelse(spec_similarity<0,0,spec_similarity)
  return(spec_similarity)
  
}



## Score_structureset_MS2 ####
Score_structureset_MS2 = function(Mset, StructureSet_df, NodeSet, Related_files, LibrarySet,
                                  # only spectra matching/similarity score > cutoff, then nodes get bonus
                                  MS2_match = 1, MS2_match_cutoff = 0.5, 
                                  MS2_similarity = 0.5, MS2_similarity_cutoff = 0.5 
){ 
  
  
  node_MS2_valid = sapply(NodeSet, function(x){!is.null(x$MS2)})
  if(!any(node_MS2_valid)){return(NULL)}
  node_MS2 = (1:length(NodeSet))[node_MS2_valid]
  pos_MS2_library = Related_files$pos_MS2_library
  neg_MS2_library = Related_files$neg_MS2_library
  MS2_library = c(pos_MS2_library, neg_MS2_library)
  pos_MS2_library_external_id = sapply(pos_MS2_library, "[[", 'external_id')
  neg_MS2_library_external_id = sapply(neg_MS2_library, "[[", 'external_id')
  MS2_library_external_id = sapply(MS2_library, "[[", 'external_id')
  ion_mode = Related_files$global_parameter$mode
 if(ion_mode==1){
    pos_count = 1:nrow(Mset$pos_raw_data)
    neg_count = NULL
    both_count =  NULL
  }else{
    neg_count = 1:nrow(Mset$neg_raw_data)
    pos_count = NULL
    both_count =  NULL
  }
  
  StructureSet_df_MS2 = StructureSet_df %>%
    filter(steps == 0) %>%
    # filter(steps == 0, transform == "") %>%
    merge(LibrarySet %>%
            dplyr::select(library_id, note, mass), 
          by.x = "parent_id",
          by.y = "library_id") %>%
    mutate(contain_msr_MS2 = node_id %in% node_MS2) %>%
    mutate(contain_lib_MS2 = ((note %in% pos_MS2_library_external_id) & (node_id %in% c(pos_count, both_count)))|
             ((note %in% neg_MS2_library_external_id) & (node_id %in% c(both_count, neg_count)))) %>%
    filter(contain_msr_MS2, contain_lib_MS2)
  
  fwd_dp = rev_dp = rep(0, nrow(StructureSet_df_MS2))
  
  for(i in 1:nrow(StructureSet_df_MS2)){
    
    H_mass = 1.007825032
    e_mass = 0.000548579
    
    node_id = StructureSet_df_MS2$node_id[i]
    if(!is.null(NodeSet[[node_id]]$MS2$pos_MS2) & is.null(NodeSet[[node_id]]$MS2$neg_MS2)){
      pos_MS2_msr = NodeSet[[node_id]]$MS2$pos_MS2$spec %>% mutate(mz = (mz- (H_mass-e_mass)*1))
      neg_MS2_msr = NULL
      if(nrow(pos_MS2_msr)==1){
        next
      }
      MS2_msr = list(pos_MS2_msr = pos_MS2_msr, neg_MS2_msr = neg_MS2_msr)
    }else if(is.null(NodeSet[[node_id]]$MS2$pos_MS2) & !is.null(NodeSet[[node_id]]$MS2$neg_MS2)){
      neg_MS2_msr = NodeSet[[node_id]]$MS2$neg_MS2$spec %>% mutate(mz = (mz- (H_mass-e_mass)*(-1)))
      pos_MS2_msr = NULL
      if(nrow(neg_MS2_msr)==1){
        next
      }
      MS2_msr = list(pos_MS2_msr = pos_MS2_msr, neg_MS2_msr = neg_MS2_msr)
    }else if(!is.null(NodeSet[[node_id]]$MS2$pos_MS2) & !is.null(NodeSet[[node_id]]$MS2$neg_MS2)){
      pos_MS2_msr = NodeSet[[node_id]]$MS2$pos_MS2$spec %>% mutate(mz = (mz- (H_mass-e_mass)*1))
      neg_MS2_msr = NodeSet[[node_id]]$MS2$neg_MS2$spec %>% mutate(mz = (mz- (H_mass-e_mass)*(-1)))
      if(nrow(pos_MS2_msr)==1 | nrow(neg_MS2_msr)==1){
        next
      }
      MS2_msr = list(pos_MS2_msr = pos_MS2_msr, neg_MS2_msr = neg_MS2_msr)
    }
    
    mz1 = StructureSet_df_MS2$mass.x[i] 
    
    
    
    HMDB_id = StructureSet_df_MS2$note[i]
    mz2 = StructureSet_df_MS2$mass.y[i] 
    lib_MS2_id = which(MS2_library_external_id == HMDB_id)
    MS2_lib= MS2_library[lib_MS2_id]
    MS2_lib = lapply(MS2_lib, function(x){
      if(x[["polarity"]] == "positive"){
        x[["spectrum"]][,1] = x[["spectrum"]][,1] - (H_mass-e_mass)*(1)
      }else{
        x[["spectrum"]][,1] = x[["spectrum"]][,1] - (H_mass-e_mass)*(-1)
      }
      return(x)
    })
    
    
    # skip if all MS2_lib spectra are one-row
    if(all(sapply(MS2_lib, function(x){nrow(x$spectrum)}) == 1)){
      next
    }
    
    if(abs(mz1 - mz2) < max(1e-3, mz2*(Related_files$global_parameter$instrument_parameter$ms1_ms2_match_ppm))){
      temp_mz_parent = (mz1 + mz2)/2
    } else {
      temp_mz_parent = -Inf
    }
    # Fwd dot product
    spec_score = 0
    for(j in 1:length(MS2_lib)){
      if(MS2_lib[[j]]$polarity == "positive"){
        spec1_df = as.data.frame(MS2_msr$pos_MS2_msr)
        if(nrow(spec1_df)<=1){next}
      }else{
        spec1_df = as.data.frame(MS2_msr$neg_MS2_msr)
        if(nrow(spec1_df)<=1){next}
      }
      spec2_df = MS2_lib[[j]][["spectrum"]][,c(1,2)]
      if(nrow(spec2_df)<=1){
        next
      }
      
      spec_merge_df = try(mergeMzIntensity(spec1_df, spec2_df, ppmTol = Related_files$global_parameter$instrument_parameter$ms1_ms2_match_ppm, absTol = 1e-3), silent = TRUE)
      if(inherits(spec_merge_df, "try-error")){
        spec_merge_df = mergeMzIntensity_backup(spec1_df, spec2_df, ppmTol = Related_files$global_parameter$instrument_parameter$ms1_ms2_match_ppm, absTol = 1e-3)
      }
      # spec_score[j] = Score_merge_MS2(spec_merge_df, mz_parent = temp_mz_parent)
      spec_score[j] = Score_merge_MS2_entropy(spec_merge_df, weight = 0.25)
      ## warning or even error occur here if one spectrum have same mz or close mz
      if(is.na(spec_score[j])){
        spec_merge_df = mergeMzIntensity_backup(spec1_df, spec2_df, ppmTol = Related_files$global_parameter$instrument_parameter$ms1_ms2_match_ppm, absTol = 1e-3)
        # spec_score[j] = Score_merge_MS2(spec_merge_df, mz_parent = temp_mz_parent)
        spec_score[j] = Score_merge_MS2_entropy(spec_merge_df, weight = 0.25)
      }
    }
    # print(spec_score)
    fwd_dp[i] = max(spec_score, na.rm = TRUE)
    
    # Rev dot product
    spec_score = 0
    for(j in 1:length(MS2_lib)){
      if(MS2_lib[[j]]$polarity == "positive"){
        spec1_df = as.data.frame(MS2_msr$pos_MS2_msr)
        if(nrow(spec1_df)<=1){next}
      }else{
        spec1_df = as.data.frame(MS2_msr$neg_MS2_msr)
        if(nrow(spec1_df)<=1){next}
      }
      spec2_df = MS2_lib[[j]][["spectrum"]][,c(1,2)]
      if(nrow(spec2_df)<=1){
        next
      }
      
      spec1_df[,1] = mz1 - spec1_df[,1]
      spec2_df[,1] = mz2 - spec2_df[,1]
      
      spec1_df = spec1_df[nrow(spec1_df):1,]
      spec2_df = spec2_df[nrow(spec2_df):1,]
      
      spec_merge_df = try(mergeMzIntensity(spec1_df, spec2_df, ppmTol = Related_files$global_parameter$instrument_parameter$ms1_ms2_match_ppm, absTol = 1e-3), silent = TRUE)
      if(inherits(spec_merge_df, "try-error")){
        spec_merge_df = mergeMzIntensity_backup(spec1_df, spec2_df, ppmTol = Related_files$global_parameter$instrument_parameter$ms1_ms2_match_ppm, absTol = 1e-3)
      }
      # spec_score[j] = Score_merge_MS2(spec_merge_df, mz_parent = 0)
      spec_score[j] = Score_merge_MS2_entropy(spec_merge_df, weight = 0.25)
      
      ## warning or even error occur here if one spectrum have same mz or close mz
      if(is.na(spec_score[j])){
        spec_merge_df = mergeMzIntensity_backup(spec1_df, spec2_df, ppmTol = Related_files$global_parameter$instrument_parameter$ms1_ms2_match_ppm, absTol = 1e-3)
        # spec_score[j] = Score_merge_MS2(spec_merge_df, mz_parent = 0)
        spec_score[j] = Score_merge_MS2_entropy(spec_merge_df, weight = 0.25)
      }
    }
    rev_dp[i] = max(spec_score, na.rm = TRUE)
  }
  
  # How to combine fwd score and rev score into a MS2 score?
  # For exact matched MS2, it takes only fwd spectrum score
  # For transformed MS2, it takes the higher one of the fwd and rev spectrum score
  StructureSet_df_MS2_ = StructureSet_df_MS2 %>%
    mutate(fwd_score = fwd_dp,
           rev_score = rev_dp,
           higher_fwd_rev = pmax(fwd_score, rev_score)) %>%
    mutate(score_MS2 = case_when(
      transform=="" & fwd_score>MS2_match_cutoff ~ fwd_score*MS2_match,
      transform!="" & higher_fwd_rev>MS2_similarity_cutoff ~ higher_fwd_rev*MS2_similarity,
      TRUE ~ 0 # The rest
    ))%>%
    dplyr::select(struct_set_id, score_MS2)
}

## Score_structureset_rdbe ####
Score_structureset_rdbe = function(StructureSet_df){
  # Penalize formula based on RDBE rule
  structureset_rdbe = StructureSet_df %>%
    mutate(score_RDBE = ifelse(rdbe >= 0 | category != "Metabolite", 0, -10)) %>%
    dplyr::select(struct_set_id, score_RDBE)
}

## Score_structureset_element_ratio ####
Score_structureset_element_ratio = function(StructureSet_df){
  # Penalize formula based on P/O and Si/O ratio
  temp_formula = StructureSet_df$formula
  
  formula_O = str_extract_all(temp_formula, "(?<=O)([:digit:]|\\.)+")
  formula_O_num = sapply(formula_O, function(x){sum(as.numeric(x))})
  formula_P = str_extract_all(temp_formula, "(?<=P)([:digit:]|\\.)+")
  formula_P_num = sapply(formula_P, function(x){sum(as.numeric(x))})
  formula_Si = str_extract_all(temp_formula, "(?<=Si)([:digit:]|\\.)+")
  formula_Si_num = sapply(formula_Si, function(x){sum(as.numeric(x))})
  
  structureset_element_ratio = StructureSet_df %>%
    mutate(score_element_ratio = case_when(
      formula_O_num < 3*formula_P_num ~ -10,
      formula_O_num < 2*formula_Si_num ~ -10,
      TRUE ~ 0
    )) %>%
    dplyr::select(struct_set_id, score_element_ratio)
}

## Score_structureset_missing_isotope ####
Score_structureset_missing_isotope = function(Related_files, Mset, StructureSet_df, NodeSet, EdgeSet_all,
                                              isotope = c("Cl"),
                                              isotope_penalty = -1,
                                              isotope_inten_cutoff = 5e4){
  input_mode = Related_files$global_parameter$mode
  # isotope = c("Cl","C")
  if(length(isotope_penalty) != length(isotope)){
    if(length(isotope_penalty) != 1){
      warning("incosistnet length for 'isotope' and 'isotope_penalty'")
    } else {
      isotope_penalty = rep(isotope_penalty, length(isotope))
    }
  }
  if(length(isotope_inten_cutoff) != length(isotope)){
    if(length(isotope_inten_cutoff) != 1){
      warning("incosistnet length for 'isotope' and 'isotope_inten_cutoff'")
    } else {
      isotope_inten_cutoff = rep(isotope_inten_cutoff, length(isotope))
    }
  }
  
  temp_missing_isotope_ls = list()
  
 
  node_inten = sapply(NodeSet, get_inten)
  for(i in 1:length(isotope)){
      iso_edge = EdgeSet_all %>%
        filter(grepl(paste0(isotope[i],"-1"), linktype))
      temp_missing_isotope = StructureSet_df %>%
        filter(str_detect(formula, paste0("(?<!\\])",isotope[i]))) %>% # Find element but not isotope formulas
        filter(node_inten[node_id] > log10(isotope_inten_cutoff[i])) %>%
        mutate(score_missing_isotope = ifelse(node_id %in% iso_edge$node1, 0, isotope_penalty[i])) %>%
        dplyr::select(struct_set_id, score_missing_isotope)
      colnames(temp_missing_isotope)[2] = paste0(colnames(temp_missing_isotope)[2],"_",isotope[i])
      temp_missing_isotope_ls[[length(temp_missing_isotope_ls)+1]] = temp_missing_isotope
  }
  
  structureset_missing_isotope = suppressMessages(
    Reduce(full_join, temp_missing_isotope_ls))
  
}


# Score_structure ####
Score_structureset = function(Related_files, Mset, StructureSet,NodeSet, LibrarySet, EdgeSet_all){
  StructureSet_df = bind_rows(StructureSet) %>%
    mutate(struct_set_id = 1:nrow(.)) %>%
    mutate(class = case_when(
      category %in% c("Unknown") ~ "Unknown",
      category %in% c("Metabolite") & steps == 0 & transform == "" ~ "Metabolite",
      category %in% c("Metabolite") & transform != "" ~ "Putative metabolite",
      TRUE ~ "Artifact"
    ))
  
  structureset_database_origin = Score_structureset_database_origin(StructureSet_df, LibrarySet,
                                                                   manual_match = 1)
  structureset_mz = Score_structureset_mz(StructureSet_df, NodeSet, 
                                          score_ppm_error_rate = Related_files$global_parameter$instrument_parameter$score_mz_ppm_error_rate)
  
  
  structureset_RT = Score_structureset_RT(StructureSet_df, NodeSet, LibrarySet,
                                          rt_match = 1, known_rt_tol = 0.5)
  
  structureset_MS2 = Score_structureset_MS2(Mset, StructureSet_df, NodeSet,Related_files, LibrarySet,
                                            # only spectra matching/similarity score > cutoff, then nodes get bonus
                                            MS2_match = 1, MS2_match_cutoff = 0.5, 
                                            MS2_similarity = 0.5, MS2_similarity_cutoff = 0.5)
  # empirical rules
  structureset_rdbe = Score_structureset_rdbe(StructureSet_df)
  structureset_element_ratio = Score_structureset_element_ratio(StructureSet_df)
  structureset_missing_isotope = Score_structureset_missing_isotope(Related_files, Mset, StructureSet_df, NodeSet, EdgeSet_all,
                                                                    isotope = c("Cl"), isotope_penalty = c(-0.5), 
                                                                    isotope_inten_cutoff = 5*Related_files$global_parameter$instrument_parameter$edge_expand_inten_cutoff)
  
  structureset_list = list(StructureSet_df,
                           structureset_database_origin,
                           structureset_mz, 
                           structureset_RT,
                           structureset_MS2,
                           structureset_rdbe,
                           structureset_element_ratio,
                           structureset_missing_isotope)
  
  structureset_list = structureset_list[!sapply(structureset_list, is.null)]
  
  StructureSet_df2 = suppressMessages(
    Reduce(left_join, structureset_list))
  
}

# Score_structure_propagation ####
Score_structure_propagation = function(StructureSet_df,
                                       # bio_decay = -1,
                                       artifact_decay = -0.5){
  
  # summing all scores as propagation_prior
  {
    score_sum = StructureSet_df %>% 
      dplyr::select(starts_with("score")) %>% 
      rowSums(na.rm = TRUE)
    
    structure_propagation = StructureSet_df %>%
      mutate(propagation_prior = score_sum) 
  }
  
  # propagation
  {
    structure_propagation = structure_propagation %>%
      mutate(score_propagation = 0) 
    
    unique_step = sort(unique(structure_propagation$steps))[2]
    
    ## Only artifact is propagated here 
    ## If nodes derived from biotransform needs scores, 
    ## dplyr::select metabolites and change step from 0.01 to 1
    for(unique_step in sort(unique(structure_propagation$steps))){
      # dplyr::select the parent and the parent_propagation_score
      structure_propagation_parent = structure_propagation %>%
        filter(propagation_prior + score_propagation + artifact_decay > 0) %>%
        filter(abs(steps - (unique_step - 0.01)) < 1e-6) %>% # 1e-6 as floating error
        mutate(parent_propagation_score = score_propagation+propagation_prior) %>%
        dplyr::select(node_id, formula, parent_propagation_score) %>%
        dplyr::rename(parent_id = node_id, 
                      parent_formula = formula)
      
      # Propagates parent_propagation_score to annotations derived from parents
      structure_propagation_target = suppressMessages(structure_propagation %>%
                                                        filter(steps == unique_step) %>%
                                                        inner_join(structure_propagation_parent, relationship="many-to-many") %>%
                                                        mutate(score_propagation = case_when(
                                                          category %in% c("Ring_artifact") ~ 0,
                                                          TRUE ~ pmax(parent_propagation_score + artifact_decay, 0)
                                                        )) %>%
                                                        dplyr::select(-parent_propagation_score) %>%
                                                        arrange(-score_propagation) %>%
                                                        distinct(struct_set_id, .keep_all = TRUE))
      
      # Update structure_propagation
      structure_propagation = suppressMessages(structure_propagation %>%
                                                 full_join(structure_propagation_target) %>%
                                                 arrange(-score_propagation) %>%
                                                 distinct(struct_set_id, .keep_all = TRUE))
    }
  }
  
  structure_propagation = structure_propagation %>% 
    dplyr::select(struct_set_id, score_propagation)
  StructureSet_df2 = suppressMessages(Reduce(left_join, list(StructureSet_df, 
                                                             structure_propagation)))
}

# initiate_ilp_nodes ####
initiate_ilp_nodes = function(StructureSet_df, NodeSet){
  # summing all scores
  {
    score_sum = StructureSet_df %>% 
      dplyr::select(starts_with("score")) %>% 
      rowSums(na.rm = TRUE)
    
    StructureSet_df = StructureSet_df %>%
      mutate(sum_score = score_sum) 
  }
  node_mz = sapply(NodeSet, get_mz)
  node_rt = sapply(NodeSet, get_rt)
  node_inten = sapply(NodeSet, get_inten)
  node_data = data.frame(node_id = as.numeric(1:length(node_mz)),
                         log10_inten = node_inten,
                         medMz = node_mz,
                         medRt = node_rt)
  ilp_nodes = suppressMessages(StructureSet_df %>%
                                 arrange(node_id, -sum_score, steps, category) %>%
                                 distinct(node_id, formula, class, .keep_all=TRUE) %>%
                                 #filter(sum_score >= 0) %>%
                                 arrange(node_id) %>%
                                 mutate(ilp_node_id = 1:nrow(.)) %>%
                                 inner_join(node_data %>%
                                              dplyr::select(node_id, log10_inten, medMz, medRt)) %>%
                                 dplyr::select(ilp_node_id, everything()))
  
  for(i in 1:ncol(ilp_nodes)){
    newval = ilp_nodes[,i]
    if(length(class(newval))!=1){
      stop("ilp_nodes format error")
      # print(i)
      # ilp_nodes[,i] = newval
    }
  }
  
  ilp_nodes
  
  # hist(ilp_nodes$score_sum)
}

# initiate_ilp_edges ####
initiate_ilp_edges = function(EdgeSet_all, CplexSet, Exclude = ""){
  ilp_nodes = CplexSet$ilp_nodes %>%
    arrange(ilp_node_id) 
  
  EdgeSet_df = EdgeSet_all %>%
    filter(category != "Heterodimer") %>%
    filter(!category %in% Exclude)
  
  ilp_nodes_assigned = ilp_nodes %>% filter(category != "Unknown")
  node_id_ilp_node_mapping = split(ilp_nodes_assigned$ilp_node_id, ilp_nodes_assigned$node_id)
  ilp_nodes_formula = ilp_nodes %>% pull(formula)
  
  ## Identify connections that reflects atomic differences
  {
    match_matrix_index_ls = list()
    
    for(i in 1:nrow(EdgeSet_df)){
      # if(i%%1000==0) print(i)
      edge_id = EdgeSet_df$edge_id[i]
      category = EdgeSet_df$category[i]
      
      node1 = EdgeSet_df$node1[i]
      ilp_nodes1 = node_id_ilp_node_mapping[[as.character(node1)]]
      formula1 = ilp_nodes_formula[ilp_nodes1]
      if(length(formula1) == 0){next}
      node2 = EdgeSet_df$node2[i]
      ilp_nodes2 = node_id_ilp_node_mapping[[as.character(node2)]]
      formula2 = ilp_nodes_formula[ilp_nodes2]
      if(length(formula2) == 0){next}
      
      ## Calculating formula_transform based on formula 1 and different transformation
      # fast_calculate_formula vs. my_calculate_formula
      if(category %in% c("Biotransform", "Adduct", "Fragment", 
                         "Natural_abundance", "Radical","Multicharge_isotope")){
        transform = EdgeSet_df$linktype[i]
        formula_transform = fast_calculate_formula(formula1, transform, 1)
      } else if(category == "Oligomer") {
        fold = as.numeric(EdgeSet_df$linktype[i])
        formula_transform = mapply(fast_calculate_formula, formula1, formula1, fold-1)
      } else if(category == "Multicharge") {
        fold = as.numeric(EdgeSet_df$linktype[i])
        formula_transform = mapply(fast_calculate_formula, formula1, formula1, fold-1)
      } else if(category == "Library_MS2_fragment"){
        formula1_position = which(formula1 == EdgeSet_df$formula1[i])
        formula_transform = rep("", length(formula1))
        formula_transform[formula1_position] = EdgeSet_df$formula2[i]
      } else if(category == "Experiment_MS2_fragment"){
        formula_transform = sub("MS2_fragment_", "", formula1)
      } else if(category == c("Ring_artifact")){
        formula_transform = paste0("Ring_artifact_", formula1)
      } else (
        stop("Unknown category in Edge list")
      )
      
      ## Matching formula_transform and formula2
      formula_match_matrix = matrix(TRUE, length(formula_transform), length(formula2))
      for(j in 1:length(formula_transform)){
        formula_match_matrix[j,] = formula_transform[j] == formula2
      }
      
      formula_match_matrix_index = which(formula_match_matrix == TRUE, arr.ind = TRUE) 
      if(dim(formula_match_matrix_index)[1] == 0){next}
      
      ## Recording
      match_matrix_index_ls[[length(match_matrix_index_ls)+1]] = list(edge_id = rep(edge_id, dim(formula_match_matrix_index)[1]),
                                                                      ilp_nodes1 = ilp_nodes1[formula_match_matrix_index[, 1]],
                                                                      ilp_nodes2 = ilp_nodes2[formula_match_matrix_index[, 2]])
    }
  }
  
  ## Filter by empirical rules that constraning node classes for connections
  {
    
    met_id = ilp_nodes %>%
      filter(class %in% c("Metabolite", "Putative metabolite")) %>% 
      pull(ilp_node_id)
    
    nonmet_id = ilp_nodes %>%
      filter(class %in% c("Artifact")) %>% 
      pull(ilp_node_id)
    
    library_met_id = ilp_nodes %>%
      filter(class %in% c("Metabolite")) %>% 
      pull(ilp_node_id)
    
    radical_id = ilp_nodes %>%
      filter(category == "Radical") %>%
      pull(ilp_node_id)
    
    ilp_edges = suppressMessages(bind_rows(match_matrix_index_ls) %>%
                                   left_join(EdgeSet_df)) %>%
      left_join(ilp_nodes %>% dplyr::select(ilp_node_id, steps), by = c("ilp_nodes1"="ilp_node_id"))
    
    ilp_edges_filter = ilp_edges %>%
      filter(case_when(
        category == "Biotransform" & ilp_nodes1 %in% met_id & ilp_nodes2 %in% met_id ~ TRUE,
        category %in% c("Fragment","Experiment_MS2_fragment") & ilp_nodes1 %in% nonmet_id ~ TRUE,
        category == "Oligomer" & steps==0 & ilp_nodes2 %in% nonmet_id  ~ TRUE, 
        category == "Multicharge" & ilp_nodes1 %in% nonmet_id & ilp_nodes2 %in% met_id ~ TRUE, 
        category == "Library_MS2_fragment" & ilp_nodes1 %in% nonmet_id & ilp_nodes2 %in% library_met_id ~ TRUE, # node2 has to be metabolite because of HMDB library matching
        category == "Radical" & ilp_nodes1 %in% radical_id ~ TRUE, # specific radical id
        category %in% c("Adduct", "Natural_abundance", "Multicharge_isotope", "Ring_artifact") & ilp_nodes2 %in% nonmet_id ~ TRUE,
        TRUE ~ FALSE
      )) %>%
      dplyr::select(-steps) %>%
      mutate(formula1 = ilp_nodes_formula[ilp_nodes1],
             formula2 = ilp_nodes_formula[ilp_nodes2]) %>%
      mutate(ilp_edge_id = 1:nrow(.))
  }
  
  ilp_edges_filter
  
}

# score_ilp_nodes ####
score_ilp_nodes = function(CplexSet, 
                           metabolite_score = 0.1,
                           putative_metabolite_score = 0.1, 
                           artifact_score = 0.1, 
                           unknown_score = -0.5)
{
  
  ilp_nodes = CplexSet$ilp_nodes
  
  ilp_nodes = ilp_nodes %>%
    mutate(score_class = case_when(
      class == "Metabolite" ~ metabolite_score,
      class == "Putative metabolite" ~ putative_metabolite_score,
      class == "Artifact" ~ artifact_score,
      class == "Unknown" ~ unknown_score
      # class == "Unknown" ~ unknown_score-(log10_inten-6)*0.2,
    ))
  
  cplex_score_node = ilp_nodes %>% 
    dplyr::select(starts_with("score_")) %>% 
    rowSums(na.rm = TRUE)
  
  ilp_nodes = ilp_nodes %>%
    mutate(cplex_score = round(cplex_score_node,7))
  
  return(ilp_nodes)
}






## Score_type_category ####
Score_type_category = function(ilp_edges,
                               NodeSet,
                               Mset, 
                               type_score_biotransform = 0, type_score_adduct = 0.5, 
                               type_score_natural_abundance = 1, type_score_fragment = 0.3, 
                               type_score_radical = 0.2, type_score_oligomer = 0.5, 
                               type_score_multicharge = 0.5, type_score_multicharge_isotope = 1, 
                               type_score_ring_artifact = 2,
                               type_score_experiment_MS2_fragment = 1, 
                               type_score_library_MS2_fragment = 0.3){
  # Score rule category 
  input_mode = Related_files$global_parameter$mode
  
  pos_count = 1:length(NodeSet)
  neg_count = 1:length(NodeSet)
  
  ilp_edges = ilp_edges %>%
    mutate(score_category = case_when(
      category == "Biotransform" ~ type_score_biotransform,
      category == "Adduct" ~ type_score_adduct, 
      # category == "Natural_abundance" ~ type_score_natural_abundance,
      category == "Natural_abundance" & ((node1 %in% pos_count & node2 %in% pos_count) | (node1 %in% neg_count & node2 %in% neg_count)) ~ type_score_natural_abundance,
      category == "Fragment" ~ type_score_fragment,
      category == "Radical" ~ type_score_radical,
      category == "Oligomer" ~ type_score_oligomer,
      category == "Multicharge" ~ type_score_multicharge, 
      # category == "Multicharge_isotope" ~ type_score_multicharge_isotope,
      category == "Multicharge_isotope" & ((node1 %in% pos_count & node2 %in% pos_count) | (node1 %in% neg_count & node2 %in% neg_count)) ~ type_score_multicharge_isotope,
      category == "Ring_artifact" ~ type_score_ring_artifact,
      category == "Experiment_MS2_fragment" ~ type_score_experiment_MS2_fragment,
      category == "Library_MS2_fragment" ~ 0, # special handling below
      TRUE ~ 0
    )) %>%
    dplyr::select(ilp_edge_id, score_category) %>%
    filter(TRUE)
}

#### Get_Ints obtain vector for Chro_Similarity ####
Get_Ints = function(rts, rt_len, rt_list, ints_list ){
  
  rt_list = round(rt_list,5)
  # rts = seq(from = rt_mean - delt, by = (rt_list[301]-rt_list[1])/300, length = rt_len)
  rts = round(rts,5)
  rt = rts[rts>=min(rt_list) & rts<=max(rt_list)]
  if(length(rt)!=0){
    if(rt[1] %in% rt_list){
      index1 = which(rt_list == rt[1])
      # index2 = length(rt_list)
      # ints = ifelse(rt[length(rt)]<=max(rt_list), ints_list[index1:(index1+rt_len-1)], ints_list[index1:index2] )
      if(rt[length(rt)]<=max(rt_list)){
        ints =  ints_list[index1:(index1+length(rt)-1)]
        a = rep(0, rt_len-length(ints))
      }else{
        index2 = length(rt_list)
        ints = ints_list[index1:index2]
        a = rep(0,length(rt)-length(ints))
      }
      
      ints = c(ints, a)
    }else{
      index1 = (which(rt_list>rt[1])-1)[1]
      index2 = (which(rt_list>rt[length(rt)]))[1]
      rate = (rt[1]-rt_list[index1])/(rt_list[index1+1]-rt_list[index1])
      ints = ints_list[index1:index2]
      ints_diff = diff(ints)*rate
      ints = ints[-length(ints)]+ints_diff
      a = rep(0, length(which(rts<min(rt_list))))
      b = rep(0, length(which(rts>max(rt_list))))
      ints = c(a, ints, b)
    }
    return(ints)
  }else{
    ints = rep(0, rt_len)
    return(ints)
  }
  
}

#### Cos_Calculate calculate cosine ####
Cos_Calculate = function(a, b){
  return(a%*%b / (sqrt(a%*%a) * sqrt(b%*%b)))
}

### Chro_Similarity calculate similarity of eligible peaks####
Chro_Similarity = function(x, NodeSet, delt, rt_len){
  node1 = x[["node1"]]
  node2 = x[["node2"]]
  rt1 = x[["rt_node1"]]
  rt2 = x[["rt_node2"]]
  
  if(! is.null(NodeSet[[as.numeric(node1)]]$pos_info)){
    chro1_pos_rt = t(as.data.frame(NodeSet[[as.numeric(node1)]]$pos_info$pos_ms$rtime_list))
    chro1_pos_ints = t(as.data.frame(NodeSet[[as.numeric(node1)]]$pos_info$pos_ms$intensity_list))
  }else{
    chro1_pos_rt = NULL
    chro1_pos_ints = NULL
  }
  if(! is.null(NodeSet[[as.numeric(node1)]]$neg_info)){
    chro1_neg_rt = t(as.data.frame(NodeSet[[as.numeric(node1)]]$neg_info$neg_ms$rtime_list))
    chro1_neg_ints = t(as.data.frame(NodeSet[[as.numeric(node1)]]$neg_info$neg_ms$intensity_list))
  }else{
    chro1_neg_rt = NULL
    chro1_neg_ints = NULL
  }
  if(! is.null(NodeSet[[as.numeric(node2)]]$pos_info)){
    chro2_pos_rt = t(as.data.frame(NodeSet[[as.numeric(node2)]]$pos_info$pos_ms$rtime_list))
    chro2_pos_ints = t(as.data.frame(NodeSet[[as.numeric(node2)]]$pos_info$pos_ms$intensity_list))
  }else{
    chro2_pos_rt = NULL
    chro2_pos_ints = NULL
  }
  if(! is.null(NodeSet[[as.numeric(node2)]]$neg_info)){
    chro2_neg_rt = t(as.data.frame(NodeSet[[as.numeric(node2)]]$neg_info$neg_ms$rtime_list))
    chro2_neg_ints = t(as.data.frame(NodeSet[[as.numeric(node2)]]$neg_info$neg_ms$intensity_list))
  }else{
    chro2_neg_rt = NULL
    chro2_neg_ints = NULL
  }
  chro1_rt = rbind(chro1_pos_rt, chro1_neg_rt)
  chro1_ints = rbind(chro1_pos_ints, chro1_neg_ints)
  chro2_rt = rbind(chro2_pos_rt, chro2_neg_rt)
  chro2_ints = rbind(chro2_pos_ints, chro2_neg_ints)
  
  res = vector()
  for(i in 1:nrow(chro1_rt)){
    for(j in 1:nrow(chro2_rt)){
      
      rt_mean = (rt1 + rt2)*60/2
      rt_range = seq(rt_mean - delt, rt_mean + delt, length = rt_len)
      
      ints1 = Get_Ints(rt_range, rt_len, chro1_rt[i,], chro1_ints[i,])
      ints2 = Get_Ints(rt_range, rt_len, chro2_rt[j,], chro2_ints[j,])
      cos = Cos_Calculate(ints1, ints2)
      res = c(res, cos)
    }
  }
  cos_median = median(res, na.rm = TRUE)
  return(cos_median)
}

## Score_Chromatogram_Similarity ####
Score_Chromatogram_Similarity = function(ilp_edges, EdgeSet_all, NodeSet,
                                         chro_cutoff_artifact = 0.8,
                                         chro_score_artifact_multiplier = -1){
  node_rt = sapply(NodeSet, get_rt)
  node = EdgeSet_all %>% 
    mutate(rt_node1 = node_rt[node1],
           rt_node2 = node_rt[node2]) %>%
    filter(category != "Biotransform") %>%
    # filter(abs(rt_node1-rt_node2)<=0.2) %>%
    dplyr::select(edge_id, node1, node2, rt_node1, rt_node2) 
  score_chro = apply(node, 1, Chro_Similarity, NodeSet, 15, 151)
  # no_cores = detectCores()-2
  # cl = makeCluster(no_cores)
  # score_chro = parApply(cl, node, 1, Chro_Similarity, NodeSet, 15, 151)
  # stopCluster(cl)
  node = node %>%
    mutate(score_chromatogram = score_chro) %>%
    dplyr::select(edge_id, score_chromatogram)
  
  ilp_edge = EdgeSet_all %>%
    dplyr::select(edge_id) %>%
    left_join(node, by = 'edge_id') %>%
    right_join(ilp_edges %>%
                 dplyr::select(edge_id, category,ilp_edge_id), by = 'edge_id') 
  ilp_edge$score_chromatogram[which(is.na(ilp_edge$score_chromatogram))] = 0
  ilp_edge_chro = ilp_edge %>%
    mutate(score_chromatogram_artifact = case_when(
      category %in% c("Biotransform") ~ 0,
      score_chromatogram > chro_cutoff_artifact ~ 0,
      TRUE ~ (1-score_chromatogram) * chro_score_artifact_multiplier)) %>%
    dplyr::select(ilp_edge_id, score_chromatogram_artifact) %>%
    filter(TRUE)
  
}

## Score_inten_isotope ####
Score_inten_isotope = function(Related_files, ilp_edges, Mset, NodeSet, 
                               inten_cutoff_isotope = 3,
                               score_sigma_isotope = 0.2
){
  
  input_mode = Related_files$global_parameter$mode
  
  node_inten = sapply(NodeSet, get_inten) %>% as.numeric()
  ilp_edges_isotope = ilp_edges %>%
    filter(category == "Natural_abundance") %>%
    mutate(inten1 = node_inten[node1],
           inten2 = node_inten[node2]) %>%
    mutate(inten_ratio_measured = inten2 - inten1, 
           inten_ratio_calculated = log10(mapply(isotopic_abundance, .$formula1, .$linktype)),
           measured_calculated_ratio = 10^(inten_ratio_measured-inten_ratio_calculated)) %>%
    mutate(p_obs = dnorm(measured_calculated_ratio, 1, score_sigma_isotope+10^(inten_cutoff_isotope-pmin(inten1, inten2))),
           p_theory = dnorm(1, 1, score_sigma_isotope+10^(inten_cutoff_isotope-pmin(inten1, inten2)))) %>%
    mutate(score_inten_isotope = log10(p_obs/p_theory+1e-10)) %>%
    dplyr::select(ilp_edge_id, score_inten_isotope) %>%
    mutate(score_inten_isotope = as.numeric(score_inten_isotope)) # remove names for score_inten_isotope
  
  return(ilp_edges_isotope)
}

## Score_MS2_library_fragment ####
Score_MS2_library_fragment = function(ilp_edges, 
                                      MS2_score_library_fragment = 0.5){
  # Score library_MS2_fragment
  # Describe if given any peak related to a HMDB metabolite, 
  # look for the metabolite's MS2, to see if any fragment is presented in MS1
  # If so, strengthen the formula pair of the fragment and parent
  if(nrow(ilp_edges %>% filter(category == "Library_MS2_fragment"))==0){
    return(NULL)
  }
  ilp_edges_MS2_library_fragment = ilp_edges %>%
    filter(category == "Library_MS2_fragment") %>%
    group_by(ilp_nodes2) %>% mutate(n = n()) %>% ungroup() %>%
    mutate(score_MS2_library_fragment = MS2_score_library_fragment/sqrt(n)) %>%
    dplyr::select(ilp_edge_id, score_MS2_library_fragment) %>%
    filter(TRUE)
  
}
### MS2_similarity edge ms2 score function####
MS2_similarity = function(spec1_df, spec2_df, ppmTol, absTol, temp_mz_parent, mz1, mz2){
  colnames(spec1_df)=colnames(spec2_df)=c("mz", "inten")
  if(nrow(spec1_df)==1 | nrow(spec2_df)==1){return(NULL)}
  
  # Fwd dot product
  spec_merge_df = try(mergeMzIntensity(spec1_df, spec2_df, ppmTol = ppmTol, absTol = absTol), silent = TRUE)
  if(inherits(spec_merge_df, "try-error")){
    spec_merge_df = mergeMzIntensity_backup(spec1_df, spec2_df, ppmTol = ppmTol, absTol = absTol)
  }
  # fwd_dp = Score_merge_MS2(spec_merge_df, mz_parent = temp_mz_parent)
  fwd_dp = Score_merge_MS2_entropy(spec_merge_df, weight = 0.25)
  ## warning or even error occur here if one spectrum have same mz or close mz
  if(is.na(fwd_dp)){
    spec_merge_df = mergeMzIntensity_backup(spec1_df, spec2_df, ppmTol = ppmTol, absTol = absTol)
    # fwd_dp = Score_merge_MS2(spec_merge_df, mz_parent = temp_mz_parent)
    fwd_dp = Score_merge_MS2_entropy(spec_merge_df, weight = 0.25)
  }
  
  # Rev dot product
  spec1_df[,1] = mz1 - spec1_df[,1]
  spec2_df[,1] = mz2 - spec2_df[,1]
  spec1_df = spec1_df[nrow(spec1_df):1,]
  spec2_df = spec2_df[nrow(spec2_df):1,]
  spec_merge_df = try(mergeMzIntensity(spec1_df, spec2_df, ppmTol = ppmTol, absTol = absTol), silent = TRUE)
  if(inherits(spec_merge_df, "try-error")){
    spec_merge_df = mergeMzIntensity_backup(spec1_df, spec2_df, ppmTol = ppmTol, absTol = absTol)
  }
  # rev_dp = Score_merge_MS2(spec_merge_df, mz_parent = 0)
  rev_dp = Score_merge_MS2_entropy(spec_merge_df, weight = 0.25)
  if(is.na(rev_dp)){
    spec_merge_df = mergeMzIntensity_backup(spec1_df, spec2_df, ppmTol = ppmTol, absTol = absTol)
    # rev_dp = Score_merge_MS2(spec_merge_df, mz_parent = 0)
    rev_dp = Score_merge_MS2_entropy(spec_merge_df, weight = 0.25)
  }
  return(c(fwd_dp,rev_dp))
}
## Score_MS2_experiment_biotransform ####
Score_MS2_experiment_biotransform = function(ilp_edges, NodeSet, 
                                             MS2_score_similarity = 1, 
                                             MS2_similarity_cutoff = 0.3,
                                             ppmTol = 10E-6, 
                                             absTol = 1e-3) 
{
  # Score Experiment MS2 similarity
  # It compares Biotransform connection in edge list where both experiment MS2 exist
  # If their MS2 are similar, then strengthen such connection
  
  node_MS2_valid = sapply(NodeSet, function(x){!is.null(x$MS2)})
  if(!any(node_MS2_valid)){return (NULL)}
  
  node_MS2 = (1:length(NodeSet))[node_MS2_valid]
  
  ilp_edges_MS2_similarity = ilp_edges %>%
    distinct(node1, node2, .keep_all = TRUE) %>%
    filter(category == "Biotransform") %>% 
    filter(node1 %in% node_MS2, node2 %in% node_MS2)
  
  fwd_dp = rev_dp = rep(0, nrow(ilp_edges_MS2_similarity))
  for(i in 1:nrow(ilp_edges_MS2_similarity)){
    
    node1 = ilp_edges_MS2_similarity$node1[i]
    node2 = ilp_edges_MS2_similarity$node2[i]
    
    # Adjust parent ion because mz in NodeSet is neutral mass
    
    
    mz = sapply(NodeSet, get_mz)
    mz1 = mz[node1]
    mz2 = mz[node2]
    if(abs(mz1 - mz2) < max(absTol, mz2*ppmTol)){
      temp_mz_parent = (mz1 + mz2)/2
    } else {
      temp_mz_parent = -Inf
    }
    
    # Read node1 and node2 MS2
    H_mass = 1.007825032
    e_mass = 0.000548579
    
    spec1 = NodeSet[[node1]]$MS2 
    spec2 = NodeSet[[node2]]$MS2 
    spec1_pos_exist = !is.null(spec1$pos_MS2)
    spec1_neg_exist = !is.null(spec1$neg_MS2)
    spec2_pos_exist = !is.null(spec2$pos_MS2)
    spec2_neg_exist = !is.null(spec2$neg_MS2)
    full_exist = spec1_pos_exist & spec1_neg_exist & spec2_pos_exist & spec2_neg_exist
    if(spec1_pos_exist & spec2_pos_exist & !full_exist){
      spec1$pos_MS2$spec[,1] = spec1$pos_MS2$spec[,1] - (H_mass-e_mass)*1
      spec2$pos_MS2$spec[,1] = spec2$pos_MS2$spec[,1] - (H_mass-e_mass)*1
      spec1_df = spec1$pos_MS2$spec
      spec2_df = spec2$pos_MS2$spec
    }else if(spec1_neg_exist & spec2_neg_exist & !full_exist){
      spec1$neg_MS2$spec[,1] = spec1$neg_MS2$spec[,1] - (H_mass-e_mass)*(-1)
      spec2$neg_MS2$spec[,1] = spec2$neg_MS2$spec[,1] - (H_mass-e_mass)*(-1)
      spec1_df = spec1$neg_MS2$spec
      spec2_df = spec2$neg_MS2$spec
    }else if(full_exist){
      spec1$pos_MS2$spec[,1] = spec1$pos_MS2$spec[,1] - (H_mass-e_mass)*1
      spec2$pos_MS2$spec[,1] = spec2$pos_MS2$spec[,1] - (H_mass-e_mass)*1
      spec1$neg_MS2$spec[,1] = spec1$neg_MS2$spec[,1] - (H_mass-e_mass)*(-1)
      spec2$neg_MS2$spec[,1] = spec2$neg_MS2$spec[,1] - (H_mass-e_mass)*(-1)
      spec1_df_pos = spec1$pos_MS2$spec
      spec1_df_neg = spec1$neg_MS2$spec
      spec2_df_pos = spec2$pos_MS2$spec
      spec2_df_neg = spec2$neg_MS2$spec
    }else{
      next
    }
    if(!full_exist){
      ms2_score = MS2_similarity(spec1_df, spec2_df, ppmTol, absTol, temp_mz_parent, mz1, mz2)
      if(is.null(ms2_score)){next}
      fwd_dp[i] = ms2_score[1]
      rev_dp[i] = ms2_score[2]
    }else{
      ms2_score_pos = MS2_similarity(spec1_df_pos, spec2_df_pos, ppmTol, absTol, temp_mz_parent, mz1, mz2)
      ms2_score_neg = MS2_similarity(spec1_df_neg, spec2_df_neg, ppmTol, absTol, temp_mz_parent, mz1, mz2)
      if(is.null(ms2_score_pos) & is.null(ms2_score_pos)){next}
      if(max(ms2_score_pos) > max(ms2_score_neg)){
        fwd_dp[i] = ms2_score_pos[1]
        rev_dp[i] = ms2_score_pos[2]
      }else{
        fwd_dp[i] = ms2_score_neg[1]
        rev_dp[i] = ms2_score_neg[2]
      }
    }
  }
  if(is.null(fwd_dp) | is.null(rev_dp)){return(NULL)}
  ilp_edges_MS2_similarity_ = ilp_edges_MS2_similarity %>%
    mutate(fwd_score = fwd_dp,
           rev_score = rev_dp,
           higher_fwd_rev = pmax(fwd_score, rev_score)) %>%
    filter(higher_fwd_rev>MS2_similarity_cutoff) %>%
    group_by(ilp_nodes1) %>% mutate(n1 = n()) %>% ungroup() %>%
    group_by(ilp_nodes2) %>% mutate(n2 = n()) %>% ungroup() %>%
    mutate(higher_n = pmax(n1, n2)) %>%
    mutate(score_MS2_similarity =  higher_fwd_rev*MS2_score_similarity/sqrt(higher_n)) %>%
    dplyr::select(ilp_edge_id, score_MS2_similarity) %>%
    filter(TRUE)
  
}
## Score_MS2_abiotic_mz_appearance ####
Score_MS2_abiotic_mz_appearance = function(ilp_edges, NodeSet, 
                                           MS2_score_abiotic_mz_appearance = 0.3){
  # Using MS2 to strengthen artifact connections such as "Adduct", "Fragment", "Multicharge", "Oligomer", "Radical"
  # where node1's mz appears in node2's MS2 spectra
  MS2_experiment_fragment = ilp_edges %>%
    filter(category == "Experiment_MS2_fragment") %>%
    distinct(node1, node2)
  
  ilp_edges_MS2_similarity_artifact = suppressMessages(ilp_edges %>%
                                                         right_join(MS2_experiment_fragment) %>%
                                                         filter(category %in% c("Adduct", "Fragment", "Multicharge", "Oligomer", "Radical")) %>% 
                                                         mutate(score_MS2_abiotic_mz_appearance = MS2_score_abiotic_mz_appearance) %>%
                                                         dplyr::select(ilp_edge_id, score_MS2_abiotic_mz_appearance) %>%
                                                         filter(TRUE))
}

# score_ilp_edges ####
score_ilp_edges = function(Related_files, EdgeSet_all, CplexSet, NodeSet, Mset){
  ilp_edges = CplexSet$ilp_edges
  
  ilp_edges_type_category = Score_type_category(ilp_edges,
                                                NodeSet,
                                                Mset,
                                                type_score_biotransform = 0, type_score_adduct = 0.5, 
                                                type_score_natural_abundance = 1, type_score_fragment = 0.3, 
                                                type_score_radical = 0.2, type_score_oligomer = 0.5, 
                                                type_score_multicharge = 0.5, type_score_multicharge_isotope = 1, 
                                                type_score_ring_artifact = 2,
                                                type_score_experiment_MS2_fragment = 0.5, 
                                                type_score_library_MS2_fragment = 0.3)
  
  # ilp_edges_rt_penalty = Score_rt_penalty(ilp_edges, NodeSet,
  #                                         rt_cutoff_artifact = 0.05,
  #                                         rt_score_artifact_multiplier = -5)
  
  ilp_edges_chro_penalty = Score_Chromatogram_Similarity(ilp_edges, EdgeSet_all, NodeSet,
                                                         chro_cutoff_artifact = 0.8,
                                                         chro_score_artifact_multiplier = -0.5)
  
  ilp_edges_inten_isotope = Score_inten_isotope(Related_files, ilp_edges, Mset, NodeSet, 
                                                inten_cutoff_isotope = 3,
                                                score_sigma = 0.2)
  
  ilp_edges_MS2_library_fragment = Score_MS2_library_fragment(ilp_edges, 
                                                              MS2_score_library_fragment = 0.5)
  ilp_edges_MS2_experiment_biotransform = Score_MS2_experiment_biotransform(ilp_edges, NodeSet, 
                                                                            MS2_score_similarity = 1, 
                                                                            MS2_similarity_cutoff = 0.3,
                                                                            ppmTol = Related_files$global_parameter$instrument_parameter$ms1_ms2_match_ppm, 
                                                                            absTol = 1e-3)
  
  ilp_edges_MS2_abiotic_mz_appearance = Score_MS2_abiotic_mz_appearance(ilp_edges,
                                                                        MS2_score_abiotic_mz_appearance = 0.3)
  
  ilp_edges_list = list(ilp_edges, 
                        ilp_edges_type_category,
                        ilp_edges_chro_penalty,
                        ilp_edges_inten_isotope,
                        ilp_edges_MS2_library_fragment,
                        ilp_edges_MS2_experiment_biotransform,
                        ilp_edges_MS2_abiotic_mz_appearance)
  
  ilp_edges_list = ilp_edges_list[!sapply(ilp_edges_list, is.null)]
  
  ilp_edges2 = suppressMessages(
    Reduce(left_join, ilp_edges_list))
  
  
  cplex_score_edge = ilp_edges2 %>% 
    dplyr::select(starts_with("score_")) %>% 
    rowSums(na.rm=TRUE)
  
  ilp_edges2 = ilp_edges2 %>%
    mutate(cplex_score = round(cplex_score_edge,7))
  
  return(ilp_edges2)
}
# initiate_heterodimer_ilp_edges ####
initiate_heterodimer_ilp_edges = function(EdgeSet_all, CplexSet, NodeSet){
  
  node_inten = sapply(NodeSet, get_inten)
  EdgeSet_df = EdgeSet_all %>%
    filter(category == "Heterodimer") %>%
    mutate(linktype = as.numeric(linktype)) %>%
    filter(node_inten[node1] > node_inten[linktype]) # node1 and linktype are interchangeable in heterodimerconnection
  
  if(nrow(EdgeSet_df) == 0){
    return(NULL)
  }
  
  ilp_nodes = CplexSet$ilp_nodes %>%
    arrange(ilp_node_id) 
  
  ilp_nodes_assigned = ilp_nodes %>% filter(category != "Unknown")
  node_id_ilp_node_mapping = split(ilp_nodes_assigned$ilp_node_id, ilp_nodes_assigned$node_id)
  
  ilp_nodes_nonmet = ilp_nodes %>% 
    filter(!category %in% c("Metabolite","Unknown")) %>%
    filter(rdbe %% 1 == 0) %>% # Remove heterodimer which may be radicals
    filter(!grepl("\\[", formula)) # Remove natural abundance in heterodimer
  
  node_id_ilp_node_mapping_nonmet = split(ilp_nodes_nonmet$ilp_node_id, ilp_nodes_nonmet$node_id)
  
  ilp_nodes_formula = ilp_nodes %>%
    pull(formula)
  
  match_matrix_index_ls = list()
  
  for(i in 1:nrow(EdgeSet_df)){
    edge_id = EdgeSet_df$edge_id[i]
    
    # node1 and node_link can be anything, node2 has to be artifact
    node1 = EdgeSet_df$node1[i]
    ilp_nodes1 = node_id_ilp_node_mapping[[as.character(node1)]]
    formula1 = ilp_nodes_formula[ilp_nodes1]
    if(length(formula1) == 0){next}
    
    node2 = EdgeSet_df$node2[i]
    ilp_nodes2 = node_id_ilp_node_mapping_nonmet[[as.character(node2)]]
    formula2 = ilp_nodes_formula[ilp_nodes2]
    if(length(formula2) == 0){next}
    
    node_link = EdgeSet_df$linktype[i]
    ilp_nodes_link = node_id_ilp_node_mapping[[as.character(node_link)]]
    formula_link = ilp_nodes_formula[ilp_nodes_link]
    if(length(formula_link) == 0){next}
    
    # if(length(formula_link)>1 & length(formula1)>1)break
    transform = formula_link
    formula_transform_all = my_calculate_formula(formula1, transform)
    
    for(k in 1:length(transform)){
      if(length(transform) == 1){
        formula_transform = formula_transform_all
      } else {
        formula_transform = formula_transform_all[, k]      
      }
      
      formula_match_matrix = matrix(TRUE, length(formula_transform), length(formula2))
      for(j in 1:length(formula_transform)){
        formula_match_matrix[j,] = formula_transform[j] == formula2
      }
      
      formula_match_matrix_index = which(formula_match_matrix == TRUE, arr.ind = TRUE) 
      
      if(dim(formula_match_matrix_index)[1] == 0){next}
      
      match_matrix_index_ls[[length(match_matrix_index_ls)+1]] = list(edge_id = rep(edge_id, dim(formula_match_matrix_index)[1]),
                                                                      ilp_nodes1 = ilp_nodes1[formula_match_matrix_index[, 1]],
                                                                      ilp_nodes2 = ilp_nodes2[formula_match_matrix_index[, 2]],
                                                                      ilp_nodes_link = rep(ilp_nodes_link[k], dim(formula_match_matrix_index)[1]))
    }
  }
  heterodimer_ilp_edges = bind_rows(match_matrix_index_ls) %>%
    merge(EdgeSet_df, all.x = TRUE) %>%
    mutate(formula1 = ilp_nodes_formula[ilp_nodes1],
           formula2 = ilp_nodes_formula[ilp_nodes2],
           formula_link = ilp_nodes_formula[ilp_nodes_link])
  return(heterodimer_ilp_edges)
}



# score_heterodimer_ilp_edges ####
score_heterodimer_ilp_edges = function(CplexSet, type_score_heterodimer = 1,
                                       MS2_score_experiment_fragment = 0.5){
  
  if(is.null(CplexSet$heterodimer_ilp_edges )){
    return (NULL)
  }
  
  heterodimer_ilp_edges = CplexSet$heterodimer_ilp_edges %>%
    mutate(ilp_edge_id = 1:nrow(.) + nrow(CplexSet$ilp_edges)) %>%
    mutate(score_category = case_when(
      category == "Heterodimer" ~ type_score_heterodimer,
      TRUE ~ 0
    )) %>%
    filter(TRUE)
  
  if(MS2_score_experiment_fragment != 0){
    heterodimer_ilp_edges_experiment_MS2_fragment = CplexSet$ilp_edges %>%
      filter(category == "Experiment_MS2_fragment") %>%
      distinct(node1, node2) %>%
      merge(heterodimer_ilp_edges) %>%
      mutate(score_experiment_MS2_fragment = MS2_score_experiment_fragment) %>%
      dplyr::select(ilp_edge_id, score_experiment_MS2_fragment) %>%
      filter(TRUE)
    
    heterodimer_ilp_edges = heterodimer_ilp_edges %>%
      merge(heterodimer_ilp_edges_experiment_MS2_fragment, all.x = TRUE) %>%
      mutate(score_experiment_MS2_fragment = ifelse(is.na(score_experiment_MS2_fragment), 
                                                    0, 
                                                    score_experiment_MS2_fragment))
    
  }
  
  cplex_score_edge = heterodimer_ilp_edges %>% 
    dplyr::select(starts_with("score_")) %>% 
    rowSums(na.rm=TRUE)
  
  heterodimer_ilp_edges = heterodimer_ilp_edges %>%
    mutate(cplex_score = cplex_score_edge)
  
  return(heterodimer_ilp_edges)
}

# Initiate_cplexset ####
Initiate_cplexset = function(CplexSet, type){
  
  ilp_nodes = CplexSet$ilp_nodes %>%
    arrange(ilp_node_id)
  ilp_edges = CplexSet$ilp_edges %>%
    arrange(ilp_edge_id)
  
  heterodimer_ilp_edges = CplexSet$heterodimer_ilp_edges
  exist_heterodimer = TRUE
  
  if(is.null(heterodimer_ilp_edges)){
    exist_heterodimer = FALSE
    heterodimer_ilp_edges = NULL
    nrow_heterodimer_ilp_edges = 0
  } else {
    heterodimer_ilp_edges = CplexSet$heterodimer_ilp_edges %>%
      arrange(ilp_edge_id) %>%
      mutate(heterodimer_ilp_edge_id = 1:nrow(.))
    nrow_heterodimer_ilp_edges = nrow(heterodimer_ilp_edges)
  }
  
  # Basic constraint ####
  {
    # triplet_nodes
    # for one node_id, sum(ilp_node_id) = 1
    {
      ilp_rows = rep(1,length = nrow(ilp_nodes))
      ilp_nodes.node_id = ilp_nodes$node_id
      for(i in 2:length(ilp_nodes.node_id)){
        if(ilp_nodes.node_id[i] == ilp_nodes.node_id[i-1]){
          ilp_rows[i] = ilp_rows[i-1]
        } else{
          ilp_rows[i] = ilp_rows[i-1] + 1
        }
      }
      
      triplet_node = ilp_nodes %>%
        mutate(ilp_row_id = ilp_rows) %>%
        dplyr::select(ilp_row_id, ilp_node_id) %>%
        transmute(i = ilp_row_id,  
                  j = ilp_node_id, 
                  v = 1)
      
      if(length(unique(triplet_node$i))!=max(ilp_rows)){
        stop("Error in triplet_node")
      }
    }
    
    # triplet_edges
    # e12 - n1 <= 0, e12 - n2 <= 0
    # Because triplet_node take max(triplet_node$i) rows, and max(triplet_node$j) columns
    {
      triplet_edge_edge1 = ilp_edges %>%
        transmute(i = ilp_edge_id * 2 - 1 + max(triplet_node$i), 
                  j = ilp_edge_id + max(triplet_node$j),
                  v = 1)
      
      triplet_edge_edge2 = ilp_edges %>%
        transmute(i = ilp_edge_id * 2 + max(triplet_node$i), 
                  j = ilp_edge_id + max(triplet_node$j),
                  v = 1)
      
      triplet_edge_node1 = ilp_edges %>%
        transmute(i = ilp_edge_id * 2 - 1 + max(triplet_node$i),
                  j = ilp_nodes1,
                  v = -1)
      
      triplet_edge_node2 = ilp_edges %>%
        transmute(i = ilp_edge_id * 2 + max(triplet_node$i),
                  j = ilp_nodes2,
                  v = -1)
      
      triplet_edge = bind_rows(triplet_edge_edge1, 
                               triplet_edge_edge2,
                               triplet_edge_node1,
                               triplet_edge_node2)
      
      if(length(unique(triplet_edge$i))!=nrow(ilp_edges) * 2){
        stop("Error in triplet_edge")
      }
    }
    
    triplet_basic = bind_rows(triplet_node, 
                              triplet_edge)
    
    # triplet_heterodimer
    # e123 <= n1, e123 <= n2, e123 <= n3
    # Because triplet_edge take max(triplet_edge$i) rows, and max(triplet_edge$j) columns
    if(exist_heterodimer){
      triplet_edge_edge1 = heterodimer_ilp_edges %>%
        transmute(i = heterodimer_ilp_edge_id * 3 - 2 + max(triplet_edge$i), 
                  j = heterodimer_ilp_edge_id + max(triplet_edge$j),
                  v = 1)
      
      triplet_edge_edge2 = heterodimer_ilp_edges %>%
        transmute(i = heterodimer_ilp_edge_id * 3 - 1 + max(triplet_edge$i), 
                  j = heterodimer_ilp_edge_id + max(triplet_edge$j),
                  v = 1)
      
      triplet_edge_edge3 = heterodimer_ilp_edges %>%
        transmute(i = heterodimer_ilp_edge_id * 3 + max(triplet_edge$i), 
                  j = heterodimer_ilp_edge_id + max(triplet_edge$j),
                  v = 1)
      
      triplet_edge_node1 = heterodimer_ilp_edges %>%
        transmute(i = heterodimer_ilp_edge_id * 3 - 2 + max(triplet_edge$i),
                  j = ilp_nodes1,
                  v = -1)
      
      triplet_edge_node2 = heterodimer_ilp_edges %>%
        transmute(i = heterodimer_ilp_edge_id * 3 - 1 + max(triplet_edge$i),
                  j = ilp_nodes2,
                  v = -1)
      
      triplet_edge_node_link = heterodimer_ilp_edges %>%
        transmute(i = heterodimer_ilp_edge_id * 3 + max(triplet_edge$i),
                  j = ilp_nodes_link,
                  v = -1)
      
      triplet_edge_heterodimer = bind_rows(triplet_edge_edge1, 
                                           triplet_edge_edge2,
                                           triplet_edge_edge3,
                                           triplet_edge_node1,
                                           triplet_edge_node2,
                                           triplet_edge_node_link) 
      
      if(length(unique(triplet_edge_heterodimer$i))!=nrow_heterodimer_ilp_edges * 3){
        stop("Error in triplet_edge_heterodimer")
      }
      
      triplet_basic = bind_rows(triplet_basic,
                                triplet_edge_heterodimer)
    }
  }
  triplet_extended = triplet_basic
  # Extended constraint ####
  {
    # triplet_isotope ####
    ## constrain an isotope formula must come with an isotope edge
    {
      
      if(nrow(ilp_edges %>% filter(category == "Natural_abundance"))==0){
        triplet_isotope = NULL
        nrow_triplet_isotope = 0
        
      }else{
        ilp_edges_isotope = ilp_edges %>%
          filter(category == "Natural_abundance") %>%
          group_by(edge_id, ilp_nodes1) %>%
          mutate(n1 = n()) %>%
          ungroup() %>%
          group_by(edge_id, ilp_nodes2) %>%
          mutate(n2 = n()) %>%
          ungroup() %>%
          group_by(ilp_nodes2) %>%
          mutate(n3 = n()) %>%
          ungroup()
        ## Simple case - one parent ilp_node comes with one isotope ilp_node
        {
          ilp_edges_isotope_1_1 = ilp_edges_isotope %>%
            filter(n1 == 1 & n2 == 1 & n3==1)
          ilp_edges_isotope_1_2 = ilp_edges_isotope %>%
            filter(n1 == 1 & n2 == 1 & n3!=1)
          if(nrow(ilp_edges_isotope_1_1)==0){
            triplet_isotope_edge_1_1 = NULL
            triplet_isotope_node_1_1 = NULL
          }else{
            triplet_isotope_edge_1_1 = ilp_edges_isotope_1_1 %>%
              transmute(i = 1:nrow(.) + max(triplet_extended$i),
                        j = ilp_edge_id + max(triplet_node$j), # locate to the isotope edge in triplet_edge
                        v = 1)
            triplet_isotope_node_1_1 = ilp_edges_isotope_1_1 %>%
              transmute(i = 1:nrow(.) + max(triplet_extended$i),
                        j = ifelse(direction == 1, ilp_nodes2, ilp_nodes1),
                        v = -1)
          }
          if(nrow(ilp_edges_isotope_1_2)==0){
            triplet_isotope_edge_1_2 = NULL
            triplet_isotope_node_1_2 = NULL
          }else{
            ilp_edges_isotope_1_2 = ilp_edges_isotope_1_2 %>%
              filter(n3 != 1) %>%
              group_by(ilp_nodes2) %>%
              group_split()
            ilp_edges_isotope_1_2 = bind_rows(ilp_edges_isotope_1_2,
                                              .id = "isotope_id") %>%
              mutate(isotope_id = as.numeric(isotope_id))
            triplet_isotope_edge_1_2 = ilp_edges_isotope_1_2 %>%
              transmute(i = isotope_id + ifelse(is.null(triplet_isotope_edge_1_1), max(triplet_extended$i), max(triplet_isotope_edge_1_1$i)),
                        j = ilp_edge_id + max(triplet_node$j), # locate to the isotope edge in triplet_edge
                        v = 1)
            triplet_isotope_node_1_2 = ilp_edges_isotope_1_2 %>%
              transmute(i = isotope_id + ifelse(is.null(triplet_isotope_edge_1_1), max(triplet_extended$i), max(triplet_isotope_edge_1_1$i)),
                        j = ifelse(direction == 1, ilp_nodes2, ilp_nodes1),
                        v = -1) %>% distinct()
          }
          
          
          }
        
        ## When two ilp_node1 point to same isotope ilp_node2, two lines needs to be combined
        {
          ilp_edges_isotope_2_1 = ilp_edges_isotope %>%
            filter((n1 == 1 & n2 != 1 & direction == 1 & n3==1))
          
          ilp_edges_isotope_2_2 = ilp_edges_isotope %>%
            filter(n1 == 1 & n2 != 1 & direction == 1 & n3 != 1)
          
          # ilp_edges_isotope_2 = ilp_edges_isotope %>%
          #   filter((n1 == 1 & n2 != 1 & direction == 1) |
          #            (n2 == 1 & n1 != 1 & direction == -1))
          
          ## Try catching bug here
          # If ilp_edges_isotope_1 + ilp_edges_isotope_2 does not represent all ilp_edges_isotope
          if((nrow(ilp_edges_isotope_1_1) + nrow(ilp_edges_isotope_1_2) + nrow(ilp_edges_isotope_2_1) + nrow(ilp_edges_isotope_2_2)) != nrow(ilp_edges_isotope)){
            stop("Bugs in isotope triplex.")
          }
          
          if(nrow(ilp_edges_isotope_2_1)!=0){
            ilp_edges_isotope_2_1_1 = ilp_edges_isotope %>%
              filter(n1 != 1) %>%
              group_by(edge_id, ilp_nodes1) %>%
              group_split()
            ilp_edges_isotope_2_1_2 = ilp_edges_isotope %>%
              filter(n2 != 1) %>%
              group_by(edge_id, ilp_nodes2) %>%
              group_split()
            
            # It is assumed that isotope_id is counting the number of list
            ilp_edges_isotope_2_1 = bind_rows(ilp_edges_isotope_2_1_1,
                                              ilp_edges_isotope_2_1_2,
                                              .id = "isotope_id") %>%
              mutate(isotope_id = as.numeric(isotope_id))
            
            # isotope_id specify which row the triplex should be in
            triplet_isotope_edge_2_1 = ilp_edges_isotope_2_1 %>%
              transmute(i = isotope_id + ifelse((is.null(triplet_isotope_edge_1_1)&is.null(triplet_isotope_edge_1_2)),max(triplet_extended$i),max(c(triplet_isotope_edge_1_1$i, triplet_isotope_edge_1_2$i))),
                        j = ilp_edge_id + max(triplet_node$j),
                        v = 1)
            triplet_isotope_node_2_1 = ilp_edges_isotope_2_1 %>%
              transmute(i = isotope_id + ifelse((is.null(triplet_isotope_edge_1_1)&is.null(triplet_isotope_edge_1_2)),max(triplet_extended$i),max(c(triplet_isotope_edge_1_1$i, triplet_isotope_edge_1_2$i))),
                        j = ifelse(direction == 1, ilp_nodes2, ilp_nodes1),
                        v = -1) %>%
              distinct()
          }else{
            triplet_isotope_edge_2_1 = NULL
            triplet_isotope_node_2_1 = NULL
          }
          
          if(nrow(ilp_edges_isotope_2_2)!=0){
            ilp_edges_isotope_2_2 = ilp_edges_isotope_2_2 %>%
              filter(n3 != 1) %>%
              group_by(ilp_nodes2) %>%
              group_split()
            ilp_edges_isotope_2_2 = bind_rows(ilp_edges_isotope_2_2,
                                              .id = "isotope_id") %>%
              mutate(isotope_id = as.numeric(isotope_id))
            # isotope_id specify which row the triplex should be in
            triplet_isotope_edge_2_2 = ilp_edges_isotope_2_2 %>%
              transmute(i = isotope_id + ifelse((is.null(triplet_isotope_edge_1_1)&is.null(triplet_isotope_edge_1_2)&is.null(triplet_isotope_edge_2_1)),max(triplet_extended$i),max(c(triplet_isotope_edge_1_1$i, triplet_isotope_edge_1_2$i, triplet_isotope_edge_2_1$i))),
                        j = ilp_edge_id + max(triplet_node$j),
                        v = 1)
            triplet_isotope_node_2_2 = ilp_edges_isotope_2_2 %>%
              transmute(i = isotope_id + ifelse((is.null(triplet_isotope_edge_1_1)&is.null(triplet_isotope_edge_1_2)&is.null(triplet_isotope_edge_2_1)),max(triplet_extended$i),max(c(triplet_isotope_edge_1_1$i, triplet_isotope_edge_1_2$i, triplet_isotope_edge_2_1$i))),
                        j = ifelse(direction == 1, ilp_nodes2, ilp_nodes1),
                        v = -1) %>%
              distinct()
          }else{
            triplet_isotope_edge_2_2 = NULL
            triplet_isotope_node_2_2 = NULL
          }
          
        }
        
        ## Output
        {
          triplet_isotope = bind_rows(triplet_isotope_edge_1_1,
                                      triplet_isotope_node_1_1,
                                      triplet_isotope_edge_1_2,
                                      triplet_isotope_node_1_2,
                                      triplet_isotope_edge_2_1,
                                      triplet_isotope_node_2_1,
                                      triplet_isotope_edge_2_2,
                                      triplet_isotope_node_2_2)
          
          # Sanity check #
          nrow_triplet_isotope_1_1 = nrow(ilp_edges_isotope_1_1)
          if(nrow(ilp_edges_isotope_1_2)!=0){
            nrow_triplet_isotope_1_2 = max(ilp_edges_isotope_1_2$isotope_id)
          }else{
            nrow_triplet_isotope_1_2 = 0
          }
          
          if(nrow(ilp_edges_isotope_2_1)!=0){
            nrow_triplet_isotope_2_1 = max(ilp_edges_isotope_2_1$isotope_id)
          }else{
            nrow_triplet_isotope_2_1 = 0
          }
          if(nrow(ilp_edges_isotope_2_2)!=0){
            nrow_triplet_isotope_2_2 = max(ilp_edges_isotope_2_2$isotope_id)
          }else{
            nrow_triplet_isotope_2_2 = 0
          }
          
          nrow_triplet_isotope = nrow_triplet_isotope_1_1 + nrow_triplet_isotope_1_2 + nrow_triplet_isotope_2_1 + nrow_triplet_isotope_2_2
          if(nrow_triplet_isotope != length(unique(triplet_isotope$i))){
            stop("error in triplet_isotope")
          }
          
        }
        triplet_extended = bind_rows(triplet_extended, triplet_isotope)
      }
      
      
      
    }
    
    # triplet multiple_edges ####
    ## triplet that force only one edge exist between two ilp_nodes
    ## E.g. adduct vs fragment of NH3, H3PO4 adduct and heterodimer, etc
    ## e12 + e12' + ... <= n1, e12 + e12' +... <= n2
    {
      if(exist_heterodimer){
        multiple_edges = bind_rows(ilp_edges, heterodimer_ilp_edges %>% mutate(linktype = as.character(linktype)))
      } else {
        multiple_edges = ilp_edges
      }
      
      multiple_edges_exist = multiple_edges %>%
        filter(category != "Library_MS2_fragment") %>% # allow Library_MS2_fragment to form double edge
        group_by(ilp_nodes1, ilp_nodes2) %>%
        mutate(n_twonodes = n()) %>%
        filter(n_twonodes > 1)
      if(nrow(multiple_edges_exist)==0){
        triplet_edge_multiple_edges = NULL
        nrow_triplet_multiple_edges = 0
      }else{
        multiple_edges = multiple_edges_exist %>%
          arrange(ilp_nodes1, ilp_nodes2) %>%
          group_split() %>%
          bind_rows(.id = "multiple_edges_id") %>%
          mutate(multiple_edges_id = as.numeric(multiple_edges_id))
        
        
        triplet_multiple_edges_regular_edge1 = multiple_edges %>%
          filter(category != "Heterodimer") %>%
          transmute(i = multiple_edges_id * 2 - 1 + max(triplet_extended$i),
                    j = ilp_edge_id + max(triplet_node$j),
                    v = 1)
        
        triplet_multiple_edges_node1 = multiple_edges %>%
          transmute(i = multiple_edges_id * 2 - 1 + max(triplet_extended$i),
                    j = ilp_nodes1,
                    v = -1) %>%
          distinct() # Duplicate entries
        
        triplet_multiple_edges_regular_edge2 = multiple_edges %>%
          filter(category != "Heterodimer") %>%
          transmute(i = multiple_edges_id * 2 + max(triplet_extended$i),
                    j = ilp_edge_id + max(triplet_node$j),
                    v = 1)
        
        triplet_multiple_edges_node2 = multiple_edges %>%
          transmute(i = multiple_edges_id * 2 + max(triplet_extended$i),
                    j = ilp_nodes2,
                    v = -1) %>%
          distinct() # Duplicate entries
        # Aviod duplicated triplet is important. Otherwise, a bomb is given by Rstudio.
        
        triplet_edge_multiple_edges = bind_rows(triplet_multiple_edges_regular_edge1,
                                                triplet_multiple_edges_node1,
                                                triplet_multiple_edges_regular_edge2,
                                                triplet_multiple_edges_node2)
        if(exist_heterodimer){
          triplet_multiple_edges_heterodimer_edge1 = multiple_edges %>%
            filter(category == "Heterodimer") %>%
            transmute(i = multiple_edges_id * 2 - 1 + max(triplet_extended$i),
                      j = heterodimer_ilp_edge_id + max(triplet_edge$j),
                      v = 1)
          
          triplet_multiple_edges_heterodimer_edge2 = multiple_edges %>%
            filter(category == "Heterodimer") %>%
            transmute(i = multiple_edges_id * 2 + max(triplet_extended$i),
                      j = heterodimer_ilp_edge_id + max(triplet_edge$j),
                      v = 1)
          
          triplet_edge_multiple_edges = bind_rows(triplet_edge_multiple_edges,
                                                  triplet_multiple_edges_heterodimer_edge1,
                                                  triplet_multiple_edges_heterodimer_edge2)
        }
        
        
        
        nrow_triplet_multiple_edges = max(multiple_edges$multiple_edges_id)
        
        # Sanity check #
        if(nrow_triplet_multiple_edges * 2 != length(unique(triplet_edge_multiple_edges$i))){
          stop("error in triplet_multiple_edges")
        }
        
        
      }
      triplet_extended = bind_rows(triplet_extended,
                                   triplet_edge_multiple_edges)
    }
    
    
    # 
    # triplet_multiple_ilp_edges ####
    # triplet that there is one and only one edge between two ilp_nodes
    # a+b-e1-e2<=1
    {
      if(nrow(multiple_edges_exist)==0){
        triplet_multiple_ilp_edges = NULL
      }else{
        triplet_multiple_edges_node1 = multiple_edges %>%
          transmute(i = multiple_edges_id * 2 - 1 + max(triplet_extended$i),
                    j = ilp_nodes1,
                    v = 1) %>%
          distinct() # Duplicate entries
        
        
        
        triplet_multiple_edges_node2 = multiple_edges %>%
          transmute(i = multiple_edges_id * 2 - 1 + max(triplet_extended$i),
                    j = ilp_nodes2,
                    v = 1) %>%
          distinct() # Duplicate entries
        
        triplet_multiple_edges_regular_edge1 = multiple_edges %>%
          filter(category != "Heterodimer") %>%
          transmute(i = multiple_edges_id * 2 - 1 + max(triplet_extended$i),
                    j = ilp_edge_id + max(triplet_node$j),
                    v = -1)
        triplet_multiple_edges_regular_edge2 = multiple_edges %>%
          filter(category != "Heterodimer") %>%
          transmute(i = multiple_edges_id * 2 + max(triplet_extended$i),
                    j = ilp_edge_id + max(triplet_node$j),
                    v = -1)
        
        
        
        # Aviod duplicated triplet is important. Otherwise, a bomb is given by Rstudio.
        
        triplet_multiple_ilp_edges = bind_rows(triplet_multiple_edges_node1,
                                               triplet_multiple_edges_node2,
                                               triplet_multiple_edges_regular_edge1,
                                               triplet_multiple_edges_regular_edge2)
        
        
        
      }
      triplet_extended = bind_rows(triplet_extended,
                                   triplet_multiple_ilp_edges)
    }
  }
  
  # CPLEX solver parameter ####
  {
    # Generate sparse matrix on left hand side
    triplet_df = triplet_extended %>% distinct()
    if(nrow(triplet_df)!=nrow(triplet_extended)){
      stop("Duplicates in triplet_df")
    }
    
    # converts the triplet into matrix
    mat = slam::simple_triplet_matrix(i=triplet_df$i,
                                      j=triplet_df$j,
                                      v=triplet_df$v)
    
    # CPLEX parameter
    
    nc <- max(mat$j)
    obj <- c(ilp_nodes$cplex_score, 
             ilp_edges$cplex_score,
             heterodimer_ilp_edges$cplex_score)
    # Sanity check
    if(nc != length(obj)){
      stop("objective length does not match number of column")
    }
    lb <- rep(0, nc)
    ub <- rep(1, nc)
    ctype <- rep(type,nc)
    
    nr <- max(mat$i)
    # nrow_triplet_isotope = 0
    # nrow_triplet_double_edges = 0
    rhs = c(rep(1, max(ilp_rows)), # sum(ilp_node_id)=1
            rep(0, nrow(ilp_edges) * 2), # e12-n1<=0;e12-n2<=0
            rep(0, nrow_heterodimer_ilp_edges * 3), # e123-n1<=0;e123-n2<=0;e123-n3<=0
            rep(0, nrow_triplet_isotope), # e12-n_isotope=0
            rep(0, nrow_triplet_multiple_edges * 2), # e12+e12'+...-n1<=0;e12+e12'+...-n2<=0
            rep(1, nrow_triplet_multiple_edges * 2))
    
    
    sense <- c(rep("E", max(ilp_rows)), 
               rep("L", nrow(ilp_edges) * 2),
               rep("L", nrow_heterodimer_ilp_edges * 3),
               rep("E", nrow_triplet_isotope),
               rep("L", nrow_triplet_multiple_edges * 2),
               rep("L", nrow_triplet_multiple_edges * 2))
    
    
    test = c(rep("E", 2), rep("L",0))
    # Sanity check
    if(nr != length(rhs) | nr != length(sense)){
      stop("constraint length does not match number of row")
    }
    
    triplet_df = triplet_df %>% arrange(j)
    cnt=as.vector(table(triplet_df$j))
    beg=vector()
    beg[1]=0
    for(i in 2:length(cnt)){beg[i]=beg[i-1]+cnt[i-1]}
    ind=triplet_df$i-1
    val = triplet_df$v
    
    CPX_MAX = -1
    CPLEX_para = list(nc = nc,
                      nr = nr,
                      CPX_MAX = CPX_MAX,
                      obj = obj,
                      rhs = rhs,
                      sense = sense,
                      beg = beg,
                      cnt = cnt,
                      ind = ind, 
                      val = val,
                      lb = lb,
                      ub = ub,
                      ctype = ctype,
                      mat = mat
                      # Without ctype the MIP may still run, but throw out error randomly
    )
  }
  
  return(CPLEX_para)
}


## Add_constraint_CPLEX ####
## Not functioning ##
Add_constraint_CPLEX = function(CplexSet, obj){
  
  # obj = obj_cplex
  env <- openEnvCPLEX()
  prob <- initProbCPLEX(env)
  
  nc = CplexSet$para$nc
  nr = CplexSet$para$nr
  CPX_MAX = CplexSet$para$CPX_MAX
  rhs = CplexSet$para$rhs
  sense = CplexSet$para$sense
  beg = CplexSet$para$beg
  cnt = CplexSet$para$cnt
  ind = CplexSet$para$ind
  val = CplexSet$para$val
  lb = CplexSet$para$lb
  ub = CplexSet$para$ub
  ctype = CplexSet$para$ctype
  
  copyLpwNamesCPLEX(env, prob, nc, nr, CPX_MAX, obj = obj, rhs, sense,
                    beg, cnt, ind, val, lb, ub, NULL, NULL, NULL)
  copyColTypeCPLEX(env, prob, ctype)
  
  
  # nnz = sum((unknown_formula$ILP_result!=0)==TRUE)
  # matbeg = 0
  # matval = rep(1,nnz)
  # matind = which(unknown_formula$ILP_result!=0)-1
  # addRowsCPLEX(env, prob, ncols=0, nrows=1, nnz=nnz, matbeg=matbeg, matind=matind, matval=matval,
  #              rhs = base::floor(nnz*.99), sense = "L",
  #              cnames = NULL, rnames = NULL)
  # addRowsCPLEX(env, prob, ncols=0, nrows=1, nnz=1, matbeg=0, matind=3909, matval=1,
  #              rhs = 1, sense = "E",
  #              cnames = NULL, rnames = NULL)
  # delRowsCPLEX(env, prob, begin = nr, end = getNumRowsCPLEX(env, prob)-1)
  # getNumRowsCPLEX(env, prob)
  
  # addMIPstartsCPLEX(env, prob, mcnt = 1, nzcnt = nc, beg = 0, varindices = 1:nc,
  #                   values = CplexSet$Init_solution2$CPLEX_x, effortlevel = 1, mipstartname = NULL)
  # 
  
  
  tictoc::tic()
  return_code = mipoptCPLEX(env, prob)
  result_solution=solutionCPLEX(env, prob)
  
  print(paste(return_codeCPLEX(return_code),"-",
              status_codeCPLEX(env, getStatCPLEX(env, prob)),
              " - OBJ_value =", result_solution$objval))
  
  tictoc::toc()
  
  # writeProbCPLEX(env, prob, "prob.lp")
  delProbCPLEX(env, prob)
  closeEnvCPLEX(env)
  return(list(obj = obj, result_solution = result_solution))
}
# Run_cplex ####
Run_cplex = function(CplexSet, obj_cplex, para_option = "para",
                     type = "C", # "ilp" enforces integer in assignment, "lp" does not.
                     relative_gap = 1e-4, total_run_time = 3000){
  optimization = ifelse(type=="B", "ilp", "lp")
  # obj_cplex = CplexSet$para$obj
  env <- openEnvCPLEX()
  prob <- initProbCPLEX(env)
  
  para = CplexSet[[para_option]]
  nc = para$nc
  nr = para$nr
  CPX_MAX = para$CPX_MAX
  rhs = para$rhs
  sense = para$sense
  beg = para$beg
  cnt = para$cnt
  ind = para$ind
  val = para$val
  lb = para$lb
  ub = para$ub
  ctype = para$ctype
  
  copyLpwNamesCPLEX(env, prob, nc, nr, CPX_MAX, obj = obj_cplex, rhs, sense,
                    beg, cnt, ind, val, lb, ub, NULL, NULL, NULL)
  setDblParmCPLEX(env, CPX_PARAM_TILIM, total_run_time) # total run time
  
  
  tictoc::tic()
  if(optimization[1] == "ilp"){
    # setting ctype decides if it is MIP. lp needs to without ctype 
    copyColTypeCPLEX(env, prob, ctype)
    
    # Conserve memory true
    setIntParmCPLEX(env, CPX_PARAM_MEMORYEMPHASIS, CPX_ON)
    # Increase probing
    setIntParmCPLEX(env, CPX_PARAM_PROBE, 3)
    # Sets which continuous optimizer will be used to solve the initial relaxation of a MIP.
    setIntParmCPLEX(env, 2025, 3) # 0 is default. one test shows 3 may works better and 4 is ok.
    # Sets an absolute tolerance on the gap between the best integer objective and the objective of the best node remaining
    setDblParmCPLEX(env, 2009, relative_gap) # Setting relative seems better
    
    # Sets an absolute tolerance on the gap between the best integer objective and the objective of the best node remaining
    # setDblParmCPLEX(env, 2008, 1000) 
    
    # setIntParmCPLEX(env, CPX_PARAM_INTSOLLIM, 2)
    # setDefaultParmCPLEX(env)
    # getChgParmCPLEX(env)
    # Assess parameters
    # getParmNameCPLEX(env, 2025)
    
    return_code = mipoptCPLEX(env, prob)
    result_solution=solutionCPLEX(env, prob)
    # Access Relative Objective Gap for a MIP Optimization Description
    rel_gap_postrun = getMIPrelGapCPLEX(env, prob)
    # result_solution_info = solnInfoCPLEX(env, prob)
    
    print(paste(return_codeCPLEX(return_code),"-",
                status_codeCPLEX(env, getStatCPLEX(env, prob)),
                " - OBJ_value =", round(result_solution$objval, 2),
                "(bestobjective - bestinteger) / (1e-10 + |bestinteger|) =",
                signif(rel_gap_postrun,5)
    )
    )
  } else if(optimization[1] == "lp"){
    
    return_code = lpoptCPLEX(env, prob)
    result_solution=solutionCPLEX(env, prob)
    print(paste(return_codeCPLEX(return_code),"-",
                status_codeCPLEX(env, getStatCPLEX(env, prob)),
                " - OBJ_value =", round(result_solution$objval, 2))
    )
  }
  
  tictoc::toc()
  
  # writeProbCPLEX(env, prob, "prob.lp")
  delProbCPLEX(env, prob)
  closeEnvCPLEX(env)
  
  return(list(obj = obj_cplex, result_solution = result_solution))
}

# Add lp_solution to CplexSet ####
add_CPLEX_solution = function(CplexSet, solution = "ilp_solution", type = "C"){
  if(type=="C"){
    res0 = CplexSet$ilp_solution$result_solution$x
    cplex_score = c(CplexSet$ilp_nodes$cplex_score, CplexSet$ilp_edges$cplex_score)
    if(!is.null(CplexSet$heterodimer_ilp_edges)){
      cplex_score = c(cplex_score, CplexSet$heterodimer_ilp_edges$cplex_score)
    }
    nodes <- CplexSet[["ilp_nodes"]][["node_id"]]
    res <- res0[c(1:length(nodes))]
    node_score = cplex_score[c(1:length(nodes))]
    
    res1 <- unlist(sapply(unique(nodes), FUN = function(x) {
      thisnode <- round(res[nodes == x],3)
      thisnode_score <- node_score[nodes == x]
      sol_node <- which(thisnode==max(thisnode))
      if(length(sol_node)>1){
        sol_node <- sol_node[which.max(thisnode_score[sol_node])]
      }
      vec <- rep(0, length(thisnode))
      vec[sol_node] <- 1
      return(vec)
    }))
    res0[c(1:length(nodes))] <- res1
    res0[res0 < 0.5] <- 0
    res0[res0 >= 0.5] <- 1
    CPLEX_x <- res0
  }else if(type=="B"){
    CPLEX_x = CplexSet[[solution]]$result_solution$x
  }else{
    warning("type must be B or C")
    return(CplexSet)
  }
  
  
  if(is.null(CPLEX_x)){
    warning(paste(solution, "does not exist."))
    return(CplexSet)
  }
  
  CplexSet$ilp_nodes[solution] = CPLEX_x[1:nrow(CplexSet$ilp_nodes)]
  CplexSet$ilp_nodes = CplexSet$ilp_nodes %>%
    dplyr::select(ilp_node_id, everything())
  
  CplexSet$ilp_edges[solution] = CPLEX_x[1:nrow(CplexSet$ilp_edges) + nrow(CplexSet$ilp_nodes)]
  CplexSet$ilp_edges = CplexSet$ilp_edges %>%
    dplyr::select(ilp_nodes1, ilp_nodes2, everything())
  
  if(!is.null(CplexSet$heterodimer_ilp_edges)){
    CplexSet$heterodimer_ilp_edges[solution] = CPLEX_x[1:nrow(CplexSet$heterodimer_ilp_edges) + 
                                                         nrow(CplexSet$ilp_nodes) + 
                                                         nrow(CplexSet$ilp_edges)]
  }
  
  
  
  return(CplexSet)
}

Run_SYM = function(CplexSet, total_run_time){
  mat <- CplexSet$para$mat;
  max <- TRUE;
  obj <- CplexSet[["para"]][["obj"]];
  sense <- CplexSet[["para"]][["sense"]];
  sense[sense == "E"] <- "==";
  sense[sense == "L"] <- "<=";
  dir <- sense;
  rhs <- CplexSet[["para"]][["rhs"]];
  
  # NOTE: "I" is integer, "B" is binary,
  # both of them are running too slowly (>1h)
  # for huge constraints matrix
  types <- "C"
  OptiSolution <- lpsymphony::lpsymphony_solve_LP(obj,
                                                  mat,
                                                  dir,
                                                  rhs,
                                                  types = types,
                                                  max = max,
                                                  verbosity = -1, gap_limit = 1e-3, time_limit = total_run_time);
}

## LPsymphony optimization
add_SYM_solution <- function(CplexSet, solution) {
  solution = CplexSet$ilp_solution
  nodes <- CplexSet[["ilp_nodes"]][["node_id"]]
  res0 <- solution[["solution"]]
  res <- res0[c(1:length(nodes))]
  
  res1 <- unlist(sapply(unique(nodes), FUN = function(x) {
    thisnode <- round(res[nodes == x],3)
    sol_node <- which.max(thisnode)
    vec <- rep(0, length(thisnode))
    vec[sol_node] <- 1
    return(vec)
  }))
  res0[c(1:length(nodes))] <- res1
  res0[res0 < 0.5] <- 0
  res0[res0 >= 0.5] <- 1
  CPLEX_x <- res0
  #CPLEX_x <- solution$solution
  #CPLEX_x[CPLEX_x < 0.5] <- 0
  #CPLEX_x[CPLEX_x >= 0.5] <- 1
  
  #CPLEX_x = CplexSet[[solution]]$result_solution$x
  if(is.null(CPLEX_x)){
    warning(paste("solution is empty."))
    return(CplexSet)
  }
  
  CplexSet$ilp_nodes["ilp_solution"] = CPLEX_x[1:nrow(CplexSet$ilp_nodes)]
  CplexSet$ilp_nodes = CplexSet$ilp_nodes %>%
    dplyr::select(ilp_node_id, everything())
  
  CplexSet$ilp_edges["ilp_solution"] = CPLEX_x[1:nrow(CplexSet$ilp_edges) + nrow(CplexSet$ilp_nodes)]
  CplexSet$ilp_edges = CplexSet$ilp_edges %>%
    dplyr::select(ilp_nodes1, ilp_nodes2, everything())
  
  if(!is.null(CplexSet$heterodimer_ilp_edges)){
    CplexSet$heterodimer_ilp_edges["ilp_solution"] = CPLEX_x[1:nrow(CplexSet$heterodimer_ilp_edges) + 
                                                               nrow(CplexSet$ilp_nodes) + 
                                                               nrow(CplexSet$ilp_edges)]
  }
  CplexSet$CPLEX_x <- CPLEX_x
  return(CplexSet)
}

# ------- Annotation functions -------- ####
## core_annotate ####
core_annotate = function(CplexSet, StructureSet_df, LibrarySet, solution = solution){
  
  # Make annotation to core peaks
  
  core_annotation = suppressMessages(StructureSet_df %>%
                                       left_join(LibrarySet %>%
                                                   dplyr::select(library_id, name, note, origin, SMILES) %>% 
                                                   dplyr::rename(parent_id = library_id)) %>%
                                       left_join(CplexSet$ilp_nodes %>%
                                                   dplyr::select(node_id, formula, class, ilp_node_id, eval(solution))) %>%
                                       mutate(annotation = case_when(
                                         category == "Unknown" ~ "Unknown",
                                         steps == 0 & transform == "" ~ paste(note, name, parent_formula),
                                         steps == 0 & direction == 1 ~ paste(name, parent_formula, "+", transform),
                                         steps == 0 & direction == -1 ~ paste(name, parent_formula, "-", transform),
                                         steps %% 1 == 0 ~ paste("Peak", node_id, formula)
                                       )))
  
  # Score different annotations and rank for same ilp_nodes
  core_annotation_score = core_annotation %>%
    mutate(score_source = case_when(
      origin == "manual_library" ~ 0.5,
      origin == "known_library" ~ 0.5,
      status == "quantified" ~ 0.6,
      status == "detected" ~ 0.3,
      status == "expected" ~ 0,
      status == "predicted" ~ -0.5
    )) %>% 
    mutate(score_transform = ifelse(transform == "", 0, -0.7)) %>%
    mutate(rank_score = core_annotation %>% 
             dplyr::select(starts_with("score")) %>% 
             rowSums(na.rm = TRUE)) %>%
    arrange(ilp_node_id, -!!as.symbol(solution), -rank_score) 
  
  core_annotation_score %>% 
    dplyr::select(ilp_node_id, eval(solution), annotation, everything())
}

## initiate_g_met ####
initiate_g_met = function(CplexSet){
  ilp_nodes_met = CplexSet$ilp_nodes %>%
    filter(class %in% c("Metabolite", "Putative metabolite")) %>%
    dplyr::select(ilp_node_id, everything())
  
  ilp_edges_met = CplexSet$ilp_edges %>%
    filter(category == "Biotransform") %>%
    mutate(direction = ifelse(direction >= 0, 1, -1)) %>%
    dplyr::select(ilp_nodes1, ilp_nodes2, everything())
  
  g_met = graph_from_data_frame(ilp_edges_met, 
                                directed = TRUE, 
                                vertices = ilp_nodes_met)
  
}
## initiate_g_nonmet ####
initiate_g_nonmet = function(CplexSet, solution){
  ilp_edges_merge = CplexSet$ilp_edges
  if(!is.null(CplexSet$heterodimer_ilp_edges)){
    ilp_edges_merge = merge(CplexSet$ilp_edges, CplexSet$heterodimer_ilp_edges, all = TRUE) 
  }
  ilp_edges_nonmet = ilp_edges_merge %>%
    filter(category != "Biotransform") %>%
    dplyr::select(ilp_nodes1, ilp_nodes2, everything())
  
  ilp_edges_nonmet_reorder1 = ilp_edges_nonmet %>%
    filter(direction == 1) %>%
    dplyr::rename(from = ilp_nodes1, to = ilp_nodes2)
  ilp_edges_nonmet_reorder2 = ilp_edges_nonmet %>%
    filter(direction == -1) %>%
    dplyr::rename(from = ilp_nodes2, to = ilp_nodes1)
  ilp_edges_nonmet_directed = bind_rows(ilp_edges_nonmet_reorder1, ilp_edges_nonmet_reorder2)
  
  ilp_edges_weights = ilp_edges_nonmet_directed %>%
    mutate(edge_weight = 2 - !!as.symbol(solution) - 0.1*cplex_score) %>%
    mutate(edge_weight = ifelse(category == "Heterodimer", edge_weight + 1.5, edge_weight)) %>%
    arrange(edge_weight)
  
  g_nonmet = graph_from_data_frame(ilp_edges_weights, 
                                   directed = TRUE, 
                                   vertices = CplexSet$ilp_nodes)
  
  return(g_nonmet)
}

## initiate_met_dist_mat ####
initiate_met_dist_mat = function(g_met, CplexSet, core_annotation){
  
  met_annotation_unique = core_annotation %>%
    filter(steps == 0) %>%
    arrange(ilp_node_id, -rank_score) %>%
    filter(class %in% c("Metabolite")) %>%
    distinct(ilp_node_id, .keep_all = TRUE)
  
  # 2273 core_met * 15342 ilp_metabolites means 267 Mb ~ 35MB
  core_met_names = met_annotation_unique %>% pull(ilp_node_id) %>% as.character()
  
  # Keep order of from and to but treated as undirected graph
  distMatrix <- shortest.paths(g_met, 
                               v=core_met_names,
                               to=V(g_met),
                               mode = "all") 
  
  distMatrix[distMatrix==0] = Inf # So that shortest distance prevent finding itself
  
  return(distMatrix)
}

## initiate_nonmet_dist_mat ####
initiate_nonmet_dist_mat = function(g_nonmet, CplexSet, core_annotation, 
                                    solution,
                                    only_solution_result = TRUE){
  
  ilp_nodes_nonmet = CplexSet$ilp_nodes %>%
    filter(class == "Artifact") %>%
    arrange(-!!as.symbol(solution))
  
  if(only_solution_result){
    ilp_nodes_nonmet = ilp_nodes_nonmet %>%
      filter(!!as.symbol(solution) > 0.01)
  }
  
  target_names = ilp_nodes_nonmet %>% pull(ilp_node_id) %>% as.character()
  
  core_names = core_annotation %>%
    filter(category != "Unknown") %>%
    filter(steps %% 1 == 0) %>% 
    distinct(ilp_node_id, .keep_all=TRUE) %>%
    pull(ilp_node_id) %>% as.character()
  
  g_edges = igraph::as_data_frame(g_nonmet, "edges") 
  
  # core_names = "61508"
  # target_names = "28"
  distMatrix <- shortest.paths(g_nonmet, 
                              v=core_names,
                              to=target_names, 
                              weights = g_edges$edge_weight,
                              mode = c("out"))
  distMatrix[distMatrix==0] = Inf 

  # An idea of compression 
  # Calculate sparse matrix on the run (Compression sparse column/row format)
  # Make sure each column can be extracted conviniently 
  
  return(distMatrix)
}
## all_network ####
all_network = function(g_met, g_nonmet){
  
  g_met2_node = igraph::as_data_frame(g_met, "vertices") 
  
  g_met2_edges = igraph::as_data_frame(g_met, "edges") 
  
  g_met2 = graph_from_data_frame(g_met2_edges,
                                 vertices = g_met2_node,
                                 directed = FALSE)
  
  if(!is.null(g_nonmet)){
    g_nonmet2_node = igraph::as_data_frame(g_nonmet, "vertices") 
    
    g_nonmet2_edges = igraph::as_data_frame(g_nonmet, "edges") 
    
    g_nonmet2 = graph_from_data_frame(g_nonmet2_edges,
                                      vertices = g_nonmet2_node,
                                      directed = FALSE)
  } else {
    g_nonmet2_node = NULL
    g_nonmet2_edges = NULL
  }
    g_all_nodes = bind_rows(g_met2_node, g_nonmet2_node) %>%
      distinct(name, .keep_all = TRUE) %>%
      filter(class != "Unknown")
    
    g_all_edges = bind_rows(g_met2_edges, g_nonmet2_edges) %>%
      distinct() %>%
      filter(from %in% g_all_nodes$name & to %in% g_all_nodes$name)
    
    g_all = graph_from_data_frame(g_all_edges,
                                  directed = FALSE,
                                  vertices = g_all_nodes)
  
}
## valid_network ####
valid_network = function(g_all){
  
  g_all_nodes = igraph::as_data_frame(g_all, "vertices")
  g_all_edges = igraph::as_data_frame(g_all, "edges")
  
  # Retain only ilp_solution nodes
  {
    g_all_nodes = g_all_nodes %>%
      filter(ilp_solution > 0.01)
    
    g_all_edges = g_all_edges %>%
      filter(from %in% g_all_nodes$name, to %in% g_all_nodes$name)
  }
  
  # Problem 1 - many single node in the output
  # Reason - Heterodimer edge score enables nodes without regular connections
  # Solutions - Remove nodes do not have regular connections
  {
    g_all_nodes = g_all_nodes %>%
      filter(name %in% g_all_edges$from | name %in% g_all_edges$to)
    
    g_all = graph_from_data_frame(g_all_edges,
                                  directed = FALSE,
                                  vertices = g_all_nodes)
  }
  
  # Problem 2 - Putative metabolites in network without connecting to knwon metabolites
  # Reason - Edge score by connection with isotope/adduct etc.
  # Solutions - Only retain subnetwork where step=0 annotation exist
  {
    
    clu=components(g_all)
    #subnetwork criteria 
    mainnetwork = igraph::groups(clu)[table(clu$membership) == max(table(clu$membership))]
    subnetwork = igraph::groups(clu)[table(clu$membership) < length(mainnetwork[[1]])]
    
    step0 = g_all_nodes %>%
      filter(steps == 0) %>%
      pull(name)
    subnetwork_valid = sapply(subnetwork, function(x){
      any(x %in% step0)
    })
    nodes_invalid = unlist(subnetwork[!subnetwork_valid])
    
    g_all_nodes_valid = g_all_nodes %>%
      filter(!name %in% nodes_invalid)
    g_all_edges_valid = g_all_edges %>%
      filter(from %in% g_all_nodes_valid$name & to %in% g_all_nodes_valid$name)
    g_all_valid = graph_from_data_frame(g_all_edges_valid,
                                        directed = FALSE,
                                        vertices = g_all_nodes_valid)
  }
  
  g_all_valid
}
# Initiate_networkset ####
Initiate_networkset = function(CplexSet, StructureSet_df, LibrarySet, 
                               solution = "ilp_solution"){
  core_annotation = core_annotate(CplexSet, StructureSet_df, LibrarySet, solution = solution)
  # prepare tracking graph for biotransform network and abiotic network
  # pre-calculate distance
  
  g_met = initiate_g_met(CplexSet)
  met_dist_mat = tryCatch(
    initiate_met_dist_mat(g_met, CplexSet, core_annotation), 
    error = function(e) { e })
  if (inherits(met_dist_mat, "simpleError")) {
    print("Insufficient memory for too many peaks. Only validated metabolites are annotated.")
    return(list())
  } 
  
  g_nonmet = initiate_g_nonmet(CplexSet, solution = solution)
  if(vcount(g_nonmet)<800000){
    nonmet_dist_mat = tryCatch(
      initiate_nonmet_dist_mat(g_nonmet, CplexSet, core_annotation,
                               solution = solution, only_solution_result = TRUE),
      error = function(e) { e })
    if (inherits(nonmet_dist_mat, "simpleError")) {
      print("Insufficient memory for too many peaks. Only metabolites are annotated.")
      g_nonmet = NULL
      nonmet_dist_mat = NULL
    } 
  }else{
    print("Insufficient memory for too many peaks. Only metabolites are annotated.")
    g_nonmet = NULL
    nonmet_dist_mat = NULL
  }
  
  
  g_all = all_network(g_met, g_nonmet)
  g_all_valid = valid_network(g_all)
  
  return(list(core_annotation = core_annotation,
              g_met = g_met,
              met_dist_mat = met_dist_mat,
              g_nonmet = g_nonmet, 
              nonmet_dist_mat = nonmet_dist_mat,
              g_all = g_all,
              g_all_valid = g_all_valid))
}

# path_annotation ####
path_annotation = function(CplexSet, NetworkSet, solution = "ilp_solution"){
  
  if(length(NetworkSet)==0){
    CplexSet$ilp_nodes = CplexSet$ilp_nodes %>%
      left_join(LibrarySet %>%
                  mutate(accession = note) %>%
                  filter(grepl('HMDB|PBCM', accession)) %>%
                  select(accession, name)) %>%
      mutate(path = ifelse(class=='Metabolite', paste(accession, name), NA))
    return(CplexSet)
  }
  
  # prepare tracking graph for biotransform network and abiotic network
  # pre-calculate distance
  {
    core_annotation = NetworkSet$core_annotation
    g_met = NetworkSet$g_met
  
    met_dist_mat = NetworkSet$met_dist_mat
    
    ilp_edges_annotate_met = igraph::as_data_frame(g_met)
    
    core_annotation_unique_met = core_annotation %>%
      filter(steps == 0) %>%
      arrange(ilp_node_id, -rank_score) %>%
      distinct(ilp_node_id, .keep_all = TRUE)
    
    # Annotate each ilp_nodes in for loop
    # Roughly 10s ~ 1000 entries(including skipped ones)
    ilp_edges_annotate_met$from <- as.integer(ilp_edges_annotate_met$from)
    ilp_edges_annotate_met$to <- as.integer(ilp_edges_annotate_met$to)
    
    ilp_nodes = CplexSet$ilp_nodes
    solution_index = as.integer(which(colnames(ilp_nodes)==solution)-1)
    class_index = as.integer(which(colnames(ilp_nodes)=="class")-1)
    
  }
  
  if(is.null(NetworkSet$g_nonmet)){
    path_annotation <- path_annotate_met_only(ilp_nodes = ilp_nodes,
                                              solution_index = solution_index,
                                              class_index = class_index,
                                              canu_met = core_annotation_unique_met,
                                              ilp_edges_anno_met = ilp_edges_annotate_met,
                                              dis_mat_met = met_dist_mat,
                                              g_annotation = g_met)
  } else {
    g_nonmet = NetworkSet$g_nonmet
    nonmet_dist_mat = NetworkSet$nonmet_dist_mat
    
    ilp_edges_annotate_nonmet = igraph::as_data_frame(g_nonmet)
    core_annotation_unique_nonmet = core_annotation %>%
      filter(steps %% 1 == 0) %>%
      arrange(ilp_node_id, -rank_score) %>%
      distinct(ilp_node_id, .keep_all = TRUE)
    
    ilp_edges_annotate_nonmet$from <- as.integer(ilp_edges_annotate_nonmet$from)
    ilp_edges_annotate_nonmet$to <- as.integer(ilp_edges_annotate_nonmet$to)
    
    path_annotation <- path_annotate(ilp_nodes = ilp_nodes,
                                     solution_index = solution_index,
                                     class_index = class_index,
                                     canu_met = core_annotation_unique_met,
                                     ilp_edges_anno_met = ilp_edges_annotate_met,
                                     dis_mat_met = met_dist_mat,
                                     g_annotation = g_met,
                                     canu_nonmet = core_annotation_unique_nonmet,
                                     ilp_edges_anno_nonmet = ilp_edges_annotate_nonmet,
                                     dis_mat_nonmet = nonmet_dist_mat,
                                     g_anno_non = g_nonmet)
  }
    
  CplexSet$ilp_nodes = ilp_nodes %>%
    mutate(path = path_annotation) %>%
    filter(T)

  return(CplexSet)
}



## get groupId for output ####
get_groupId = function(result,NodeSet){
  node_id = result$peak_id
  groupId = lapply(node_id, function(x){
    node = NodeSet[[as.numeric(x)]]
    pos_id = node$pos_info$pos_raw$id
    neg_id = node$neg_info$neg_raw$id
    if(is.null(pos_id)){
      pos_id=NA
    }
    if(is.null(neg_id)){
      neg_id = NA
    }
    id = data.frame(pos_id=pos_id,
                    neg_id=neg_id)
  })
  groupId = bind_rows(groupId)
}

# Output ####
get_NetID_output = function(CplexSet, NodeSet, mode, simplified = TRUE){
  H_mass = 1.00782503224
  e_mass = 0.00054857990943
  
  NetID_output = CplexSet$ilp_nodes %>% 
    filter(ilp_solution > 1e-6) 
  if(is.null(NetID_output$score_MS2)){
    NetID_output = mutate(NetID_output, score_MS2 = rep(NA, nrow(NetID_output)))
  }
  if(is.null(NetID_output$score_known_rt)){
    NetID_output = mutate(NetID_output, score_known_rt = rep(NA, nrow(NetID_output)))
  }
  NetID_output = NetID_output %>%
    mutate(ms2 = ifelse(is.na(score_MS2),0,score_MS2),
           rt = ifelse(is.na(score_known_rt),0,score_known_rt)) %>%
    mutate(level = case_when(
      class=="Metabolite" & rt * ms2 > 0 ~ 1,
      class=="Metabolite" & rt==0 & ms2!=0 ~ 2.1,
      class=="Metabolite" & rt!=0 & ms2==0 ~ 2.2,
      class=="Metabolite" & (rt + ms2) == 0 ~ 3,
      class=="Putative metabolite" ~ 4
    )) %>% dplyr::select(-c(rt,ms2))
  NetID_output = NetID_output %>% 
    mutate(mass = signif(mass, 7),
           medMz = signif(medMz, 7),
           medRt = round(medRt, 2),
           log10_inten = round(log10_inten, 2),
           ppm_error = (mass-medMz)/medMz*1e6,
           ppm_error = round(ppm_error, 2)) %>%
    dplyr::rename(peak_id = node_id,
                  annotation = path)
  id = get_groupId(NetID_output, NodeSet)
  if(mode==1){
    NetID_output = NetID_output %>% mutate(groupId=id$pos_id)
  }else if(mode==-1){
    NetID_output = NetID_output %>% mutate(groupId=id$neg_id)
  }
  if(simplified){
   
      NetID_output = NetID_output %>%
        dplyr::select(peak_id, groupId, medMz, medRt, log10_inten, class, formula, ppm_error, score_known_rt, score_MS2, level, accession, annotation) %>%
        mutate(medMz = medMz + (H_mass-e_mass)*mode)
    
  }else{
   
      NetID_output = NetID_output %>%
        dplyr::select(peak_id, groupId, medMz, medRt, log10_inten, class, formula, ppm_error, score_known_rt, score_MS2, level, accession, annotation, 
                      everything()) %>%
        mutate(medMz = medMz + (H_mass-e_mass)*mode)
  }
  
  NetID_output = NetID_output  %>%
    dplyr::rename('id' = 'groupId')
  
  NetID_output
}

# NetID_interface ####
NetID_interface = function(Related_files, Mset){
  sink(file = paste(MASTER_ADDRESS, "/log/NetID_log.txt", sep = ""), split=TRUE)
  printtime = Sys.time()
  timestamp = paste(unlist(regmatches(printtime, gregexpr("[[:digit:]]+", printtime))),collapse = '')
  
  cat(crayon::green("Setting up NodeSet and EdgeSet...\n"))
  # Setting up NodeSet and EdgeSet ####
  {
    NodeSet = Initiate_Nodeset(Mset)
    
    EdgeSet = Initiate_edgeset(Mset, Related_files, NodeSet, 
                               mz_tol_abs = Related_files$global_parameter$instrument_parameter$ppm * 20, 
                               mz_tol_ppm = Related_files$global_parameter$instrument_parameter$ppm, 
                               rt_tol_bio = Inf, rt_tol_nonbio = 0.2)
  } 
  
  cat(crayon::green("Setting up LibrarySet and StructureSet...\n"))
  # Setting up LibrarySet and StructureSet ####
  {
    LibrarySet = Initiate_Libraryset(Mset, Related_files)
    LibrarySet_expand = Expand_Libraryset(LibrarySet)
    
    StructureSet = Initilize_Empty_Structureset(NodeSet)
    StructureSet = Match_Library_Structureset(LibrarySet,
                                              LibrarySet_expand,
                                              StructureSet, 
                                              NodeSet,
                                              ppm_tol = Related_files$global_parameter$instrument_parameter$ppm)
  }
  
  cat(crayon::green("start EdgeSet_expand...\n"))
  # EdgeSet_expand ####
  {
    EdgeSet_expand = Expand_edgeset(Mset,
                                    Related_files,
                                    NodeSet,
                                    LibrarySet,
                                    StructureSet,
                                    EdgeSet,
                                    RT_cutoff = 0.2, inten_cutoff = Related_files$global_parameter$instrument_parameter$edge_expand_inten_cutoff,
                                    types = c("ring_artifact",
                                              "oligomer_multicharge",
                                              # "heterodimer",
                                              "experiment_MS2_fragment",
                                              "library_MS2_fragment"
                                    ))
    
    EdgeSet_all = merge_edgeset(EdgeSet, EdgeSet_expand)
  }
  
  cat(crayon::green("Candidate formula pool propagation and scoring...\n"))
  # Candidate formula pool propagation and scoring ####
  {
    StructureSet = Propagate_structureset(Mset, 
                                          NodeSet,
                                          StructureSet,
                                          EdgeSet_all,
                                          biotransform_step = 2,
                                          artifact_step = 3,
                                          propagation_ppm_threshold = Related_files$global_parameter$instrument_parameter$propagation_ppm,
                                          propagation_abs_threshold = 0e-4,
                                          record_RT_tol = 0.2,
                                          record_ppm_tol = Related_files$global_parameter$instrument_parameter$record_ppm)
    
    StructureSet_df = Score_structureset(Related_files, Mset, StructureSet, NodeSet, LibrarySet, EdgeSet_all)
    StructureSet_df = Score_structure_propagation(StructureSet_df, artifact_decay = -0.5)
  }
  
  cat(crayon::green("setting up CplexSet...\n"))
  # CplexSet and ilp_nodes and ilp_edges #### 
  {
    CplexSet = list()
    
    # Initialize
    CplexSet[["ilp_nodes"]] = initiate_ilp_nodes(StructureSet_df, NodeSet)
    CplexSet[["ilp_edges"]] = initiate_ilp_edges(EdgeSet_all, CplexSet, Exclude = "")
    CplexSet[["heterodimer_ilp_edges"]] = initiate_heterodimer_ilp_edges(EdgeSet_all, CplexSet, NodeSet)
    print("Finish CplexSet initialization")
    
    # Score
    CplexSet[["ilp_nodes"]] = score_ilp_nodes(CplexSet, 
                                              metabolite_score = 0.1,
                                              putative_metabolite_score = 0, 
                                              artifact_score = 0, 
                                              unknown_score = -0.5)
    
    CplexSet[["ilp_edges"]] = score_ilp_edges(Related_files, EdgeSet_all, CplexSet, NodeSet, Mset)
    
    
    CplexSet[["heterodimer_ilp_edges"]] = score_heterodimer_ilp_edges(CplexSet,
                                                                      type_score_heterodimer = 0,
                                                                      MS2_score_experiment_fragment = 0.5)
    print("Finish CplexSet scoring")
    
    
    
    
    CplexSet[["para"]] = Initiate_cplexset(CplexSet, type="C")
    # CplexSet[["para_reduce"]] = Initiate_cplexset_reduce(CplexSet)
    
    print(paste("Complexity is", CplexSet$para$nc, "variables and", CplexSet$para$nr, "constraints."))
    print(sapply(CplexSet$para, length))
    # print(paste("Complexity is", CplexSet$para_reduce$nc, "variables and", CplexSet$para_reduce$nr, "constraints."))
    # print(sapply(CplexSet$para_reduce, length))
  }
  
  cat(crayon::green("Run optimization...\n"))
  # Run optimization ####
  {
    
    CplexSet[["ilp_solution"]] = Run_cplex(CplexSet, obj_cplex = CplexSet$para$obj,
                                           optimization = c("lp"),
                                           relative_gap = 1e-3, total_run_time = 500000)
    CplexSet = add_CPLEX_solution(CplexSet,
                                  solution = "ilp_solution", type="C")
    # CplexSet[["ilp_solution"]] = Run_SYM(CplexSet, total_run_time = -1)
    # CplexSet = add_SYM_solution(CplexSet)
    # 
    #CplexSet$ilp_nodes = CplexSet$ilp_nodes %>% mutate(ilp_solution=ifelse(ilp_solution>=0.9999, 1,0))
    
    
  }
  
  cat(crayon::green("path annotation...\n"))
  # Path annotation ####
  {
    time5 = Sys.time()
    NetworkSet = list()
    NetworkSet = Initiate_networkset(CplexSet, StructureSet_df, LibrarySet,
                                     solution = "ilp_solution")
    
    CplexSet = path_annotation(CplexSet, NetworkSet, solution = "ilp_solution")
    
    save(NodeSet, EdgeSet_all, LibrarySet, StructureSet_df, CplexSet, NetworkSet, file = paste(MASTER_ADDRESS, "/NetID_output.RData", sep = ""))
    print(Sys.time()-time5)
  }
  
  # Output ####
  {
    NetID_output =  get_NetID_output(CplexSet, NodeSet, Related_files$global_parameter$mode, simplified = TRUE)
  }
  cat(crayon::green("NetID total run time:"))
  print(Sys.time()-printtime)
  sink()
  return(NetID_output)
}
# ---- End ---- ####