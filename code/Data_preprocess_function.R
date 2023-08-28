# Data preprocess functions
library(tidyverse)
library(openxlsx)
library(stringi)
library(Rcpp)

library(mzR)
library(xcms)
library(mclust)
library(tseries)

# mzmine format to ElMAVEN detail format  ####
mzmine2detail = function(mzmine){
  test = mzmine %>%
    rename_all(function(x){
      x %>% 
        str_remove_all('datafile:|.mzXML:height|.mzXML') %>%
        stri_replace_last_fixed('mz_range:max','mzmax') %>%
        stri_replace_last_fixed('mz_range:min','mzmin') %>%
        stri_replace_last_fixed('rt_range:max','rtmax') %>%
        stri_replace_last_fixed('rt_range:min','rtmin') %>%
        stri_replace_last_fixed('area','peakArea')
    }) %>%
    dplyr::rename(groupId = id) %>%
    select(everything(), -matches(":"), contains("mz"),contains("rt"),contains("peakArea")) %>%
    select(-c(2:10)) %>%
    gather(key="key",value="value",-groupId) %>%
    # filter(!is.na(value)) %>%
    separate(key,into=c("sample","key"),sep=":") %>%
    mutate(key = ifelse(is.na(key), "peakIntensity",key)) %>%
    spread(key, value) %>% dplyr::rename(peakMz = mz)
}
### Remove peaks that are below signal-blank ration ####
Remove_By_SBR <- function(em_info, sbr) {
  
  blk_filter = em_info[grepl("blank", em_info$sample), ]
  if(nrow(blk_filter != 0)){
    blk_site =blk_filter[1,]$groupId
    em_info_no_blk = em_info %>% filter(groupId < blk_site)
    em_info = em_info %>% filter(groupId>=blk_site)
    o_sample_names <- unique(em_info$sample)  
    blk_ind <- grep("blank", o_sample_names)
    ints <- matrix(em_info$peakIntensity, nrow = length(o_sample_names))
    blk_int <- apply(ints[blk_ind, ], 2, max)
    below_sbr <- ints[-blk_ind, ] > sbr * matrix(rep(blk_int, length(o_sample_names) - length(blk_ind)),
                                                 nrow = length(o_sample_names) - length(blk_ind),
                                                 byrow  = TRUE)
    inds <- which(!below_sbr, arr.ind = TRUE)
    index <- inds[, 1] + (inds[, 2] - 1) * length(o_sample_names)
    em_info$peakMz[index] = 0
    em_info$rt[index] = 0
    print(paste(length(index), "are removed by signal-blank ratio filter."))
    em_info = bind_rows(em_info_no_blk, em_info)
  }
  
  return(em_info)
}

### remove identical groups ####
Remove_Identical <- function(em_info, mode) {
  nrows <- sum(!grepl("blank", unique(em_info$sample)))
  samples <- length(unique(em_info$sample))
  # mz_thr <- (10 * 100 / 10^6) * sqrt(nrows)
  # rt_thr <- 0.001
  em_info = em_info %>% arrange(., sample) %>% arrange(., groupId)
  blk <- grepl("blank", em_info$sample)
  mz <- matrix(em_info$peakMz[!blk], nrow = nrows)
  rt <- matrix(em_info$rt[!blk], nrow = nrows)
  d_mz <- as.matrix(proxy::dist(mz, by_rows = FALSE)) == 0
  d_rt <- as.matrix(proxy::dist(rt, by_rows = FALSE)) == 0
  d <- d_mz & d_rt
  kill <- vector()
  # k_log <- matrix(ncol = 2)
  for( i in 1 : (ncol(d_mz) - 1)) {
    id = em_info$groupId[i * samples]
    tk <- (unname(which(d[(i + 1) : nrow(d), i])) + i) 
    if(length(tk)) {
      tk_id <- em_info$groupId[tk * samples]
      kill = unique(c(kill, tk_id))
      # k_log = rbind(k_log, matrix(c(rep(id, length(tk)), tk_id), ncol = 2))
    }
    
  }
  
  # k_log = k_log[2 : nrow(k_log), ] %>% cbind(., rep(mode, times = nrow(.)))
  # k_log = data.frame(remain = k_log[2 : nrow(k_log), 1],
  #                    killed = k_log[2 : nrow(k_log), 2],
  #                    mode = rep(mode, times = nrow(k_log) - 1))
  # base::colnames(k_log) = c("remain", "killed", "mode")
  print(paste(length(kill), "identical groups are removed."))
  # write_csv(k_log, paste(MASTER_ADDRESS, "/log/kill_log.csv", sep = ""), append = TRUE)
  # print(paste("kill_log written at ", MASTER_ADDRESS, "/log/kill_log.csv", sep = ""))
  em_info = filter(em_info, !(groupId %in% kill))
  return(em_info)
}




###### fuzzy_group ####
fuzzy_group = function(data, abs_tol = 0.002, rel_tol = 10e-6){
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
##### mz_RT_group ####
mz_RT_group = function(mzs, rts,
                       mz_abs_tol = 0.002, mz_rel_tol = 10e-6,
                       RT_abs_tol = 0.2, RT_rel_tol = 10e-6
)
{
  
  group1 = fuzzy_group(mzs,abs_tol = mz_abs_tol, rel_tol = mz_rel_tol) 
  
  group2_split = split(rts, group1)
  group2_ls = sapply(group2_split, fuzzy_group, RT_abs_tol, RT_rel_tol) %>% unlist()
  
  seq_split_by_group1 = split(1:length(mzs), group1) %>% unlist()
  group2 = numeric(length(mzs))
  group2[seq_split_by_group1] = group2_ls 
  
  result_group = split(1:length(mzs), paste(group1, group2))
  
  return(result_group)
}


#### Pair_Group find the peak group that has the same mz,rt ####
Pair_Group = function(em_info){
  nrows <- sum(!grepl("blank", unique(em_info$sample)))
  blk <- grepl("blank", em_info$sample)
  em_info = em_info[!blk,]
  mz = em_info %>% pull(peakMz)
  rt = em_info %>% pull(rt)
  
  group = mz_RT_group(mz, rt, mz_abs_tol = 0, mz_rel_tol = 0, 
                       RT_abs_tol = 0, RT_rel_tol = 0) %>% 
    lapply(.,t) %>% lapply(.,as.data.frame) %>% bind_rows() 
  if(ncol(group)==1){return(NULL)}
  else{
    group = group  %>% filter(V1<V2) %>%
      filter(mz[V1]!=0)
  }
 
  group = group[, apply(group, 2, function(y) any(!is.na(y)))]
  group = group %/% nrows + 1
  group_full = group %>% distinct(.,.keep_all = TRUE) %>%
    apply(.,1,function(x){
      x = na.omit(x)
      x = as.data.frame(t(combn(x,2)))
      return(x)
    }) %>% bind_rows()
  
  g = igraph::graph_from_data_frame(group_full, directed = F)
  clu = igraph::components(g)
  clu_group = igraph::groups(clu)
  clu_group = lapply(clu_group, as.numeric)
  return(clu_group)
}


##### calculate the sd of rnorm ####
Var_Calculator <- function(dm) {
  dm[which(!dm)] = NA
  # dst <- apply(dm, 2, proxy::dist) %>% "["(., !is.na(.))
  # dst = dst[which(dst != 0)] %>% -log10(.)
  # z <- (dst - mean(dst))/sd(dst)
  # qqnorm(z)
  # abline(0,1)
  var <- apply(dm, 2, function(x) stats::sd(x, na.rm = TRUE)) %>% 
    "["(., !is.na(.)) %>% "["(., !(!.)) %>% -log10(.)
  z<-(var - mean(var))/sd(var)
  qqnorm(z)
  abline(0,1)
  print(c(mean(var), sd(var)))
  # print(c(mean(dst), sd(dst)))
  return(c(mean(var), sd(var)))
}

##### mclust the mz,rt ####
Vote_by_Foot <- function(group, coeff_mz, coeff_rt, target_em, vis) {
  u_group <- unique(group[, c(2, 3, 5)])
  GENERATE = 50 # # how many points shall be generated for GMM
  if (GENERATE > 0) {
    # In order for the code results to be repeatable, set the random number seed
    # 2000 for mz, 1000 for rt, they can't be equal
    gmz <- sapply(u_group$mz, 
                  function(x){
                    set.seed(2000)
                    rnorm(n = GENERATE, mean = x, sd = 10^-coeff_mz)})
    
    grt <- sapply(u_group$rt, 
                  function(x){
                    set.seed(1000)
                    rnorm(n = GENERATE, mean = x, sd = 10^-coeff_rt)})
    g_group <- data.frame(mz = c(u_group$mz, as.vector(gmz)), rt = c(u_group$rt, as.vector(grt)))
  }
  mc <- Mclust(g_group)
  cls <- mc$classification[1 : nrow(u_group)]
  group$code = 0
  for (j in 1 : nrow(group)) {
    for (k in 1 : nrow(u_group)) {
      if (prod(group[j, 2:3] == u_group[k, 1 : 2]) != 0 && identical(group$sample[j], u_group$sample[k])) {
        group$group[j] = cls[k]
        group$code[j] = k
        break
      }
    }
  }
  # print(group)
  u_group = cbind(u_group, cls)
  split_code <- split(group$code, group$ids)
  split_group <- split(group$group, group$ids) %>% lapply(., function(x) {
    (table(c(x, 1 : mc$G)) - 1) / length(x)
  })
  u_group$gt <- 0
  for(i in 1 : nrow(u_group)) {
    ca <- sapply(split_code, function(x) i %in% x)
    candidates <- names(which(ca))
    rank <- sapply(candidates, function(x) {
      unname(split_group[[x]][as.character(cls[i])])
    }) %>% base::sort(., decreasing = TRUE)
    u_group$gt[i] = as.numeric(names(rank[1]))
  }
  if (vis) {
    title <- paste(unlist(strsplit(base::basename(target_em), "\\."))[1], 
                   paste0(unique(group$ids), collapse = "-"), sep = "_")
    jpeg(filename = paste(MASTER_ADDRESS, "/regroupphotos", "/", title, ".jpg", sep = ''))
    plot.Mclust(mc, "classification")
    title(title)
    gts <- unique(u_group$gt)
    smps <- unique(u_group$sample)
    nick <- sapply((1 : length(smps)) + 64, intToUtf8)
    palette("Dark 2")
    points(u_group, 
           col = sapply(u_group$gt, function(x) which(gts == x)) + 1, 
           pch = 18, cex = 1.5)
    legend("topleft", legend = gts, col = (1 : length(gts)) + 1, pch = 18, cex = 2)
    legend("bottomright", legend = paste(smps, nick, sep = "="))
    text(u_group,sapply(u_group$sample, function(x) nick[which(smps == x)]), 
         pos = sample(1 : 4, size = nrow(u_group), replace = TRUE), 
         cex = 1.5)
    dev.off()
    
    # plot.Mclust(mc, "BIC")
    # plot.Mclust(mc, "uncertainty")
    # plot.Mclust(mc, "density")
  }
  return(u_group)
}

#### regroup the peaks by mz,rt ####
Group_Builder <- function(group, mzm, rtm, target_em, samplenames, id_table, em_info, vis = TRUE) {
  if(is.null(group)){return(em_info)}
  patch <- em_info[1, ]
  coeff_mz <- Var_Calculator(mzm)
  coeff_rt <- Var_Calculator(rtm)
  c_mz <- coeff_mz[1] - 0 * coeff_mz[2]
  TRUST <- 2 # sd = miu - sigma * TRUST, miu & sigma are the gaussian with in a peakgroup among different samples, calculated from input csv data
  c_rt <- coeff_rt[1] - TRUST * coeff_rt[2]
  i <- 0
  smps <- length(unique(em_info$sample))
  for(g in group) {
    mz <- as.vector(mzm[, g])
    rt <- as.vector(rtm[, g])
    temp <- data.frame(ids = rep(id_table[g], each = nrow(mzm)),
                       mz = mz, rt = rt, group = rep(0, length(mz)), 
                       sample = rep(samplenames, times = length(g)))
    temp = temp[apply(temp[, 2:3], 1, prod) != 0, ]
    result <- Vote_by_Foot(temp, c_mz, c_rt, target_em, vis)
    sub_em <- em_info[as.vector(sapply(g, function(x) ((x-1) * smps + 1) : ((x-1) * smps + smps))), ]
    ind <- vector()
    for(j in 1 : nrow(result)) {
      ind <- c(ind,
               which((sub_em$groupId == result$gt[j]) & (sub_em$peakMz == result$mz[j]) & (sub_em$rt == result$rt[j]) & (sub_em$sample == result$sample[j])))
      
    }
    ind = sort(c(ind, which(grepl("blank", sub_em$sample))))
    sub_em$peakMz[setdiff(1 : nrow(sub_em), ind)] = 0
    sub_em$rt[setdiff(1 : nrow(sub_em), ind)] = 0
    
    i = i + 1
    patch = rbind(patch, sub_em)
  }
  returnval <- filter(em_info, !(groupId %in% id_table[(unlist(group))])) %>% 
    rbind(., patch[-1, ]) %>% arrange(., sample) %>% arrange(., groupId)
}


### regroup peaks by mz,rt using mclust ####
Regroup <- function(em_info, mode) {
  nrows <- sum(!grepl("blank", unique(em_info$sample)))
  samplenames <- unique(em_info$sample)[!grepl("blank", unique(em_info$sample))]
  em_info = em_info %>% arrange(., sample) %>% arrange(., groupId)
  unions = Pair_Group(em_info)
  blk <- grepl("blank", em_info$sample)
  mz <- matrix(em_info$peakMz[!blk], nrow = nrows)
  rt <- matrix(em_info$rt[!blk], nrow = nrows)
  # unions = Pair_Group(mz, rt)
  id_table <- unique(em_info$groupId)
  returnval <- Group_Builder(unions, mz, rt, target_em, samplenames, id_table, em_info)
  mode_label = case_when(
    mode==1 ~ "pos table",
    mode==-1 ~ "neg table"
  )
  print(paste(mode_label, "regroup finished."))
  return(returnval)
}

## data clean including Remove_By_SBR, Remove_Identical, Regroup ####
Data_Clean <- function(em_info,
                       mode = 1){
  
  o_sample_names <- unique(em_info$sample)  # this includes blank!!!
  sample_names <- o_sample_names %>% grepl("blank", .) %>% "!"() %>% "["(o_sample_names, .)
  # em_info <- Remove_Identical(o_em_info, mode) %>% Remove_By_SBR(., sbr = SBR)
  # sbr: signal-blank ratio for peak filtering
  em_info <- Remove_By_SBR(em_info, sbr = 2) 
  em_info <- Remove_Identical(em_info, mode)
  em_info <- Regroup(em_info, mode)
  return(em_info)
}

## check system error by using lm ####
Check_Sys_Error = function(peak_table, Related_files, ppm_tol = 5e-6, Rt_match = FALSE, mode){

  Metabolites_HMDB = Related_files$HMDB_library %>%
    dplyr::rename(note = accession) %>%
    mutate(origin = "HMDB_library")
  
  Metabolites_known = Related_files$known_library %>%
    mutate(category = "Metabolite") %>%
    dplyr::rename(note = accession) %>%
    mutate(origin = "known_library") 
  if(Related_files$global_parameter$LC_method %in% colnames(Related_files[["known_library"]])){
    Metabolites_known = Metabolites_known %>%
      mutate(rt = .[,eval(Related_files$global_parameter$LC_method)])
  }else{
    Rt_match = FALSE
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
  
  o_sample_names = unique(peak_table$sample)
  sample_names = o_sample_names %>% grepl("blank", .) %>% "!"() %>% "["(o_sample_names, .)
  
  
  for(i in 1:length(sample_names)){
    mass_original = peak_table[peak_table$sample == sample_names[i],]$peakMz
    H_mass = 1.00782503224
    e_mass = 0.00054857990943
    id_original = as.numeric(1 : length(mass_original))
    mass_rt_data = data.frame(mass_original, id_original) %>% arrange(mass_original)
    mass = mass_rt_data$mass_original
    id = mass_rt_data$id_original
    mass_match = lapply(1:length(mass_original), function(x){
      data.frame(
        node_id = x,
        formula = "Unknown",
        mass = mass_original[x], 
        rdbe = 0, 
        category = "Unknown", 
        parent_id = x, 
        parent_formula = "Unknown",
        transform = "",
        direction = 0,
        stringsAsFactors = FALSE
        
      )
    })
    temp_id = as.numeric(id)
    
    LibrarySet = bind_rows(Metabolites_HMDB, Metabolites_known, Manual, Adducts) %>%
      distinct(SMILES, formula, .keep_all = TRUE) %>%
      mutate(library_id = (1+length(mass_original)) : (length(mass_original)+nrow(.))) %>%
      filter(TRUE)
    
    
    lib = LibrarySet %>%
      mutate(node_id = library_id,
             formula = formula,
             mass =mass,
             parent_id = library_id,
             parent_formula = formula,
             transform = "",
             direction = 1,
             rdbe = rdbe,
             category = category) %>%
      dplyr::select(node_id, formula, mass, rdbe, parent_id, parent_formula, transform, direction, category) %>%
      filter(TRUE) 
    lib$mass = lib$mass + (H_mass-e_mass)*mode
    lib = lib %>% arrange(mass)
    lib_mass = lib$mass
    length_lib = length(lib_mass)
    
    j=j_min=j_max=1
    while(j <= length(mass)){
      
      mass_tol = mass[j] * ppm_tol
      while(lib_mass[j_min] < mass[j] - mass_tol & j_min < length_lib){
        j_min = j_min + 1
      }
      if(lib_mass[j_min] > mass[j] + mass_tol){
        j = j+1
        next
      }
      j_max = j_min
      
      while(lib_mass[j_max] - mass[j] < mass_tol & j_max < length_lib){
        j_max = j_max + 1
      }
      j_max = j_max - 1
      
      if(j_min > j_max){break}
      candidate = lib[j_min:j_max,] 
      
      candidate[["node_id"]] = temp_id[j]
      
      mass_match[[temp_id[j]]] = bind_rows(mass_match[[temp_id[j]]], candidate)
      j = j+1
    }
    if(length(mass_match) == 0){
      peak_table[peak_table$sample == sample_names[i],]$peakMz = mass_original
      
    }else{
      known_library = LibrarySet %>%
        filter(origin %in% c("known_library","manual_library"))
      if(nrow(known_library) == 0){
        known_library = LibrarySet
      }
      
      
      
      known_library_msr = bind_rows(mass_match) %>%
        filter(parent_id %in% known_library$library_id & transform == "") %>%
        mutate(msr_mass = mass_original[node_id], 
               mass_diff = mass - msr_mass,
               ppm_mass_diff = mass_diff/mass * 1e6) %>%
        filter(quantile(ppm_mass_diff, 0.1)<ppm_mass_diff,
               quantile(ppm_mass_diff, 0.9)>ppm_mass_diff) %>%
        distinct(node_id, .keep_all = TRUE) %>%
        filter(TRUE)
      
      if(Rt_match & (!is.null(LibrarySet$rt))){
        rt_original = peak_table[peak_table$sample == sample_names[i],]$rt
        library_RT = LibrarySet$rt
        names(library_RT) = LibrarySet$library_id
        
        known_library_msr = known_library_msr %>%
          mutate(msr_rt = rt_original[node_id],
                 lib_rt = library_RT[as.character(parent_id)]) %>%
            mutate(rt_match = abs(msr_rt - lib_rt) < 1) %>%
            arrange(-rt_match) %>%
            distinct(node_id, .keep_all = TRUE) %>%
            filter(rt_match)
         
      }
      # calculate the skewness of error, if skewness < 1, check sys error by using lm, else not
      # ppm_norm_test = shapiro.test(known_library_msr$ppm_mass_diff)
      
      # ppm_skewness = skewness(known_library_msr$ppm_mass_diff)
      ppm_jb_test = jarque.bera.test(known_library_msr$ppm_mass_diff)
      
      if(ppm_jb_test$p.value > 0.05 & nrow(known_library_msr) > 10){
        lsq_result = lm(known_library_msr$mass_diff~known_library_msr$msr_mass)
        
        
        ppm_adjust = as.numeric(lsq_result$coefficients[2]*10^6)
        abs_adjust = as.numeric(lsq_result$coefficients[1])
        
        fitdistData = fitdistrplus::fitdist(known_library_msr$ppm_mass_diff, "norm")
        plot(fitdistData)
        #print(fitdistData)
        mass_adjustment = abs(ppm_adjust + abs_adjust/250*1e6) > 0.2
        if(mass_adjustment){
          for(n in 1:length(mass_original)){
            if(mass_original[n] != 0){
              
              mass_original[n] = mass_original[n] * (ppm_adjust * 1e-6 + 1) + abs_adjust
            }
            
          }
          peak_table[peak_table$sample == sample_names[i],]$peakMz = mass_original
        }
        
      }
      
      
    }
    
    
  }
  return(peak_table)
  
}

# data clean interface ####
Data_Clean_Procedure = function(DataSet, 
                                Related_files){
  dir.create(paste(MASTER_ADDRESS, "/regroupphotos", sep = ""))
 
  input_mode = Related_files$global_parameter$mode # input_mode:  1(input only pos data) ;-1(input only neg data)
  if(input_mode == 1){
    DataSet[["pos_peak_table"]] = Data_Clean(Mset$pos_peak_table, mode = 1)
  }else{
    DataSet[["neg_peak_table"]] = Data_Clean(Mset$neg_peak_table, mode = -1)
  }
  return(DataSet)
}

# data correction and normalization ####
Data_Std = function(Mset, Related_files){
  H_mass = 1.00782503224
  e_mass = 0.00054857990943
  ppm_tol = Related_files$global_parameter$instrument_parameter$ppm
  input_mode = Related_files$global_parameter$mode
  if(input_mode>=0){
    Mset[["pos_peak_table"]] = Check_Sys_Error(Mset[["pos_peak_table"]], Related_files, ppm_tol = ppm_tol, Rt_match = FALSE, mode = 1) %>%
      mutate(peakMz = peakMz - (H_mass-e_mass)*1) %>% dplyr::rename(id = groupId) 
    Mset[["pos_raw_data"]] = Write_NetID_Input(Mset$pos_msdata, Mset$pos_peak_table, 1) %>% mutate(medMz = medMz - (H_mass-e_mass)*1)
    
  }
 
  if(input_mode<=0){
    Mset[["neg_peak_table"]] = Check_Sys_Error(Mset[["neg_peak_table"]], Related_files, ppm_tol = ppm_tol, Rt_match = FALSE, mode = -1) %>%
      mutate(peakMz = peakMz - (H_mass-e_mass)*(-1)) %>% dplyr::rename(id = groupId) 
    Mset[["neg_raw_data"]] = Write_NetID_Input(Mset$neg_msdata, Mset$neg_peak_table, -1) %>% mutate(medMz = medMz - (H_mass-e_mass)*(-1))
  }
  return(Mset)
}

#### find m data around max within n scans ####
Find_M_Data <- function(chr, m , n) {
  
  int <- intensity(chr)
  inds <- which(!is.na(int))
  dff <- diff(inds, m - 1, 1)
  t_dff <- inds[which(dff <= n - 1)]
  pk <- which(int == max(int, na.rm = TRUE))[1]
  if( sum((pk - t_dff) < n & (pk - t_dff) >= 0) > 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

### fill in missing values ####
Fill_Miss <- function(delt, chrs, thr, LENGTHOUT, PEAK_DATA_LENGTH, door) {
  if( is.null(dim(chrs))) {
    
    gate <- ((length(intensity(chrs)) - sum(is.na(intensity(chrs)))) >= thr) && Find_M_Data(chrs, PEAK_DATA_LENGTH[1], PEAK_DATA_LENGTH[2])
    if(gate | door) {
      # each chromatogram has 301 points, within 60s
      a = length(intensity(chrs))
      b = sum(is.na(intensity(chrs)))
      bins <- seq(min(rtime(chrs)), by = 2*delt/(LENGTHOUT-1), length.out = LENGTHOUT)
      m <- min(intensity(chrs), na.rm = TRUE)
      a_fit <- approx(intensity(chrs) ~ rtime(chrs),  xout = bins, yleft = m, yright = m, ties = "ordered")
      chrs@intensity <- a_fit$y
      chrs@rtime <- a_fit$x
    } else {
      chrs@intensity <- 0
      chrs@rtime <- 0
    }
    
    return(chrs)
  } 
  
  for( j in 1 : ncol(chrs)) {
    chr <- chrs[1, j]
    # plot(chr, col = 1, type = "o")
    # check whether the EIC is good enough
    gate <- ((length(intensity(chr)) - sum(is.na(intensity(chr)))) >= thr) && Find_M_Data(chr)
    gate=TRUE
    if(gate | door) {
      # bins <- seq(min(rtime(chr)), max(rtime(chr)), length.out = 10 * length(rtime(chr)))
      bins <- seq(min(rtime(chr)), by = 2*delt/(LENGTHOUT-1), length.out = LENGTHOUT)
      m <- min(intensity(chr), na.rm = TRUE)
      a_fit <- approx(intensity(chr) ~ rtime(chr),  xout = bins, yleft = m, yright = m, ties = "ordered")
      chrs[1, j]@intensity <- a_fit$y
      chrs[1, j]@rtime <- a_fit$x
      # lines(a_fit, col = 4, lty = 3)
    } else {
      chrs[1, j]@intensity <- 0
      chrs[1, j]@rtime <- 0
    }
    
  }
  return(chrs)
}


## get EICs details ####
Get_Raw_Peaks_N <- function(ppm, delt, em_info, mzXML_file_path, mode, abs = TRUE) {
  # neg <- Get_Raw_Peaks_N(ppm = Related_files$global_parameter$instrument_parameter$ppm, delt = 30, neg_em_info, 
  #                        target_raw = grep("^.*(neg).*(.mzXML)$", input_files, value = TRUE), -1, abs = FALSE)
  
  PEAK_DATA_LENGTH = c(5,30) # peaks with 5 data around max within 30 scans 
  if( abs) ppm = ppm * 100  # for absolute ppm
  else ppm = ppm # for relative ppm
  
  pf <- setClass("peak_info", slots = c(id = "numeric", Rt = "list", Mz = "list", 
                                        intensity_list = "list", rtime_list = "list" 
  ))
  
  
  sample_names <- names(mzXML_file_path)
  for(sample_ID_i in 1:length(sample_names)) {
    sample_ID = sample_names[sample_ID_i]
    t1 <- Sys.time()
    s_em_info <- filter(em_info, sample == sample_ID)
    s_em_info$rt <- s_em_info$rt * 60
    if(sample_ID == sample_names[1]) { beginner <- TRUE}
    else { beginner <- FALSE}
    if( beginner) {
      info <- pf(id = unique(em_info$groupId))
    }
    sample_file <- mzXML_file_path[sample_ID_i]
    rawdata <- readMSData(sample_file, msLevel = 1, mode = "onDisk")
    if( abs) mzr <- t(sapply(s_em_info$peakMz, "+", c(-ppm, ppm))) # for absolute ppm
    else mzr <- t(sapply(s_em_info$peakMz, "*", c(1 - ppm, 1 + ppm))) # for relative ppm
    
    rtr <- t(sapply(s_em_info$rt, "+", c(-delt, delt)))
    chrs <- chromatogram(rawdata, rt = rtr, mz = mzr)
    for(i in 1 : length(chrs)) {
      chr <- chrs[i]
      ints <- list() # its' component(ought to only have 1 component) is a list containing data for the peak, 
      rts <- list() # and the component is named by the sample_ID
      Mz <- vector()
      Rt <- vector()
      chr = Fill_Miss(delt, chr, thr = PEAK_DATA_LENGTH[1] - 1, LENGTHOUT = 301, PEAK_DATA_LENGTH, F)
      int <- unname(intensity(chr))
      rt <- unname(rtime(chr))
      if(!sum(int, na.rm = TRUE)) {
        ints[[sample_ID]] <- NA
        rts[[sample_ID]] <- NA
        Rt[sample_ID] <- NA
        Mz[sample_ID] <- NA
      } else {
        ints[[sample_ID]] <- int
        rts[[sample_ID]] <- rt
        Rt[sample_ID] <- s_em_info$rt[i]
        Mz[sample_ID] <- s_em_info$peakMz[i]
      }
      
      if( beginner) {
        info@intensity_list[[i]] <- ints # if the 1st sample, creat a list for this groupId
        info@rtime_list[[i]] <- rts
        info@Rt[[i]] <- Rt
        info@Mz[[i]] <- Mz
      } else {
        info@intensity_list[[i]] <- c(info@intensity_list[[i]], ints)
        info@rtime_list[[i]] <- c(info@rtime_list[[i]], rts)
        info@Rt[[i]] <- c(info@Rt[[i]], Rt)
        info@Mz[[i]] <- c(info@Mz[[i]], Mz)
      }
    }
    
    t2 <- Sys.time()
    print(paste(sample_ID, difftime(t2, t1, units = "mins"), "mins"))
  }
  # na_ind <- which(is.na(sapply(info@Mz, mean)))
  chrom_na = sapply(info@Mz,function(x){
    if(all(is.na(x))){
      return(NA)
    }else if(length(which(is.na(x)))<length(x) & length(which(is.na(x)))>0){
      return(1)
    }else{return(0)}
  })
  single_na = which(chrom_na==1)
  if(length(single_na>0)){
    for(chrom in single_na){
      na_index = which(is.na(info@Mz[[chrom]]))
      com_index = which(!is.na(info@Mz[[chrom]]))[1]
      info@Rt[[chrom]][na_index] = em_info[length(sample_names)*(chrom-1)+na_index, ]$rt * 60
      info@Mz[[chrom]][na_index] = em_info[length(sample_names)*(chrom-1)+na_index, ]$peakMz
      info@rtime_list[[chrom]][na_index] = info@rtime_list[[chrom]][com_index]
      info@intensity_list[[chrom]][na_index] = info@intensity_list[[chrom]][com_index]
    }
  }
  
  na_ind = which(is.na(chrom_na))
  miss <- data.frame(id = info@id[na_ind], mode = rep(mode, times = length(na_ind)))
  write_csv(miss, paste(MASTER_ADDRESS, "/log/miss_log.csv", sep = ""), append = TRUE)
  print(paste("miss_log file written at ", MASTER_ADDRESS, ".log/miss_log.csv", sep = ""))
  if( length(na_ind)) {
    info@id = info@id[-na_ind]
    info@Rt = info@Rt[-na_ind]
    info@Mz = info@Mz[-na_ind]
    info@intensity_list = info@intensity_list[-na_ind]
    info@rtime_list = info@rtime_list[-na_ind]
  }
  return(info)
} 

## merge pos and neg chromatogram data ####
Merge_2 <- function(pos1, pos2) {
  # d1 <- read.csv("/Users/zhoutianyu/Desktop/NetID data_N/pos/pos_1_full.csv")
  # d2 <- read.csv("/Users/zhoutianyu/Desktop/NetID data_N/pos/pos_2_full.csv")
  # d2$groupId = d2$groupId + max(d1$groupId)
  # d3 <- rbind(d1, d2)
  # write.csv(d3, "/Users/zhoutianyu/Desktop/NetID data_N/pos/pos_full.csv")
  pf <- setClass("peak_info", slots = c(id = "numeric", Rt = "list", Mz = "list", 
                                        intensity_list = "list", rtime_list = "list" 
  ))
  # id <- c(c(pos2@id, pos1@id + max(pos2@id)))
  id <- c(pos1@id, pos2@id)
  Rt <- c(pos1@Rt, pos2@Rt)
  Mz <- c(pos1@Mz, pos2@Mz)
  it <- c(pos1@intensity_list, pos2@intensity_list)
  rtime <- c(pos1@rtime_list, pos2@rtime_list)
  pos <- pf(id = id, Rt = Rt, Mz = Mz, intensity_list = it, rtime_list = rtime)
  return(pos)
}

# get EICs interface ####
Get_EICs = function(DataSet, Related_files, mzXML_file_path){
  
  
  input_mode = Related_files$global_parameter$mode
  pos_em_info = DataSet$pos_peak_table
  neg_em_info = DataSet$neg_peak_table
  if(input_mode == 1){
    pos <- Get_Raw_Peaks_N(ppm = Related_files$global_parameter$instrument_parameter$ppm, delt = 30, pos_em_info,
                           mzXML_file_path = mzXML_file_path, 1, abs = FALSE)
    DataSet[["pos_msdata"]] = pos
  }else if(input_mode == -1){
    neg <- Get_Raw_Peaks_N(ppm = Related_files$global_parameter$instrument_parameter$ppm, delt = 30, neg_em_info, 
                           mzXML_file_path = mzXML_file_path, -1, abs = FALSE)
    DataSet[["neg_msdata"]] = neg
  }
  
  return(DataSet)
}

