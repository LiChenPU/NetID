# Library ####
{
  library(tidyverse)
  library(Rcpp)
  library(rstudioapi)
  options(scipen = 999)
  
  code_dir = dirname(rstudioapi::getSourceEditorContext()$path)
  main_dir = dirname(code_dir)
  main_data_dir = paste(main_dir, "data", sep = .Platform$file.sep)
}

# Function ####
{
  source(paste0(main_dir, "/code/Data_preprocess_function.R"))
  source(paste0(main_dir, "/code/NetID_function.R"))
  Rcpp::sourceCpp(paste0(main_dir, "/code/src/FastNetID.cpp"))
}

# Parameter setup ####
{
  # put the peak table file directory here
  analyzed_folder = "QTOF_6600_demo"

  setwd(paste(main_data_dir, analyzed_folder, sep = .Platform$file.sep))
  analyzed_dir = getwd()
  
  LC_method = "Hilicon"
  
  #sample filter
  sample_filter = "293T"
  
  # ion mode
  mode = 1
  
  # ms2 option
  ms2_option = T
  
  # manual_library option
  manual_library_option = T
  
  # optimization type
  # "C": continuous linear programming - fast & reasonably accurate; 
  # "B": integer(binary) linear programming - slow & slightly better
  Ctype = "B"   
  
  # instrument_parameter
  MS_type = 'QTOF' # QTOF, Orbitrap

  instrument_parameter_type = "Default" # Custom or Default

  instrument_parameter_custom = data.frame(ppm = 10e-6,
                                    ms1_ms2_match_ppm = 20e-6,
                                    propagation_ppm = 10e-6,
                                    record_ppm = 10e-6,
                                    edge_expand_inten_cutoff = 1e3,
                                    score_mz_ppm_error_rate = -0.125)
}   

# File path setup ####
{
  mode_label = case_when(
    mode==1 ~ "pos",
    mode==-1 ~ "neg"
  )
  MASTER_ADDRESS <<- paste(analyzed_dir, mode_label, "NetID_output", sep = .Platform$file.sep)
  dir.create(MASTER_ADDRESS, showWarnings = FALSE)
}

# Loading files ####
{
  # dependency files
  {
    # Reference compound library
    str_lib_file = paste(main_dir, "dependence/hmdb_pubchemlite_merge_result_simple.rds", sep = .Platform$file.sep)
    
    # known library
    known_library_file = paste(main_dir, "dependence/known_library.xlsx", sep = .Platform$file.sep)
    
    # MS2 library
    neg_MS2_library_file = paste(main_dir, "dependence/HMDB_exp_MS2_neg.rds", sep = .Platform$file.sep)
    pos_MS2_library_file = paste(main_dir, "dependence/HMDB_exp_MS2_pos.rds", sep = .Platform$file.sep)
    
    # empirical file
    empirical_rule_file = paste(main_dir, "dependence/empirical_rules.csv", sep = .Platform$file.sep)
    
    # propagation file
    propagation_rule_file = paste(main_dir, "dependence/propagation_rule.csv", sep = .Platform$file.sep)
    
    # manual library file
    if(manual_library_option){
      manual_library_file = paste(main_dir, "dependence/manual_library.csv", sep = .Platform$file.sep)
    }else{
      manual_library_file = ""
    }
  }
  
  # raw data
  {
    # Peak tables
    
    mzmine_files = paste(analyzed_dir, mode_label, "raw_data/mzmine.csv", sep = .Platform$file.sep)
     
    peak_table_files = paste(analyzed_dir, mode_label, "peak_table.csv", sep = .Platform$file.sep)
     
    peak_tables = read_csv(peak_table_files)
    detail_peak_tables = read_csv(mzmine_files) %>% mzmine2detail() %>% filter(groupId %in% peak_tables$id)
    
    detail_peak_table_files = paste(MASTER_ADDRESS, "/", mode_label,"_detail_peak_table.csv", sep = "")
    write_csv(detail_peak_tables, file = detail_peak_table_files)
  
    # mzXML files
    untarget_raw_files = dir(path = paste(analyzed_dir, mode_label, "raw_data", sep = .Platform$file.sep), 
                             full.names = T, pattern = ".mzXML")
    untarget_raw_files = grepl(sample_filter, basename(untarget_raw_files)) %>% "["(untarget_raw_files, .)
    
    # sample names 
    sample_names = str_remove_all(basename(untarget_raw_files), '.mzXML')
    names(untarget_raw_files) = sample_names
    untarget_raw_files=untarget_raw_files[!grepl("blank|blk|bk", sample_names, ignore.case = T)]
    
    # ms2 files
    ms2_file_path = list(pos = paste(analyzed_dir, "pos/raw_data/MS2_Data", sep = .Platform$file.sep),
                     neg = paste(analyzed_dir, "neg/raw_data/MS2_Data", sep = .Platform$file.sep))
  }
}

# Main
# Read files for NetID ####
{
  # Read dependency files ####
  Related_files = list()
  Related_files = Read_Related_Files(mode, neg_MS2_library_file, pos_MS2_library_file,
                                     str_lib_file, known_library_file, LC_method,
                                     manual_library_file, empirical_rule_file, propagation_rule_file,
                                     MS_type, instrument_parameter_type, instrument_parameter_custom)
  
  # Read data files ####
  Mset = list()
  Mset = Read_Files(Mset,
                    Related_files,
                    detail_peak_table_files,
                    untarget_raw_files,
                    sample_filter)
    
  if(ms2_option){
    Mset = Read_MS2data(Mset, Related_files,
                        MS2_folder=ms2_file_path) 
  }
}

# Annotated by NetID ####
{
  sink(file = paste(MASTER_ADDRESS, "/log/NetID_log.txt", sep = ""), split=TRUE)
  printtime = Sys.time()
  timestamp = paste(unlist(regmatches(printtime, gregexpr("[[:digit:]]+", printtime))),collapse = '')
  
  # Setting up NodeSet and EdgeSet ####
  cat(crayon::green("Setting up NodeSet and EdgeSet...\n"))
  {
    NodeSet = Initiate_Nodeset(Mset)
    
    EdgeSet = Initiate_edgeset(Mset, Related_files, NodeSet, 
                               mz_tol_abs = Related_files$global_parameter$instrument_parameter$ppm * 20, 
                               mz_tol_ppm = Related_files$global_parameter$instrument_parameter$ppm, 
                               rt_tol_bio = Inf, rt_tol_nonbio = 0.2)
  } 
  
  # Setting up LibrarySet and StructureSet ####
  cat(crayon::green("Setting up LibrarySet and StructureSet...\n"))
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
  
  # EdgeSet_expand ####
  cat(crayon::green("start EdgeSet_expand...\n"))
  {
    EdgeSet_expand = Expand_edgeset(Mset,
                                    Related_files,
                                    NodeSet,
                                    LibrarySet,
                                    StructureSet,
                                    EdgeSet,
                                    RT_cutoff = 0.2, 
                                    inten_cutoff = Related_files$global_parameter$instrument_parameter$edge_expand_inten_cutoff,
                                    types = c("ring_artifact",
                                              "oligomer_multicharge",
                                              "experiment_MS2_fragment",
                                              "library_MS2_fragment"
                                    ))
    
    EdgeSet_all = merge_edgeset(EdgeSet, EdgeSet_expand)
  }

  # Candidate formula pool propagation and scoring ####
  cat(crayon::green("Candidate formula pool propagation and scoring...\n"))
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
    if(nrow(StructureSet_df)>5e6){
      warning("Too many candidates. Network optimization may be slow. 
              Consider changing parameters")
    }
  }
  
  # CplexSet and ilp_nodes and ilp_edges #### 
  cat(crayon::green("setting up CplexSet...\n"))
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
  
    CplexSet[["para"]] = Initiate_cplexset(CplexSet, type = Ctype)
    print(paste("Complexity is", CplexSet$para$nc, "variables and", CplexSet$para$nr, "constraints."))
    print(sapply(CplexSet$para, length))
  }
  
  # Run optimization ####
  cat(crayon::green("Run optimization...\n"))
  {
    CplexSet[["ilp_solution"]] = Run_cplex(CplexSet, obj_cplex = CplexSet$para$obj,
                                           type = Ctype,
                                           relative_gap = 1e-3, total_run_time = 500000)
    CplexSet = add_CPLEX_solution(CplexSet,
                                  solution = "ilp_solution", type=Ctype)
  }
  
  # Path annotation ####
  cat(crayon::green("path annotation...\n"))
  {
    time5 = Sys.time()
    NetworkSet = list()
    NetworkSet = Initiate_networkset(CplexSet, StructureSet_df, LibrarySet,
                                     solution = "ilp_solution")
    
    CplexSet = path_annotation(CplexSet, NetworkSet, solution = "ilp_solution")
    
    save(Related_files, Mset, NodeSet, EdgeSet_all, LibrarySet, StructureSet_df, CplexSet, NetworkSet, file = paste(MASTER_ADDRESS, "/NetID_output.RData", sep = ""))
    print(Sys.time()-time5)
  }
  
  # Output ####
  {
    NetID_output = get_NetID_output(CplexSet, NodeSet, 
                                    Related_files$global_parameter$mode, 
                                    simplified = TRUE)

    # Node and edge lists for cytoscape output
    cyto_nodes = CplexSet$ilp_nodes %>%
      filter(ilp_solution > 0.01)
    cyto_edges = CplexSet$ilp_edges %>%
      filter(ilp_solution > 0.01 | 
               (ilp_nodes1 %in% cyto_nodes$ilp_node_id &
                  ilp_nodes2 %in% cyto_nodes$ilp_node_id &
                  category == "Biotransform"))
  }
  cat(crayon::green("NetID total run time:"))
  print(Sys.time()-printtime)
  sink()
}

# output
write_csv(NetID_output, file = paste(MASTER_ADDRESS, "NetID_annotation.csv", sep = .Platform$file.sep), na="")
write_csv(cyto_nodes, file = paste(MASTER_ADDRESS, "cyto_nodes.csv", sep = .Platform$file.sep), na="")
write_csv(cyto_edges, file = paste(MASTER_ADDRESS, "cyto_edges.csv", sep = .Platform$file.sep), na="")


