# Main Codes ####
# Setting path ####
{
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  source("NetID_function.R")
  
  work_dir = "../Mouse_liver_neg/"
  setwd(work_dir)
  printtime = Sys.time()
  timestamp = paste(unlist(regmatches(printtime, gregexpr("[[:digit:]]+", printtime))),collapse = '')
}
# Read data and files ####
{
  Mset = list()
  # Read in files 
  Mset = read_files(filename = "raw_data.csv",
                    LC_method = "Hilic_25min_QE", # "Hilic_Rutgers_QEPlus" "Hilic_25min_QE", lipids is empty
                    ion_mode = -1 # 1 for pos mode and -1 for neg mode
                    )
  Mset = read_MS2data(Mset,
                      MS2_folder = "MS2_neg_200524") # MS2_neg_200524
}

# Data cleaning ####
{
  # Analyze cohort info, identify blank samples
  Mset[["Cohort"]]=Cohort_Info(Mset, first_sample_col_num = 15)
  print(Mset$Cohort)
  
  # Clean up duplicate peaks 
  Mset[["Data"]] = Peak_cleanup(Mset,
                                mz_tol=2/10^6, rt_tol=0.1, # merge peaks within mz_tol and rt_tol
                                inten_cutoff=0e4,
                                high_blank_cutoff = 0,
                                first_sample_col_num = 15)
  
  print(c(nrow(Mset$raw_data), nrow(Mset$Data)))
}

# Setting up NodeSet, LibrarySet and StructureSet ####
{
  NodeSet = Initiate_nodeset(Mset)
  
  LibrarySet = Initiate_libraryset(Mset)
  LibrarySet_expand = Expand_libraryset(LibrarySet)
  
  StructureSet = Initilize_empty_structureset(NodeSet)
  StructureSet = Match_library_structureset(LibrarySet,
                                            LibrarySet_expand,
                                            StructureSet, 
                                            NodeSet,
                                            ppm_tol = 10e-6)
}

# Adjust measured mass by matching to known metabolites ####
{
  measured_mz_adjust = T 
  if(measured_mz_adjust){
    Sys_msr_error = Check_sys_error(NodeSet, StructureSet, LibrarySet, 
                                    RT_match = F)
    
    mass_adjustment = abs(Sys_msr_error$ppm_adjust + Sys_msr_error$abs_adjust/250*1e6) > 0.2
    
    if(mass_adjustment){
      NodeSet = lapply(NodeSet, function(Node){
        Node$mz = Node$mz * (Sys_msr_error$ppm_adjust * 1e-6 + 1) + Sys_msr_error$abs_adjust
        return(Node)
      })
      StructureSet = Initilize_empty_structureset(NodeSet)
      StructureSet = Match_library_structureset(LibrarySet,
                                                LibrarySet_expand,
                                                StructureSet, 
                                                NodeSet,
                                                ppm_tol = 10e-6)
    }
  }
}

# Setting up EdgeSet ####
{
  EdgeSet = Initiate_edgeset(Mset, NodeSet, 
                             mz_tol_abs = 2e-4, mz_tol_ppm = 10, 
                             rt_tol_bio = Inf, rt_tol_nonbio = 0.2)
  
  EdgeSet_expand = Expand_edgeset(EdgeSet,
                                  RT_cutoff = 0.2, inten_cutoff = 1e4,
                                  types = c("ring_artifact",
                                            "oligomer_multicharge",
                                            "heterodimer",
                                            "experiment_MS2_fragment",
                                            "library_MS2_fragment"
                                  ))
  
  EdgeSet_all = merge_edgeset(EdgeSet, EdgeSet_expand)
}

# Candidate formula pool propagation and scoring ####
{
  StructureSet = Propagate_structureset(Mset, 
                                        NodeSet,
                                        StructureSet,
                                        EdgeSet_all,
                                        biotransform_step = 2,
                                        artifact_step = 3,
                                        propagation_ppm_threshold = 10e-6,
                                        propagation_abs_threshold = 0e-4,
                                        record_RT_tol = 0.2,
                                        record_ppm_tol = 5e-6)
  
  StructureSet_df = Score_structureset(StructureSet)
  StructureSet_df = Score_structure_propagation(StructureSet_df, artifact_decay = -0.5)   
}

# CplexSet and ilp_nodes and ilp_edges #### 
{
  CplexSet = list()
  
  # Initialize
  CplexSet[["ilp_nodes"]] = initiate_ilp_nodes(StructureSet_df)
  CplexSet[["ilp_edges"]] = initiate_ilp_edges(EdgeSet_all, CplexSet, Exclude = "")
  CplexSet[["heterodimer_ilp_edges"]] = initiate_heterodimer_ilp_edges(EdgeSet_all, CplexSet, NodeSet)
  print("Finish CplexSet initialization")
  
  # Score
  CplexSet[["ilp_nodes"]] = score_ilp_nodes(CplexSet, 
                                            metabolite_score = 0.1,
                                            putative_metabolite_score = 0, 
                                            artifact_score = 0, 
                                            unknown_score = -0.5)
  
  CplexSet[["ilp_edges"]] = score_ilp_edges(CplexSet, NodeSet)
  
  CplexSet[["heterodimer_ilp_edges"]] = score_heterodimer_ilp_edges(CplexSet, 
                                                                    type_score_heterodimer = 0,
                                                                    MS2_score_experiment_fragment = 0.5)
  print("Finish CplexSet scoring")
  
  CplexSet[["para"]] = Initiate_cplexset(CplexSet)
  CplexSet[["para_reduce"]] = Initiate_cplexset_reduce(CplexSet)
  
  print(paste("Complexity is", CplexSet$para$nc, "variables and", CplexSet$para$nr, "constraints."))
  print(sapply(CplexSet$para, length))
}

# Run optimization ####
{
  CplexSet[["ilp_solution"]] = Run_cplex(CplexSet, obj_cplex = CplexSet$para$obj, 
                                         optimization = c("ilp"),
                                         relative_gap = 1e-3, total_run_time = 5000)
  CplexSet = add_CPLEX_solution(CplexSet, Mset, 
                                solution = "ilp_solution")
  # 
  # CplexSet[["lp_solution"]] = Run_cplex(CplexSet, obj_cplex = CplexSet$para$obj,
  #                                       optimization = c("lp"),
  #                                       total_run_time = 5000)
  # CplexSet = add_CPLEX_solution(CplexSet, Mset,
  #                               solution = "lp_solution")

}

# Path annotation ####
{
  NetworkSet = list()
  NetworkSet = Initiate_networkset(CplexSet, StructureSet_df, LibrarySet, 
                                   solution = "ilp_solution")
  
  CplexSet = path_annotation(CplexSet, NetworkSet, solution = "ilp_solution")
  
  save(CplexSet,StructureSet_df,LibrarySet,NetworkSet,
          file = paste("NetID", "output.RData", sep="_"))
}

# Output ####
{
  NetID_output =  get_NetID_output(CplexSet$ilp_nodes, simplified = T)
  write.csv(NetID_output,"NetID_output.csv", row.names = F, na="")
}

save.image()
print("total run time")
print(Sys.time()-printtime)

# Test 
{
  test = CplexSet$ilp_nodes %>%
    filter(node_id == 5115)

  test = NetID_output %>%
    filter(grepl("Ni", formula))
}


