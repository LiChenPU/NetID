library(dplyr)
library(xlsx)
library(janitor)
WORK_DIR = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(WORK_DIR)

## statistics of global networks ####
{
  setwd(WORK_DIR)
  
  datapath = "Sc_neg_ymdb"
  setwd(datapath)
  ## Run NetID first and copy the environment .RData to the datapath first ##
  filename = ".RData" 
  load(filename)
  datapath = Sc_neg_ymdb
  
  # Global parameter in the background ####
  {
    ilp_nodes = CplexSet$ilp_nodes %>%
      mutate(medMz = signif(medMz, 7),
             medRt = round(medRt, 2),
             log10_inten = round(log10_inten, 2),
             ppm_error = (mass-medMz)/medMz*1e6,
             ppm_error = round(ppm_error, 2))
    
    ilp_edges = CplexSet$ilp_edges
    
    g_met = NetworkSet$g_met
    
    g_nonmet = NetworkSet$g_nonmet
    
  }
  
  
  # Node ####
  {
    
    test = ilp_nodes %>% 
      filter(ilp_solution > 0.01)
    test2 = test %>% 
      filter(class == "Metabolite")
    
    tabyl(test, class)
    # print(c(nrow(test3), nrow(test3)-nrow(test3_filter)))
    
  }
  # Edge ####
  {
    ilp_edges_filter = ilp_edges %>%
      arrange(-ilp_solution) %>%
      distinct(edge_id, .keep_all = T) %>%
      filter(ilp_nodes1 %in% test$ilp_node_id, 
             ilp_nodes2 %in% test$ilp_node_id) %>%
      filter(!(category != "Biotransform" & ilp_solution == 0))
    
    
  }
  # Network ####
  
  # Statistics
  {
    # Initial steps
    {
      # Formulas
      formulaset_initial = StructureSet_df %>%
        filter(class != "Unknown") %>%
        filter(transform == "") %>%
        distinct(node_id, formula, .keep_all = T) %>%
        filter(steps == 0)
      
      
      # After edge expansion
      formulaset_expand = StructureSet_df %>%
        filter(class != "Unknown") %>%
        distinct(node_id, formula, class, .keep_all=T)
      
      
      
    }
    
  }
  # Output ####
  {
    setwd(WORK_DIR)

    wb = loadWorkbook("summary.xlsx", password=NULL)
    all_sheets = getSheets(wb)
    sheet_id = which(names(all_sheets) == datapath)
    if(length(sheet_id)==0){
      sheet <- createSheet(wb, sheetName = datapath)
    } else {
      sheet = all_sheets[[sheet_id]]
    }
    
    addDataFrame(tabyl(test$class), sheet, startRow=2, startColumn=1)
    addDataFrame(data.frame(total = nrow(ilp_edges_filter),
                            biotransform = nrow(ilp_edges_filter %>% filter(category == "Biotransform")),
                            artifact = nrow(ilp_edges_filter %>% filter(category != "Biotransform"))
    ),
    sheet, startRow=7, startColumn=1)
    
    addDataFrame(data.frame(unique_nodes = nrow(distinct(formulaset_initial, node_id)),
                            initial_formula = nrow(formulaset_initial),
                            initial_edge = nrow(EdgeSet_all),
                            biochemical = nrow(EdgeSet_all %>% filter(category == "Biotransform")),
                            abiotic = nrow(EdgeSet_all %>% filter(category != "Biotransform")),
                            unique_expanded = nrow(distinct(formulaset_expand, node_id)),
                            formula_expanded = nrow(formulaset_expand)),
                 sheet, startRow=10, startColumn=1)
    
    addDataFrame(table(formulaset_expand$class),
                 sheet, startRow=13, startColumn=1)
    addDataFrame(table(table(formulaset_initial$node_id)),
                 sheet, startRow=18, startColumn=1)
    
    addDataFrame(data.frame(total = nrow(ilp_edges_filter),
                            biochemical = nrow(ilp_edges_filter %>% filter(category == "Biotransform")),
                            abiotic = nrow(ilp_edges_filter %>% filter(category != "Biotransform"))
                            ),
                 sheet, startRow=7, startColumn=6)
    
    saveWorkbook(wb, "summary.xlsx")
  }
}
