# Shiny library ####
{
  library(shiny)
  library(shinyjs)
  library(igraph)
  library(reactlog)
  library(ShinyTester)
  library(shinythemes)
  library(visNetwork)
  library(dplyr)
  library(RColorBrewer)
  library(stringr)
  library(ChemmineR)
  library(ChemmineOB)
  library(enviPat)
  library(lc8)
  data("isotopes")
  
  options(shiny.maxRequestSize=1*1024^3)
}

# function ####

## show_peak_list ####
show_peak_list = function(ilp_nodes, input_interest, ion_form, mz_ppm){
  
  if(is.null(ilp_nodes)){
    return(NULL)
  }
  ilp_nodes_filter = ilp_nodes %>%
    filter(ilp_solution > 1e-6) %>%
    dplyr::select(node_id, medMz, medRt, log10_inten, class, formula, ppm_error) %>%
    dplyr::rename(peak_id = node_id)
  
  mz_interest = 0
  if(!is.na(as.numeric(input_interest))){
    mz_interest = as.numeric(input_interest)
  } else {
    input_interest = gsub(" ", "", input_interest)
    formula_check = enviPat::check_chemform(isotopes, input_interest)
    if(!formula_check$warning){
      mz_interest = formula_mz(formula_check$new_formula)
    } 
  }
  
  if(mz_interest != 0){
    if(ion_form == "M"){mz_adjust = mz_interest}
    if(ion_form == "M+H"){mz_adjust = mz_interest - 1.007276}
    if(ion_form == "M-H"){mz_adjust = mz_interest + 1.007276}
    ilp_nodes_filter = ilp_nodes_filter %>%
      filter(abs(medMz - mz_adjust) < medMz * mz_ppm * 1e-6)
  }
  
  return(ilp_nodes_filter)
}

## network_partner_met ####  
network_partner_met = function(query_ilp_id, 
                               search_degree = 2,
                               g_annotation = g_met,
                               core_ilp_node = core_met,
                               optimized_only = T){
  
  # g_annotation = g_met
  # query_ilp_id = 1
  # core_ilp_node = core_met
  # optimized_only = T
  # search_degree = 2
  
  if(optimized_only){
    g_nodes = igraph::as_data_frame(g_annotation, "vertices") %>%
      filter(ilp_solution != 0)
    g_edges = igraph::as_data_frame(g_annotation, "edges") %>%
      filter(from %in% g_nodes$name, to %in% g_nodes$name)
    
    g_annotation = graph_from_data_frame(g_edges,
                                         directed = T,
                                         vertices = g_nodes)
    
    core_ilp_node = core_ilp_node %>%
      filter(ilp_solution != 0)
  }
  
  core_ilp_id = as.character(core_ilp_node$ilp_node_id)
  query_ilp_id = as.character(query_ilp_id)
  
  valid_ilp_ids = igraph::as_data_frame(g_annotation, "vertices") %>% pull(name)
  
  if(!any(query_ilp_id == valid_ilp_ids)){
    return(NULL)
  }
  
  query_distMatrix <- shortest.paths(g_annotation,
                                     v=core_ilp_id,
                                     to=query_ilp_id,
                                     mode = "all")
  
  
  
  query_distMatrix_min = min(query_distMatrix[query_distMatrix!=0])
  
  if(query_distMatrix_min > search_degree){
    sub_nodes = igraph::as_data_frame(g_annotation, "vertices") %>%
      filter(name == query_ilp_id)
    
    sub_edges = igraph::as_data_frame(g_annotation, "edges") %>%
      filter(FALSE)
    
    g_met_interest = graph_from_data_frame(sub_edges,
                                           directed = T,
                                           vertices = sub_nodes)
    return(g_met_interest)
  }
  
  shortest_ilp_nodes = which(query_distMatrix <= search_degree)
  parent_selected = core_ilp_id[shortest_ilp_nodes]
  paths_connect_ij_nodes = lapply(parent_selected, function(x){
    shortest_paths(g_annotation, 
                   from = x, 
                   to = query_ilp_id, 
                   mode = "all", 
                   output = "vpath")
  })
  
  ilp_node_path = lapply(paths_connect_ij_nodes, function(x){
    x[[1]] %>% unlist() %>% names() %>% as.numeric()
  })
  
  ilp_node_interest = unlist(ilp_node_path)
  
  sub_nodes = igraph::as_data_frame(g_annotation, "vertices") %>%
    filter(name %in% as.character(ilp_node_interest))
  
  sub_edges = igraph::as_data_frame(g_annotation, "edges") %>%
    filter(from %in% as.character(ilp_node_interest) & 
             to %in% as.character(ilp_node_interest))
  
  g_met_interest = graph_from_data_frame(sub_edges,
                                         directed = F,
                                         vertices = sub_nodes)
  # plot.igraph(g_met_interest)
  
  return(g_met_interest)
}

## network_partner_nonmet ####  
network_partner_nonmet = function(query_ilp_id, 
                                  search_degree = 2,
                                  g_annotation = g_nonmet,
                                  core_ilp_node = core_nonmet,
                                  weighted = F,
                                  weight_tol = 0,
                                  optimized_only = T){
  
  # query_ilp_id = 63662
  # g_annotation = g_nonmet
  # core_ilp_node = core_nonmet
  # weighted = F
  # weight_tol = 1
  # optimized_only = T
  
  
  if(optimized_only){
    g_nodes = igraph::as_data_frame(g_annotation, "vertices") %>%
      filter(ilp_solution != 0)
    # g_edges = igraph::as_data_frame(g_annotation, "edges") %>%
    #   filter(ilp_solution != 0)
    g_edges = igraph::as_data_frame(g_annotation, "edges") %>%
      filter(from %in% g_nodes$name, to %in% g_nodes$name) %>%
      arrange(-ilp_solution) %>%
      distinct(from, to, .keep_all=T) %>%
      filter(ilp_solution>1e-6)
    
    g_annotation = graph_from_data_frame(g_edges, 
                                         directed = T, 
                                         vertices = g_nodes)
    
    core_ilp_node = core_ilp_node %>%
      filter(ilp_solution != 0)
  }
  
  core_ilp_id = as.character(core_ilp_node$ilp_node_id)
  query_ilp_id = as.character(query_ilp_id)
  
  valid_ilp_ids = igraph::as_data_frame(g_annotation, "vertices") %>% pull(name)
  
  if(!any(query_ilp_id == valid_ilp_ids)){
    return(NULL)
  }
  
  
  g_edge = igraph::as_data_frame(g_annotation, "edges")
  
  
  if(weighted){
    query_distMatrix <- shortest.paths(g_annotation,
                                       v=core_ilp_id,
                                       to=query_ilp_id,
                                       mode = "out",
                                       weights = g_edge$edge_weight)
  } else {
    query_distMatrix <- shortest.paths(g_annotation,
                                       v=core_ilp_id,
                                       to=query_ilp_id,
                                       mode = "out")
  }
  
  
  query_distMatrix_min = min(query_distMatrix[query_distMatrix!=0])
  
  if(is.infinite(query_distMatrix_min)){
    sub_nodes = igraph::as_data_frame(g_annotation, "vertices") %>%
      filter(name == query_ilp_id)
    
    sub_edges = igraph::as_data_frame(g_annotation, "edges") %>%
      filter(FALSE)
    
    g_interest = graph_from_data_frame(sub_edges,
                                       directed = F,
                                       vertices = sub_nodes)
    return(g_interest)
  }
  
  shortest_ilp_nodes = which(query_distMatrix <= search_degree + weight_tol)
  # shortest_ilp_nodes = which(query_distMatrix == query_distMatrix_min)
  parent_selected = core_ilp_id[shortest_ilp_nodes]
  paths_connect_ij_nodes = lapply(parent_selected, function(x){
    all_shortest_paths(g_annotation, 
                       from = x, 
                       to = query_ilp_id, 
                       mode = "out")
  })
  
  ilp_node_path = lapply(paths_connect_ij_nodes, function(x){
    x[[1]] %>% unlist() %>% names() %>% as.numeric()
  })
  
  ilp_node_interest = unlist(ilp_node_path)
  
  sub_nodes = igraph::as_data_frame(g_annotation, "vertices") %>%
    filter(name %in% as.character(ilp_node_interest))
  
  sub_edges = igraph::as_data_frame(g_annotation, "edges") %>%
    filter(from %in% as.character(ilp_node_interest) & 
             to %in% as.character(ilp_node_interest))
  
  g_interest = graph_from_data_frame(sub_edges,
                                     directed = T,
                                     vertices = sub_nodes)
  # plot.igraph(g_interest)
  return(g_interest)
}



## network_child_met ####  
network_child_met = function(query_ilp_id, 
                                g_annotation = g_met,
                                connect_degree = 1, 
                                optimized_only = T){
  
  # query_ilp_id = 52577
  # g_annotation = g_met
  # connect_degree = 1
  # optimized_only = T
  
  query_ilp_id = as.character(query_ilp_id)
  g_nodes = igraph::as_data_frame(g_annotation, "vertices")
  g_edges = igraph::as_data_frame(g_annotation, "edges")
  
  if(optimized_only){
    g_nodes = g_nodes %>%
      filter(ilp_solution != 0)
    g_edges = g_edges %>%
      filter(from %in% g_nodes$name, to %in% g_nodes$name) %>%
      arrange(-ilp_solution) %>%
      distinct(from, to, .keep_all=T)
    
    g_annotation = graph_from_data_frame(g_edges, 
                                         directed = T, 
                                         vertices = g_nodes)
  }
  
  if(!as.numeric(query_ilp_id) %in% g_nodes$name){
    return(NULL)
  }
  
  
  g_partner = make_ego_graph(g_annotation, 
                             connect_degree,
                             nodes = query_ilp_id, 
                             mode = c("all"))[[1]]
  
  # plot.igraph(g_partner)
  return(g_partner)
}

## network_child_nonmet ####  
network_child_nonmet = function(query_ilp_id, 
                                g_annotation = g_nonmet,
                                connect_degree = 1, 
                                optimized_only = T){
  
  # query_ilp_id = 63662
  # g_annotation = g_nonmet
  # core_ilp_node = core_nonmet
  
  query_ilp_id = as.character(query_ilp_id)
  
  if(optimized_only){
    g_nodes = igraph::as_data_frame(g_annotation, "vertices") %>%
      filter(ilp_solution != 0)
    # g_edges = igraph::as_data_frame(g_annotation, "edges") %>%
    #   filter(ilp_solution != 0)
    g_edges = igraph::as_data_frame(g_annotation, "edges") %>%
      filter(from %in% g_nodes$name, to %in% g_nodes$name) %>%
      arrange(-ilp_solution) %>%
      distinct(from, to, .keep_all=T) %>%
      filter(ilp_solution > 0.01)
    
    g_annotation = graph_from_data_frame(g_edges, 
                                         directed = T, 
                                         vertices = g_nodes)
  }
  
  g_partner = make_ego_graph(g_annotation, 
                             connect_degree,
                             nodes = query_ilp_id, 
                             mode = c("out"))[[1]]
  
  # plot.igraph(g_partner)
  return(g_partner)
}

## Plot_g_interest ####
Plot_g_interest = function(g_interest, query_ilp_node, 
                           show_node_labels = T, show_edge_labels = T,
                           log_inten_cutoff = 0)
{
  
  # query_ilp_node = 1
  # show_node_labels = T
  # show_edge_labels = T
  # log_inten_cutoff = 4
  # g_interest = network_child_met(query_ilp_node, connect_degree=2)
  # temp_nodes = igraph::as_data_frame(g_interest, "vertices")
  
  if(!is.igraph(g_interest)){
    print("g_interest is not a graph")
    return(NULL)
  }
  
  my_palette = brewer.pal(5, "Set3")
  nodes = igraph::as_data_frame(g_interest, "vertices") %>%
    dplyr::rename(id = name) %>%
    filter(log10_inten > log_inten_cutoff) %>%
    mutate(size = log10_inten * 2) %>%
    mutate(label = "") %>%
    mutate(color.border = "black") %>% 
    mutate(color.background = case_when(
      # assigned to the first color, not overwitten by later assignment
      id == as.character(query_ilp_node) ~ my_palette[4],
      class == "Metabolite" ~ my_palette[1],
      class == "Putative metabolite" ~ my_palette[2],
      class == "Artifact" ~ my_palette[3],
      class == "Unknown" ~ "#666666"
    )) %>%
    # [:digit:] means 0-9, \\. means ".", + means one or more, \\1 means the content in (), <sub> is HTML language
    mutate(formula_sub = str_replace_all(formula,"([[:digit:]|\\.|-]+)","<sub>\\1</sub>")) %>%
    mutate(title = paste0("Formula:", formula_sub, "<br>",
                          "ID:", node_id, "<br>",
                          "mz:", round(medMz,4), "<br>",
                          "RT:", round(medRt,2), "<br>",
                          "TIC:", format(signif(10^log10_inten,5), scientific = T))
    )
  
  
  if(show_node_labels){
    nodes = nodes %>%
      mutate(label = formula)
  }
  
  edges = igraph::as_data_frame(g_interest, "edges") %>%
    filter(from %in% nodes$id, to %in% nodes$id) %>%
    mutate(arrows = "") %>%
    # mutate(arrows = "to") %>%
    mutate(label = "") %>%
    mutate(color = "#666666") %>%
    filter(T)
  
  if(show_edge_labels & nrow(edges)!=0){
    edges = edges %>%
      mutate(font.color = brewer.pal(9, "Blues")[7]) %>% 
      mutate(label = case_when(
        category == "Multicharge" ~ paste0(linktype, "-charge"),
        category == "Oligomer" ~ paste0(linktype, "-mer"),
        category == "Heterodimer" ~ paste0("Peak ", linktype),
        category == "Library_MS2_fragment" ~ linktype,
        # direction == -1 ~ paste0("-", linktype),
        T ~ linktype
      ))
  }
  Output = visNetwork(nodes, edges) %>%
    visLegend() %>%
    visOptions(manipulation = TRUE,
               highlightNearest = TRUE
               # height = "100%", width = "100%"
    ) %>%
    visEvents(click = "function(nodes){
                  Shiny.onInputChange('click', nodes.nodes[0]);
              ;}") %>%
    visInteraction(navigationButtons = TRUE) %>%
    visIgraphLayout(layout = "layout_with_fr",
                    randomSeed = 1234)

  
}
## download_g_interest #### 
download_g_interest = function(g_interest, query_ilp_node){
  if(!is.igraph(g_interest)){
    print("g_interest is not a graph")
    return(NULL)
  }
  
  my_palette = brewer.pal(5, "Set3")
  nodes = igraph::as_data_frame(g_interest, "vertices") %>%
    pull(name)
  
  ilp_nodes %>%
    filter(ilp_node_id %in% as.numeric(nodes)) %>%
    dplyr::rename(peak_id = node_id,
                  annotation = path) %>%
    dplyr::select(peak_id, medMz, medRt, log10_inten, class, formula, ppm_error, annotation, 
                  everything())
    
}

## my_SMILES2structure ####
my_SMILES2structure = function(SMILES){
  SDF = try(smiles2sdf(SMILES), silent=T)
  if(inherits(SDF, "try-error")){
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    # text(x=0.5, y=0.5, paste("Invalid SMILES"))
    return(0)
  }
  tryError = try(ChemmineR::plotStruc(SDF[[1]]), silent=T)
  if(inherits(tryError, "try-error")){
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    # text(x=0.5, y=0.5, paste("Invalid SDF"))
    return(0)
  }
  return(1)
  
}