## server ##
server <- function(input, output, session) {
  
  
  # top part - Peak list ####
  {
    peak_list <- reactive({
      print("enter show_peak_list")
      show_peak_list(ilp_nodes,
                     input_interest = input$input_interest, 
                     ion_form = input$ion_form, 
                     mz_ppm = input$mz_ppm)
      
    })
    
    output$peak_list <- DT::renderDataTable({
      peak_list()
    }, options = list(pageLength = 5,
                      aoColumnDefs = list(list(sClass="alignLeft"))
    ),
    
    filter = 'bottom',
    rownames = FALSE
    )
  }
  
  
  # bottom part - Network visualization ####
  {
    ## selected peak ####
    node_selected <- reactive({
      print("enter node_selected")
      req(is.numeric(as.numeric(input$peak_id)))
      node_selected = ilp_nodes %>%
        filter(node_id == input$peak_id) %>%
        arrange(-ilp_solution, -cplex_score)
      if(input$optimized_only){
        node_selected = node_selected %>%
          slice(1)
      }
      node_selected
    })
    
    observeEvent(peak_list(), {
      print("enter update selected peak_id")
      x = peak_list() 
      # updateSelectInput(session, "formula_select",
      #                   # label = paste("Formula", length(x)),
      #                   choices = x,
      #                   selected = head(x, 1)
      updateSelectInput(session, "peak_id",
                        # label = "Peak ID",
                        # choices = as.character(x$peak_id),
                        selected = as.character(x$peak_id[1]))
    })
    
    query_ilp_id = reactive({
      print("enter query_ilp_id")
      node_selected() %>%
        filter(formula == req(query_formula())) %>%
        filter(class == req(query_class())) %>%
        pull(ilp_node_id)
      
    })
    
    # Use reactiveVal to pass along the latest formula and class whenever a change happened
    # A bug will happen otherwise because input$formula and input$class are not invalidated in 
    # the update event
    query_formula = reactiveVal()
    query_class = reactiveVal()
    observe(
      query_formula(input$formula)
    )
    observe(
      query_class(input$class)
    )
    observe({
      list(input$peak_id)
      x = node_selected()
      print("enter observe")
      
      peak_formula = x$formula
      updateSelectInput(session, "formula",
                        choices = peak_formula,
                        selected = peak_formula[1]
      )
      query_formula(peak_formula[1])
      
      peak_class = node_selected() %>% filter(formula == peak_formula) %>% pull(class)
      updateSelectInput(session, "class",
                        choices = peak_class,
                        selected = peak_class[1]
      )
      query_class(peak_class[1])
      print(peak_formula)
      print(peak_class)
      
    })
    
    ## Network graph ####
    {
      g_parent <- reactive({
        print("enter g_parent")
        g_interest = NULL
        if(input$class %in% c("Metabolite", "Putative metabolite")){
          g_interest = network_partner_met(query_ilp_id = query_ilp_id(),
                                           search_degree = input$connect_degree,
                                           g_annotation = g_met,
                                           core_ilp_node = core_met,
                                           optimized_only = isolate(input$optimized_only))
        }
        if(input$class == "Artifact"){
          g_interest = network_partner_nonmet(query_ilp_id = query_ilp_id(),
                                              search_degree = input$connect_degree,
                                              g_annotation = g_nonmet,
                                              core_ilp_node = core_nonmet,
                                              weight_tol = 1,
                                              optimized_only = isolate(input$optimized_only))
        }
        g_interest
      })
      
      g_child_met <- reactive({
        print("enter g_child_met")
        g_child_met = network_child_met(query_ilp_id = query_ilp_id(),
                                        g_annotation = g_met,
                                        connect_degree = input$connect_degree,
                                        optimized_only = isolate(input$optimized_only))
      })
      g_child_nonmet <- reactive({
        print("enter g_child_nonmet")
        g_child_artifact = network_child_nonmet(query_ilp_id = query_ilp_id(),
                                                g_annotation = g_nonmet,
                                                connect_degree = input$connect_degree,
                                                optimized_only = isolate(input$optimized_only))
      })
      
      g_interest <- eventReactive(input$plot_network, {
        print("enter g_interest")
        merge_nodes = NULL
        merge_edges = NULL
        if(input$biochemical_graph & input$class %in% c("Metabolite", "Putative metabolite")){
          # req(g_parent())
          if(!is.null(g_parent())){
            g_nodes = igraph::as_data_frame(g_parent(), "vertices")
            g_edges = igraph::as_data_frame(g_parent(), "edges")
            merge_nodes = bind_rows(merge_nodes, g_nodes)
            merge_edges = bind_rows(merge_edges, g_edges)
          }
          
        }
        
        if(input$abiotic_graph & input$class == "Artifact"){
          # req(g_parent())
          if(!is.null(g_parent())){
            g_nodes = igraph::as_data_frame(g_parent(), "vertices")
            g_edges = igraph::as_data_frame(g_parent(), "edges")
            merge_nodes = bind_rows(merge_nodes, g_nodes)
            merge_edges = bind_rows(merge_edges, g_edges)
          }
          
        }
        
        if(input$biochemical_graph){
          # req(g_child_met())
          if(!is.null(g_child_met())){
            g2_nodes = igraph::as_data_frame(g_child_met(), "vertices")
            g2_edges = igraph::as_data_frame(g_child_met(), "edges")
            merge_nodes = bind_rows(merge_nodes, g2_nodes)
            merge_edges = bind_rows(merge_edges, g2_edges)
          }
          
        }
        if(input$abiotic_graph){
          # req(g_child_nonmet())
          if(!is.null(g_child_nonmet())){
            g2_nodes = igraph::as_data_frame(g_child_nonmet(), "vertices")
            g2_edges = igraph::as_data_frame(g_child_nonmet(), "edges")
            merge_nodes = bind_rows(merge_nodes, g2_nodes)
            merge_edges = bind_rows(merge_edges, g2_edges)
          }
          
        }
        
        if(is.null(merge_edges)){
          print("merge_edges is null.")
          return(NULL)
        }
        
        
        merge_nodes = merge_nodes %>%
          distinct()
        
        g_interest = graph_from_data_frame(merge_edges,
                                           directed = T,
                                           merge_nodes)
        
      })
      
      output$Network_plot <- renderVisNetwork({
        print("enter output$Network_plot renderVisNetwork")
        # req(g_interest())
        Plot_g_interest(g_interest(), 
                        query_ilp_node = isolate(query_ilp_id()), 
                        show_node_labels = input$node_labels, 
                        show_edge_labels = input$edge_labels,
                        log_inten_cutoff = 0
        )
      })
      
      # try to update the netwrok whenever peak_id is updated, make the code too fragile
      # observeEvent(input$peak_id, {
      #   click("plot_network")
      # }, ignoreInit = T)
    }
    
    
    ## show structure_table ####
    {
      structure_table_trigger <- reactiveVal()
      observeEvent(input$plot_network, {
        structure_table_trigger("plot_network")
      })
      observeEvent(input$click, {
        structure_table_trigger("click")
      })
      
      # Use reactiveVal + observeEvant to update internal value
      structure_table <- reactive({
        print("enter structure_table")
        req(structure_table_trigger())
        req(query_ilp_id())
        if(structure_table_trigger() == "plot_network"){
          print(query_ilp_id())
          core_annotation %>%
            filter(ilp_node_id == query_ilp_id())
        } else if(structure_table_trigger() == "click"){
          core_annotation %>%
            filter(ilp_node_id == input$click)
        }
      })
      
      output$structure_list <- DT::renderDataTable({
        print("enter output$structure_list")
        structure_table() %>%
          dplyr::select(class, annotation, origin, note)
      }, options = list(pageLength = 5),
      fillContainer = F,
      rownames = TRUE
      )
    }
    
    
    ## show structure_plot ####
    {
      # initialize the two parameter
      structure_plot_counter = reactiveVal(1)  
      structure_table_rows = reactive(nrow(structure_table()))
      
      # goto next smiles when click -> button
      observeEvent(input$next_struct, {
        req(structure_plot_counter() <= structure_table_rows())
        structure_plot_counter(structure_plot_counter()+1)
      })
      
      # goto previous smiles when click <- button
      observeEvent(input$previous_struct, {
        req(structure_plot_counter() > 1)
        structure_plot_counter(structure_plot_counter()-1)
      })
      
      # reset when new structure table is trigger / dblclick also resets
      observeEvent(c(structure_table(), input$struct_plot_dblclick), {
        shinyjs::show(id = "previous_struct", time = 0)
        shinyjs::show("next_struct")
        shinyjs::show("struct_num")
        structure_plot_counter(1)
      })
      
      # hide initially and show when structure_table is first called
      observe({
        shinyjs::hide("previous_struct")
        shinyjs::hide("next_struct")
        shinyjs::hide("struct_num")
      })
      
      observeEvent(input$struct_num, {
        structure_plot_counter(as.numeric(input$struct_num))
      })
      
      observeEvent(structure_plot_counter(), {
        updateSelectInput(session, "struct_num",
                          choices = 1:structure_table_rows(),
                          selected = structure_plot_counter()
        )
      })
      
      
      
      output$structure <- renderPlot({
        print("enter output$structure")
        smiles = structure_table() %>%
          pull(SMILES)
        my_SMILES2structure(smiles[structure_plot_counter()])
      })
      
      output$struct_annotation <- renderText({
        print("enter output$struct_annotation")
        struct_annotation = structure_table() %>%
          pull(annotation)
        struct_annotation[structure_plot_counter()]
      })
      
    }
  }
}





