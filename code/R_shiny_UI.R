## ui ####
ui <- tagList(
  # shinythemes::themeSelector(),
  useShinyjs(),
  navbarPage(
    theme = shinytheme("cerulean"), # other theme can be viewed from themeSelector()
    "NetID",
    ## Input file tab ####
    # tabPanel("Input file",
    #   fileInput("file", "Choose NetID_output rds File",
    #             multiple = FALSE,
    #             accept = NULL)
    #          
    # ),
    
    ## Peak list tab ####
    tabPanel("Peak list",
      ## sidebarPanel ####
      # sidebarPanel(
      #   width = 3,
      fluidRow(
        ## filter and parameter ####
        column(3, 
               wellPanel(
                 textInput(inputId = "input_interest", 
                              label = "Enter a mz or formula of interest",
                              value = 0),
                 selectInput(inputId = "ion_form", 
                             label = "Select ionization",
                             choices = c("M+H","M-H","M"),
                             selected = "M"),
                 numericInput(inputId = "mz_ppm", 
                              label = "ppm",
                              value = 3)
                 
                 
                 )
               
        ),
        
        ## peak_list ####
        column(9,
               DT::dataTableOutput("peak_list")
               )
        
      ),
      
      br(),
      hr(),
      
      
      fluidRow(
        ## Network plot and options ####
        column(6, 
               wellPanel(
                 fluidRow(
                 column(10, 
                   column(3,
                          selectInput(inputId = "peak_id",
                                       label = "Peak ID",
                                      choices = unique(ilp_nodes$node_id))),
                   column(4,
                          selectInput(inputId = "formula",
                                      label = "Formula",
                                      choices = character(0))),
                   column(3,
                          selectInput(inputId = "class",
                                      label = "Class",
                                      choices = character(0))),
                   column(2,
                          numericInput(inputId = "connect_degree", 
                                       label = "Degree",
                                       value = 1)),
                   
                   column(3,
                          checkboxInput("biochemical_graph", "Biochemical graph",
                                        value = T)),
                   column(3,
                          checkboxInput("abiotic_graph", "Abiotic graph",
                                        value = F)),
                   column(2,
                          checkboxInput("node_labels", "Node labels",
                                        value = T)),
                   column(2,
                          checkboxInput("edge_labels", "Edge labels",
                                        value = T)),
                   column(2,
                          checkboxInput("optimized_only", "Optimized only",
                                        value = T))
                 ),
                 column(2, 
                        # column(6,
                               actionButton("plot_network", label = h2("Plot"),
                                            with = "400px")
                               # )
                        
                        ),
                 hr(),
                 downloadButton("download_html", "Download plot"),
                 column(12,
                        visNetworkOutput("Network_plot", width = "100%",  height = "750px")
                 )
               ))
          
          
        ),
        
        
        
        ## Structure ####
        column(6,
               wellPanel(fluidRow(
                 plotOutput("structure", 
                              dblclick = dblclickOpts(id = "struct_plot_dblclick")),
                   
                 textOutput("struct_annotation"),
                 # Offset with previous section
                 column(2, offset = 1, 
                        actionButton("previous_struct", "<-")),
                 column(2, offset = 0, 
                        actionButton("next_struct", "->")),
                 column(2, offset = 0,
                        selectInput(inputId = "struct_num",
                                    label = NULL,
                                    choices = numeric(0))),
                 column(2, offset = 1,
                        # Download Button
                        downloadButton("download_csv", "Download csv")
                        ),
                 
                 DT::dataTableOutput("structure_list")
               )))
        ) # End of lower graph ####
    ) # End of tabPanel 
  ), # End of navbarPage
  
  tags$style(type='text/css', "#plot_network { width: 100%; color: #000000; margin-top: 25px;}")
) 
      
    
      
     





### old ######

# wellPanel(
#   sliderInput(inputId = "Peak_inten_range", 
#               label = "Peak Intensity (log10)",
#               min = 0, max = 10, step = 0.01,
#               value = c(0,10)),
#   sliderInput(inputId = "mz_range", 
#               label = "m/z range",
#               min = 0, max = 1500, step = 1,
#               value = c(0,1500)),
#   sliderInput(inputId = "rt_range", 
#               label = "RT range",
#               min = 0, max = 20, step = .01,
#               value = c(0,20))
# )
# mainPanel(
#   tabPanel("peak_table", DT::dataTableOutput("peak_list"))
# )
# ),

## Network visualization tab ####
# tabPanel("Network visualization",
#          ## First row ####
#          fluidRow(
#            column(3,
#                   numericInput(inputId = "peak_id",
#                                label = "Peak ID",
#                                value = 1
#                   )
#            ),
#            column(3,
#                   verbatimTextOutput("peak_mz")
#            ),
#            column(3,
#                   verbatimTextOutput("peak_rt")
#            ),
#            column(3,
#                   verbatimTextOutput("peak_inten")
#            )
#          ),
#          
#          ## Second row ####
#          fluidRow(
#            column(3, 
#                   selectInput(inputId = "formula",
#                               label = "Formula",
#                               choices = character(0)
#                   )
#            ),
#            column(3,
#                   selectInput(inputId = "class",
#                               label = "Class",
#                               choices = character(0)
#                   )
#            ),
#            column(2,
#                   checkboxInput("optimized_only", "optimized_only",
#                                 value = T)
#            )
#          ),
#          hr(),
#          ## Network + structure ####
#          fluidRow(
#            ## Network ####
#            column(7,
#                   visNetworkOutput("Network_plot", height = "100%", width = "100%"),
#                   column(2,
#                          checkboxInput("node_labels", "Node labels",
#                                        value = T)
#                   ),
#                   column(2,
#                          checkboxInput("edge_labels", "Edge labels",
#                                        value = T)
#                   ),
#                   column(2,
#                          checkboxInput("parent_graph", "Parent graph",
#                                        value = T)
#                   ),
#                   column(2,
#                          checkboxInput("child_graph", "Child graph",
#                                        value = F)
#                   ),
#                   column(2, offset = 1,
#                          actionButton("plot_network", "Plot network")
#                   )
#                   
#            ),
#            ## Structure ####
#            column(5,
#                   plotOutput("structure", 
#                              click = clickOpts(id = "struct_plot_click"),
#                              dblclick = dblclickOpts(id = "struct_plot_dblclick"),
#                              hover = hoverOpts(id = "struct_plot_hover"),
#                              brush = brushOpts(id = "struct_plot_brush")),
#                   DT::dataTableOutput("structure_list")
#                   
#                   ## Structure here ##
#                   )
#          )
# 