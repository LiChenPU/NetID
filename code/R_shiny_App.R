
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('R_shiny_functions.R')

# options(digits = 8)

# Read in files ####
datapath = "./Sc_neg"

setwd(datapath)
filename = "20201103114511_output.RData" # Sc_neg
load(filename)

# Global parameter in the background ####
{
  ilp_nodes = CplexSet$ilp_nodes %>%
    mutate(medMz = signif(medMz, 7),
           medRt = round(medRt, 2),
           log10_inten = round(log10_inten, 2),
           ppm_error = (mass-medMz)/medMz*1e6,
           ppm_error = round(ppm_error, 2))
  
  g_met = NetworkSet$g_met
  core_met = ilp_nodes %>%
    filter(steps == 0) %>%
    filter(class %in% c("Metabolite", "Putative metabolite"))
  
  g_nonmet = NetworkSet$g_nonmet
  core_nonmet = ilp_nodes %>%
    filter(steps %% 1 == 0) %>%
    filter(class != "Unknown")
  
  core_annotation = NetworkSet$core_annotation %>%
    arrange(ilp_node_id, -rank_score) %>%
    # distinct(ilp_node_id, .keep_all = T) %>%
    filter(T)
  
}


## Run shiny ####
{
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  # reactlogReset()
  # options(shiny.reactlog=TRUE) # Enable reactlog here
  
  source('R_shiny_functions.R')
  source('R_shiny_UI.R', local = TRUE)
  source('R_shiny_Server.R')
  app = shinyApp(ui = ui, server = server)
  runApp(app)
  # shiny::reactlogShow() # Run after closing the app
  
  ## To do ####
  # Add the output button
  
  # fix a minor bug where ilp_node_id = 3106 
  ## Column `from` can't be converted from numeric to character
  
  
}


