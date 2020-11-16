netid_packages <- function(){
  
  netid_pkgs <- c("tidyverse","enviPat","slam","readr","stringi","stringr","fitdistrplus",
                  "pracma","readxl","igraph","janitor","tictoc","rstudioapi",
                  "shiny","shinyjs","shinythemes","ShinyTester","ChemmineR","ChemmineOB")
  
  list_installed <- installed.packages()
  
  new_pkgs <- subset(netid_pkgs, !(netid_pkgs %in% list_installed[, "Package"]))
  
  if(length(new_pkgs)!=0){
    if (!requireNamespace("devtools", quietly = TRUE))
      install.packages("devtools")
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(new_pkgs)
    print(c(new_pkgs, " packages added..."))
  }
  
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
}

netid_packages()
devtools::install_github("LiChenPU/Formula_manipulation")
