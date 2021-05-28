
WORK_DIR = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(WORK_DIR)
## set WORD_DIR to the absolute directory of the FDR example folder
REP <- 10 ## REP is the time of repeats for Target-Decoy running 
RATIO <- 1 ## RATIO = # of decoy / # of target, 1 is recommended

source("Target-Decoy Generator.R")

Target_Decoy_Generator(RATIO, REP, "hmdb")
Target_Decoy_Generator(RATIO, REP, "ymdb")
Target_Decoy_Generator(RATIO, REP, "pbcm")
Target_Decoy_Generator(RATIO, REP, "pbcm_bio")


NetID_wrapper(RATIO, REP, "hmdb_decoy")
NetID_wrapper(RATIO, REP, "ymdb_decoy")
NetID_wrapper(RATIO, REP, "pbcm_decoy")
NetID_wrapper(RATIO, REP, "pbcm_bio_decoy")


hmdb_mm <- ScareCrow_Repeater("hmdb_decoy")
ymdb_mm <- ScareCrow_Repeater("ymdb_decoy")
pbcm1_mm <- ScareCrow_Repeater("pbcm_decoy")
pbcm2_mm <- ScareCrow_Repeater("pbcm_bio_decoy")
print(Sys.time())

save.image(file = "FDR.RData")






