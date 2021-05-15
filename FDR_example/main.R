WORK_DIR <- "C:/Users/liche/Desktop/revision/FDR_plot/NetID_FDR/FDR example/FDR_example/" 
## set WORD_DIR to the absolute directory of the FDR example folder
REP <- 2
## REP is the time of repeats for Target-Decoy running
RATIO <- 1
## RATIO = # of decoy / # of target, 1 is recommended

setwd(WORK_DIR)

source("Target-Decoy Generator.R")
source("plot.R")


Target_Decoy_Generator(RATIO, REP, "hmdb")
Target_Decoy_Generator(RATIO, REP, "ymdb")
# Target_Decoy_Generator(RATIO, REP, "pbcm")
Target_Decoy_Generator(RATIO, REP, "pbcm_new")
# Target_Decoy_Generator(RATIO, REP, "pbcm_bio")
Target_Decoy_Generator(RATIO, REP, "pbcm_bio_new")

NetID_wrapper(RATIO, REP, "hmdb_decoy")
NetID_wrapper(RATIO, REP, "ymdb_decoy")
# NetID_wrapper(RATIO, REP, "pbcm_decoy")
NetID_wrapper(RATIO, REP, "pbcm_new_decoy")
# NetID_wrapper(RATIO, REP, "pbcm_bio_decoy")
NetID_wrapper(RATIO, REP, "pbcm_bio_new_decoy")

hmdb_mm <- ScareCrow_Repeater("hmdb_decoy")
ymdb_mm <- ScareCrow_Repeater("ymdb_decoy")
# pbcm1_mm <- ScareCrow_Repeater("pbcm_decoy")
pbcm1_mm <- ScareCrow_Repeater("pbcm_new_decoy")
# pbcm2_mm <- ScareCrow_Repeater("pbcm_bio_decoy")
pbcm2_mm <- ScareCrow_Repeater("pbcm_bio_new_decoy")


gt_all <- rbind(MM_2_GT(hmdb_mm, "HMDB", REP), MM_2_GT(ymdb_mm, "YMDB", REP),
                MM_2_GT(pbcm1_mm, "PBCM_BIO", REP), MM_2_GT(pbcm2_mm, "PBCM", REP))
fdr_all <- rbind(MM_2_FDR(hmdb_mm, "HMDB", REP), MM_2_FDR(ymdb_mm, "YMDB", REP),
                 MM_2_FDR(pbcm1_mm, "PBCM_BIO", REP), MM_2_FDR(pbcm2_mm, "PBCM", REP))

gt_all$method = factor(gt_all$method, levels = c('mz', 'node', 'edge', 'ni'))
fdr_all$method = factor(fdr_all$method, levels = c('mz', 'node', 'edge', 'ni'))

ggplot(data = gt_all, mapping = aes(x = method, y = precision, color = method)) + 
    geom_boxplot() + #geom_point(data = stan, mapping = aes(x = method, y = precision, fill = method), shape = 22, color = "black") +
    geom_point(alpha = 1 / 5, position = "jitter") +
    ggtitle("Precision by Ground Truth") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~ library, nrow = 2) #+coord_flip()

ggplot(data = fdr_all, mapping = aes(x = method, y = FDR, color = method)) + 
    geom_boxplot() + 
    geom_point(alpha = 1 / 5, position = "jitter") +
    ggtitle("FDR by Nature Method") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~ library, nrow = 2) #+coord_flip()

ggplot(data = gt_all, mapping = aes(x = method, y = precision, color = method)) +
    stat_summary(fun.min = min, fun.max = max, fun = mean) + 
    geom_point(alpha = 1 / 5, position = "jitter") + 
    #geom_point(data = stan, mapping = aes(x = method, y = precision, fill = method), shape = 22, color = "black") + 
    facet_wrap(~ library, nrow = 2) #+coord_flip()

ggplot(data = fdr_all, mapping = aes(x = method, y = FDR, color = method)) +
    stat_summary(fun.min = min, fun.max = max, fun = mean) + 
    geom_point(alpha = 1 / 5, position = "jitter") + 
    facet_wrap(~ library, nrow = 2) #+coord_flip()