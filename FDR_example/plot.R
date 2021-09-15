library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(dplyr)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
load("FDR.RData")

# Functions ####
{
  MM_2_GT <- function(mm, library, rep) {
    ni <- sapply(mm, function(x) x[2, 1])
    mz <- sapply(mm, function(x) x[2, 2])
    node <- sapply(mm, function(x) x[2, 3])
    edge <- sapply(mm, function(x) x[2, 4])
    gt <- data.frame(precision = c(ni, mz, node, edge),
                     method = rep(c('ni', 'mz', 'node', 'edge'), each = rep),
                     library = rep(library, rep * 4))
  }
  
  MM_2_FDR <- function(mm, library, rep) {
    ni <- sapply(mm, function(x) x[1, 1])
    mz <- sapply(mm, function(x) x[1, 2])
    node <- sapply(mm, function(x) x[1, 3])
    edge <- sapply(mm, function(x) x[1, 4])
    fdr <- data.frame(FDR = c(ni, mz, node, edge),
                      method = rep(c('ni', 'mz', 'node', 'edge'), each = rep),
                      library = rep(library, rep * 4))
  }
}


## plot - Annotation precision by methods ####
{
  ## Data ####
  {
    gt_all <- rbind(MM_2_GT(hmdb_mm, "HMDB", REP), MM_2_GT(ymdb_mm, "YMDB", REP),
                    MM_2_GT(pbcm1_mm, "PBCM_BIO", REP), MM_2_GT(pbcm2_mm, "PBCM", REP))
    
    gt_all$method = factor(gt_all$method, levels = c('mz', 'node', 'edge', 'ni'))
    
    gt_all$method = as.character(gt_all$method)
    gt_all$method[gt_all$method=="mz"] = "m/z\nonly"
    gt_all$method[gt_all$method=="node"] = "node\nscores"
    gt_all$method[gt_all$method=="edge"] = "node+edge\nscores"
    gt_all$method[gt_all$method=="ni"] = "NetID\noptimization"
    gt_all$method = factor(gt_all$method, levels = c('m/z\nonly', 'node\nscores', 'node+edge\nscores', 'NetID\noptimization'))
    
    
  }
  
  ## Plot ####
  {
    p_all_precision = ggplot(data = gt_all, mapping = aes(x = method, y = precision, color = method)) + 
      geom_boxplot(outlier.shape = NA) + #geom_point(data = stan, mapping = aes(x = method, y = precision, fill = method), shape = 22, color = "black") +
      geom_point(alpha = 1 / 5, position = "jitter") +
      # ggtitle("Precision by Ground Truth") + 
      theme_classic(base_size = 12 # edit font size for all non-data text
      ) + 
      guides(color = guide_legend(
        title = "Annotation\nmethod",
        reverse = F
      )) +
      labs(x = NULL,
           # title = "U13C-glucose",
           y = "Fraction correct") +
      scale_y_continuous(limits = c(0.7, 1.05),
                         expand = c(0,0),
                         breaks = c(.6,.7,.8,.9, 1)
      ) +
      theme(plot.title = element_text(size = 12, hjust = 0.5),
            plot.margin = margin(0.2,0.2,0.2,0.2,"inch"),
            axis.text.x = element_text(angle = 0, hjust = .5, vjust = .7)
            # axis.ticks.x = element_blank()
      ) +
      # theme(plot.title = element_text(hjust = 0.5)) +
      facet_wrap(~ library, nrow = 2) #+coord_flip()
    
    p1_precision = ggplot(data = gt_all %>% filter(library=="HMDB"), mapping = aes(x = method, y = precision, color = method)) + 
      geom_boxplot(outlier.shape = NA) + #geom_point(data = stan, mapping = aes(x = method, y = precision, fill = method), shape = 22, color = "black") +
      geom_point(alpha = 1 / 5, position = "jitter") +
      # ggtitle("Precision by Ground Truth") +
      # guides(color = F) +
      guides(color = guide_legend(
        title = "Annotation\nmethod",
        reverse = F
      )) +
      labs(x = NULL,
           # title = "U13C-glucose",
           y = "Fraction correct") +
      scale_y_continuous(limits = c(0.7, 1.05),
                         expand = c(0,0),
                         breaks = c(.6,.7,.8,.9, 1)
      ) +
      theme_classic(base_size = 14 # edit font size for all non-data text
      ) + 
      theme(plot.title = element_text(size = 12, hjust = 0.5),
            plot.margin = margin(0.2,0.2,0.2,0.2,"inch"),
            axis.text.x = element_text(angle = 0, hjust = .5, vjust = .7)
            # axis.ticks.x = element_blank()
      ) 
    # theme(plot.title = element_text(hjust = 0.5)) +
    # facet_wrap(~ library, nrow = 2) #+coord_flip()
    print(p_all_precision)
    print(p1_precision)
  }
  
  ## Output ####
  {
    ggpubr::ggarrange(
      p1_precision,
      
      common.legend = T, legend = "none",
      # align = "hv",
      nrow = 1, ncol = 1
    ) %>%
      ggexport(filename = "plot_7_1.pdf", width = 4, height = 2.5)
    
    ggpubr::ggarrange(
      p_all_precision,
      
      common.legend = T, legend = "none",
      # align = "hv",
      nrow = 1, ncol = 1
    ) %>%
      ggexport(filename = "plot_7_1_all.pdf", width = 8, height = 5)
  }
}
## plot - FDR by methods ####
{
  ## Data ####
  {
    fdr_all <- rbind(MM_2_FDR(hmdb_mm, "HMDB", REP), MM_2_FDR(ymdb_mm, "YMDB", REP),
                     MM_2_FDR(pbcm1_mm, "PBCM_BIO", REP), MM_2_FDR(pbcm2_mm, "PBCM", REP))
    
    fdr_all$method = factor(fdr_all$method, levels = c('mz', 'node', 'edge', 'ni'))
    
    fdr_all$method = as.character(fdr_all$method)
    fdr_all$method[fdr_all$method=="mz"] = "m/z\nonly"
    fdr_all$method[fdr_all$method=="node"] = "node\nscores"
    fdr_all$method[fdr_all$method=="edge"] = "node+edge\nscores"
    fdr_all$method[fdr_all$method=="ni"] = "NetID\noptimization"
    fdr_all$method = factor(fdr_all$method, levels = c('m/z\nonly', 'node\nscores', 'node+edge\nscores', 'NetID\noptimization'))
    
    
  }
  
  ## Plot ####
  {
    p_all_FDR = ggplot(data = fdr_all, mapping = aes(x = method, y = FDR, color = method)) + 
      geom_boxplot(outlier.shape = NA) + #geom_point(data = stan, mapping = aes(x = method, y = precision, fill = method), shape = 22, color = "black") +
      geom_point(alpha = 1 / 5, position = "jitter") +
      # ggtitle("Precision by Ground Truth") + 
      theme_classic(base_size = 12 # edit font size for all non-data text
      ) + 
      guides(color = guide_legend(
        title = "Annotation\nmethod",
        reverse = F
      )) +
      labs(x = NULL,
           # title = "U13C-glucose",
           y = "False discovery rate ") +
      scale_y_continuous(limits = c(0, 0.2),
                         expand = c(0,0),
                         breaks = c(0, 0.05, 0.1, 0.15, 0.2)
      ) +
      theme(plot.title = element_text(size = 12, hjust = 0.5),
            plot.margin = margin(0.2,0.2,0.2,0.2,"inch"),
            axis.text.x = element_text(angle = 0, hjust = .5, vjust = .7)
            # axis.ticks.x = element_blank()
      ) +
      # theme(plot.title = element_text(hjust = 0.5)) +
      facet_wrap(~ library, nrow = 2) #+coord_flip()
    
    p1_FDR = ggplot(data = fdr_all %>% filter(library=="HMDB"), mapping = aes(x = method, y = FDR, color = method)) + 
      geom_boxplot(outlier.shape = NA) + #geom_point(data = stan, mapping = aes(x = method, y = precision, fill = method), shape = 22, color = "black") +
      geom_point(alpha = 1 / 5, position = "jitter") +
      # ggtitle("Precision by Ground Truth") + 
      guides(color = guide_legend(
        title = "Annotation\nmethod",
        reverse = F
      )) +
      labs(x = NULL,
           # title = "U13C-glucose",
           y = "False discovery rate ") +
      scale_y_continuous(limits = c(0, 0.2),
                         expand = c(0,0),
                         breaks = c(0, 0.05, 0.1, 0.15, 0.2)
      ) +
      theme_classic(base_size = 14 # edit font size for all non-data text
      ) + 
      theme(plot.title = element_text(size = 12, hjust = 0.5),
            plot.margin = margin(0.2,0.2,0.2,0.2,"inch"),
            axis.text.x = element_text(angle = 0, hjust = .5, vjust = .7)
            # axis.ticks.x = element_blank()
      ) 
    # theme(plot.title = element_text(hjust = 0.5)) +
    # facet_wrap(~ library, nrow = 2) #+coord_flip()
    print(p_all_FDR)
    print(p1_FDR)
  }
  
  ## Output ####
  {
    ggpubr::ggarrange(
      p1_FDR,
      
      common.legend = T, legend = "none",
      # align = "hv",
      nrow = 1, ncol = 1
    ) %>%
      ggexport(filename = "plot_7_2.pdf", width = 4, height = 2.5)
    
    ggpubr::ggarrange(
      p_all_FDR,
      
      common.legend = T, legend = "none",
      # align = "hv",
      nrow = 1, ncol = 1
    ) %>%
      ggexport(filename = "plot_7_2_all.pdf", width = 8, height = 5)
  }
}

## output ####
{
  ggpubr::ggarrange(
    p1_FDR,
    p1_precision,
    common.legend = F, legend = "none",
    align = "hv",
    nrow = 2, ncol = 1
  ) %>%
    ggexport(filename = "plot_merge.pdf", width = 4, height = 5)
}
