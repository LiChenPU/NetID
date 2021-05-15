library(ggplot2)

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


