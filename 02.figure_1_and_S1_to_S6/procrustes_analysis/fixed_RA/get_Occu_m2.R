#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(parallelDist)
library(vegan)

args = commandArgs(trailingOnly = TRUE)
all_bc <- args[1]
all_asv <- args[2]
m2_file <- args[3]

all_bc <- readRDS(all_bc)
asv <- readRDS(all_asv)

log <- c()

for (j in 1:20){
        ## top j% Occu
        top_j <- round(nrow(asv) * j / 100)

        occu <- rowSums(asv > 0) / ncol(asv)
        occu <- occu[order(-occu)]
        occu_top <- names(occu[1 : top_j])

        this_rep <- occu_top
        this_asv <- asv[rownames(asv) %in% this_rep, ]
        ## re-normalization
        this_asv <- apply(this_asv, 2, function(x) x/sum(x))
        this_bc <- as.matrix(parDist(t(this_asv), method = "bray", threads = 80))

        this_pro <- procrustes(all_bc, this_bc, symmetric = TRUE, digits = 3)

        this_file <- paste0("Occu_", j, ".rds")

        saveRDS(this_pro, this_file)
        log <- rbind(log, c(j, this_pro$ss))
}


colnames(log) <- c("Occu", "M2")
write.table(log, m2_file, row.names = F, sep = "\t", quote = F)
