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
for (i in 1 : 20) {
    ## top i% RA
    top_i <- round(nrow(asv) * i / 100)

        ra <- rowSums(asv) / ncol(asv)
        ra <- ra[order(-ra)]
        ra_top <- names(ra[1 : top_i])

        this_rep <- ra_top

        this_asv <- asv[rownames(asv) %in% this_rep, ]
        ## re-normalization
        this_asv <- apply(this_asv, 2, function(x) x/sum(x))
        this_bc <- as.matrix(parDist(t(this_asv), method = "bray", threads = 80))

        this_pro <- procrustes(all_bc, this_bc, symmetric = TRUE, digits = 3)

        this_file <- paste0("RA_", i, ".rds")

        saveRDS(this_pro, this_file)
        print(paste(i, this_pro$ss))
        log <- rbind(log, c(i, this_pro$ss))

}

colnames(log) <- c("RA", "M2")
write.table(log, m2_file, row.names = F, sep = "\t", quote = F)
