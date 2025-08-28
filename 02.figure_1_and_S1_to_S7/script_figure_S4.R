
library(parallelDist)
library(devtools)
library(ggridges)
library(cowplot)
source("../01.common_scripts/plot_setting.R")
source("../01.common_scripts/pcoa_function.R")
source("../01.common_scripts/dual_plot_functions.R")
source("../01.common_scripts/sig_test.R")

asv_file <- "../00.data/Bac_ASV_raref.rds"
alpha_file <- "../00.data/Bac_alpha_mean.txt"
tax_file <- "../00.data/Bac_ASV_tax.txt"

asv <- readRDS(asv_file)
design <- read.table("../00.data/design_854.txt", header = T, sep = "\t")
alpha <- read.table(alpha_file, header = T, sep = "\t")

asv <- asv[, colnames(asv) %in% design$Sample_ID]
asv <- asv[rowSums(asv) > 0, ]

des_soil <- design[design$Compartment == "soil", ]
des_soil$Group <- "soil"

des_nsoil <- design[design$Compartment != "soil", ]
des_nsoil$Group <- paste0(des_nsoil$Host.Species, "_", des_nsoil$Compartment)

des <- rbind(des_soil, des_nsoil)

### for panel C
tax <- read.table(tax_file, header = T, sep = "\t")

p0 <- plot_occu_RA("soil")
p1 <- plot_occu_RA("chlamydomonas_reinhardtii_phycosphere")

p2a <- plot_occu_RA("lotus_japonicus_rhizosphere")
p2b <- plot_occu_RA("lotus_japonicus_root")

p4a <- plot_occu_RA("arabidopsis_thaliana_rhizosphere")
p4b <- plot_occu_RA("arabidopsis_thaliana_rhizoplane")
p4c <- plot_occu_RA("arabidopsis_thaliana_root")

p_c <- plot_grid(p0, p1, p2a, p2b,
                 p4a, p4b, p4c, ncol = 1, align = "v", axis = "l")

ggsave("Figure_S4.pdf", p_c, height = 7, width = 7)
