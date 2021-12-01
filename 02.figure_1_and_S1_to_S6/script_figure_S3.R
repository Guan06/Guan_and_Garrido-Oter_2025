#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(cowplot)

source("../01.common_scripts/plot_setting.R")
source("../01.common_scripts/dual_plot_functions.R")

### for bacteria
asv_file <- "../00.data/Bac_ASV_raref.rds"
design_file <- "../00.data/Bac_design_3809.txt"
tax_file <- "../00.data/Bac_ASV_tax.txt"

asv <- readRDS(asv_file)
des <- read.table(design_file, header = T, sep = "\t")
tax <- read.table(tax_file, header = T, sep = "\t")

des$Group <- des$Compartment

p1 <- plot_occu_RA_tax("soil", t = "Class")
p2 <- plot_occu_RA_tax("rhizosphere", t = "Class")
p3 <- plot_occu_RA_tax("rhizoplane", t = "Class")
p4 <- plot_occu_RA_tax("root", t = "Class")

p_bac <- plot_grid(p1, p2, p3, p4, ncol = 1, align = "vl",
                   labels = c('a', '', '', ''))

### for fungi
asv_file <- "../00.data/Fun_ASV_raref.rds"
design_file <- "../00.data/Fun_design_2232.txt"
tax_file <- "../00.data/Fun_ASV_tax.txt"

asv <- readRDS(asv_file)
des <- read.table(design_file, header = T, sep = "\t")
tax <- read.table(tax_file, header = T, sep = "\t")

des$Group <- des$Compartment

p1 <- plot_occu_RA_tax("soil", t = "Class")
p2 <- plot_occu_RA_tax("rhizosphere", t = "Class")
p3 <- plot_occu_RA_tax("rhizoplane", t = "Class")
p4 <- plot_occu_RA_tax("root", t = "Class")

p_fun <- plot_grid(p1, p2, p3, p4, ncol = 1, align = "vl",
                   labels = c('b', '', '', ''))

p <- plot_grid(p_bac, p_fun, nrow = 1)
ggsave("Figure_S3.pdf", p, height = 4, width = 8)
