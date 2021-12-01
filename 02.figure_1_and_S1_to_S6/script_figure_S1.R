#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

source("../01.common_scripts/plot_setting.R")
source("../01.common_scripts/pcoa_function.R")

## for figure panel A
bc_file <- "../00.data/Bac_ASV_bc.rds"
design_file <- "../00.data/Bac_design_3809.txt"

bc <- readRDS(bc_file)
design <- read.table(design_file, header = T, sep = "\t")
dmr <- cmdscale(bc, k = 4, eig = T)

p_a <- pcoa(dmr, design, 12, "Compartment", "Host.Species")

## for figure panel B
bc_file <- "../00.data/Fun_ASV_bc.rds"
design_file <- "../00.data/Fun_design_2232.txt"

bc <- readRDS(bc_file)
design <- read.table(design_file, header = T, sep = "\t")
dmr <- cmdscale(bc, k = 4, eig = T)
p_b <- pcoa(dmr, design, 12, "Compartment", "Host.Species")
p_b <- p_b + theme(legend.position = "none")

## organize the panels with cowplot
library(cowplot)
p_a_legend <- get_legend(p_a)
p_a <- p_a + theme(legend.position = "none")

p <- plot_grid(p_a, p_b, p_a_legend, nrow = 1,
                   labels = c('a', 'b', ''),
                   rel_widths = c(1, 1, 0.3),
                   rel_heights = c(1, 1, 1))

ggsave("Figure_S1.pdf", p, height = 4.3, width = 10)
