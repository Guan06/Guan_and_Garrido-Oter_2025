## define the R used in the scripts
#### #!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(cowplot)
library(ggridges)

source("../01.common_scripts/plot_setting.R")
source("../01.common_scripts/dual_plot_functions.R")
source("../01.common_scripts/order.R")
source("../01.common_scripts/pcoa_function.R")

#### Plotting panel a and b in Figure 1
################################################################################
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
p_a_legend <- get_legend(p_a)
p_a <- p_a + theme(legend.position = "none")

p_ab <- plot_grid(p_a, p_b, p_a_legend, nrow = 1,
               labels = c('a', 'b', ''),
               rel_widths = c(1, 1, 0.3),
               rel_heights = c(1, 1, 1))
################################################################################


#### Plotting panel c and d in Figure 1
################################################################################
### for panel c, bacteria
asv_file <- "../00.data/Bac_ASV_raref.rds"
design_file <- "../00.data/Bac_design_3809.txt"

asv <- readRDS(asv_file)
des <- read.table(design_file, header = T, sep = "\t")
des$Group <- des$Compartment

p_a1 <- plot_occu_RA_ASV("soil")
p_a2 <- plot_occu_RA_ASV("rhizosphere")
p_a3 <- plot_occu_RA_ASV("rhizoplane")
p_a4 <- plot_occu_RA_ASV("root")

p_a <- plot_grid(p_a1, p_a2, p_a3, p_a4, ncol = 1, align = "vl")

### for panel d, fungi
asv_file <- "../00.data/Fun_ASV_raref.rds"
design_file <- "../00.data/Fun_design_2232.txt"

asv <- readRDS(asv_file)
des <- read.table(design_file, header = T, sep = "\t")
des$Group <- des$Compartment

p_b1 <- plot_occu_RA_ASV("soil")
p_b2 <- plot_occu_RA_ASV("rhizosphere")
p_b3 <- plot_occu_RA_ASV("rhizoplane")
p_b4 <- plot_occu_RA_ASV("root")

p_b <- plot_grid(p_b1, p_b2, p_b3, p_b4, ncol = 1, align = "vl")

p_cd <- plot_grid(p_a, p_b, nrow = 1, labels = 'auto')
################################################################################

#### Plotting panel e and f in Figure 1
################################################################################
asv_file <- "../00.data/Bac_ASV_raref.rds"
rep_file <- "../00.data/Bac_ASV_rep.rds"
design_file <- "../00.data/Bac_design_3809.txt"

asv <- readRDS(asv_file)
rep <- readRDS(rep_file)
design <- read.table(design_file, header = T, sep = "\t")

asv <- asv[rownames(asv) %in% rownames(rep), ]

RA0 <- data.frame(Sample_ID = colnames(asv),
                 RA_rep = colSums(asv))

design <- design[, colnames(design) %in% c("Sample_ID", "Host.Species",
                                           "Compartment", "Soil")]

RA <- merge(RA0, design)
ra_file <- "Figure_1e.txt"
write.table(RA, ra_file, quote = F, sep = "\t", row.names = F)

RA_p1 <- RA[RA$Compartment %in%
            c("soil", "rhizosphere", "rhizoplane", "root"), ]

p_e <- ggplot(RA_p1, aes(x = RA_rep, y = Compartment, fill = Compartment),
             color = "gray48") +
    geom_jitter(aes(shape = Host.Species),
                size = 1.4, color = "gray48", alpha = 0.6,
                position = position_jitter(width = 0.1)) +
    geom_density_ridges(scale = 0.5, quantile_lines = TRUE, quantiles = 2,
                      color = "gray48", alpha = 0.6) +
    scale_fill_manual(values = c_Com) +
    scale_color_manual(values = c_Com) +
    scale_shape_manual(values = s_Host) +
    scale_y_discrete(limits = c("soil", "rhizosphere",
                                "rhizoplane", "root")) +
    xlim(c(0,1)) + main_theme +
    theme(legend.position = "none") +
    coord_flip()

## for fungi
asv_file <- "../00.data/Fun_ASV_raref.rds"
rep_file <- "../00.data/Fun_ASV_rep.rds"
design_file <- "../00.data/Fun_design_2232.txt"

asv <- readRDS(asv_file)
rep <- readRDS(rep_file)
design <- read.table(design_file, header = T, sep = "\t")

asv <- asv[rownames(asv) %in% rownames(rep), ]

RA0 <- data.frame(Sample_ID = colnames(asv),
                 RA_rep = colSums(asv))

design <- design[, colnames(design) %in% c("Sample_ID", "Host.Species",
                                           "Compartment", "Soil")]

RA <- merge(RA0, design)
ra_file <- "Figure_1f.txt"
write.table(RA, ra_file, quote = F, sep = "\t", row.names = F)

RA_p2 <- RA[RA$Compartment %in%
            c("soil", "rhizosphere", "rhizoplane", "root"), ]

p_f <- ggplot(RA_p2, aes(x = RA_rep, y = Compartment, fill = Compartment),
             color = "gray48") +
    geom_jitter(aes(shape = Host.Species),
                size = 1.4, color = "gray48", alpha = 0.7,
                position = position_jitter(width = 0.1)) +
    geom_density_ridges(scale = 0.5, quantile_lines = TRUE, quantiles = 2,
                      color = "gray48", alpha = 0.6) +
    scale_fill_manual(values = c_Com) +
    scale_color_manual(values = c_Com) +
    scale_shape_manual(values = s_Host) +
    scale_y_discrete(limits = c("soil", "rhizosphere",
                                "rhizoplane", "root")) +
    xlim(c(0,1)) + main_theme +
    theme(legend.position = "none") +
    coord_flip()

p_ef <- plot_grid(p_e, p_f, nrow = 1, labels = c('e', 'f'))
################################################################################


## Merge the panels
p <- plot_grid(p_ab, p_cd, p_ef, ncol = 1,
               rel_heights = c(3.5, 4, 2))

ggsave("Figure_1.pdf", p, width = 12, height = 12.5)
