#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(cowplot)
library(ggridges)
library(parallelDist)

source("../01.common_scripts/plot_setting.R")
source("../01.common_scripts/dual_plot_functions.R")
source("../01.common_scripts/order.R")
source("../01.common_scripts/pcoa_function.R")

### for panel a, bacteria
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

### for panel b, fungi
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

p_ab <- plot_grid(p_a, p_b, nrow = 1, labels = 'auto')

### for panel c and d
ra_file <- "./procrustes_analysis/fixed_Occu/Bac_RA_m2.txt"
occu_file <- "./procrustes_analysis/fixed_RA/Bac_Occu_m2.txt"

ra <- read.table(ra_file, header = T, sep = "\t")
occu <- read.table(occu_file, header = T, sep = "\t")

colnames(ra)[1] <- colnames(occu)[1] <- "Percent"

ra$Factor <- "Relative abundance"
occu$Factor <- "Occupancy"

dat <- rbind(ra, occu)

p_c <- ggplot(dat, aes(Percent, M2, color = Factor)) +
    geom_line() + geom_point() +
    scale_color_manual(values = c("#FFD662FF", "#00539CFF")) +
    main_theme +
    theme(legend.position = "none") +
    geom_segment(x = 6, y = 0.005, xend = 6, yend = 0.035,
                 color = "#00539CFF", linetype = "dashed") +
    geom_segment(x = 6, y = 0.015, xend = 6, yend = 0.045,
                 color = "#FFD662FF", linetype = "dashed")

###
ra_file <- "./procrustes_analysis/fixed_Occu/Fun_RA_m2.txt"
occu_file <- "./procrustes_analysis/fixed_RA/Fun_Occu_m2.txt"

ra <- read.table(ra_file, header = T, sep = "\t")
occu <- read.table(occu_file, header = T, sep = "\t")

colnames(ra)[1] <- colnames(occu)[1] <- "Percent"

ra$Factor <- "RA"
occu$Factor <- "Occupancy"

dat <- rbind(ra, occu)

p_d <- ggplot(dat, aes(Percent, M2, color = Factor)) +
    geom_line() + geom_point() +
    scale_color_manual(values = c("#FFD662FF", "#00539CFF")) +
    main_theme +
    theme(legend.position = "none") +
    geom_segment(x = 6, xend = 6, y = 0, yend = 0.03,
                 color = "#00539CFF", linetype = "dashed") +
    geom_segment(x = 7, xend = 7, y = 0.025, yend = 0.055,
                 color = "#FFD662FF", linetype = "dashed")

p_cd <- plot_grid(p_c, p_d, nrow = 1, labels = c('c', 'd'))

#### for panel e, f
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
    geom_density_ridges(scale = 0.5, quantile_lines = TRUE, quantiles = 2,
                        color = "gray48") +
    geom_jitter(aes(shape = Host.Species),
                size = 1.4, color = "gray48", alpha = 0.6,
                position = position_jitter(width = 0.1)) +
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

p_f <- ggplot(RA_p1, aes(x = RA_rep, y = Compartment, fill = Compartment),
             color = "gray48") +
    geom_density_ridges(scale = 0.5, quantile_lines = TRUE, quantiles = 2,
                        color = "gray48") +
    geom_jitter(aes(shape = Host.Species),
                size = 1.4, color = "gray48", alpha = 0.7,
                position = position_jitter(width = 0.1)) +
    scale_fill_manual(values = c_Com) +
    scale_color_manual(values = c_Com) +
    scale_shape_manual(values = s_Host) +
    scale_y_discrete(limits = c("soil", "rhizosphere",
                                "rhizoplane", "root")) +
    xlim(c(0,1)) + main_theme +
    theme(legend.position = "none") +
    coord_flip()

p_ef <- plot_grid(p_e, p_f, nrow = 1, labels = c('e', 'f'))

### for panel g, h, i
asv_file <- "../00.data/Bac_ASV_rep.rds"
design_file <- "../00.data/Bac_design_3809.txt"

asv <- readRDS(asv_file)
bc <- parDist(t(asv), method = "bray", threads = 40)
bc <- as.matrix(bc)

design <- read.table(design_file, header = T, sep = "\t")
dmr <- cmdscale(bc, k = 4, eig = T)
p_g <- pcoa(dmr, design, 12, "Compartment", "Host.Species")
p_g <- p_g + theme(legend.position = "none")
bac_rep_bc <- bc

asv_file <- "../00.data/Fun_ASV_rep.rds"
design_file <- "../00.data/Fun_design_2232.txt"

asv <- readRDS(asv_file)
bc <- parDist(t(asv), method = "bray", threads = 40)
bc <- as.matrix(bc)

design <- read.table(design_file, header = T, sep = "\t")
dmr <- cmdscale(bc, k = 4, eig = T)
p_h <- pcoa(dmr, design, 12, "Compartment", "Host.Species")
p_h <- p_h + theme(legend.position = "none")
fun_rep_bc <- bc

bac_bc_file <- "../00.data/Bac_ASV_bc.rds"
fun_bc_file <- "../00.data/Fun_ASV_bc.rds"

bac_bc <- readRDS(bac_bc_file)
fun_bc <- readRDS(fun_bc_file)
c_bac <- bac_rep_bc / bac_bc
c_fun <- fun_rep_bc / fun_bc

flatten_mat <- function(x) {
    ut <- upper.tri(x)
    data.frame(
        sample_1 = rownames(x)[row(x)[ut]],
        sample_2 = rownames(x)[col(x)[ut]],
        C = x[ut]
    )
}

c_bac_df <- flatten_mat(c_bac)
c_fun_df <- flatten_mat(c_fun)

c_bac_df$Kingdom <- "Bacteria"
c_fun_df$Kingdom <- "Fungi"
c_df <- rbind(c_bac_df, c_fun_df)

p_i <- ggplot(c_df, aes(Kingdom, C)) +
    geom_violin(color = "gray64", fill = "gray48") +
    ylim(0, 1.2) +
    theme(legend.position = "none") + labs(x = '') +
    main_theme

p_ghi <- plot_grid(p_g, p_h, p_i, nrow = 1, rel_widths = c(1, 1, 0.45),
                   labels = c('g', 'h', 'i'), rel_heights = c(1, 1, 1))

p <- plot_grid(p_ab, p_cd, p_ef, p_ghi, ncol = 1,
               rel_heights = c(4, 2, 3, 3.5))

ggsave("Figure_1.pdf", p, width = 12, height = 12.5)
