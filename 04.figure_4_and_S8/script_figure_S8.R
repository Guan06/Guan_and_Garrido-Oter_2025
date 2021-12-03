#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(mina)
library(cowplot)
library(ggridges)
library(parallelDist)
source("../01.common_scripts/plot_setting.R")
source("../01.common_scripts/pcoa_function.R")
source("../01.common_scripts/sig_test.R")

design_file <- "../00.data/design_854.txt"
ra_file <- "RepASVs_RA.txt"
alpha_file <- "../00.data/Bac_alpha_mean_repASVs.txt"
asv_file <- "../00.data/Bac_ASV_rep.rds"

### setting the orders for panel a and b
this_order <- c("soil",
                "lotus_japonicus_rhizosphere",
                "lotus_japonicus_root",
                "chlamydomonas_reinhardtii_phycosphere",
                "arabidopsis_thaliana_rhizosphere",
                "arabidopsis_thaliana_rhizoplane",
                "arabidopsis_thaliana_root")

this_label <- c("CAS", "Lj_rhizos", "Lj_root",
                "Chlamy_phyco",
                "At_rhizos",
                "At_rhizop", "At_root")
## r2
bio_factors <- c("Compartment", "Soil.Batch",
                  "Host.Species", "Host.Genotype", "Conditions")
tech_factors <- c("Study", "Run")

###

## read in design and get Group info
design <- read.table(design_file, header = T, sep = "\t")
des_soil <- design[design$Compartment == "soil", ]
des_soil$Group <- "soil"

des_nsoil <- design[design$Compartment != "soil", ]
des_nsoil$Group <- paste0(des_nsoil$Host.Species, "_", des_nsoil$Compartment)
des <- rbind(des_soil, des_nsoil)

ra <- read.table(ra_file, header = T, sep = "\t")
ra <- ra[, 1:2]
ra <- merge(ra, des)

p_a <- ggplot(ra, aes(x = RA_rep, y = Group, fill = Compartment,
                     shape = Host.Species), color = "gray48") +
    geom_density_ridges(scale = 0.5, quantile_lines = TRUE, quantiles = 2,
                        color = "gray48") +
    geom_jitter(aes(shape = Host.Species, color = Compartment),
                size = 1, alpha = 0.7,
                position = position_jitter(width = 0.1)) +
    scale_color_manual(values = c_Com) +
    scale_fill_manual(values = c_Com) +
    scale_shape_manual(values = s_Host) +
    scale_y_discrete(limits = this_order, labels = this_label) +
    main_theme +
    labs(y = "", x = "aRA fof repASVs") +
    theme(legend.position = "none",
          axis.text.x = element_blank()) +
    coord_flip()

### for panel b
alpha <- read.table(alpha_file, header = T, sep = "\t")
alpha <- merge(alpha, des)
p_b <- ggplot(alpha, aes(x = Shannon, y = Group, fill = Compartment,
                     shape = Host.Species), color = "gray48") +
    geom_density_ridges(scale = 0.5, quantile_lines = TRUE, quantiles = 2,
                        color = "gray48") +
    geom_jitter(aes(shape = Host.Species, color = Compartment),
                size = 1, alpha = 0.7,
                position = position_jitter(width = 0.1)) +
    scale_color_manual(values = c_Com) +
    scale_fill_manual(values = c_Com) +
    scale_shape_manual(values = s_Host) +
    scale_y_discrete(limits = this_order) +
    main_theme +
    labs(y = "") +
    theme(legend.position = "none",
          axis.text.x = element_blank()) +
    coord_flip()

## for panel c
asv <- readRDS(asv_file)
asv <- asv[, colnames(asv) %in% des$Sample_ID]
bc <- parDist(t(asv), method = "bray", threads = 40)
bc <- as.matrix(bc)

print(get_r2(bc, design, group = bio_factors))
print(get_r2(bc, design, group = c(bio_factors, tech_factors)))

dmr <- cmdscale(bc, k = 4, eig = T)
p_c <- pcoa(dmr, design, 12, "Compartment", "Host.Species", size = 2.6)
p_c <- p_c + theme(legend.position = "none")

## putting together panels
p_ab <- plot_grid(p_a, p_b, labels = 'auto',
                  ncol = 1,
                  align = 'v', axis = 'l')

p <- plot_grid(p_ab, p_c, labels = c('', 'c'), nrow = 1,
               rel_widths = c(1.1, 1))

ggsave("Figure_S8.pdf", p, height = 5, width = 12)
