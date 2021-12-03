#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(mina)
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

bc <- parDist(t(asv), method = "bray", threads = 40)
bc <- as.matrix(bc)

## r2
bio_factors <- c("Compartment", "Soil.Batch",
                  "Host.Species", "Host.Genotype", "Conditions")
tech_factors <- c("Study", "Run")


print(get_r2(bc, design, group = bio_factors))
print(get_r2(bc, design, group = c(bio_factors, tech_factors)))

dmr <- cmdscale(bc, k = 4, eig = T)
p_a <- pcoa(dmr, design, 12, "Compartment", "Host.Species", size = 2.4)
p_a <- p_a + theme(legend.position = "none")

### for alpha diversity
des_soil <- design[design$Compartment == "soil", ]
des_soil$Group <- "soil"

des_nsoil <- design[design$Compartment != "soil", ]
des_nsoil$Group <- paste0(des_nsoil$Host.Species, "_", des_nsoil$Compartment)

des <- rbind(des_soil, des_nsoil)
alpha <- merge(alpha, des, by.x = "Sample_ID", by.y = "Original.SampleID")

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
    scale_y_discrete(limits = this_order, labels = this_label) +
    main_theme +
    labs(y = "") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust=1, size = 8)) +
    coord_flip()

## function for significant test
sig <- box_sig(alpha, "Group", "Shannon")
sig$FDR <- p.adjust(sig$Significance, method = "fdr")
write.table(sig, "Figure_4b.txt", row.names = F, sep = "\t", quote = F)

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

all_left <- plot_grid(p_a, p_b, labels = c('a', 'b'), ncol = 1, align = "v",
                      axis = "l", rel_heights = c(1, 0.7))
all <- plot_grid(all_left, p_c, nrow = 1, labels = c('', "c"), align = "h",
                 rel_widths = c(1, 1.3))

ggsave("Figure_4.pdf", all, height = 7, width = 12)
