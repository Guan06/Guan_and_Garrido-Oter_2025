source("../01.common_scripts/plot_setting.R")
source("../01.common_scripts/pcoa_function.R")
source("../01.common_scripts/order.R")

### Panel a, alpha-diversity of bacterial community calculated from repASVs
################################################################################
design_file <- "../00.data/Bac_design_3809.txt"
alpha_rep_file <- "../00.data/Bac_alpha_mean_repASVs.txt"

alpha_rep <- read.table(alpha_rep_file, header = T, sep = "\t")

design <- read.table(design_file, header = T, sep = "\t")

host_lst <- c("zea_mays_l", "lotus_japonicus",
              "arabidopsis_thaliana", "chlamydomonas_reinhardtii")

design <- design[design$Host.Species %in% host_lst, ]
design$Group <- paste0(design$Host.Species, "_", design$Soil, "_",
                       design$Compartment)
design$Com_ST <- paste0(design$Compartment, "_", design$Soil)


alpha_rep_design <- merge(alpha_rep, design)
alpha_rep_design$Shannon <- as.numeric(alpha_rep_design$Shannon)

p_a <- ggplot(alpha_rep_design, aes(x = Group, y = Shannon, color = Compartment,
                                    fill = Compartment, shape = Host.Species)) +
  geom_boxplot(outlier.size = 0.3, outlier.shape = NULL) +
  scale_color_manual(values = c_Com) +
  scale_fill_manual(values = c_Com) +
  scale_shape_manual(values = s_Host) +
  scale_x_discrete(limits = plot_order, labels = bac_labels) +
  main_theme +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8))
################################################################################

### Panel b, alpha-diversity of fungal community calculated from repASVs
################################################################################
alpha_rep_file <- "../00.data/Fun_alpha_mean_rep.txt"
design_file <- "../00.data/Fun_design_2232.txt"

alpha_rep <- read.table(alpha_rep_file, header = T, sep = "\t")
design <- read.table(design_file, header = T, sep = "\t")

host_lst <- c("zea_mays_l", "lotus_japonicus",
              "arabidopsis_thaliana")

design <- design[design$Host.Species %in% host_lst, ]
design$Group <- paste0(design$Host.Species, "_", design$Soil, "_",
                       design$Compartment)

design$Com_ST <- paste0(design$Compartment, "_", design$Soil)

alpha_rep_design <- merge(alpha_rep, design)
alpha_rep_design$Shannon <- as.numeric(alpha_rep_design$Shannon)
p_b <- ggplot(alpha_rep_design, aes(x = Group, y = Shannon, color = Compartment,
                                    fill = Compartment, shape = Host.Species)) +
  geom_boxplot(outlier.size = 0.3, outlier.shape = NULL) +
  scale_color_manual(values = c_Com) +
  scale_fill_manual(values = c_Com) +
  scale_shape_manual(values = s_Host) +
  scale_x_discrete(limits = fun_plot_order, labels = fun_labels) +
  main_theme +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8))
################################################################################

### Panel c-d, Beta-diversity calculated from repASVs
################################################################################
asv_file <- "../00.data/Bac_ASV_rep.rds"
design_file <- "../00.data/Bac_design_3809.txt"

asv <- readRDS(asv_file)
#bc <- parDist(t(asv), method = "bray", threads = 40)
bc <- vegan::vegdist(t(asv), method = "bray")
bc <- as.matrix(bc)

design <- read.table(design_file, header = T, sep = "\t")
dmr <- cmdscale(bc, k = 4, eig = T)
p_g <- pcoa(dmr, design, 12, "Compartment", "Host.Species")
p_S4c <- p_g + theme(legend.position = "none")

asv_file <- "../00.data/Fun_ASV_rep.rds"
design_file <- "../00.data/Fun_design_2232.txt"

asv <- readRDS(asv_file)
#bc <- parDist(t(asv), method = "bray", threads = 40)
bc <- vegan::vegdist(t(asv), method = "bray")
bc <- as.matrix(bc)

design <- read.table(design_file, header = T, sep = "\t")
dmr <- cmdscale(bc, k = 4, eig = T)
p_h <- pcoa(dmr, design, 12, "Compartment", "Host.Species")
p_S4d <- p_h + theme(legend.position = "none")

## merge panels 
################################################################################
p <- plot_grid(p_a, p_b, p_S4c, p_S4d, ncol = 2, labels = 'auto',
          align = "hv", axis = "lb", rel_heights = c(1, 1.5))

ggsave("Figure_S4.pdf", p, width = 10, height = 8)
