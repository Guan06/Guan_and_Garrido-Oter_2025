### #!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

source("../01.common_scripts/plot_setting.R")
source("../01.common_scripts/pcoa_function.R")
source("../01.common_scripts/order.R")

## for panel a, alpha diversity of bacterial community
alpha_file <- "../00.data/Bac_alpha_mean.txt"
design_file <- "../00.data/Bac_design_3809.txt"

alpha <- read.table(alpha_file, header = T, sep = "\t")

design <- read.table(design_file, header = T, sep = "\t")

host_lst <- c("zea_mays_l", "lotus_japonicus",
              "arabidopsis_thaliana", "chlamydomonas_reinhardtii")

design <- design[design$Host.Species %in% host_lst, ]
design$Group <- paste0(design$Host.Species, "_", design$Soil, "_",
                       design$Compartment)

inter <- intersect(alpha$Sample_ID, design$Original.SampleID)

alpha <- alpha[alpha$Sample_ID %in% inter, ]

design <- design[design$Original.SampleID %in% inter, ]
design <- design[match(alpha$Sample_ID, design$Original.SampleID), ]
design$Com_ST <- paste0(design$Compartment, "_", design$Soil)

alpha$Sample_ID <- design$Sample_ID
alpha_design <- merge(alpha, design)

alpha_design$Shannon <- as.numeric(alpha_design$Shannon)


p_a <- ggplot(alpha_design, aes(x = Group, y = Shannon, color = Compartment,
                                fill = Compartment, shape = Host.Species)) +
  geom_boxplot(outlier.size = 0.3, outlier.shape = NULL) +
  scale_color_manual(values = c_Com) +
  scale_fill_manual(values = c_Com) +
  scale_shape_manual(values = s_Host) +
  scale_x_discrete(limits = plot_order, labels = bac_labels) +
  main_theme +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8))

### for panel b, alpha diversity of fungal community

alpha_file <- "../00.data/Fun_alpha_mean.txt"
design_file <- "../00.data/Fun_design_2232.txt"

alpha <- read.table(alpha_file, header = T, sep = "\t")
design <- read.table(design_file, header = T, sep = "\t")

host_lst <- c("zea_mays_l", "lotus_japonicus",
              "arabidopsis_thaliana")

design <- design[design$Host.Species %in% host_lst, ]
design$Group <- paste0(design$Host.Species, "_", design$Soil, "_",
                       design$Compartment)

design$Com_ST <- paste0(design$Compartment, "_", design$Soil)
alpha_design <- merge(alpha, design)

alpha_design$Shannon <- as.numeric(alpha_design$Shannon)

p_b <- ggplot(alpha_design, aes(x = Group, y = Shannon, color = Compartment,
                                fill = Compartment, shape = Host.Species)) +
  geom_boxplot(outlier.size = 0.3, outlier.shape = NULL) +
  scale_color_manual(values = c_Com) +
  scale_fill_manual(values = c_Com) +
  scale_shape_manual(values = s_Host) +
  scale_x_discrete(limits = fun_plot_order, labels = fun_labels) +
  main_theme +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8))

p <- plot_grid(p_a, p_b, nrow = 1, labels = 'auto',
               align = "h", axis = "b")

ggsave("Figure_S1.pdf", p, height = 3, width = 12)
