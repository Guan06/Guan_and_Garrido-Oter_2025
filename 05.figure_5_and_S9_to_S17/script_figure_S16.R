#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(ggridges)
library(cowplot)
source("../01.common_scripts/plot_setting.R")

dis_folder1 <- "./Figure_S14/Network_cluster/spectra_dis/"
dis_folder2 <- "./Figure_S14/Family/spectra_dis/"
dis_folder3 <- "./Figure_S16/netcls_egv/"
dis_folder4 <- "./Figure_S16/family_egv/"

## for CAS vs At root, with both network cluster and family group, panel a
#x_order <- c("All_BS", "After_Cluster_67", "Cluster_67",
#             "After_Nitrospiraceae", "Nitrospiraceae", "All_PM")

x_order <- c("All_BS", "Cluster_67", "Nitrospiraceae", "All_PM")
x <- "arabidopsis_thaliana_root_vs_CAS"

## for permutating groups *Before* changing point
this_cmp1 <- read.table(paste0(dis_folder1, "accu_dis_",
                               x, ".txt"), header = T)
this_cmp2 <- read.table(paste0(dis_folder2, "accu_dis_",
                               x, ".txt"), header = T)

this_cmp1 <- this_cmp1[this_cmp1$Family %in% c("All_BS", "Cluster_67"), ]

this_cmp2 <- this_cmp2[this_cmp2$Family %in% c("Nitrospiraceae", "All_PM"), ]

## for permutating groups *After* changing point
this_cmp3 <- read.table(paste0(dis_folder3, "CAS_vs_At_root/dis_spectra_pm_",
                               x, ".txt"), header = T)
this_cmp4 <- read.table(paste0(dis_folder4, "CAS_vs_At_root/dis_spectra_pm_",
                               x, ".txt"), header = T)

this_cmp3$Family <- "After_Cluster_67"
this_cmp3$Mean <- mean(this_cmp3$Distance)
this_cmp3$Compare <- paste0(this_cmp3$Group1, "_vs_", this_cmp4$Group2)

this_cmp4$Family <- "After_Nitrospiraceae"
this_cmp4$Mean <- mean(this_cmp4$Distance)
this_cmp4$Compare <- paste0(this_cmp4$Group1, "_vs_", this_cmp4$Group2)

this_cmp12 <- rbind(this_cmp1, this_cmp2)
this_cmp34 <- rbind(this_cmp3, this_cmp4)

this_cmp <- rbind(this_cmp12, this_cmp34)

p_a <- ggplot(this_cmp12, aes(x = Distance, y = Family)) +
        geom_jitter(shape = 1, color = "gray", alpha = 0.7, size = 0.4,
                    position = position_jitter(width = 0.05)) +
        geom_density_ridges(color = "gray48", fill = "#1B9E77", alpha = 0.7,
                            scale = 0.5, quantile_lines = TRUE, quantiles = 2) +
          scale_y_discrete(limits = x_order) +
         labs(y = "") +
            main_theme +
            theme(plot.background = element_blank(),
                  legend.position = "none",
                  plot.title = element_text(size = 8),
                axis.text.x = element_text(colour = "black", angle = 60,
                size = 8, hjust = 1)) +
        coord_flip()

## for CAS vs Lj root, with both network cluster and family group, panel b
#x_order <- c("All_BS", "After_Cluster_121", "Cluster_121",
#             "After_Myxococcaceae", "Myxococcaceae", "All_PM")

x_order <- c("All_BS", "Cluster_121", "Myxococcaceae", "All_PM")

x <- "CAS_vs_lotus_japonicus_root"

this_cmp1 <- read.table(paste0(dis_folder1, "accu_dis_",
                               x, ".txt"), header = T)
this_cmp2 <- read.table(paste0(dis_folder2, "accu_dis_",
                               x, ".txt"), header = T)

this_cmp1 <- this_cmp1[this_cmp1$Family %in% c("All_BS", "Cluster_121"), ]

this_cmp2 <- this_cmp2[this_cmp2$Family %in% c("Myxococcaceae", "All_PM"), ]

## for permutating groups *After* changing point
this_cmp3 <- read.table(paste0(dis_folder3, "CAS_vs_Lj_root/dis_spectra_pm_",
                               x, ".txt"), header = T)
this_cmp4 <- read.table(paste0(dis_folder4, "CAS_vs_Lj_root/dis_spectra_pm_",
                               x, ".txt"), header = T)

this_cmp3$Family <- "After_Cluster_121"
this_cmp3$Mean <- mean(this_cmp3$Distance)
this_cmp3$Compare <- paste0(this_cmp3$Group1, "_vs_", this_cmp4$Group2)

this_cmp4$Family <- "After_Myxococcaceae"
this_cmp4$Mean <- mean(this_cmp4$Distance)
this_cmp4$Compare <- paste0(this_cmp4$Group1, "_vs_", this_cmp4$Group2)

this_cmp12 <- rbind(this_cmp1, this_cmp2)
this_cmp34 <- rbind(this_cmp3, this_cmp4)

this_cmp <- rbind(this_cmp12, this_cmp34)

p_b <- ggplot(this_cmp12, aes(x = Distance, y = Family)) +
        geom_jitter(shape = 1, color = "gray", alpha = 0.7, size = 0.4,
                    position = position_jitter(width = 0.05)) +
        geom_density_ridges(color = "gray48", fill = "#D95F02", alpha = 0.7,
                            scale = 0.5, quantile_lines = TRUE, quantiles = 2) +
        scale_y_discrete(limits = x_order) +
        labs(y = "") +
        main_theme +
        theme(plot.background = element_blank(),
                  legend.position = "none",
                  plot.title = element_text(size = 8),
                axis.text.x = element_text(colour = "black", angle = 60,
                size = 8, hjust = 1)) +
        coord_flip()

p <- plot_grid(p_a, p_b, nrow = 1, labels = "auto")

#ggsave("Figure_S16.pdf", p, height = 3, width = 5)
ggsave("Figure_S16.pdf", p, height = 3, width = 6)
