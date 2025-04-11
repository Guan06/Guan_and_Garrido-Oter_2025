### #!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(cowplot)
source("../01.common_scripts/plot_setting.R")

#### for family ###############################################################
## permutation excluding "Unassigned" ASVs
dis_file <- "./Figure_S7/Family_single_pm_dis.txt"
dis_mean_file <- "./Figure_S7/Family_single_pm_dis_mean.txt"

dis <- read.table(dis_file, header = T, sep = "\t")
dis_mean <- read.table(dis_mean_file, header = T, sep = "\t")
dis_mean <- dis_mean[order(dis_mean$Decrease_Mean, decreasing = T), ]
dis_mean_order <- as.character(dis_mean$Family)

x_order <- dis_mean_order

dis$Color <- dis$Comparison
dis$Color[!(dis$Comparison %in% c("arabidopsis_thaliana_root_vs_CAS",
                             "CAS_vs_lotus_japonicus_root"))] <- "Others"
dis_p1 <- dis[, c("Family", "Decrease", "Comparison", "Color")]

dis_mean$Comparison <- dis_mean$Color <- "Mean"
dis_p2 <- dis_mean
colnames(dis_p2) <- c("Family", "Decrease", "Comparison", "Color")

dis_merge <- rbind(dis_p1, dis_p2)

p_b <- ggplot(dis_merge) +
    geom_line(aes(Family, Decrease, group = Comparison, color = Color)) +
    geom_point(aes(Family, Decrease, color = Color))+
    scale_x_discrete(limits = x_order) +
    scale_color_manual(values = c_cmp) +
    main_theme +
    geom_hline(yintercept = 0, color = "navyblue", linetype = "dashed") +
    theme(legend.position = "none",
              plot.title = element_text(size = 8),
            axis.text.x = element_text(colour = "black", angle = 90,
            size = 6, hjust = 1))

write.table(dis_merge, "Figure_S7b.txt",
            quote = F, sep = "\t", row.names = F)
## calculate correlation between CAS_vs_At_root and CAS_vs_Lj_root
d_at <- dis[dis$Comparison == "arabidopsis_thaliana_root_vs_CAS", 1:2]
d_lj <- dis[dis$Comparison == "CAS_vs_lotus_japonicus_root", 1:2]
colnames(d_at)[2] <- "Decrease_At"
colnames(d_lj)[2] <- "Decrease_Lj"
d12 <- merge(d_at, d_lj)
cor.test(d12$Decrease_At, d12$Decrease_Lj)

######## for network cluster ##################################################
dis_file <- "./Figure_S7/Cluster_single_pm_dis.txt"
dis_mean_file <- "./Figure_S7/Cluster_single_pm_dis_mean.txt"

dis <- read.table(dis_file, header = T, sep = "\t")
dis_mean <- read.table(dis_mean_file, header = T, sep = "\t")
dis_mean <- dis_mean[order(dis_mean$Decrease_Mean, decreasing = T), ]
dis_mean_order <- as.character(dis_mean$Cluster)

x_order <- dis_mean_order

dis$Color <- dis$Comparison

dis$Color[!(dis$Comparison %in% c("arabidopsis_thaliana_root_vs_CAS",
                             "CAS_vs_lotus_japonicus_root"))] <- "Others"

dis_p1 <- dis[, c("Cluster", "Decrease", "Comparison", "Color")]

dis_mean$Comparison <- dis_mean$Color <- "Mean"
dis_p2 <- dis_mean
colnames(dis_p2) <- c("Cluster", "Decrease", "Comparison", "Color")

dis_merge <- rbind(dis_p1, dis_p2)

p_a <- ggplot(dis_merge) +
    geom_line(aes(Cluster, Decrease, group = Comparison, color = Color)) +
    geom_point(aes(Cluster, Decrease, color = Color))+
    scale_x_discrete(limits = x_order) +
    scale_color_manual(values = c_cmp) +
    main_theme +
    geom_hline(yintercept = 0, color = "navyblue", linetype = "dashed") +
    theme(legend.position = "none",
              plot.title = element_text(size = 8),
            axis.text.x = element_text(colour = "black", angle = 90,
            size = 6, hjust = 1))

write.table(dis_merge, "Figure_S7a.txt",
            quote = F, sep = "\t", row.names = F)

d_at <- dis[dis$Comparison == "arabidopsis_thaliana_root_vs_CAS", 1:2]
d_lj <- dis[dis$Comparison == "CAS_vs_lotus_japonicus_root", 1:2]
colnames(d_at)[2] <- "Decrease_At"
colnames(d_lj)[2] <- "Decrease_Lj"
d12 <- merge(d_at, d_lj)
cor.test(d12$Decrease_At, d12$Decrease_Lj)

p <- plot_grid(p_a, p_b, nrow = 2, align = "v", axis = "l", labels = 'auto',
               rel_heights = c(1, 1.2))
ggsave("Figure_S7.pdf", p, width = 10, height = 7)
