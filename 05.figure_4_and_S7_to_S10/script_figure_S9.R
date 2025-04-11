### #!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(cowplot)
library(ggridges)
source("../01.common_scripts/plot_setting.R")

#### for family ###############################################################
dis_mean_file <- "./Figure_S7/Family_single_pm_dis_mean.txt"

dis_mean <- read.table(dis_mean_file, header = T, sep = "\t")
dis_mean <- dis_mean[order(dis_mean$Decrease_Mean, decreasing = T), ]
family <- as.character(dis_mean$Family)
x_order <- c("All_BS", family, "All_PM")

dis_folder <- "./Figure_S9/Family/spectra_dis/"

dis_files <- list.files(path = dis_folder, pattern = "_mean.txt",
                        full.names = TRUE, recursive = FALSE)

DIS <- c()

for (dis_file in dis_files) {
    dis <- read.table(dis_file, header = T)
    dis <- dis[match(x_order, dis$Family), ]

    ## get All_BS distance and remove this row
    dis_bs <- dis$Mean[1]
    dis_pm <- dis$Mean[length(dis$Mean)]

    dis <- dis[-1, ]

    dis$Decrease <- dis_bs - dis$Mean
    dis$Decrease_ratio <- dis$Decrease / dis_bs

    this_cmp <- unlist(strsplit(dis_file, "accu_dis_"))[2]
    this_cmp <- unlist(strsplit(this_cmp, "_mean.txt"))[1]

    dis$Comparison <- this_cmp

    DIS <- rbind(DIS, dis)
}

DIS$Color <- DIS$Comparison

DIS$Color[!(DIS$Comparison %in% c("arabidopsis_thaliana_root_vs_CAS",
                             "CAS_vs_lotus_japonicus_root"))] <- "Others"
x_order <- x_order[-1]

p1 <- ggplot(DIS, aes(x = Family, y = Decrease_ratio)) +
    geom_boxplot(color = "steelblue3", fill = NA, outlier.shape = NA) +
    scale_x_discrete(limits = x_order) +
    scale_color_manual(values = c_cmp) +
    main_theme +
    labs(x = "", y = "Distance decrease ratio") +
    theme(plot.background = element_blank(),
          legend.position = "none",
              plot.title = element_text(size = 8),
            axis.text.x = element_text(colour = "black", angle = 90,
            size = 6, hjust = 1))

DIS$Sig <- ifelse(DIS$P_value_W < 0.05, "sig", "non-sig")
p2 <- ggplot(DIS, aes(x = Family, y = P_value_F)) +
    geom_line(aes(group = Comparison, color = Color), alpha = 0.8, size = 0.3) +
    geom_point(aes(color = Color, shape = Sig), alpha = 0.8, size = 1) +
    scale_shape_manual(values = c(1, 16)) +
    scale_x_discrete(limits = x_order) +
    scale_color_manual(values = c_cmp) +
    main_theme +
    labs(x = "") +
    scale_y_continuous(position = "right", limits = c(0, 1)) +
    labs(y = "P value") +
    geom_hline(yintercept = 0.05, color = "gold3", linetype = "dashed") +
    geom_hline(yintercept = 0.5, color = "gray", alpha = 0.7) +
    geom_hline(yintercept = 0.4, color = "gray", alpha = 0.7,
                linetype = "dashed") +
    geom_hline(yintercept = 0.6, color = "gray", alpha = 0.7,
                linetype = "dashed") +
    theme(plot.background = element_blank(),
          legend.position = "none",
              plot.title = element_text(size = 8),
            axis.text.x = element_text(colour = "black", angle = 90,
            size = 6, hjust = 1))

aligned_plots <- align_plots(p1, p2, align = "hv", axis = "tblr")
p <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

write.table(DIS, "Figure_S9.txt", quote = F, sep = "\t", row.names = F)
ggsave("Figure_S9.pdf", p, width = 10, height = 3.5)
