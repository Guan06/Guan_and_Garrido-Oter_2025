#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(ggplot2)
library(cowplot)

source("../01.common_scripts/plot_setting.R")
c_density <- c("Original dataset" = "#4682B4",
            "Permutated dataset" = "#B47846")

plot_density <- function(bs_file, pm_file){
    dis_bs <- read.table(bs_file, header = T, sep = "\t")
    dis_pm <- read.table(pm_file, header = T, sep = "\t")

    dis_bs$Group <- paste0(dis_bs$Group1, "_", dis_bs$Group2)
    dis_bs$Factor <- "Original dataset"

    dis_pm$Group <- paste0(dis_pm$Group1, "_", dis_pm$Group2)
    dis_pm$Factor <- "Permutated dataset"

    dis <- rbind(dis_bs, dis_pm)

    p <- ggplot(dis, aes(Distance, fill = Factor, color = Factor)) +
        geom_density(aes(y= ..scaled..), alpha = 0.7) +
        main_theme +
        scale_fill_manual(values = c_density) +
        scale_color_manual(values = c_density) +
        theme(legend.position = "none")
}

p3 <- plot_density("./Figure_S11/Jaccard_dis_bs_arabidopsis_thaliana_root_vs_CAS.txt",
                  "./Figure_S11/Jaccard_dis_pm_arabidopsis_thaliana_root_vs_CAS.txt")
p4 <- plot_density("./Figure_S11/Jaccard_dis_bs_CAS_vs_lotus_japonicus_root.txt",
                   "./Figure_S11/Jaccard_dis_pm_CAS_vs_lotus_japonicus_root.txt")

p1 <- plot_density("./Figure_S11/Spectra_dis_bs_arabidopsis_thaliana_root_vs_CAS.txt",
                   "./Figure_S11/Sepctra_dis_pm_arabidopsis_thaliana_root_vs_CAS.txt")
p2 <- plot_density("./Figure_S11/Spectra_dis_bs_CAS_vs_lotus_japonicus_root.txt",
                   "./Figure_S11/Sepctra_dis_pm_CAS_vs_lotus_japonicus_root.txt")

p <- plot_grid(p1, p2, p3, p4, nrow = 2, labels = 'auto')
ggsave("Figure_S11.pdf", height = 3)
