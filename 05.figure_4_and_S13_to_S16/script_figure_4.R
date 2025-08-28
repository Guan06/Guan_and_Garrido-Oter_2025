### #!/biodata/dep_psl/grp_rgo/guan/bin//Rscript

library(vegan)
library(cowplot)
library(dplyr)

source("../01.common_scripts/plot_setting.R")
source("../01.common_scripts/plot_CAS_setting.R")

######### Figure 4 panel a 
################################################################################
## run the script ./script_figure_S7.R and ./get_grp_size.R before plotting the
## panel a and figure S8. Also run ./get_group_aRA_occupancy.R to get the input
## file aRA_occupancy_network_clusters.txt and aRA_occupancy_families.txt

dis_file <- "./Figure_S7a.txt"
size_file <- "./Figure_4a.txt"

ra_occu_file <- "./aRA_occupancy_network_clusters.txt"

#### function to plot the scatters from data frame
plot_mean <- function(x, cmp){
    dat <- merge(x, decrease)
    t <- cor.test(dat[[cmp]], dat$Decrease, method = "spearman", exact = FALSE)
    edge <- t$estimate
    sig <- t$p.value

    edge <- round(edge, 2)
#    sig <- round(sig, 4)
    print(sig)

    p <- ggplot(dat, aes_string(cmp, "Decrease")) +
        scale_x_log10() +
        geom_smooth(method = "lm", se = FALSE,
                    color = "#017F97") +
        geom_point(color = "#7570B3", alpha = 0.8) + main_theme +
        labs(title = paste0("r = ", edge))

}

## plot mean of all groups
dis <- read.table(dis_file, header = T, sep = "\t")
dis_mean <- dis[dis$Comparison == "Mean", ]
colnames(dis_mean)[1] <- "Group"
decrease <- dis_mean[, c("Group", "Decrease")]

size <- read.table(size_file, header = T, sep = "\t")
size_mean <- as.data.frame(size %>% group_by(Group) %>%
                           summarise(Size = mean(Number_of_ASVs)))

size_mean$Group <- paste0("Cluster_", size_mean$Group)
p1 <- plot_mean(size_mean, "Size") +
    labs(x = "Mean number of ASVs")

ra_occu <- read.table(ra_occu_file, header = T, sep = "\t")
colnames(ra_occu)[1] <- "Group"

ra <- ra_occu[, c("Group", "RA", "Compartment")]
ra_mean <- as.data.frame(ra %>% group_by(Group) %>%
                        summarise(Mean_aRA = mean(RA)))

ra_mean$Group <- paste0("Cluster_", ra_mean$Group)
p2 <- plot_mean(ra_mean, "Mean_aRA") +
    labs(x = "Mean aRA")

### get the normalized aRA for each family under each condition
size2 <- size[, c("Group", "Compartment", "Number_of_ASVs")]

ra_norm <- ra %>% inner_join(size2, by = c("Group", "Compartment"))
ra_norm$aRA_norm <- ra_norm$RA / ra_norm$Number_of_ASVs

ra_norm <- as.data.frame(ra_norm %>% group_by(Group) %>%
                         summarise(Mean_aRA_norm = mean(aRA_norm)))

ra_norm$Group <- paste0("Cluster_", ra_norm$Group)
p3 <- plot_mean(ra_norm, "Mean_aRA_norm") +
    labs(x = "Mean normalized aRA")

occu <- ra_occu[, c("Group", "Occupancy", "Compartment")]
occu_mean <- as.data.frame(occu %>% group_by(Group) %>%
                           summarise(Mean_Occupancy = mean(Occupancy)))

occu_mean$Group <- paste0("Cluster_", occu_mean$Group)
p4 <- plot_mean(occu_mean, "Mean_Occupancy") +
    labs(x = "Mean occupancy")

p_a <- plot_grid(p1, p2, p3, p4, nrow = 1, align = "h", axis = "b",
                 labels = c('a', '', '', ''))

## width = 10, height = 2.8
################################################################################

######## Figure 4 panel b
################################################################################
dis_mean_file <- "./Figure_S7/Cluster_single_pm_dis_mean.txt"
dis_mean <- read.table(dis_mean_file, header = T, sep = "\t")
dis_mean <- dis_mean[order(dis_mean$Decrease_Mean, decreasing = T), ]
family <- as.character(dis_mean$Cluster)

x_order <- c("All_BS", family, "All_PM")

dis_folder <- "./Figure_S9/Network_cluster/spectra_dis/"
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
p_b <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
write.table(DIS, "Figure_4b.txt", quote = F, sep = "\t", row.names = F)
################################################################################

p <- plot_grid(p_a, p_b, ncol = 1, labels = c('', 'b'),
                   rel_heights = c(2.5, 3))

ggsave("Figure_4.pdf", p, width = 10, height = 5.5)
