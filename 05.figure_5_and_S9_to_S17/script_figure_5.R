#!/biodata/dep_psl/grp_rgo/guan/bin//Rscript

library(parallelDist)
library(vegan)
library(cowplot)
library(dplyr)

source("../01.common_scripts/plot_setting.R")
source("../01.common_scripts/plot_CAS_setting.R")

egv_folder <- "./CAS_spectra_egv/"

bs_files <- list.files(egv_folder, pattern = "^spectra_bs_", full.names = TRUE)
pm_files <- list.files(egv_folder, pattern = "^spectra_pm_", full.names = TRUE)

## for bootstrap results and plot
bs_vectors <- c()
for (bs in bs_files) {
    this <- readRDS(bs)
    bs_vectors <- rbind(bs_vectors, this)
}

dis_bs <- as.matrix(parDist(bs_vectors, method = "euclidean"))

rownames(dis_bs) <- colnames(dis_bs) <- paste0(rownames(dis_bs), "_", seq(1 : nrow(dis_bs)))

## get the group information
idx <- seq(1,length(colnames(dis_bs)) * 2, by = 2)
design <- data.frame(Net = colnames(dis_bs),
    Group = unlist(strsplit(colnames(dis_bs), "_bs"))[idx])

ad <- adonis(formula = dis_bs ~ Group,
             data = design, by = "margin", add = F, parallel = 40)
print(ad)

dis_bs_dmr <- cmdscale(dis_bs, k = 4, eig = T)
eig <- dis_bs_dmr$eig
eig[eig < 0] <- 0
eig1 <- 100 * eig[1] / sum(eig)
eig2 <- 100 * eig[2] / sum(eig)

points <- as.data.frame(dis_bs_dmr$points[, 1:2])
colnames(points) <- c("x", "y")
points$Net <- rownames(points)
points <- merge(points, design)

p_a <- ggplot(points, aes(x, y, color = Group, shape = Group)) +
    geom_point(alpha = 0.6) +
    main_theme +
    scale_color_manual(values = c_group, labels = tags) +
    scale_shape_manual(values = s_group, labels = tags) +
    coord_fixed(ratio = 1) +
    theme(legend.position = "top") +
    labs (x = paste0("PC1 (", format(eig1, digits = 4), "%)"),
          y = paste0("PC2 (", format(eig2, digits = 4), "%)"))

plot_tile <- function(x){
    dis <- read.table(x, header = T, sep = "\t")

	dis$Sig[dis$P < 0.05] <- "sig"
	dis$Sig[dis$P >=0.05] <- "non-sig"

	c_sig <- c("sig" = "gold",
	           "non-sig" = "gray")

	## set d to 2 for jaccard and d = 0 for spectra distance
	max_dis <- max(dis$Distance_Mean)
	if (max_dis <= 1) {
	    d <- 2
	} else {
	    d <- 0
	}

	p <- ggplot(dis, aes(Group1, Group2)) +
	    geom_tile(aes(fill = Distance_Mean, color = Sig), size = 1.2) +
	    geom_text(aes(label = round(Distance_Mean, d), color = Sig), size = 4) +
	    scale_color_manual(values = c_sig) +
	    main_theme +
	    theme(axis.text.x = element_text(colour = "black", angle = 90,
	                                     size = 4, hjust = 1),
	          axis.text.y = element_text(size = 4),
	          aspect.ratio=1)
}

dis_b <- "Figure_5b.txt"
p_b <- plot_tile(dis_b)

p_ab <- plot_grid(p_a, p_b, nrow = 1, labels = c('a', 'b'), rel_widths = c(1, 1.1))
#ggsave("test.pdf", p_ab, width = 10, height = 4)

######### for panel c ##########################################################
## run the script ./script_figure_S12.R and ./get_grp_size.R before plotting the
## panel c and figure S13. Also run ./get_group_aRA_occupancy.R to get the input
## file aRA_occupancy_network_clusters.txt and aRA_occupancy_families.txt

dis_file <- "./Figure_S12a.txt"
size_file <- "./Figure_5c.txt"

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

p_c <- plot_grid(p1, p2, p3, p4, nrow = 1, align = "h", axis = "b",
                 labels = c('c', '', '', ''))

## width = 10, height = 2.8
######## for panel d ##################################################
dis_mean_file <- "./Figure_S12/Cluster_single_pm_dis_mean.txt"
dis_mean <- read.table(dis_mean_file, header = T, sep = "\t")
dis_mean <- dis_mean[order(dis_mean$Decrease_Mean, decreasing = T), ]
family <- as.character(dis_mean$Cluster)

x_order <- c("All_BS", family, "All_PM")

dis_folder <- "./Figure_S14/Network_cluster/spectra_dis/"
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
p_d <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
write.table(DIS, "Figure_5d.txt", quote = F, sep = "\t", row.names = F)

### height = 3.5, width = 10
p_all <- plot_grid(p_ab, p_c, p_d, ncol = 1, labels = c('', '', 'd'),
                   rel_heights = c(4, 2.5, 3))

ggsave("Figure_5.pdf", p_all, width = 10, height = 9.5)
