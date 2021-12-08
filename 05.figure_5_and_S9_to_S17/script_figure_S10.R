#!/biodata/dep_psl/grp_rgo/guan/bin//Rscript

library(parallelDist)
library(vegan)
library(cowplot)
library(dplyr)

source("../01.common_scripts/plot_setting.R")
source("../01.common_scripts/plot_CAS_setting.R")

egv_folder <- "./spectra_egv_all/"

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

## for permutation dataset, panel b
pm_vectors <- c()

for (pm in pm_files){
    cmp <- unlist(strsplit(pm, "spectra_pm_"))[2]
    cmp <- unlist(strsplit(cmp, ".rds"))[1]
    g1 <- unlist(strsplit(cmp, "_vs_"))[1]
    g2 <- unlist(strsplit(cmp, "_vs_"))[2]
    if (g1 == g2) next
    this <- readRDS(pm)
    pm_vectors <- rbind(pm_vectors, this)
}

dis_pm <- as.matrix(dist(pm_vectors))
rownames(dis_pm) <- colnames(dis_pm) <- paste0(rownames(dis_pm), "_",
                                               seq(1 : nrow(dis_pm)))
idx <- seq(1,length(colnames(dis_pm)) * 2, by = 2)

design <- data.frame(Net = colnames(dis_pm),
    Group = unlist(strsplit(colnames(dis_pm), "_pm"))[idx])

## error here because of phycosphere if only use "_p"!!!

ad <- adonis(formula = dis_pm ~ Group,
             data = design, by = "margin", add = F, parallel = 40)
print(ad)

dis_pm_dmr <- cmdscale(dis_pm, k = 4, eig = T)
eig <- dis_pm_dmr$eig
eig[eig < 0] <- 0
eig1 <- 100 * eig[1] / sum(eig)
eig2 <- 100 * eig[2] / sum(eig)

points <- as.data.frame(dis_pm_dmr$points[, 1:2])
colnames(points) <- c("x", "y")
points$Net <- rownames(points)
points <- merge(points, design)

p_b <- ggplot(points, aes(x, y, color = Group, shape = Group)) +
    geom_point(alpha = 0.6) +
    main_theme +
    scale_color_manual(values = c_group, labels = tags) +
    scale_shape_manual(values = s_group, labels = tags) +
    coord_fixed(ratio = 1) +
    theme(legend.position = "none") +
    labs (x = paste0("PC1 (", format(eig1, digits = 4), "%)"),
          y = paste0("PC2 (", format(eig2, digits = 4), "%)"))

p <- plot_grid(p_a, p_b, ncol = 1, labels = 'auto')
ggsave("Figure_S10.pdf", p, height = 7, width = 5)
