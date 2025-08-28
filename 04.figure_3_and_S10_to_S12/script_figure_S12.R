### #!/biodata/dep_psl/grp_rgo/guan/bin//Rscript

library(parallelDist)
library(cowplot)
library(vegan)

source("../01.common_scripts/plot_setting.R")
source("../01.common_scripts/pcoa_function.R")

args = commandArgs(trailingOnly = TRUE)
rds_folder <- args[1]
plot_file <- args[2]

bs_files <- list.files(rds_folder, pattern = "spectra_bs_", full.names = TRUE)
pm_files <- list.files(rds_folder, pattern = "spectra_pm_", full.names = TRUE)
bs_vectors <- c()
pm_vectors <- c()

for (bs in bs_files) {
    this <- readRDS(bs)
    bs_vectors <- rbind(bs_vectors, this)
}

for (pm in pm_files) {
    cmp <- unlist(strsplit(pm, "spectra_pm_"))[2]
    cmp <- unlist(strsplit(cmp, ".rds"))[1]
    g1 <- unlist(strsplit(cmp, "_vs_"))[1]
    g2 <- unlist(strsplit(cmp, "_vs_"))[2]
    if (g1 == g2) next
    this <- readRDS(pm)
    pm_vectors <- rbind(pm_vectors, this)
}

dis_bs <- as.matrix(parDist(bs_vectors, method = "euclidean"))
dis_pm <- as.matrix(parDist(pm_vectors, method = "euclidean"))

### plot for original dataset
rc_n <- paste0(rownames(dis_bs), "_", seq(1 : nrow(dis_bs)))
rownames(dis_bs) <- colnames(dis_bs) <- rc_n

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

p1 <- ggplot(points, aes(x, y, color = Group)) +
    geom_point(alpha = 0.7) +
    scale_colour_manual(values = c_Com) +
    main_theme +
    theme(legend.position = "none") +
    coord_fixed(ratio = 1) +
    labs (x = paste0("PC1 (", format(eig1, digits = 4), "%)"),
          y = paste0("PC2 (", format(eig2, digits = 4), "%)"))

### plot for permutated dataset
rc_n <- paste0(rownames(dis_pm), "_", seq(1 : nrow(dis_pm)))
rownames(dis_pm) <- colnames(dis_pm) <- rc_n

## get the group information
idx <- seq(1,length(colnames(dis_pm)) * 2, by = 2)
design <- data.frame(Net = colnames(dis_pm),
    Group = unlist(strsplit(colnames(dis_pm), "_pm"))[idx])
design$Compartment <- design$Group

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

p2 <- ggplot(points, aes(x, y, color = Group)) +
    geom_point(alpha = 0.7) +
    scale_colour_manual(values = c_Com) +
    main_theme +
    theme(legend.position = "none") +
    coord_fixed(ratio = 1) +
    labs (x = paste0("PC1 (", format(eig1, digits = 4), "%)"),
          y = paste0("PC2 (", format(eig2, digits = 4), "%)"))

p <- plot_grid(p1, p2, ncol = 1)
ggsave(plot_file, p, width = 5, height = 7)
