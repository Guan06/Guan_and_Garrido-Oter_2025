#!/biodata/dep_psl/grp_rgo/guan/bin//Rscript

library(parallelDist)
library(vegan)

source("../01.common_scripts/plot_setting.R")
source("../01.common_scripts/plot_CAS_setting.R")

egv_folder <- "./CAS_spectra_egv/"

pm_files <- list.files(egv_folder, pattern = "^spectra_pm_", full.names = TRUE)
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

p <- ggplot(points, aes(x, y, color = Group, shape = Group)) +
    geom_point(alpha = 0.6) +
    main_theme +
    scale_color_manual(values = c_group, labels = tags) +
    scale_shape_manual(values = s_group, labels = tags) +
    coord_fixed(ratio = 1) +
    theme(legend.position = "top") +
    labs (x = paste0("PC1 (", format(eig1, digits = 4), "%)"),
          y = paste0("PC2 (", format(eig2, digits = 4), "%)"))
ggsave("Figure_S9.pdf", p, height = 5, width = 5)
