### #!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(cowplot)
source("../01.common_scripts/plot_setting.R")

## for panel A --> bacteria
ra_file <- "./fixed_Occu/Bac_RA_m2.txt"
occu_file <- "./fixed_RA/Bac_Occu_m2.txt"

ra <- read.table(ra_file, header = T, sep = "\t")
occu <- read.table(occu_file, header = T, sep = "\t")

colnames(ra)[1] <- colnames(occu)[1] <- "Percent"

ra$Factor <- "Relative abundance"
occu$Factor <- "Occupancy"

dat <- rbind(ra, occu)

p_a <- ggplot(dat, aes(Percent, M2, color = Factor)) +
    geom_line() + geom_point() +
    scale_color_manual(values = c("#FFD662FF", "#00539CFF")) +
    main_theme +
    theme(legend.position = "top") +
    geom_segment(x = 6, y = 0.005, xend = 6, yend = 0.035,
                 color = "#00539CFF", linetype = "dashed") +
    geom_segment(x = 6, y = 0.015, xend = 6, yend = 0.045,
                 color = "#FFD662FF", linetype = "dashed")

## for panel B --> fungi
ra_file <- "./fixed_Occu/Fun_RA_m2.txt"
occu_file <- "./fixed_RA/Fun_Occu_m2.txt"

ra <- read.table(ra_file, header = T, sep = "\t")
occu <- read.table(occu_file, header = T, sep = "\t")

colnames(ra)[1] <- colnames(occu)[1] <- "Percent"

ra$Factor <- "RA"
occu$Factor <- "Occupancy"

dat <- rbind(ra, occu)

p_b <- ggplot(dat, aes(Percent, M2, color = Factor)) +
    geom_line() + geom_point() +
    scale_color_manual(values = c("#FFD662FF", "#00539CFF")) +
    main_theme +
    theme(legend.position = "none") +
    geom_segment(x = 6, xend = 6, y = 0, yend = 0.03,
                 color = "#00539CFF", linetype = "dashed") +
    geom_segment(x = 7, xend = 7, y = 0.025, yend = 0.055,
                 color = "#FFD662FF", linetype = "dashed")

## organize plots
p_legend <- get_legend(p_a)
p_a <- p_a + theme(legend.position = "none")
p <- plot_grid(p_a, p_b, nrow = 1, labels = 'auto',
               align = "h", axis = "b")

p_merge <- plot_grid(p_legend, p, ncol = 1,
                     rel_heights = c(0.1, 1))

ggsave("Figure_S5.pdf", p_merge, height = 2.5)
