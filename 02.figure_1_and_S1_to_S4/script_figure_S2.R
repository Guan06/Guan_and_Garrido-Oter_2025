#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

source("../01.common_scripts/plot_setting.R")
source("../01.common_scripts/order.R")

## for panel A
bac_file <- "../00.data/Bac_tax_sum.txt"
sum <- read.table(bac_file, header = T, sep = "\t")
cc <- length(unique(sum$Phylum))

p_a <- ggplot(sum, aes(x = Group, y = Sum, fill = Phylum)) +
    geom_bar(stat = 'identity', position = 'stack') +
    main_theme +
    scale_x_discrete(limits = plot_order, labels = bac_labels) +
    scale_fill_manual(values = getOI(cc)) +
    theme(axis.text.x = element_text(colour = "black", angle = 90,
                                     size = 8, hjust = 1)) +
    labs(x = "", y = "")

## for panel B
fun_file <- "../00.data/Fun_tax_sum.txt"
sum <- read.table(fun_file, header = T, sep = "\t")
cc <- length(unique(sum$Phylum))

p_b <- ggplot(sum, aes(x = Group, y = Sum, fill = Phylum)) +
    geom_bar(stat = 'identity', position = 'stack') +
    main_theme +
    scale_x_discrete(limits = fun_plot_order, labels = fun_labels) +
    scale_fill_manual(values = getDark2(cc)) +
    theme(axis.text.x = element_text(colour = "black", angle = 90,
                                     size = 8, hjust = 1)) +
    labs(x = "", y = "")

## organize the panels with cowplot
library(cowplot)

p_a_legend <- get_legend(p_a)
p_b_legend <- get_legend(p_b)
p_a <- p_a + theme(legend.position = "none")
p_b <- p_b + theme(legend.position = "none")
p_aligned <- plot_grid(p_a, p_b,
                       plot_grid(p_a_legend, p_b_legend, nrow = 1),
                       ncol = 1, align = "v", labels = c('a', 'b', ''))

ggsave("Figure_S2.pdf", p_aligned, width = 7, height = 10)
