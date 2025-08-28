### #!/biodata/dep_psl/grp_rgo/guan/bin//Rscript

source("../01.common_scripts/plot_setting.R")

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
	                                     size = 6, hjust = 1),
	          axis.text.y = element_text(size = 6),
	          aspect.ratio=1)
}

dis_b <- "Figure_S11a.txt"
dis_c <- "Figure_S11b.txt"

p_b <- plot_tile(dis_b)
p_c <- plot_tile(dis_c)

library(cowplot)
p_bc <- plot_grid(p_b, p_c, nrow = 1)
ggsave("Figure_S11.pdf", width = 12, height = 6)
