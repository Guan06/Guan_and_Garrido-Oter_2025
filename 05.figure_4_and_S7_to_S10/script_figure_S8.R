### #!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(cowplot)
library(dplyr)

source("../01.common_scripts/plot_setting.R")

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
        labs(title = paste0("r = ", edge),
             y = "Family permutation decrease")

}
######### for family ##########################################################
## If any file not exist, run the corresponding scripts before plotting for S8.
## script_figure_S7.R -> Figure_S7b.txt
## get_grp_size.R -> Figure_S8.txt

dis_file <- "Figure_S7b.txt"
size_file <- "Figure_S8.txt"

ra_occu_file <- "./aRA_occupancy_families.txt"

## plot mean of all groups
dis <- read.table(dis_file, header = T, sep = "\t")
dis_mean <- dis[dis$Comparison == "Mean", ]
colnames(dis_mean)[1] <- "Group"
decrease <- dis_mean[, c("Group", "Decrease")]
decrease <- decrease[decrease$Group != "Unassigned", ]

size <- read.table(size_file, header = T, sep = "\t")
size_mean <- as.data.frame(size %>% group_by(Group) %>%
                           summarise(Size = mean(Number_of_ASVs)))

p1 <- plot_mean(size_mean, "Size") +
    labs(x = "Mean number of ASVs")

ra_occu <- read.table(ra_occu_file, header = T, sep = "\t")
colnames(ra_occu)[1] <- "Group"

ra <- ra_occu[, c("Group", "RA", "Compartment")]
ra_mean <- as.data.frame(ra %>% group_by(Group) %>%
                        summarise(Mean_aRA = mean(RA)))

p2 <- plot_mean(ra_mean, "Mean_aRA") +
    labs(x = "Mean aRA")

### get the normalized aRA for each family under each condition
size2 <- size[, c("Group", "Compartment", "Number_of_ASVs")]

ra_norm <- ra %>% inner_join(size2, by = c("Group", "Compartment"))
ra_norm$aRA_norm <- ra_norm$RA / ra_norm$Number_of_ASVs

ra_norm <- as.data.frame(ra_norm %>% group_by(Group) %>%
                         summarise(Mean_aRA_norm = mean(aRA_norm)))
p3 <- plot_mean(ra_norm, "Mean_aRA_norm") +
    labs(x = "Mean normalized aRA")

occu <- ra_occu[, c("Group", "Occupancy", "Compartment")]
occu_mean <- as.data.frame(occu %>% group_by(Group) %>%
                           summarise(Mean_Occupancy = mean(Occupancy)))
p4 <- plot_mean(occu_mean, "Mean_Occupancy") +
    labs(x = "Mean occupancy")

p <- plot_grid(p1, p2, p3, p4, nrow = 1, align = "h", axis = "b",
                 labels = 'auto')

ggsave("Figure_S8.pdf", p, height = 2.7, width = 10)
