library(ggridges)
library(tidyr)

source("../01.common_scripts/plot_setting.R")
source("../01.common_scripts/plot_CAS_setting.R")
source("../01.common_scripts/network_function.R")

asv_file <- "../00.data/CAS_rep_ASV_reNorm.rds"
des_file <- "../00.data/design_854.txt"

asv <- readRDS(asv_file)
des <- read.table(des_file, header = T, sep = "\t")

des <- des[des$Host.Species != "chlamydomonas_reinhardtii", ]
des <- des[des$Host.Species != "lotus_japonicus" | des$Compartment == "soil", ]
dim(des)

asv <- asv[rowSums(asv) > 0, colnames(asv) %in% des$Sample_ID]
dim(asv)

c_lst <- unique(des$Compartment)

all_node <- c()
for (c in c_lst) {
  this_des <- des[des$Compartment == c, ]
  this_asv <- asv[, colnames(asv) %in% this_des$Sample_ID]
  this_asv <- this_asv[rowSums(this_asv) > 0, ]
  print(dim(this_des))
  print(dim(this_asv))
  
  this_net <- rcorr(t(this_asv), type = "spearman")
  mat_cor <- this_net$r
  mat_p <- this_net$P
  
  vect_q <- p.adjust(mat_p, method = "fdr")
  mat_q <- matrix(vect_q, nrow = nrow(mat_p), ncol = ncol(mat_p))
  
  mat_cor[mat_q > 0.05] <- 0
  
  this_edge <- flattenCorrMatrix(mat_cor, mat_q)
  this_edge <- this_edge[this_edge$edge != 0, ]
  
  no_node <- length(unique(c(this_edge$node1, this_edge$node2)))
  print(paste(c, no_node, nrow(this_edge)))
  
  write.table(this_edge, paste0("edge_", c, ".txt"), 
              quote = F, sep = "\t", row.names = F)
  
  print(paste("Number of positive edges:", sum(this_edge$edge > 0)))
  print(paste("Number of negative edges:", sum(this_edge$edge < 0)))
  
  ## density
  this_den <- nrow(this_edge) / (no_node * (no_node - 1) / 2)
  print(paste("Density: ", this_den))
  
  this_node <- get_node_features(this_edge)
  this_node$Compartment <- c
  all_node <- rbind(all_node, this_node)
}

all_node_long <- all_node %>% pivot_longer(cols = Degree:Betweenness,
                              names_to = "Feature", 
                              values_to = "Value")
all_node_long <- all_node_long[all_node_long$Value != 0, ]
#all_node_long <- all_node_long[!is.na(all_node_long$Value), ]

ggplot(all_node_long, aes(log(Value), Compartment, fill = Compartment)) +
  geom_density_ridges(color = "gray47", quantile_lines = TRUE, 
                      quantiles = 2, scale = 1.8) + #,
                      # jittered_points = TRUE,
                      # position = position_points_jitter(width = 0.05, height = 0),
                      # point_shape = '|', point_size = 3, 
                      # point_alpha = 1, alpha = 0.7) + 
  theme_ridges() +
  coord_cartesian(clip = "off") +
  facet_grid(~Feature, scales = "free", space = "free_y") +
  scale_y_discrete(limits = c("root","rhizoplane", "rhizosphere", 
                              "soil")) +
  scale_fill_manual(values = c_Com) +
  theme(legend.position = "NA", 
        strip.background = element_blank(),
        strip.clip = "off") +
  labs(x = "Log-transformed feature value", y = "")

ggsave("network_node_feature_dis.pdf", width = 8, height = 4)

## Go to Cytoscape to generate the network plots.

tax <- read.table("../00.data/Bac_ASV_tax.txt", header = T, sep = "\t")
tax <- tax[rownames(tax) %in% rownames(asv), ]
tax$node <- rownames(tax)
write.table(tax, "CAS_node.txt", row.names = F, quote = F, sep = "\t")
