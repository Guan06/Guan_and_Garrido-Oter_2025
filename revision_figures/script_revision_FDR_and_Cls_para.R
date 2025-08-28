### #!/biodata/dep_psl/grp_rgo/guan/bin/Rscript
#load_all("/biodata/dep_psl/grp_rgo/guan/Rpackage_Mina/mina/")
#devtools::install_github("Guan06/mina", dependencies = TRUE)

## function def, re_format_AP
################################################################################
re_format_AP <- function(x) {
  exemplars <- names(ap_exemplars(x))
  clusterids <- which(names(ap_exemplars(x)) %in% exemplars)
  
  clusters <- ap_clusters(x)[clusterids]
  cls <- sapply(clusters, length)
  ulc <- unlist(clusters)
  y <- data.frame(ID = names(ulc),
                  Exemplar = factor(rep(exemplars, cls)),
                  Cluster = rep(clusterids, cls),
                  Cluster_size = rep(cls, cls))
  return(y)
}

ap_exemplars <- function(x) x@exemplars
ap_clusters <- function(x) x@clusters
################################################################################

if (!requireNamespace("mina", quietly = TRUE)) {
  BiocManager::install("mina", force = TRUE)
}

library(mina)
library(Hmisc)
library(apcluster)

source("../01.common_scripts/plot_setting.R")
source("../01.common_scripts/pcoa_function.R")

asv_file <- "../00.data/Bac_ASV_rep.rds"
design_file <- "../00.data/Bac_design_3809.txt"

## factors for R2 calculation
bio_factors <- c("Compartment", "Soil", "Soil.Batch",
                  "Host.Species", "Host.Genotype", "Conditions")
tech_factors <- c("Study", "Run")
###

asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")

net <- rcorr(t(asv), type = "spearman")
mat_cor <- net$r
mat_p <- net$P

diag(mat_p) <- 1
vect_q <- p.adjust(mat_p, method = "fdr")
mat_q <- matrix(vect_q, nrow = nrow(mat_p), ncol = ncol(mat_p))

mat_cor[mat_q > 0.05] <- 0
rm(mat_p)
rm(mat_q)

################################################################################
### 2025.06.26 Test multiple threshold and how they affect cluster number
## FDR
## mat_q <- p.adjust(mat_p, method = "fdr")
lst <- c(0.05, 0.01, 0.001, 0.0001)

for (t in lst) {
  print(t)
  mat_cor_cp1 <- mat_cor_cp2 <- mat_cor
  print(sum(mat_p < t))
  print(sum(mat_q < t))
  
  no_node <- nrow(mat_cor)
  denp <- sum(mat_p < t) / 2 / (no_node * (no_node - 1) / 2) 
  denq <- sum(mat_q < t) / 2 / (no_node * (no_node - 1) / 2) 
  print(denp)
  print(denq)
  
  mat_cor_cp1[mat_p > t] <- 0
  this_ap1 <- apcluster(mat_cor_cp1, p = 0)
  print(length(this_ap1))
  
  mat_cor_cp2[mat_q > t] <- 0
  this_ap2 <- apcluster(mat_cor_cp2, p = 0)
  print(length(this_ap2))
  
}
################################################################################

#cls_ap <- net_cls(mat_cor, method = "ap", cutoff = 0)

ap_default <- apcluster(mat_cor)
cls_ap_default <- re_format_AP(ap_default)

## test different parameters for p and q
all <- c()
for (p in seq(0, 1, 0.1)) {
  for (q in seq(0, 1, 0.1)){
    this_ap <- apcluster(mat_cor, p = p, q = q)
    
    this_cls <- re_format_AP(this_ap)
    
    print(paste(p, q, length(unique(this_cls$Exemplar))))
    
    this_cls$p <- p
    this_cls$q <- q
    
    all <- rbind(all, this_cls)
  }
}

saveRDS(all, "20250709_pq_test_AP.rds")

################################################################################
#### test a range of inflation values 
library(MCL)

# Run MCL with each inflation value
all <- c()

# Scale the matrix to a range of 1 to 2 as MCL only takes positive values
#scaled_mat_cor <- mat_cor + 1 
#diag(scaled_mat_cor) <- 0

mat_cor2 <- mat_cor
## use only positive edges
mat_cor2[mat_cor2 < 0] <- 0

exp_range <- seq(2, 5, 0.5)
inf_range <- seq(2, 5, 0.5)
for (expansion in exp_range) {
  for (inflation in inf_range) {
    # Run MCL with the current inflation value
    mcl_result <- mcl(mat_cor2, addLoops = T,
                      expansion = expansion,
                      inflation = inflation)
    
    # Evaluate the resulting clusters
    num_clusters <- length(unique(mcl_result$Cluster))
    
    print(paste(expansion, inflation, num_clusters))
  }
}

# Plot the results
ggplot(results, aes(x = inflation, y = num_clusters)) +
  geom_point() +
  geom_line() +
  labs(x = "Inflation", y = "Number of Clusters")


# ap <- apcluster(mat_cor, p = 0)
# cls_ap <- re_format_AP(ap)
# 
# cls_ap_tab <- get_net_cls_tab(asv, cls_ap)

################################################################################
### Plot cluster size pf AP
t <- unique(cls_ap[, 3:4])
ggplot(t, aes(Cluster_size)) + 
  geom_histogram(binwidth = 1) + 
  theme_bw() + 
  labs(y = "Number of cluster", x = "Cluster size")
ggsave("revision_cluster_size_AP.pdf", width = 4, height = 3)
################################################################################

dis_ap <- com_dis(cls_ap_tab, method = "bray")

print("R2 of AP cluster(bio_factors) for bacterial community: ")
print(get_r2(dis_ap, design, group = bio_factors))
print("R2 of AP cluster(bio and tech facotrs) for bacterial community: ")
print(get_r2(dis_ap, design, group = c(bio_factors, tech_factors)))

ap_file <- "Figure_S5a_spearman_ap.txt"
write.table(cls_ap, ap_file, quote = F, row.names = F, sep = "\t")

## diversity analysis based on network clusters
dmr_ap <- dmr(dis_ap)

p1 <- pcoa(dmr_ap, design, 12, "Compartment", "Host.Species")

## for Fungi
asv_file <- "../00.data/Fun_ASV_rep.rds"
design_file <- "../00.data/Fun_design_2232.txt"

asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")

net <- rcorr(t(asv), type = "spearman")
mat_cor <- net$r
mat_p <- net$P
mat_cor[mat_p > 0.05] <- 0

#cls_ap <- net_cls(mat_cor, method = "ap", cutoff = 0)
ap <- apcluster(mat_cor, p = 0)
cls_ap <- re_format_AP(ap)

cls_ap_tab <- get_net_cls_tab(asv, cls_ap)
dis_ap <- com_dis(cls_ap_tab, method = "bray")

print("R2 of AP cluster(bio_factors) for fungal community: ")
print(get_r2(dis_ap, design, group = bio_factors))
print("R2 of AP cluster(bio and tech facotrs) for fungal community: ")
print(get_r2(dis_ap, design, group = c(bio_factors, tech_factors)))

ap_file <- "Figure_S5b_spearman_ap.txt"
write.table(cls_ap, ap_file, quote = F, row.names = F, sep = "\t")

## diversity analysis based on network clusters
dmr_ap <- dmr(dis_ap)

p2 <- pcoa(dmr_ap, design, 12, "Compartment", "Host.Species")

## put panels together
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")

p <- plot_grid(p1, p2, nrow = 1)
ggsave("Figure_S5.pdf", p, width = 9)
