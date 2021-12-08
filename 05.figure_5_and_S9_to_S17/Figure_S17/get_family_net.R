#!/biodata/dep_psl/grp_rgo/guan/bin//Rscript

source("../../01.common_scripts/network_function.R")

asv_file <- "../../00.data/CAS_rep_ASV_reNorm.rds"
design_file <- "../../00.data/design_854.txt"
tax_file <- "../../00.data/Bac_ASV_tax.txt"

asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")
tax <- read.table(tax_file, header = T, sep = "\t")
tax <- tax[rownames(tax) %in% rownames(asv), ]

des_soil <- design[design$Compartment == "soil", ]
des_soil$Group <- des_soil$Soil

des_nsoil <- design[design$Compartment != "soil", ]
des_nsoil$Group <- paste0(des_nsoil$Host.Species, "_", des_nsoil$Compartment)
des <- rbind(des_soil, des_nsoil)

dis_mean_file <- "../../15.figure_7.17/Family_single_pm_dis_mean.txt"
dis_mean <- read.table(dis_mean_file, header = T, sep = "\t")
dis_mean <- dis_mean[order(dis_mean$Decrease_Mean, decreasing = T), ]
family <- as.character(dis_mean$Family)

family_at <- family[1:16]

### get node and edge file
get_net_feature <- function(this_asv, group, prefix) {
    this_des <- des[des$Group == group, ]
    this_asv <- this_asv[, colnames(this_asv) %in% this_des$Sample_ID]

    ## re-normalization
    this_asv <- apply(this_asv, 2, function(x) x/sum(x))

    this_sp <- rcorr(t(this_asv), type = "spearman")
    net <- flattenCorrMatrix(this_sp$r, this_sp$P)
    net[is.na(net)] <- 0
 #   net <- net[net$edge != 0, ]
    net <- net[net$p < 0.05, ]

    node <- get_node_features(net)
    tax$ASV_ID <- rownames(tax)
    node <- merge(node, tax)

    ## add the RA and Occupancy information
    ra <- data.frame(ASV_ID = rownames(this_asv),
                     RA = rowSums(this_asv) / ncol(this_asv),
                     Occupancy = rowSums(this_asv>0) / ncol(this_asv))

    node <- merge(node, ra)

    colnames(net)[1 : 3] <- c("source", "target", "weight")
    colnames(node)[1] <- "id"
    edge_file <- paste0("Family_", prefix, "_", group, "_edge.txt")
    node_file <- paste0("Family_", prefix, "_", group, "_node.txt")
    write.table(net, edge_file, quote = F, sep = "\t", row.names = F)
    write.table(node, node_file, quote = F, sep = "\t", row.names = F)
}

#### for CAS vs At root
this_tax <- tax[tax$Family %in% family_at, ]
this_asv <- asv[rownames(asv) %in% rownames(this_tax), ]
get_net_feature(this_asv, "CAS", "CAS_vs_At")
get_net_feature(this_asv, "arabidopsis_thaliana_root", "CAS_vs_At")

