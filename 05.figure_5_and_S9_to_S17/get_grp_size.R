#!/biodata/dep_psl/grp_rgo/guan/bin//Rscript

library(dplyr)
rep_file <- "../00.data/Bac_ASV_rep.rds"
design_file <- "../00.data/design_854.txt"
rep <- readRDS(rep_file)
design <- read.table(design_file, header = T, sep = "\t")
des_soil <- design[design$Compartment == "soil", ]
des_soil$Group <- "soil"

des_nsoil <- design[design$Compartment != "soil", ]
des_nsoil$Group <- paste0(des_nsoil$Host.Species, "_", des_nsoil$Compartment)

des <- rbind(des_soil, des_nsoil)
g_lst <- unique(des$Group)

### for family features
grp_file <- "../00.data/CAS_ASV_1774_tax.txt"
group <- read.table(grp_file, header = T, sep = "\t")
group$Group <- group$Family

df_size <- c()
for (g in g_lst) {
    this_des <- des[des$Group == g, ]
    this_asv <- rep[rownames(rep) %in% rownames(group),
                    colnames(rep) %in% this_des$Sample_ID]

    this_asv <- this_asv[rowSums(this_asv) > 0, colSums(this_asv) > 0]
    this_grp <- group[rownames(group) %in% rownames(this_asv), ]

    this_size <- as.data.frame(this_grp %>% group_by(Group) %>% tally)
    this_size$Compartment <- g

    df_size <- rbind(df_size, this_size)
}
colnames(df_size) <- c("Group", "Number_of_ASVs", "Compartment")
write.table(df_size, "Figure_S13.txt", quote = F, sep = "\t", row.names = F)

### for network clusters
group_file <- "../00.data/CAS_spearman_ap_tab.txt"
group <- read.table(group_file, header = T, sep = "\t")
group$Group <- group$Cluster <- as.character(group$Cluster)

df_size <- c()
for (g in g_lst) {
    this_des <- des[des$Group == g, ]
    this_asv <- rep[rownames(rep) %in% group$ID,
                    colnames(rep) %in% this_des$Sample_ID]

    this_asv <- this_asv[rowSums(this_asv) > 0, colSums(this_asv) > 0]
    this_grp <- group[group$ID %in% rownames(this_asv), ]

    this_size <- as.data.frame(this_grp %>% group_by(Group) %>% tally)
    this_size$Compartment <- g

    df_size <- rbind(df_size, this_size)
}

colnames(df_size) <- c("Group", "Number_of_ASVs", "Compartment")
write.table(df_size, "Figure_5c.txt", quote = F, sep = "\t", row.names = F)
