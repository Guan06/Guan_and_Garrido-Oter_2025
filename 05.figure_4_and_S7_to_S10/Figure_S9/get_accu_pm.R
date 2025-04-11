#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(mina)
source("../Figure_S12/function_get_pm.R")

asv_file <- "../../00.data/CAS_rep_ASV_reNorm.rds"
design_file <- "../../00.data/design_854.txt"

asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")

des_soil <- design[design$Compartment == "soil", ]
des_soil$Group <- des_soil$Soil

des_nsoil <- design[design$Compartment != "soil", ]
des_nsoil$Group <- paste0(des_nsoil$Host.Species, "_", des_nsoil$Compartment)
des <- rbind(des_soil, des_nsoil)

######### for family ##########################################################
grp_file <- "../Figure_S12/Family_single_pm_dis_mean.txt"
asv2grp_file <- "../../00.data/CAS_ASV_1774_tax.txt"

out_dir <- "./Family/"

grp <- read.table(grp_file, header = T, sep = "\t")
grp <- grp[order(grp$Decrease_Mean, decreasing = T), ]
group <- as.character(grp$Family)

asv2grp <- read.table(asv2grp_file, header = T, sep = "\t")

out_dir1 <- paste0(out_dir, "/accu_grp_pm/")
out_dir2 <- paste0(out_dir, "/accu_grp_pm_egv/")

if (!dir.exists(out_dir)) dir.create(out_dir)
if (!dir.exists(out_dir1)) dir.create(out_dir1)
if (!dir.exists(out_dir2)) dir.create(out_dir2)

len_f <- length(group)

for (i in 1:len_f) {
    f_name <- group[i]
    grp_lst <- group[1:i]

    this_asv <- asv2grp[asv2grp$Family %in% grp_lst, ]

    asv_lst <- rownames(this_asv)

    this_out_dir <- paste0(out_dir1, f_name, "/")
    if (!dir.exists(this_out_dir)) dir.create(this_out_dir)

    this_pm <- get_pm_group(asv, des, asv_lst = asv_lst,
                            group = "Group", g_size = 40, s_size = 20,
                            pm = 33,
                            dir = this_out_dir)

    this_out_dir_egv <- paste0(out_dir2, f_name, "/")
    if (!dir.exists(this_out_dir_egv)) dir.create(this_out_dir_egv)

    get_egv(this_out_dir, out_dir = this_out_dir_egv)
}

######### for network cluster #################################################
grp_file <- "../Figure_S12/Cluster_single_pm_dis_mean.txt"
asv2grp_file <- "../../00.data/CAS_spearman_ap_tab.txt"
out_dir <- "./Network_cluster/"

grp <- read.table(grp_file, header = T, sep = "\t")
grp <- grp[order(grp$Decrease_Mean, decreasing = T), ]
group <- as.character(grp$Cluster)

asv2grp <- read.table(asv2grp_file, header = T, sep = "\t")
asv2grp$Cluster <- paste0("Cluster_", asv2grp$Cluster)

out_dir1 <- paste0(out_dir, "/accu_grp_pm/")
out_dir2 <- paste0(out_dir, "/accu_grp_pm_egv/")

if (!dir.exists(out_dir)) dir.create(out_dir)
if (!dir.exists(out_dir1)) dir.create(out_dir1)
if (!dir.exists(out_dir2)) dir.create(out_dir2)

len_f <- length(group)

for (i in 1:len_f) {
    f_name <- group[i]
    grp_lst <- group[1:i]

    this_asv <- asv2grp[asv2grp$Cluster %in% grp_lst, ]

    asv_lst <- this_asv$ID

    this_out_dir <- paste0(out_dir1, f_name, "/")
    if (!dir.exists(this_out_dir)) dir.create(this_out_dir)

    this_pm <- get_pm_group(asv, des, asv_lst = asv_lst,
                             group = "Group", g_size = 40, s_size = 20,
                            pm = 33,
                            dir = this_out_dir)

    this_out_dir_egv <- paste0(out_dir2, f_name, "/")
    if (!dir.exists(this_out_dir_egv)) dir.create(this_out_dir_egv)

    get_egv(this_out_dir, out_dir = this_out_dir_egv)
}
