#!/biodata/dep_psl/grp_rgo/guan/bin//Rscript

## This script is to get the networks of permutation dataset, when all the
## non-distinctive groups are permutated, to test if the distinctive groups
## only can already recapitualte the netowrk distance.

library(parallelDist)
library(vegan)
library(mina)

source("../../01.common_scripts/plot_setting.R")
source("../../01.common_scripts/plot_CAS_setting.R")
source("../Figure_S12/function_get_pm.R")

asv_file <- "../../00.data/CAS_rep_ASV_reNorm.rds"
design_file <- "../../00.data/design_854.txt"
tax_file <- "../../00.data/Bac_ASV_tax.txt"
nc_file <- "../../00.data/CAS_spearman_ap_tab.txt"

asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")
tax <- read.table(tax_file, header = T, sep = "\t")
tax <- tax[rownames(tax) %in% rownames(asv), ]

nc <- read.table(nc_file, header = T, sep = "\t")

des_soil <- design[design$Compartment == "soil", ]
des_soil$Group <- des_soil$Soil

des_nsoil <- design[design$Compartment != "soil", ]
des_nsoil$Group <- paste0(des_nsoil$Host.Species, "_", des_nsoil$Compartment)
des <- rbind(des_soil, des_nsoil)

nc$Cluster <- paste0("Cluster_", nc$Cluster)

######### for family ##########################################################
bs_pm_dir <- "./family_bs_pm/"
egv_dir <- "./family_egv/"

if (!dir.exists(bs_pm_dir)) dir.create(bs_pm_dir)
if (!dir.exists(egv_dir)) dir.create(egv_dir)

dis_mean_file <- "../Figure_S12/Family_single_pm_dis_mean.txt"
dis_mean <- read.table(dis_mean_file, header = T, sep = "\t")
dis_mean <- dis_mean[order(dis_mean$Decrease_Mean, decreasing = T), ]
family <- as.character(dis_mean$Family)

### for CAS vs At root
#family <- family[1:26]
this_family <- family[17:length(family)]
print(this_family)

this_tax <- tax[tax$Family %in% this_family, ]
asv_lst <- rownames(this_tax)
print(length(asv_lst))

bs_pm_dir1 <- paste0(bs_pm_dir, "/CAS_vs_At_root/")
egv_dir1 <- paste0(egv_dir, "/CAS_vs_At_root/")
if (!dir.exists(bs_pm_dir1)) dir.create(bs_pm_dir1)
if (!dir.exists(egv_dir1)) dir.create(egv_dir1)

this_pm <- get_pm_group(asv, des, asv_lst = asv_lst,
                        group = "Group", g_size = 40, s_size = 20,
                        pm = 33, dir = bs_pm_dir1)

get_egv(bs_pm_dir1, out_dir = egv_dir1)

### for CAS vs Lj root
this_family <- family[68:length(family)]
print(this_family)

this_tax <- tax[tax$Family %in% this_family, ]
asv_lst <- rownames(this_tax)
print(length(asv_lst))

bs_pm_dir2 <- paste0(bs_pm_dir, "/CAS_vs_Lj_root/")
egv_dir2 <- paste0(egv_dir, "/CAS_vs_Lj_root/")

if (!dir.exists(bs_pm_dir2)) dir.create(bs_pm_dir2)
if (!dir.exists(egv_dir2)) dir.create(egv_dir2)

this_pm <- get_pm_group(asv, des, asv_lst = asv_lst,
                        group = "Group", g_size = 40, s_size = 20,
                        pm = 33, dir = bs_pm_dir2)

get_egv(bs_pm_dir2, out_dir = egv_dir2)

######### for cluster ##########################################################
bs_pm_dir <- "./netcls_bs_pm/"
egv_dir <- "./netcls_egv/"

if (!dir.exists(bs_pm_dir)) dir.create(bs_pm_dir)
if (!dir.exists(egv_dir)) dir.create(egv_dir)

dis_mean_file <- "../Figure_S12/Cluster_single_pm_dis_mean.txt"
dis_mean <- read.table(dis_mean_file, header = T, sep = "\t")
dis_mean <- dis_mean[order(dis_mean$Decrease_Mean, decreasing = T), ]
cluster <- as.character(dis_mean$Cluster)

### for CAS vs At root
this_cluster <- cluster[4:length(cluster)]
print(this_cluster)

this_nc <- nc[nc$Cluster %in% this_cluster, ]
asv_lst <- this_nc$ID
print(length(asv_lst))

bs_pm_dir1 <- paste0(bs_pm_dir, "/CAS_vs_At_root/")
egv_dir1 <- paste0(egv_dir, "/CAS_vs_At_root/")
if (!dir.exists(bs_pm_dir1)) dir.create(bs_pm_dir1)
if (!dir.exists(egv_dir1)) dir.create(egv_dir1)

this_pm <- get_pm_group(asv, des, asv_lst = asv_lst,
                        group = "Group", g_size = 40, s_size = 20,
                        pm = 33, dir = bs_pm_dir1)

get_egv(bs_pm_dir1, out_dir = egv_dir1)

### for CAS vs Lj root
this_cluster <- cluster[18:length(cluster)]
print(this_cluster)

this_nc <- nc[nc$Cluster %in% this_cluster, ]
asv_lst <- this_nc$ID
print(length(asv_lst))

bs_pm_dir2 <- paste0(bs_pm_dir, "/CAS_vs_Lj_root/")
egv_dir2 <- paste0(egv_dir, "/CAS_vs_Lj_root/")
if (!dir.exists(bs_pm_dir2)) dir.create(bs_pm_dir2)
if (!dir.exists(egv_dir2)) dir.create(egv_dir2)

this_pm <- get_pm_group(asv, des, asv_lst = asv_lst,
                        group = "Group", g_size = 40, s_size = 20,
                        pm = 33, dir = bs_pm_dir2)

get_egv(bs_pm_dir2, out_dir = egv_dir2)
