#!/biodata/dep_psl/grp_rgo/guan/bin//Rscript

library(parallelDist)
library(vegan)
library(mina)

source("../../01.common_scripts/plot_setting.R")
source("../../01.common_scripts/plot_CAS_setting.R")

asv_file <- "../../00.data/CAS_rep_ASV_reNorm.rds"
design_file <- "../../00.data/design_854.txt"
tax_file <- "../../00.data/Bac_ASV_tax.txt"
nc_file <- "../../00.data/CAS_spearman_ap_tab.txt"

asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")
tax <- read.table(tax_file, header = T, sep = "\t")
nc <- read.table(nc_file, header = T, sep = "\t")

des_soil <- design[design$Compartment == "soil", ]
des_soil$Group <- des_soil$Soil

des_nsoil <- design[design$Compartment != "soil", ]
des_nsoil$Group <- paste0(des_nsoil$Host.Species, "_", des_nsoil$Compartment)
des <- rbind(des_soil, des_nsoil)

nc$Cluster <- paste0("Cluster_", nc$Cluster)

######### for family ##########################################################
bs_pm_dir <- "./distinctive_family_bs_pm/"
egv_dir <- "./distinctive_family_egv/"

if (!dir.exists(bs_pm_dir)) dir.create(bs_pm_dir)
if (!dir.exists(egv_dir)) dir.create(egv_dir)

dis_mean_file <- "../Figure_S12/Family_single_pm_dis_mean.txt"
dis_mean <- read.table(dis_mean_file, header = T, sep = "\t")
dis_mean <- dis_mean[order(dis_mean$Decrease_Mean, decreasing = T), ]
family <- as.character(dis_mean$Family)

family <- family[1:26]

print(family)

this_tax <- tax[tax$Family %in% family, ]
this_asv <- asv[rownames(asv) %in% rownames(this_tax), ]
print(dim(this_asv))

dat <- new("mina", norm = this_asv, des = des)

dat <- bs_pm(dat, group = "Group", g_size = 40, s_size = 20,
            rm = FALSE, sig = TRUE, bs = 33, pm = 33,
            individual = TRUE, out_dir = bs_pm_dir)

dat <- net_dis_indi(bs_pm_dir, method = "spectra", egv = TRUE, dir = egv_dir)

######### for cluster ##########################################################
bs_pm_dir <- "./distinctive_netcls_bs_pm/"
egv_dir <- "./distinctive_netcls_egv/"

if (!dir.exists(bs_pm_dir)) dir.create(bs_pm_dir)
if (!dir.exists(egv_dir)) dir.create(egv_dir)

dis_mean_file <- "../Figure_S12/Cluster_single_pm_dis_mean.txt"
dis_mean <- read.table(dis_mean_file, header = T, sep = "\t")
dis_mean <- dis_mean[order(dis_mean$Decrease_Mean, decreasing = T), ]
cluster <- as.character(dis_mean$Cluster)

cluster <- cluster[1:26]
print(cluster)

this_nc <- nc[nc$Cluster %in% cluster, ]
this_asv <- asv[rownames(asv) %in% this_nc$ID, ]
print(dim(this_asv))

dat <- new("mina", norm = this_asv, des = des)

dat <- bs_pm(dat, group = "Group", g_size = 40, s_size = 20,
            rm = FALSE, sig = TRUE, bs = 33, pm = 33,
            individual = TRUE, out_dir = bs_pm_dir)

dat <- net_dis_indi(bs_pm_dir, method = "spectra", egv = TRUE, dir = egv_dir)
