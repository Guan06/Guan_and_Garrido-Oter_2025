#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

source("function_get_accu_dis.R")
source("../Figure_S12/cmp.lst")

egv_folder <- "../CAS_spectra_egv/"

######### for family ##########################################################
accu_egv_folder <- "./Family/accu_grp_pm_egv/"
dis_mean_file <- "../Figure_S12/Family_single_pm_dis_mean.txt"

dis_mean <- read.table(dis_mean_file, header = T, sep = "\t")
dis_mean <- dis_mean[order(dis_mean$Decrease_Mean, decreasing = T), ]
family <- as.character(dis_mean$Family)

out_dir <- "./Family/spectra_dis/"
if (!dir.exists(out_dir)) dir.create(out_dir)

lapply(cmp_lst1, FUN=function(x) get_accu_dis(x, egv_folder, accu_egv_folder,
                                         out_dir))
lapply(cmp_lst2, FUN=function(x) get_accu_dis(x, egv_folder, accu_egv_folder,
                                         out_dir))
lapply(cmp_lst3, FUN=function(x) get_accu_dis(x, egv_folder, accu_egv_folder,
                                         out_dir))
lapply(cmp_lst4, FUN=function(x) get_accu_dis(x, egv_folder, accu_egv_folder,
                                         out_dir))

######### for network cluster ################################################
accu_egv_folder <- "./Network_cluster/accu_grp_pm_egv/"
dis_mean_file <- "../Figure_S12/Cluster_single_pm_dis_mean.txt"

dis_mean <- read.table(dis_mean_file, header = T, sep = "\t")
dis_mean <- dis_mean[order(dis_mean$Decrease_Mean, decreasing = T), ]
family <- as.character(dis_mean$Cluster)

out_dir <- "./Network_cluster/spectra_dis/"
if (!dir.exists(out_dir)) dir.create(out_dir)

lapply(cmp_lst1, FUN=function(x) get_accu_dis(x, egv_folder, accu_egv_folder,
                                         out_dir))
lapply(cmp_lst2, FUN=function(x) get_accu_dis(x, egv_folder, accu_egv_folder,
                                         out_dir))
lapply(cmp_lst3, FUN=function(x) get_accu_dis(x, egv_folder, accu_egv_folder,
                                         out_dir))
lapply(cmp_lst4, FUN=function(x) get_accu_dis(x, egv_folder, accu_egv_folder,
                                         out_dir))
