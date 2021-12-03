#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(dplyr)

########################### Before plotting
## Run get_pm_dis.R to get the partially permutated networks and then run the
## get_dis_and_mean.R [this script] to get the mean distance change.
################################

## all inter-condition comparison
cmp_file <- "cmp.lst"
bs_pm_folder <- "../CAS_spectra_egv/"

######### for family ##########################################################
egv_folder <- "./Family/single_grp_pm_egv/"
out_dir <- "./Family/spectra_dis/"
if (!dir.exists(out_dir)) dir.create(out_dir)

cmp_lst <- read.table(cmp_file)$V1
family <- list.dirs(path = egv_folder, full.names = FALSE, recursive = FALSE)
len_f <- length(family)

for (cmp in cmp_lst){
    this_bs_dis <- paste0(bs_pm_folder, "/dis_spectra_bs_", cmp, ".txt")
    this_bs_dis <- read.table(this_bs_dis, header = T, sep = "\t")
  
    this_cmp_test <- c()

    for (i in 1:len_f){
        this_dis <- paste0(egv_folder, "/", family[i], "/dis_spectra_pm_",
                           cmp, ".txt")
        this_dis <- read.table(this_dis, header = T, sep = "\t")

        ## test if two distance are significantly different
        this_test <- t.test(this_bs_dis$Distance,
                            this_dis$Distance, alternative = "greater")

#        if (this_test$p.value > 0.05) next
        m1 <- mean(this_bs_dis$Distance)
        m2 <- mean(this_dis$Distance)
        m12 <- m1 - m2

        this_cmp_test <- rbind(this_cmp_test,
                               c(family[i], this_test$p.value, m1, m2, m12))
    }

    colnames(this_cmp_test) <-c("Family", "P_value",
                            "BS_Mean", "PM_fam_Mean", "Decrease")
    write.table(this_cmp_test, paste0(out_dir, "/", cmp, ".txt"),
            quote = F, sep = "\t", row.names = F)
}

dis_folder <- out_dir
dis_files <- list.files(dis_folder, pattern = ".txt$", full.names = TRUE)

all_dis <- c()

for (d in dis_files) {
    this_cmp <- unlist(strsplit(d, "/"))[length(unlist(strsplit(d, "/")))]
    this_cmp <- unlist(strsplit(this_cmp, ".txt"))[1]

    g1 <- unlist(strsplit(this_cmp, "_vs_"))[1]
    g2 <- unlist(strsplit(this_cmp, "_vs_"))[2]

    if (g1 == g2) next

    this_dis <- read.table(d, header = T, sep = "\t")
    this_dis <- this_dis[, c("Family", "Decrease", "P_value")]

    this_dis$Comparison <- this_cmp

    all_dis <- rbind(this_dis, all_dis)
}

write.table(all_dis, "Family_single_pm_dis.txt",
            quote = F, sep = "\t", row.names = F)

all_dis <- all_dis[, c("Family", "Decrease")]

all_dis_mean <- all_dis %>% group_by(Family) %>%
    summarise(Decrease_Mean = mean(Decrease))

all_dis_mean <- as.data.frame(all_dis_mean)

write.table(all_dis_mean, "Family_single_pm_dis_mean.txt",
            quote = F, sep = "\t", row.names = F)

######### for network cluster #################################################
egv_folder <- "./Network_cluster/single_grp_pm_egv/"
out_dir <- "./Network_cluster/spectra_dis/"
if (!dir.exists(out_dir)) dir.create(out_dir)

cluster <- list.dirs(path = egv_folder, full.names = FALSE, recursive = FALSE)
len_f <- length(cluster)

for (cmp in cmp_lst){

    this_bs_dis <- paste0(bs_pm_folder, "/dis_spectra_bs_", cmp, ".txt")
    this_bs_dis <- read.table(this_bs_dis, header = T, sep = "\t")
    m1 <- mean(this_bs_dis$Distance)

    this_cmp_test <- c()

    for (i in 1:len_f){
        this_dis <- paste0(egv_folder, "/", cluster[i], "/dis_spectra_pm_",
                           cmp, ".txt")
        this_dis <- read.table(this_dis, header = T, sep = "\t")

        this_test <- t.test(this_bs_dis$Distance,
                            this_dis$Distance, alternative = "greater")

        m2 <- mean(this_dis$Distance)
        m12 <- m1 - m2

        this_cmp_test <- rbind(this_cmp_test,
                               c(cluster[i], this_test$p.value, m1, m2, m12))
    }

    colnames(this_cmp_test) <-c("Cluster", "P_value",
                            "BS_Mean", "PM_fam_Mean", "Decrease")
    write.table(this_cmp_test, paste0(out_dir, "/", cmp, ".txt"),
            quote = F, sep = "\t", row.names = F)
}

dis_folder <- out_dir
dis_files <- list.files(dis_folder, pattern = ".txt$", full.names = TRUE)

all_dis <- c()

for (d in dis_files) {
    this_cmp <- unlist(strsplit(d, "/"))[length(unlist(strsplit(d, "/")))]
    this_cmp <- unlist(strsplit(this_cmp, ".txt"))[1]

    g1 <- unlist(strsplit(this_cmp, "_vs_"))[1]
    g2 <- unlist(strsplit(this_cmp, "_vs_"))[2]

    if (g1 == g2) next

    this_dis <- read.table(d, header = T, sep = "\t")
    this_dis <- this_dis[, c("Cluster", "Decrease", "P_value")]

    this_dis$Comparison <- this_cmp

    all_dis <- rbind(this_dis, all_dis)
}

write.table(all_dis, "Cluster_single_pm_dis.txt",
            quote = F, sep = "\t", row.names = F)

all_dis <- all_dis[, c("Cluster", "Decrease")]

all_dis_mean <- all_dis %>% group_by(Cluster) %>%
    summarise(Decrease_Mean = mean(Decrease))

all_dis_mean <- as.data.frame(all_dis_mean)

write.table(all_dis_mean, "Cluster_single_pm_dis_mean.txt",
            quote = F, sep = "\t", row.names = F)
