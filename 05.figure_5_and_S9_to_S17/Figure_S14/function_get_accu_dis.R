get_accu_dis <- function(x, egv_folder, accu_egv_folder, out_dir) {

    bs_dis_tab <- read.table(paste0(egv_folder, "/dis_spectra_bs_", x, ".txt"),
                             header = T, sep = "\t")
    pm_dis_tab <- read.table(paste0(egv_folder, "/dis_spectra_pm_", x, ".txt"),
                             header = T, sep = "\t")

    bs_dis_tab$Family <- "All_BS"
    pm_dis_tab$Family <- "All_PM"

    ## filter the intra comparison of All_BS and All_pm
    bs_dis_tab <- bs_dis_tab[bs_dis_tab$Group1 != bs_dis_tab$Group2, ]
    pm_dis_tab <- pm_dis_tab[pm_dis_tab$Group1 != pm_dis_tab$Group2, ]

    bs_dis_tab$Mean <- mean(bs_dis_tab$Distance)
    pm_dis_tab$Mean <- mean(pm_dis_tab$Distance)

    bs_dis_tab$Compare <- paste0(bs_dis_tab$Group1, "_vs_", bs_dis_tab$Group2)
    pm_dis_tab$Compare <- paste0(pm_dis_tab$Group1, "_vs_", pm_dis_tab$Group2)

    ### F-test for All_BS and All_PM
    cmp_bs <- bs_dis_tab[, c("Compare", "Distance")]
    cmp_pm <- pm_dis_tab[, c("Compare", "Distance")]
    colnames(cmp_pm)[2] <- "Distance_PM"

    tmp <- merge(cmp_bs, cmp_pm)
    all_bs_sig <- (sum(tmp$Distance_PM > tmp$Distance) + 1) / (nrow(tmp) + 1)

    cmp_bs <- pm_dis_tab[, c("Compare", "Distance")]
    tmp <- merge(cmp_bs, cmp_pm)
    all_pm_sig <- (sum(tmp$Distance_PM > tmp$Distance) + 1) / (nrow(tmp) + 1)
    ####

    this_cmp <- c()
    this_F_test <- c()
    this_W_test <- c()
    this_T_test <- c()
    
#    family <- list.dirs(accu_egv_folder, full.names = F, recursive = FALSE)

    len_f <- length(family)

    for (i in 1:len_f){
            this_dis <- paste0(accu_egv_folder, "/", family[i],
                               "/dis_spectra_pm_", x, ".txt")
            this_dis <- read.table(this_dis, header = T, sep = "\t")
            this_dis$Family <- family[i]
            this_dis$Mean <- mean(this_dis$Distance)

            ## Add F-test to compare this_dis with the bs_dis_tab
            this_dis$Compare <- paste0(this_dis$Group1, "_vs_", this_dis$Group2)

            cmp_bs <- this_dis[, c("Compare", "Distance")]

            tmp <- merge(cmp_bs, cmp_pm)
            sig <- (sum(tmp$Distance_PM > tmp$Distance) + 1) / (nrow(tmp) + 1)
            this_F_test <- rbind(this_F_test, c(family[i], sig))

            ## add other nonparametric tests
            sig2 <- wilcox.test(this_dis$Distance, pm_dis_tab$Distance,
                                alternative = "greater")$p.value
            this_W_test <- rbind(this_W_test, c(family[i], sig2))

            this_cmp <- rbind(this_cmp, this_dis)

            ## add t test
            sig3 <- t.test(this_dis$Distance, pm_dis_tab$Distance,
                           alternative = "greater")$p.value
            this_T_test <- rbind(this_T_test, c(family[i], sig3))

    }

    this_F_test <- rbind(this_F_test, c("All_BS", all_bs_sig))
    this_F_test <- rbind(this_F_test, c("All_PM", all_pm_sig))
    this_F_test <- as.data.frame(this_F_test)
    colnames(this_F_test) <- c("Family", "P_value_F")

    all_bs_sig2 <- wilcox.test(bs_dis_tab$Distance, pm_dis_tab$Distance,
                               alternative = "greater")$p.value
    all_pm_sig2 <- wilcox.test(pm_dis_tab$Distance, pm_dis_tab$Distance,
                               alternative = "greater")$p.value

    this_W_test <- rbind(this_W_test, c("All_BS", all_bs_sig2))
    this_W_test <- rbind(this_W_test, c("All_PM", all_pm_sig2))
    this_W_test <- as.data.frame(this_W_test)
    colnames(this_W_test) <- c("Family", "P_value_W")
    this_W_test$P_adjust_W <- p.adjust(this_W_test$P_value_W)

    all_bs_sig3 <- t.test(bs_dis_tab$Distance, pm_dis_tab$Distance,
                          alternative = "greater")$p.value
    all_pm_sig3 <- t.test(pm_dis_tab$Distance, pm_dis_tab$Distance,
                          alternative = "greater")$p.value

    this_T_test <- rbind(this_T_test, c("All_BS", all_bs_sig3))
    this_T_test <- rbind(this_T_test, c("All_PM", all_pm_sig3))
    this_T_test <- as.data.frame(this_T_test)
    colnames(this_T_test) <- c("Family", "P_value_T")
    this_T_test$P_adjust_T <- p.adjust(this_T_test$P_value_T)

    both_tests <- merge(this_F_test, this_W_test)
    all_tests <- merge(both_tests, this_T_test)

    class(all_tests$P_value_F) <- class(all_tests$P_adjust_W) <- "numeric"
    class(all_tests$P_value_T) <- "numeric"

    this_cmp <- rbind(bs_dis_tab, this_cmp)
    this_cmp <- rbind(this_cmp, pm_dis_tab)

    write.table(this_cmp, paste0(out_dir, "accu_dis_", x, ".txt"),
                quote = F, sep = "\t", row.names = F)

    this_mean <- unique(this_cmp[, c("Family", "Mean")])
    this_mean <- merge(this_mean, all_tests, all = T)

#   this_mean <- this_mean[match(x_order, this_mean$Family), ]

    write.table(this_mean, paste0(out_dir, "accu_dis_", x, "_mean.txt"),
                quote = F, sep = "\t", row.names = F)
}
