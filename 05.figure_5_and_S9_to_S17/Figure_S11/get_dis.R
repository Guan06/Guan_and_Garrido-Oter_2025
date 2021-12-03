#!/biodata/dep_psl/grp_rgo/guan/bin//Rscript

library(RSpectra)
library(reshape2)
library(stringr)
source("../../01.common_scripts/plot_setting.R")

get_both_dis <- function(x, evk = 100) {
    bs1_files <- sort(list.files(x, pattern = "_bs1.rds",
                                 full.names = TRUE))
    bs2_files <- sort(list.files(x, pattern = "_bs2.rds",
                                 full.names = TRUE))
    pm1_files <- sort(list.files(x, pattern = "_pm1.rds",
                                 full.names = TRUE))
    pm2_files <- sort(list.files(x, pattern = "_pm2.rds",
                                full.names = TRUE))
    len <- length(bs1_files)
    log <- c()

    for (i in 1:len) {
        bs1 <- readRDS(bs1_files[i])
        bs2 <- readRDS(bs2_files[i])
        group_mn <- strsplit(basename(bs1_files[i]), "_bs1.rds")[[1]][1]
        group_m <- strsplit(group_mn, "_vs_")[[1]][1]
        group_n <- strsplit(group_mn, "_vs_")[[1]][2]

        y_bs <- list()
        y_bs[1] <- bs1
        y_bs[2] <- bs2

        this_m <- y_bs[[1]]
        this_n <- y_bs[[2]]

        ## calculate bootstrap distance
        bs_len <- length(this_m)

        ## for Spectra distance
        spectra_m <- spectra_n <- matrix(nrow = bs_len, ncol = evk)
        flag <- 0
        for (j in 1 : bs_len) {
            adj_m <- unlist(this_m[[j]])
            adj_n <- unlist(this_n[[j]])
            adj_m[is.na(adj_m)] <- 0
            adj_n[is.na(adj_n)] <- 0

            log <- rbind(log, paste0(group_mn, " bs_", j, ": ", nrow(adj_m)))
            #skip this comparison if the adj_m and adj_n is smaller than evk
            if (nrow(adj_m) < evk ||nrow(adj_n) < evk) {
                  flag <- 1
                  break
              }
            spectra_m[j, ] <- get_spectra(adj_m, k = evk)
            spectra_n[j, ] <- get_spectra(adj_n, k = evk)
        }
        if (flag) next
        seqs <- seq(1 : bs_len)
        rownames(spectra_m) <- paste0(group_m, "_bs", seqs)
        rownames(spectra_n) <- paste0(group_n, "_bs", seqs)
        spectra_mn <- rbind(spectra_m, spectra_n)
        this_dis_bs <- get_dis_df(dist(spectra_mn))
        # filter intra group network comparison when comparing
        # networks from different environments
        if (group_m != group_n) {
            r <- this_dis_bs$Group1 != this_dis_bs$Group2
            this_dis_bs <- this_dis_bs[r, ]
        }
        write.table(this_dis_bs, paste0("Spectra_dis_bs_", group_mn, ".txt"),
                    quote = F, sep = "\t", row.names = F)

        ## for Jaccard distance
        jaccard_mn <- c()
        for (j1 in 1 : bs_len) {
            adj_m <- unlist(this_m[[j1]])
            adj_m[is.na(adj_m)] <- 0

            log <- rbind(log,
                         paste0(group_mn, " bs_", j1, ": ", nrow(adj_m)))

            m_j1 <- paste0(group_m, "_b", j1)

            for (j2 in 1 : bs_len) {
                adj_n <- unlist(this_n[[j2]])
                adj_n[is.na(adj_n)] <- 0

                n_j2 <- paste0(group_n, "_b", j2)

                contrast <- sum(abs(adj_m - adj_n))
                max <- sum(pmax(abs(adj_m), abs(adj_n)))

                dis <- contrast / max
                this <- data.frame(C1 = m_j1,
                                   C2 = n_j2,
                                   Distance = dis,
                                   Group1 = group_m,
                                   Group2 = group_n)
                jaccard_mn <- rbind(jaccard_mn, this)
            }
        }
        dis_bs_file <- paste0("Jaccard_dis_bs_", group_mn, ".txt")
        write.table(jaccard_mn, dis_bs_file, quote = F, sep = "\t", row.names = F)

        ## the permutation results
        pm1 <- readRDS(pm1_files[i])
        pm2 <- readRDS(pm2_files[i])

        y_pm <- list()
        y_pm[1] <- pm1
        y_pm[2] <- pm2

        this_mp <- y_pm[[1]]
        this_np <- y_pm[[2]]

        pm_len <- length(this_mp)

        ## for Spectra distance

        spectra_mp <- spectra_np <- matrix(nrow = pm_len, ncol = evk)

        flag <- 0
        for (k in 1 : pm_len) {
            adj_mp <- unlist(this_mp[[k]])
            adj_np <- unlist(this_np[[k]])
            adj_mp[is.na(adj_mp)] <- 0
            adj_np[is.na(adj_np)] <- 0

            if (nrow(adj_mp) < evk ||nrow(adj_np) < evk) {
                flag <- 1
                break
            }

            spectra_mp[k, ] <- get_spectra(adj_mp, k = evk)
            spectra_np[k, ] <- get_spectra(adj_np, k = evk)
        }

        if (flag) next
        seqs <- seq(1 : pm_len)
        rownames(spectra_mp) <- paste0(group_m, "_pm", seqs)
        rownames(spectra_np) <- paste0(group_n, "_pm", seqs)
        spectra_mnp <- rbind(spectra_mp, spectra_np)

        this_dis_pm <- get_dis_df(dist(spectra_mnp))

        # filter intra group network comparison when
        # comparing networks from different environments
        if (group_m != group_n) {
            r <- this_dis_pm$Group1 != this_dis_pm$Group2
            this_dis_pm <- this_dis_pm[r, ]
        }
        write.table(this_dis_pm, paste0("Sepctra_dis_pm_", group_mn, ".txt"),
                    quote = F, sep = "\t", row.names = F)

        ## for Jaccard distance
         jaccard_mnp <- c()
        for (k1 in 1 : pm_len) {
            adj_mp <- unlist(this_mp[[k1]])
            adj_mp[is.na(adj_mp)] <- 0

            log <- rbind(log, paste0(group_mn, " pm_", k1,
                              ": ", nrow(adj_mp)))

            mp_k1 <- paste0(group_m, "_p", k1)

            for (k2 in 1 : pm_len) {
                adj_np <- unlist(this_np[[k2]])
                adj_np[is.na(adj_np)] <- 0

                np_k2 <- paste0(group_n, "_p", k2)

                contrast <- sum(abs(adj_mp - adj_np))
                max <- sum(pmax(abs(adj_mp), abs(adj_np)))

                dis <- contrast / max
                this <- data.frame(C1 = mp_k1,
                                   C2 = np_k2,
                                   Distance = dis,
                                   Group1 = group_m,
                                   Group2 = group_n)

                jaccard_mnp <- rbind(jaccard_mnp, this)
            }
        }
        dis_pm_file <- paste0("Jaccard_dis_pm_", group_mn, ".txt")
        write.table(jaccard_mnp, dis_pm_file, quote = F, sep = "\t", row.names = F)
    }
}

## functions needed
get_spectra <- function(x,  k = 100){
    x <- as.matrix(x)
    x[is.na(x)] <- 0
    spectra <- eigs_sym(x, k, opts = list(retvec = FALSE))
    y <- spectra$values
}

get_dis_df <- function(x) {
    x <- as.matrix(x)
    x[upper.tri(x)] <- NA
    diag(x) <- NA
    x <- melt(x)

    colnames(x) <- c("C1", "C2", "Distance")
    x <- x[!is.na(x$Distance), ]

    if (str_detect(x$C1, "_bs")) {
        x <- cbind(x, do.call("rbind",
                              strsplit(as.character(x$C1), "_bs")))
        x <- x[, 1 : 4]

        x <- cbind(x, do.call("rbind",
                              strsplit(as.character(x$C2), "_bs")))
        x <- x[, 1 : 5]
    } else if (str_detect(x$C1, "_pm")) {
        x <- cbind(x, do.call("rbind",
                              strsplit(as.character(x$C1), "_pm")))
        x <- x[, 1 : 4]

        x <- cbind(x, do.call("rbind",
                              strsplit(as.character(x$C2), "_pm")))
        x <- x[, 1 : 5]
    }
    colnames(x)[4:5] <- c("Group1", "Group2")
    return(x)
}

## change this to the path with the folder where bs and pm files (rds) were
## stored
get_both_dis("./CAS_vs_Root/")
