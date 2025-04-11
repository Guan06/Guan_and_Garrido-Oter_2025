#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(RSpectra)
library(reshape2)
library(stringr)

get_pm_group <- function(asv, des, asv_lst, group = "Group",
                          g_size = 40, s_size = 20,
                          pm = 33, dir = "./"){

    if (s_size >= g_size) {
        stop("`s_size` can not be larger than `g_size`!")
    }

    mat <- asv
    des <- des

    mat <- mat[rowSums(mat) > 0, ]
    message(nrow(mat),
        " components are used for pm before filtering.")

    # fit the quantitative table with descriptive table
    fit <- intersect(colnames(mat), des$Sample_ID)
    mat <- mat[, colnames(mat) %in% fit]
    des <- des[des$Sample_ID %in% fit, ]
    message(nrow(des), " samples are used for pm before filtering.")

    # filter the group does not have enough samples (i.e. < g_size)
    des[[group]] <- as.factor(des[[group]])
    lst <- levels(des[[group]])
    len <- length(lst)

    for (i in 1 : len) {
        g <- lst[i]
        size <- sum(des[[group]] == g)
        if (size < g_size) lst[i] <- NA
    }

    lst <- lst[!is.na(lst)]
    len <- length(lst)

    des <- des[des[[group]] %in% lst, ]
    mat <- mat[, colnames(mat) %in% des$Sample_ID]
    message(len, " groups with ", ncol(mat), " samples used for permutation.")

    for (m in 1 : len) {
        group_m <- lst[m]
        des_m <- des[des[[group]] == group_m, ]
        mat_m <- mat[, colnames(mat) %in% des_m$Sample_ID]
        num_m <- nrow(des_m)
 
        for (n in m : len) {
            if (n == m) {
                group_n <- group_m
                this_mat_m <- this_mat_n <- mat_m
                mat_mn <- mat_n <- mat_m
                num_mn <- num_n <- num_m
 
            } else {
                group_n <- lst[n]
                des_n <- des[des[[group]] == group_n, ]
                mat_n <- mat[, colnames(mat) %in% des_n$Sample_ID]
                num_n <- nrow(des_n)
                mat_mn <- cbind(mat_m, mat_n)
                num_mn <- num_m + num_n
                this_mat_m <- mat_m
                this_mat_n <- mat_n
            }
            # start permutation
            MPLST <- list()
            NPLST <- list()

            for (p in 1 : pm) {
                bs_m <- sample.int(num_m, s_size)
                bs_n <- sample.int(num_n, s_size)
                
                mat_bs_m <- this_mat_m[, bs_m]
                mat_bs_n <- this_mat_n[, bs_n]

                pm_mn <- sample.int(num_mn, s_size * 2)
                pm_m <- pm_mn[1 : s_size]
                pm_n <- pm_mn[(s_size + 1) : (s_size * 2)]
                mat_pm_m <- mat_mn[, pm_m]
                mat_pm_n <- mat_mn[, pm_n]

                index <- rownames(mat_mn) %in% asv_lst
                mat_bs_m[index, ] <- mat_pm_m[index, ]
                mat_bs_n[index, ] <- mat_pm_n[index, ]
    
                cor_pm <- adj(mat_bs_m, method = "spearman")
                cor_pn <- adj(mat_bs_n, method = "spearman")
    
                cor_pm[is.na(cor_pm)] <- 0
                cor_pn[is.na(cor_pn)] <- 0
    
                MPLST[[p]] <- cor_pm
                NPLST[[p]] <- cor_pn
                names(MPLST)[p] <- paste0(group_m, "_", p)
                names(NPLST)[p] <- paste0(group_n, "_", p)
            }
            prefix <- paste0(dir, group_m, "_vs_", group_n)
            saveRDS(list(MPLST), paste0(prefix, "_pm1.rds"))
            saveRDS(list(NPLST), paste0(prefix, "_pm2.rds"))
            rm(MPLST, NPLST)
            gc(reset = TRUE)
        }
    }
}

get_egv <- function(x, evk = 100, out_dir) {
    pm1_files <- sort(list.files(x, pattern = "_pm1.rds",
                                 full.names = TRUE))
    pm2_files <- sort(list.files(x, pattern = "_pm2.rds",
                                full.names = TRUE))
    len <- length(pm1_files)
    log <- c()

    for (i in 1:len) {
        group_mn <- strsplit(basename(pm1_files[i]), "_pm1.rds")[[1]][1]
        group_m <- strsplit(group_mn, "_vs_")[[1]][1]
        group_n <- strsplit(group_mn, "_vs_")[[1]][2]

        ## the permutation results
        pm1 <- readRDS(pm1_files[i])
        pm2 <- readRDS(pm2_files[i])

        y_pm <- list()
        y_pm[1] <- pm1
        y_pm[2] <- pm2

        this_mp <- y_pm[[1]]
        this_np <- y_pm[[2]]

        pm_len <- length(this_mp)

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

        saveRDS(spectra_mnp,
                file = paste0(out_dir, "/spectra_pm_", group_m, "_vs_",
                        group_n, ".rds"))

        this_dis_pm <- get_dis_df(dist(spectra_mnp))

        # filter intra group network comparison when
        # comparing networks from different environments
        if (group_m != group_n) {
            r <- this_dis_pm$Group1 != this_dis_pm$Group2
            this_dis_pm <- this_dis_pm[r, ]
        }
        write.table(this_dis_pm, paste0(out_dir, "/dis_spectra_pm_", group_m,
                                        "_vs_", group_n, ".txt"),
                    quote = F, sep = "\t", row.names = F)
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

    if (str_detect(x$C1[1], "_bs")) {
        x <- cbind(x, do.call("rbind",
                              strsplit(as.character(x$C1), "_bs")))
        x <- x[, 1 : 4]

        x <- cbind(x, do.call("rbind",
                              strsplit(as.character(x$C2), "_bs")))
        x <- x[, 1 : 5]
    } else if (str_detect(x$C1[1], "_pm")) {
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
