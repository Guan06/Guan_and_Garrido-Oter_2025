#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

rarefaction_subsample <- function(x, depth = 1000, replace = TRUE){
    # initial rarefied vector
    rare <- numeric(length(x))
    if (replace) {
        suppressWarnings(subsample <- sample(1:length(x), depth, replace = TRUE,
                                             prob = x))
    } else {
        df <- data.frame(sample = (1:length(x)), time = as.numeric(x))
        obs <- apply(df, 1, function(x) {
                         rep_len(x["sample"], x["time"])
        })
        obs <- unlist(obs, use.names = FALSE)
        suppressWarnings(subsample <- sample(obs, depth, replace = FALSE))
    }

    subsample_tab <- table(subsample)
    rare[as(names(subsample_tab), "integer")] <- subsample_tab
    return(rare)
}

norm_by_raref <- function(x, depth = 1000, replace = TRUE) {

    # Remove samples contain few reads than `depth`
    if (min(colSums(x)) < depth) {
        rmsamples <- colnames(x)[colSums(x) < depth]
        message(length(rmsamples), " samples removed for low depth")
        x <- x[, colSums(x) >= depth]
    }

    x_norm <- apply(x, 2, rarefaction_subsample,
                    depth = depth, replace = replace)

    rownames(x_norm) <- rownames(x)
    x_norm <- x_norm[rowSums(x_norm) > 0, ]
    return(x_norm)
}
