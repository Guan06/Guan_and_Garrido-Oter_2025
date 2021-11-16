box_sig <- function(x, group, compare) {
    lst <- as.character(unique(x[[group]]))
    len <- length(lst)
    sig <- c()
    for (i in 1 : (len - 1 )) {
        for (j in (i + 1) : len) {
            sig_i <- x[x[[group]] == lst[i], ]
            sig_j <- x[x[[group]] == lst[j], ]
            sig_ij <- wilcox.test(sig_i[[compare]], sig_j[[compare]],
                                  alternative = "two.sided")
            this_sig <- c(lst[i], lst[j], sig_ij$p.value)
            sig <- rbind(sig, this_sig)
        }
    }
    colnames(sig) <- c("Group1", "Group2", "Significance")
    rownames(sig) <- c()
    sig <- as.data.frame(sig)
    return(sig)
}
