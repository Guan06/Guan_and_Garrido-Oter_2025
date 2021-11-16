pcoa <- function(x, des, dim, color, shape, size = 1.2) {
    eig <- x$eig
    eig[eig < 0] <- 0

    if (dim == "12") {
        points <- x$points[, 1 : 2]
        eig1 <- 100 * eig[1] / sum(eig)
        eig2 <- 100 * eig[2] / sum(eig)
    } else if (dim == "34") {
        points <- x$points[, 3 : 4]
        eig3 <- 100 * eig[3] / sum(eig)
        eig4 <- 100 * eig[4] / sum(eig)
    } else if (dim == "13") {
        points <- x$points[, c(1, 3)]
        eig1 <- 100 * eig[1] / sum(eig)
        eig3 <- 100 * eig[3] / sum(eig)
    }

    points <- as.data.frame(points)
    colnames(points) <- c("x", "y")

    des <- des[match(rownames(points), des[["Sample_ID"]]), ]

    points <- cbind(points, des[[color]])
    colnames(points)[3] <- color
    points <- cbind(points, des[[shape]])
    colnames(points)[4] <- shape

    color_df <- get_color_df(color)
    shape_df <- get_shape_df(shape)

    p <- ggplot(points,
                aes_string(x = "x", y = "y",
                           color = color, shape = shape)) +
         geom_point(size = size, alpha = 0.7) +
         coord_fixed(ratio = 1) +
         scale_colour_manual(values = color_df) +
         scale_shape_manual(values = shape_df) +
         main_theme

     if (dim == "12") {
        p <- p + labs (x = paste0("PCo1 (", format(eig1, digits = 4), "%)"),
                       y = paste0("PCo2 (", format(eig2, digits = 4), "%)"))
     } else if (dim == "34") {
        p <- p + labs (x = paste0("PCo3 (", format(eig3, digits = 4), "%)"),
                       y = paste0("PCo4 (", format(eig4, digits = 4), "%)"))
     } else if (dim == "13") {
         p <- p + labs (x = paste0("PCo1 (", format(eig1, digits = 4), "%)"),
                        y = paste0("PCo3 (", format(eig3, digits = 4), "%)"))
     }
     return(p)
}

