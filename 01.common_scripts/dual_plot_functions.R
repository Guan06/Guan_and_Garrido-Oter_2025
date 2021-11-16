
plot_dual <- function(x) {
     p1 <- ggplot(x, aes(color = Compartment, fill = Compartment)) +
        geom_histogram(aes(x = Occu), binwidth = 1) +
        main_theme +
        scale_color_manual(values = c_Com) +
        scale_fill_manual(values = c_Com) +
        labs(x = "Occupancy", y = "Number of taxa") +
        theme(legend.position = "none",
              plot.background = element_blank(),
              axis.ticks = element_blank(),
              axis.line.x = element_blank(),
              axis.text = element_text(size = 6),
              axis.title = element_text(size = 6),
              plot.title = element_text(size = 6))

    p2 <- ggplot(x) +
        geom_line(aes(x = Occu, y = cumsum(RA)), color = "gray23", alpha = 0.8) +
        scale_y_continuous(limits = c(0,1), position = "right") +
        main_theme +
        theme_minimal_hgrid(color = "gray88", line_size = 0.3) +
        labs(y = "aRA of taxa") +
        theme(legend.position = "none",
              plot.background = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.line.x = element_line(color = "black"),
              axis.line.y = element_line(color = "black"),
              axis.text.y = element_text(size = 6),
              axis.title.y = element_text(size = 6),
              axis.title.x = element_blank())

    aligned_plots <- align_plots(p2, p1, align="hv", axis="tblr")
    ggdraw(aligned_plots[[2]]) + draw_plot(aligned_plots[[1]])
}

## define the plotting function
plot_occu_RA <- function(g, t = "Family"){
    this_des <- des[des$Group == g, ]
    this_asv <- asv[, colnames(asv) %in% this_des$Sample_ID]
    this_asv <- this_asv[rowSums(this_asv) > 0, colSums(this_asv) > 0]

    this_fill <- unique(this_des$Compartment)
    occu <- data.frame(Taxonomy = rownames(this_asv),
                       Occu = rowSums(this_asv > 0),
                       RA = rowSums(this_asv) / ncol(this_asv),
                       Compartment = this_fill,
                       Rank = "ASV",
                       Tax_num = nrow(this_asv))
    occu <- occu[order(occu$Occu),]

    ## plot at the ASV level
    p12 <- plot_dual(occu)


    this_tax <- tax[rownames(tax) %in% rownames(this_asv), ]
    this_tax[[t]][is.na(this_tax[[t]])] <- "Unassigned"

    TAX <- as.factor(this_tax[[t]])
    tab_TAX <- apply(this_asv, 2, function(x) rowsum(x, TAX))
    rownames(tab_TAX) <- levels(TAX)

    un_ratio <- sum(tab_TAX[rownames(tab_TAX) == "Unassigned", ]) /
            ncol(tab_TAX)
    print(paste(g, t, un_ratio))

    this_tax_occu <- data.frame(Taxonomy = rownames(tab_TAX),
                            Occu = rowSums(tab_TAX > 0),
                            RA = rowSums(tab_TAX) / ncol(tab_TAX),
                            Compartment = this_fill,
                            Rank = t,
                            Tax_num = nrow(tab_TAX))

    this_tax_occu <- this_tax_occu[order(this_tax_occu$Occu),]

    ## filter unassigned
    this_tax_occu <- this_tax_occu[this_tax_occu$Taxonomy != "Unassigned", ]
    ## plot at the Family level
    p34 <- plot_dual(this_tax_occu)

    p <- plot_grid(p12, p34, nrow = 1, align = "h", axis = "b")
}

