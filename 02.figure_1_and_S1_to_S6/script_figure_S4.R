#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

source("../01.common_scripts/plot_setting.R")

## define the plotting function
plot_occu <- function(g){
    this_des <- des[des$Group == g, ]
    this_asv <- asv[, colnames(asv) %in% this_des$Sample_ID]
    this_asv <- this_asv[rowSums(this_asv) > 0, colSums(this_asv) > 0]

    this_fill <- unique(this_des$Compartment)
    occu <- data.frame(Taxonomy = rownames(this_asv),
                       Occu = rowSums(this_asv > 0),
                       Compartment = this_fill,
                       Rank = "ASV",
                       Tax_num = nrow(this_asv))

    for (t in tax_levels){
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
                            Compartment = this_fill,
                            Rank = t,
                            Tax_num = nrow(tab_TAX))
        occu <- rbind(occu, this_tax_occu)
    }

    this_p <- ggplot(occu, aes(x = Occu, y = ..scaled.., color = Rank)) +
        geom_density(alpha = 0.8) +
        main_theme +
        scale_color_manual(values = c_tax) +
        theme(legend.position = "none",
              axis.text = element_text(size = 6),
              axis.title = element_text(size = 6),
              plot.title = element_text(size = 6)) +
        ggtitle(g)
}

## for figure panel A
asv_file <- "../00.data/Bac_ASV_raref.rds"
design_file <- "../00.data/Bac_design_3809.txt"
tax_file <- "../00.data/Bac_ASV_tax.txt"

asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")
tax <- read.table(tax_file, header = T, sep = "\t")

des_soil <- design[design$Compartment == "soil", ]
des_soil$Group <- "soil"

host_lst <- c("zea_mays_l", "lotus_japonicus",
              "arabidopsis_thaliana")
des_nsoil <- design[design$Host.Species %in% host_lst, ]
des_nsoil <- des_nsoil[des_nsoil$Compartment != "soil", ]
des_nsoil$Group <- paste0(des_nsoil$Host.Species, "_", des_nsoil$Compartment)

des <- rbind(des_soil, des_nsoil)
tax_levels <- c("Phylum", "Class", "Order", "Family", "Genus")

p0 <- plot_occu("soil")

p1a <- plot_occu("zea_mays_l_rhizosphere")
p1b <- plot_occu("zea_mays_l_root")

p2a <- plot_occu("lotus_japonicus_rhizosphere")
p2b <- plot_occu("lotus_japonicus_root")

p4a <- plot_occu("arabidopsis_thaliana_rhizosphere")
p4b <- plot_occu("arabidopsis_thaliana_rhizoplane")
p4c <- plot_occu("arabidopsis_thaliana_root")

p_bac <- plot_grid(p0, p4b, p4a, p4c,
                   p1a, p1b, p2a, p2b,
                   ncol = 2, align = "v", axis = "l")

## for panel B
asv_file <- "../00.data/Fun_ASV_raref.rds"
design_file <- "../00.data/Fun_design_2232.txt"
tax_file <- "../00.data/Fun_ASV_tax.txt"

asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")
tax <- read.table(tax_file, header = T, sep = "\t")

des_soil <- design[design$Compartment == "soil", ]
des_soil$Group <- "soil"

des_nsoil <- design[design$Host.Species %in% host_lst, ]
des_nsoil <- des_nsoil[des_nsoil$Compartment != "soil", ]
des_nsoil$Group <- paste0(des_nsoil$Host.Species, "_", des_nsoil$Compartment)

des <- rbind(des_soil, des_nsoil)

p0 <- plot_occu("soil")
p1a <- plot_occu("zea_mays_l_rhizosphere")
p1b <- plot_occu("zea_mays_l_root")

p2a <- plot_occu("lotus_japonicus_rhizosphere")
p2b <- plot_occu("lotus_japonicus_root")

p4a <- plot_occu("arabidopsis_thaliana_rhizosphere")
p4b <- plot_occu("arabidopsis_thaliana_rhizoplane")
p4c <- plot_occu("arabidopsis_thaliana_root")

p_fun <- plot_grid(p0, p4b, p4a, p4c,
                   p1a, p1b, p2a, p2b,
                   ncol = 2, align = "v", axis = "l")

p_aligned <- plot_grid(p_bac, p_fun, ncol = 2, labels = 'auto')
ggsave("Figure_S4.pdf", height = 6, width = 14)
