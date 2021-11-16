
library(ggplot2)
library(RColorBrewer)
library(cowplot)

getDark2 = colorRampPalette(brewer.pal(8, "Dark2"))
getSet1 = colorRampPalette(brewer.pal(9, "Set1"))

BrBG <- brewer.pal(n = 11, name = "BrBG")

palette_OkabeIto <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                      "#0072B2", "#D55E00", "#CC79A7", "#999999")
getOI = colorRampPalette(palette_OkabeIto)

alpha <- .7

c_rhizos <- rgb(255 / 255, 130 / 255,   0 / 255, alpha)
c_root   <- rgb( 50 / 255, 150 / 255, 100 / 255, alpha)
c_soil   <- rgb(101 / 255,  67 / 255,  33 / 255, alpha)
c_phyco  <- rgb(139 / 255, 139 / 255,   0 / 255, alpha)
c_rhizop <- "lightgreen"

c_Com <- c("rhizoplane" = c_rhizop,
            "rhizosphere" = c_rhizos, "root" = c_root, "soil" = c_soil,
            "phycosphere" = c_phyco)

c_phyla <- c("Proteobacteria" = "#E69F00",
            "Actinobacteriota" = "#BAA546",
            "Chloroflexi" = "#8EAB8D",
            "Bacteroidota" = "#62B2D4",
            "Myxococcota" = "#43AFCF",
            "Acidobacteriota" = "#29A8AB",
            "Firmicutes" = "#0EA187",
            "Nitrospirota" = "#1FA76C",
            "Gemmatimonadota" = "#68BC5D",
            "Bdellovibrionota" = "#B1D14E",
            "RCP2-54" = "#E5DF46",
            "Planctomycetota" = "#9CBC68",
            "Fibrobacterota" = "#53998B",
            "NB1-j" = "#0A76AD",
            "Patescibacteria" = "#376C83",
            "Armatimonadota" = "#78664D",
            "Verrucomicrobiota" = "#B96017",
            "Spirochaetota" = "#D3621D",
            "Desulfobacterota" ="#D06A4F",
            "Crenarchaeota" = "#CD7382",
            "Unassigned" = "#999999",
            "Deinococcota" = "#C77BA5",
            "Latescibacterota" = "#B885A1",
            "Entotheonellaeota" = "#A88F9D")

c_tax <- c("Phylum" = BrBG[11], "Class" = BrBG[9],
            "Order" = BrBG[7], "Family" = BrBG[5],
            "Genus" = BrBG[3], "ASV" = BrBG[1])

## colors from Dark2
my_gray <- rgb(220/255, 220/255, 220/255, alpha = 0.5)

c_cmp <- c("arabidopsis_thaliana_root_vs_CAS" = "#1B9E77",
           "CAS_vs_lotus_japonicus_root" = "#D95F02",
           "Mean" = "#7570B3",
           "Others" = my_gray)

s_Host <- c("arabidopsis_thaliana" = 16, "neighbor" = 1,
            "capsella_rubella" = 5, "cardamine" = 3, "wild_grass" = 4,
            "zea_mays_l" = 2, "lotus_japonicus" = 9, "host_soil" = 8,
            "chlamydomonas_reinhardtii" = 10,
            "medicago_truncatula" = 11)

s_Soil <- c("DEMO" = 2, "DOK" = 17, "CAS" = 16,
            "Spain" = 5, "Germany" = 1, "France" = 8, "Sweden" = 9, "Italy" = 10,
            "HAS" = 4)
s_Study <- c("Reconstruct_pilot_2016" = 3, "Reconstruct_BK2017" = 2,
            "Reconstruct_VG2017" = 1, "Reconstruct_VG2018" = 1,
            "Reconstruct_RP2019" = 16)
## for size
size_Host <- c("arabidopsis_thaliana" = 1.4, "lotus_japonicus" = 2.6,
            "host_soil" = 1.4)

get_color_df <- function(x) {
    if (x == "Compartment") return (c_Com)
    if (x == "Plot") return (c_Plo)
}
get_shape_df <- function(x) {
    if (x == "Host.Species") return (s_Host)
    if (x == "Soil") return(s_Soil)
    if (x == "Study") return(s_Study)
}

main_theme <- theme(panel.background = element_blank(),
                    plot.background = element_blank(),
                    panel.grid = element_blank(),
                    axis.line.x = element_line(color = "black"),
                    axis.line.y = element_line(color = "black"),
                    axis.ticks = element_blank(),
                    axis.text = element_text(colour = "black", size = 10),
                    text = element_text(family="sans"))

