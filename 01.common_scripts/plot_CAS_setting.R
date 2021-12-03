lst <- c("arabidopsis_thaliana_rhizoplane",
    "arabidopsis_thaliana_rhizosphere",
    "arabidopsis_thaliana_root",
    "CAS",
    "chlamydomonas_reinhardtii_phycosphere",
    "lotus_japonicus_rhizosphere",
    "lotus_japonicus_root")

tags <- c("At_rhizop", "At_rhizos", "At_root", "CAS",
          "Chlamy", "Lj_rhizos", "Lj_root")

c_group <- c("lightgreen", c_rhizos, c_root, c_soil,
             c_phyco, c_rhizos, c_root)

s_group <- c(rep(16, 3), 8, 10, rep(9, 2))

names(tags) <- names(c_group) <- names(s_group) <- lst
