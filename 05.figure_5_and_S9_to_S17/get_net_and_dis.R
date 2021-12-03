#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(mina)

asv_file <- "../00.data/CAS_rep_ASV_reNorm.rds"
design_file <- "../00.data/design_854.txt"
out_dir <- "./CAS_spectra_egv/"

asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")

des_soil <- design[design$Compartment == "soil", ]
des_soil$Group <- des_soil$Soil

des_nsoil <- design[design$Compartment != "soil", ]
des_nsoil$Group <- paste0(des_nsoil$Host.Species, "_", des_nsoil$Compartment)
des <- rbind(des_soil, des_nsoil)

dat <- new("mina", norm = asv, des = des)

if (!dir.exists(bs_pm_dir)) dir.create(bs_pm_dir)
if (!dir.exists(out_dir)) dir.create(out_dir)

## for the stat calculation
dat <- bs_pm(dat, group = "Group", g_size = 40, s_size = 20,
            rm = FALSE, sig = TRUE, bs = 33, pm = 33,
            individual = FALSE)

dat <- net_dis(dat, method = "spectra", egv = TRUE, dir = out_dir)
dis <- dis_stat(dat)
write.table(dis, "Figure_5b.txt", quote = F, sep = "\t", row.names = F)
