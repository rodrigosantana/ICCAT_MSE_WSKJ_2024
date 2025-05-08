########################################################################
## Description: Projecting Operating Models based on the new Hists
## Files...
##
## Maintainer: Datenkraft - ICCAT/SCRS - TT MSE
## Author: Rodrigo Sant'Ana
## Created: seg set 02 14:20:53 2024 (-0300)
## Version: 0.0.1
##
## URL:
## Doc URL:
##
## Database info:
##
### Commentary:
##
### Code:
########################################################################

########################################################################
######@> Setup R...

######@> Update openMSE R packages...

#####@> MSEtool (3.7.2)...
if(packageVersion("MSEtool") == "3.7.9999") {
    print("MSEtool R package already installed is the correct version")
} else {
    pak::pkg_install("blue-matter/MSEtool")
}

#####@> SAMtool (1.7.0)...
if(packageVersion("SAMtool") == "1.8.1") {
    print("SAMtool R package already installed is the correct version")
} else {
    pak::pkg_install("blue-matter/SAMtool")
}

#####@> openMSE (1.1.1)...
if(packageVersion("openMSE") == "1.0.1") {
    print("openMSE R package already installed is the correct version")
} else {
    pak::pkg_install("blue-matter/openMSE")
}

# #####@> nswo-mse (0.29.0)...
# if(packageVersion("SWOMSE") == "0.29.0") {
#     print("SWOMSE R package already installed is the correct version")
# } else {
#     pak::pkg_install("ICCAT/nswo-mse")
# }

######@> Package list...
library(dplyr)
library(openMSE)
# library(ggmse)

######@> ggplot theme...
extrafont::loadfonts(device = "postscript")
rgb01 <- "black"
rgb02 <- "black"
seta <- grid::arrow(length = grid::unit(0.2, "cm"), type = "open")
seta2 <- grid::arrow(length = grid::unit(0.2, "cm"), type = "open",
                     ends = "both")
my_theme <- function(base_size = 18, base_family = "Helvetica") {
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
        theme(axis.ticks = element_line(colour = rgb01),
              axis.line = element_line(colour = rgb01, linewidth = 0.2),
              axis.text = element_text(colour = rgb02, size = 14),
              axis.title = element_text(size = 18),
              legend.background = element_blank(),
              legend.key = element_blank(),
              panel.background = element_blank(),
              panel.grid = element_line(linetype = "solid",
                                        linewidth = 0.1,
                                        colour = "gray90"),
              plot.background = element_blank(),
              complete = TRUE)
}

########################################################################
######@> Setup MSE...

######@> Creating a folder to receive the MSE objects...
if(!dir.exists("04_MSEs"))
    dir.create("04_MSEs")

######@> OM folder...
path01 <- "03_Hists"

######@> MPs folder...
path02 <- "../R-Work/"

########################################################################
######@> Loading data...

# ######@> Loading history objects...
# Hists <- sapply(dir("03_Hists",
#                     pattern = "IVInds_CORRECTED_ver03",
#                     full.names = TRUE), readRDS)
# 
# ####@> Setting names for Hists...
# names(Hists) <- paste0("OM", sprintf("%03d", 1:9))
# 
# ######@> Loading MPs...
# source("03_script_prepare_MPs_ver00.R")
# MPs <- c("FMSYref", "FMSYref75", "FMSYref110", ## Reference FMSY
#          "CE_01", "CE_02", "CE_03", ## Index-based: Exploitation rate
#          "IR_01", "IR_02", "IR_03", ## Index-based: I_ratio
#          "IS_01", "IS_02", "IS_03", ## Index-based: I_slope
#          "SP_01", "SP_02", "SP_03", "SP_04") ## Model-based: SP

########################################################################
######@> Simulating Projections Data...


##### MADE DIFFERENT VERSION IN 05A_Project_Hists.R #####

# Hists <- readRDS("03_Hists/HistList.rda")
# 
# ManagementOption <- 'DataLag_1_Interval_3'
# 
# tunedMPs <- list.files(file.path('TunedMPs', ManagementOption), full.names = TRUE)
# 
# MPnames <- gsub('.mp', '', basename(tunedMPs))
# 
# for (i in seq_along(tunedMPs)) {
#   mp <- readRDS(tunedMPs[i])
#   assign(MPnames[i], mp, envir=.GlobalEnv)
# }
# 
# 
# 
# ######@> Looping to run the MSE for all histories data...
# for(i in seq_along(Hists)) {
#     Hist <- Hists[[i]]
#     nm <- gsub('OM', '', Hist@OM@Name) |> trimws() 
#     nm <- paste(sprintf("%03d", i), nm, sep="_")
#     nm <- paste0(nm, '.mse')
#     if (!dir.exists(file.path("04_MSEs", ManagementOption)))
#       dir.create(file.path("04_MSEs", ManagementOption))
#     
#     MSE <- Project(Hist, MPs = MPnames, parallel = FALSE)
#     saveRDS(MSE, file.path("04_MSEs", ManagementOption, nm))
# }

########################################################################
##
##                  Creative Commons License 4.0
##                       (CC BY-NC-SA 4.0)
##
##  This is a humam-readable summary of (and not a substitute for) the
##  license (https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
##
##  You are free to:
##
##  Share - copy and redistribute the material in any medium or format.
##
##  The licensor cannot revoke these freedoms as long as you follow the
##  license terms.
##
##  Under the following terms:
##
##  Attribution - You must give appropriate credit, provide a link to
##  license, and indicate if changes were made. You may do so in any
##  reasonable manner, but not in any way that suggests the licensor
##  endorses you or your use.
##
##  NonCommercial - You may not use the material for commercial
##  purposes.
##
##  ShareAlike - If you remix, transform, or build upon the material,
##  you must distributive your contributions under the same license
##  as the  original.
##
##  No additional restrictions — You may not apply legal terms or
##  technological measures that legally restrict others from doing
##  anything the license permits.
##
########################################################################
