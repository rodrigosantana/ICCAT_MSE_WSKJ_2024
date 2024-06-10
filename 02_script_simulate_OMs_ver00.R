########################################################################
## Description: Simulating Operating Models based on the new OMs...
##
## Maintainer: Datenkraft - ICCAT/SCRS - TT MSE
## Author: Rodrigo Sant'Ana
## Created: dom set 10 14:20:53 2023 (-0300)
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
######@> Loading R packages...

######@> Package list...
library(dplyr)
library(openMSE)
## library(SWOMSE)
library(ggmse)

######@> Packages versions...
packageVersion("openMSE")
## [1] ‘1.1.1’
packageVersion("SAMtool")
## [1] ‘1.6.5’
packageVersion("MSEtool")
## [1] ‘3.7.9999’
packageVersion("DLMtool")
## [1] ‘6.0.6’
packageVersion("SWOMSE")
## [1] ‘0.26.0’

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

######@> Creating a folder to receive the hist objects...
if(!dir.exists("02_Hists"))
    dir.create("02_Hists")

######@> Creating a folder to receive the MSE objects...
if(!dir.exists("03_MSEs"))
    dir.create("03_MSEs")

######@> OM folder...
path01 <- "01_OMs/"

######@> MPs folder...
path02 <- "../R-Work/"

######@> Catch for recent years...
Catchdf <- data.frame(
    Year = 2021:2022,
    Catch = c(20048.21, 21377.24))

########################################################################
######@> Loading data...

######@> Operating models...
OMs <- readRDS(paste0(path01, "SS_Operating_models_IVInds_ver01.rds"),
               refhook = NULL)

#####@> Setting names for OMs...
names(OMs) <- paste0("OM", sprintf("%03d", 1:27))

######@> Management procedures...
source("03_script_prepare_MPs_ver00.R")

#####@> Defining MPs...
MPs <- c("curE", "CC_20kt", "CC_30kt", "CC_40kt",
         "GBslope_MOD", "Iratio_MOD", "Islope1_MOD",
         ## "SCA_00", "SCA_01", "SCA_02",
         "SP_01", "SP_02", "SP_03", "SP_04")

########################################################################
######@> Simulating Historical Data...

######@> List of objects...
OM_Objects <- names(OMs)

######@> Looping to create historical data...
for(i in seq_along(OM_Objects)) {
    OM <- OMs[[i]]
    OM@interval <- 3
    Hist <- Simulate(OM, parallel = FALSE, silent = FALSE)
    nm <- paste0(OM_Objects[i], "_IVInds", ".hist")
    saveRDS(Hist, file.path("02_Hists", nm))
}

########################################################################
######@> Simulating Projections Data...

######@> Loading history objects...
Hists <- sapply(dir("02_Hists", full.names = TRUE), readRDS)

#####@> Setting names for Hists...
names(Hists) <- paste0("OM", sprintf("%03d", 1:27))

######@> Looping to run the MSE for all histories data...
for(i in seq_along(Hists)) {
    Hist <- Hists[[i]]
    MSE <- Project(Hist, MPs = MPs, parallel = FALSE)
    nm <- paste0("MSE", sprintf("%03d", i), "_IVInds", "_00.mse")
    saveRDS(MSE, file.path("03_MSEs", nm))
}

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
