########################################################################
## Description: Importing SS3 Models and Data to Operating Models
##
## Maintainer: DatenKraft - SCRS/ICCAT (TT Species Group)
## Author: Rodrigo Sant'Ana
## Created: ter ago 23 17:10:34 2023 (-0300)
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
if(packageVersion("MSEtool") == "3.7.2") {
    print("MSEtool R package already installed is the correct version")
} else {
    pak::pkg_install("blue-matter/MSEtool")
}

#####@> SAMtool (1.7.0)...
if(packageVersion("SAMtool") == "1.7.0") {
    print("SAMtool R package already installed is the correct version")
} else {
    pak::pkg_install("blue-matter/SAMtool")
}

#####@> openMSE (1.1.1)...
if(packageVersion("openMSE") == "1.1.1") {
    print("openMSE R package already installed is the correct version")
} else {
    pak::pkg_install("blue-matter/openMSE")
}

#####@> nswo-mse (0.29.0)...
if(packageVersion("SWOMSE") == "0.29.0") {
    print("SWOMSE R package already installed is the correct version")
} else {
    pak::pkg_install("ICCAT/nswo-mse")
}

######@> Loading R packages...
library(openMSE)
library(dplyr)
library(reshape2)
library(readr)
library(ggplot2)
library(tidyverse)
library(viridis)
library(janitor)
library(readxl)
library(ggmse)

######@> Instaling MSEextra...
## MSEextra(force = TRUE)
library(MSEextra)
## library(SWOMSE)

######@> Defining output folder to the OMs...
if(dir.exists("02_OMs")) {
    print("OK, 02_OMs directory was already created!");
    output.dir <- "02_OMs/"
} else {
    dir.create("02_OMs");
    output.dir <- "02_OMs/"
}

######@> Path to the OMs...
path01 <- "01_SS/"

######@> Setting up ggplot theme...

#####@> Importing fonts...
extrafont::loadfonts(device = "postscript")

#####@> ggplot theme...
rgb01 <- "black"
rgb02 <- "black"
seta <- grid::arrow(length = grid::unit(0.2, "cm"), type = "open")
seta2 <- grid::arrow(length = grid::unit(0.2, "cm"), type = "open",
                     ends = "both")
my_theme <- function(base_size = 18, base_family = "Helvetica") {
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
        theme(axis.ticks = element_line(colour = rgb01),
              axis.line = element_line(colour = rgb01, linewidth = 0.2),
              axis.text = element_text(colour = rgb02, size = 18),
              axis.title = element_text(size = 18),
              strip.text = element_text(size = 22,
                                        margin = ggplot2::margin(0.3,
                                                                 0.3,
                                                                 0.3,
                                                                 0.3,
                                                                 "cm"),
                                        face = "bold"),
              legend.background = element_blank(),
              legend.key = element_blank(),
              legend.text = element_text(size = 18),
              legend.title = element_text(size = 20),
              panel.background = element_blank(),
              panel.grid = element_line(linetype = "solid",
                                        linewidth = 0.2,
                                        colour = "gray90"),
              plot.background = element_blank(),
              complete = TRUE)
}

#####@> Testing theme...
df <- data.frame(x = rnorm(10), y = rnorm(10),
                 z = rep(c("A", "B"), each = 5))
ggplot(data = df, aes(x = x, y = y, fill = z)) +
    geom_point(pch = 21, size = 5) +
    facet_grid(~z) +
    my_theme()

########################################################################
######@> Loading prepared dataset...

######@> Loading indices scenarios...
load("05_Results/tsIndex_ver02.RData")

######@> Loading extracted quantities (ie. MSY, SSB_MSY, F_MSY) from the
######@> Stock Assessment and catches from T1NC for the recent year...
load("05_Results/Reference_Quantities_ver02.RData")

########################################################################
######@> Importing OMs...

######@> Set Some Global Objects for All OMs...

#####@> Set interval to 1
interval <- 1

#####@> Create Data object to be used in all OMs...
SSdir <- paste0(path01, "WSKJ_EstRec93_Qnt25_h6")
OM1.Data <- SS2Data(SSdir,
                    Name = "OM1 Data WSKJ_EstRec93_Qnt25_h6",
                    Common_Name = "Skipjack",
                    Species = "Katsuwonus pelamis",
                    Region = "Western Atlantic Ocean")

####@> Historical Data should be identical across all OMs...
WSKJ_Data <- new('Data')
## slots_to_copy <- c("Year", "Cat", "CV_Cat")
slots_to_copy <- slotNames(OM1.Data)
for (sl in slots_to_copy) {
    slot(WSKJ_Data, sl) <- slot(OM1.Data, sl)
}

####@> Inverse-variance weighted index manually added...
index <- c(rep(NA, 29),
           tsIndex$Obs[tsIndex$Fleet ==
                       "Inverse variance weighted average"][1:40])
index <- index/mean(index, na.rm = TRUE)
names(index) <- OM1.Data@Year
cv_index <- c(rep(NA, 29), rep(0.2, 40))
names(cv_index) <- OM1.Data@Year
WSKJ_Data@Ind <- array(index, dim = c(1, length(index)))
WSKJ_Data@CV_Ind <- array(cv_index, dim = c(1, length(index)))
## WSKJ_Data@Year <- 1952:2022

## ####@> Catches for the recent period manually added...
## catches <- c(WSKJ_Data@Cat[1, ], 20048.21, 21377.24)
## cv_catches <- c(WSKJ_Data@CV_Cat, 0.2, 0.2)
## WSKJ_Data@Cat <- array(catches, dim = c(1, length(catches)))
## WSKJ_Data@CV_Cat <- array(cv_catches, dim = c(1, length(cv_catches)))

####@> Function to add WSJK_Data and I_beta = 1 to cpars...
update_cpars <- function(OM) {
    OM@cpars$Data <- WSKJ_Data
    OM@cpars$I_beta <- rep(1, OM@nsim)
    OM
}

######@>================================================================
######@> WSKJ_EstRec93_Qnt25_h6...

######@> Path to the file...
SSdir <- paste0(path01, "WSKJ_EstRec93_Qnt25_h6")

######@> W-SKJ EstRec93 Qnt25_h6 - [Method SS2OM]...
OM1 <- SS2OM(SSdir,
             nsim = 100,
             proyears = 30,
             reps = 100,
             maxF = 3,
             seed = 1,
             interval = interval,
             pstar = 0.5,
             Obs = MSEtool::Precise_Unbiased,
             Imp = MSEtool::Perfect_Imp,
             gender = 1,
             silent = FALSE,
             Name = "OM WSKJ_EstRec93_Qnt25_h6",
             Source = "ICCAT 2022 SKJ assessment",
             Author = "Rodrigo Sant'Ana & Bruno Leite Mourato",
             report = FALSE,
             filename = "OM_WSKJ_EstRec93_Qnt25_h6",
             dir = output.dir,
             open_file = TRUE)

#####@> Update OM1 cpars data compartment...
OM1 <- update_cpars(OM1)

#####@> Looking to the OM...
plot_SS2OM(OM1,
           SSdir,
           gender = 1,
           filename = "WSKJ_EstRec93_Qnt25_h6",
           dir = output.dir,
           open_file = TRUE,
           silent = FALSE)

######@>================================================================
######@> WSKJ_EstRec93_Qnt50_h6

######@> Path to the file...
SSdir <- paste0(path01, "WSKJ_EstRec93_Qnt50_h6")

######@> W-SKJ EstRec93 Qnt50_h6 - [Method SS2OM]...
OM2 <- SS2OM(SSdir,
             nsim = 100,
             proyears = 30,
             reps = 100,
             maxF = 3,
             seed = 1,
             interval = interval,
             pstar = 0.5,
             Obs = MSEtool::Precise_Unbiased,
             Imp = MSEtool::Perfect_Imp,
             gender = 1,
             silent = FALSE,
             Name = "OM WSKJ_EstRec93_Qnt50_h6",
             Source = "ICCAT 2022 SKJ assessment",
             Author = "Rodrigo Sant'Ana & Bruno Leite Mourato",
             report = FALSE,
             filename = "OM_WSKJ_EstRec93_Qnt50_h6",
             dir = output.dir,
             open_file = TRUE)

#####@> Update OM2 cpars data compartment...
OM2 <- update_cpars(OM2)

#####@> Looking to the OM...
plot_SS2OM(OM2,
           SSdir,
           gender = 1,
           filename = "WSKJ_EstRec93_Qnt50_h6",
           dir = output.dir,
           open_file = TRUE,
           silent = FALSE)

######@>================================================================
######@> WSKJ_EstRec93_Qnt75_h6

######@> Path to the file...
SSdir <- paste0(path01, "WSKJ_EstRec93_Qnt75_h6")

######@> W-SKJ EstRec93 Qnt75_h6 - [Method SS2OM]...
OM3 <- SS2OM(SSdir,
             nsim = 100,
             proyears = 30,
             reps = 100,
             maxF = 3,
             seed = 1,
             interval = interval,
             pstar = 0.5,
             Obs = MSEtool::Precise_Unbiased,
             Imp = MSEtool::Perfect_Imp,
             gender = 1,
             silent = FALSE,
             Name = "OM WSKJ_EstRec93_Qnt75_h6",
             Source = "ICCAT 2022 SKJ assessment",
             Author = "Rodrigo Sant'Ana & Bruno Leite Mourato",
             report = FALSE,
             filename = "OM_WSKJ_EstRec93_Qnt75_h6",
             dir = output.dir,
             open_file = TRUE)

#####@> Update OM3 cpars data compartment...
OM3 <- update_cpars(OM3)

#####@> Looking to the OM...
plot_SS2OM(OM3,
           SSdir,
           gender = 1,
           filename = "WSKJ_EstRec93_Qnt75_h6",
           dir = output.dir,
           open_file = TRUE,
           silent = FALSE)

######@>================================================================
######@> WSKJ_EstRec93_Qnt25_h7

######@> Path to the file...
SSdir <- paste0(path01, "WSKJ_EstRec93_Qnt25_h7")

######@> W-SKJ EstRec93 Qnt25_h7 - [Method SS2OM]...
OM4 <- SS2OM(SSdir,
             nsim = 100,
             proyears = 30,
             reps = 100,
             maxF = 3,
             seed = 1,
             interval = interval,
             pstar = 0.5,
             Obs = MSEtool::Precise_Unbiased,
             Imp = MSEtool::Perfect_Imp,
             gender = 1,
             silent = FALSE,
             Name = "OM WSKJ_EstRec93_Qnt25_h7",
             Source = "ICCAT 2022 SKJ assessment",
             Author = "Rodrigo Sant'Ana & Bruno Leite Mourato",
             report = FALSE,
             filename = "OM_WSKJ_EstRec93_Qnt25_h7",
             dir = output.dir,
             open_file = TRUE)

#####@> Update OM4 cpars data compartment...
OM4 <- update_cpars(OM4)

#####@> Looking to the OM...
plot_SS2OM(OM4,
           SSdir,
           gender = 1,
           filename = "WSKJ_EstRec93_Qnt25_h7",
           dir = output.dir,
           open_file = TRUE,
           silent = FALSE)

######@>================================================================
######@> WSKJ_EstRec93_Qnt50_h7

######@> Path to the file...
SSdir <- paste0(path01, "WSKJ_EstRec93_Qnt50_h7")

######@> W-SKJ EstRec93 Qnt50_h7 - [Method SS2OM]...
OM5 <- SS2OM(SSdir,
             nsim = 100,
             proyears = 30,
             reps = 100,
             maxF = 3,
             seed = 1,
             interval = interval,
             pstar = 0.5,
             Obs = MSEtool::Precise_Unbiased,
             Imp = MSEtool::Perfect_Imp,
             gender = 1,
             silent = FALSE,
             Name = "OM WSKJ_EstRec93_Qnt50_h7",
             Source = "ICCAT 2022 SKJ assessment",
             Author = "Rodrigo Sant'Ana & Bruno Leite Mourato",
             report = FALSE,
             filename = "OM_WSKJ_EstRec93_Qnt50_h7",
             dir = output.dir,
             open_file = TRUE)

#####@> Update OM5 cpars data compartment...
OM5 <- update_cpars(OM5)

#####@> Looking to the OM...
plot_SS2OM(OM5,
           SSdir,
           gender = 1,
           filename = "WSKJ_EstRec93_Qnt50_h7",
           dir = output.dir,
           open_file = TRUE,
           silent = FALSE)

######@>================================================================
######@> WSKJ_EstRec93_Qnt75_h7

######@> Path to the file...
SSdir <- paste0(path01, "WSKJ_EstRec93_Qnt75_h7")

######@> W-SKJ EstRec93 Qnt75_h7 - [Method SS2OM]...
OM6 <- SS2OM(SSdir,
             nsim = 100,
             proyears = 30,
             reps = 100,
             maxF = 3,
             seed = 1,
             interval = interval,
             pstar = 0.5,
             Obs = MSEtool::Precise_Unbiased,
             Imp = MSEtool::Perfect_Imp,
             gender = 1,
             silent = FALSE,
             Name = "OM WSKJ_EstRec93_Qnt75_h7",
             Source = "ICCAT 2022 SKJ assessment",
             Author = "Rodrigo Sant'Ana & Bruno Leite Mourato",
             report = FALSE,
             filename = "OM_WSKJ_EstRec93_Qnt75_h7",
             dir = output.dir,
             open_file = TRUE)

#####@> Update OM6 cpars data compartment...
OM6 <- update_cpars(OM6)

#####@> Looking to the OM...
plot_SS2OM(OM6,
           SSdir,
           gender = 1,
           filename = "WSKJ_EstRec93_Qnt75_h7",
           dir = output.dir,
           open_file = TRUE,
           silent = FALSE)

######@>================================================================
######@> WSKJ_EstRec93_Qnt25_h8

######@> Path to the file...
SSdir <- paste0(path01, "WSKJ_EstRec93_Qnt25_h8")

######@> W-SKJ EstRec93 Qnt25_h8 - [Method SS2OM]...
OM7 <- SS2OM(SSdir,
             nsim = 100,
             proyears = 30,
             reps = 100,
             maxF = 3,
             seed = 1,
             interval = interval,
             pstar = 0.5,
             Obs = MSEtool::Precise_Unbiased,
             Imp = MSEtool::Perfect_Imp,
             gender = 1,
             silent = FALSE,
             Name = "OM WSKJ_EstRec93_Qnt25_h8",
             Source = "ICCAT 2022 SKJ assessment",
             Author = "Rodrigo Sant'Ana & Bruno Leite Mourato",
             report = FALSE,
             filename = "OM_WSKJ_EstRec93_Qnt25_h8",
             dir = output.dir,
             open_file = TRUE)

#####@> Update OM7 cpars data compartment...
OM7 <- update_cpars(OM7)

#####@> Looking to the OM...
plot_SS2OM(OM7,
           SSdir,
           gender = 1,
           filename = "WSKJ_EstRec93_Qnt25_h8",
           dir = output.dir,
           open_file = TRUE,
           silent = FALSE)

######@>================================================================
######@> WSKJ_EstRec93_Qnt50_h8

######@> Path to the file...
SSdir <- paste0(path01, "WSKJ_EstRec93_Qnt50_h8")

######@> W-SKJ EstRec93 Qnt50_h8 - [Method SS2OM]...
OM8 <- SS2OM(SSdir,
             nsim = 100,
             proyears = 30,
             reps = 100,
             maxF = 3,
             seed = 1,
             interval = interval,
             pstar = 0.5,
             Obs = MSEtool::Precise_Unbiased,
             Imp = MSEtool::Perfect_Imp,
             gender = 1,
             silent = FALSE,
             Name = "OM WSKJ_EstRec93_Qnt50_h8",
             Source = "ICCAT 2022 SKJ assessment",
             Author = "Rodrigo Sant'Ana & Bruno Leite Mourato",
             report = FALSE,
             filename = "OM_WSKJ_EstRec93_Qnt50_h8",
             dir = output.dir,
             open_file = TRUE)

#####@> Update OM8 cpars data compartment...
OM8 <- update_cpars(OM8)

#####@> Looking to the OM...
plot_SS2OM(OM8,
           SSdir,
           gender = 1,
           filename = "WSKJ_EstRec93_Qnt50_h8",
           dir = output.dir,
           open_file = TRUE,
           silent = FALSE)

######@>================================================================
######@> WSKJ_EstRec93_Qnt75_h8

######@> Path to the file...
SSdir <- paste0(path01, "WSKJ_EstRec93_Qnt75_h8")

######@> W-SKJ EstRec93 Qnt75_h8 - [Method SS2OM]...
OM9 <- SS2OM(SSdir,
             nsim = 100,
             proyears = 30,
             reps = 100,
             maxF = 3,
             seed = 1,
             interval = interval,
             pstar = 0.5,
             Obs = MSEtool::Precise_Unbiased,
             Imp = MSEtool::Perfect_Imp,
             gender = 1,
             silent = FALSE,
             Name = "OM WSKJ_EstRec93_Qnt75_h8",
             Source = "ICCAT 2022 SKJ assessment",
             Author = "Rodrigo Sant'Ana & Bruno Leite Mourato",
             report = FALSE,
             filename = "OM_WSKJ_EstRec93_Qnt75_h8",
             dir = output.dir,
             open_file = TRUE)

#####@> Update OM9 cpars data compartment...
OM9 <- update_cpars(OM9)

#####@> Looking to the OM...
plot_SS2OM(OM9,
           SSdir,
           gender = 1,
           filename = "WSKJ_EstRec93_Qnt75_h8",
           dir = output.dir,
           open_file = TRUE,
           silent = FALSE)

######@> Listing OMs...
OMs <- list()
OMs[[1]] <- OM1
OMs[[2]] <- OM2
OMs[[3]] <- OM3
OMs[[4]] <- OM4
OMs[[5]] <- OM5
OMs[[6]] <- OM6
OMs[[7]] <- OM7
OMs[[8]] <- OM8
OMs[[9]] <- OM9

######@> Image file...
saveRDS(OMs,
        file = "02_OMs/SS_Operating_models_IVInds_ver02.rds")

######@> OMs Descriptions...
OM.Description <- data.frame(
    OMs = paste0("OM", 1:9),
    Name = rep(c("Reference"), each = 9),
    Description = sapply(OMs, function(x) x@Name),
    Origem = rep(c("SS3 Uncertainty Grid"), each = 9),
    Index = rep("Inverse-Variance Weighting", 9)
)

#####@> Exporting OM Description table...
write.table(OM.Description, "02_OMs/OM_Description_ver00.csv",
            row.names = TRUE, sep = ";", dec = ",")

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
