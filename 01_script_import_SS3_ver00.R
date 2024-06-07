########################################################################
## Description: Importing SS3 Models and Data to Operating Models
##
## Maintainer: DatenKraft - SCRS/ICCAT (TT Species Group)
## Author: Rodrigo Sant'Ana
## Created: ter abr 30 22:01:34 2023 (-0300)
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

######@> Instaling the required versions of the main packages...

#####@> MSEtool (3.7.2) however, Adrian had answered that version
#####@> 3.7.9999 is OK!
if(packageVersion("MSEtool") == "3.7.9999") {
    print("MSEtool R package already installed is the correct version")
} else {
    pak::pkg_install("blue-matter/MSEtool")
}

#####@> SAMtool (1.6.4)...
if(packageVersion("SAMtool") == "1.6.4") {
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

#####@> nswo-mse (0.25.1)...
if(packageVersion("SWOMSE") == "0.25.1") {
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
library(SWOMSE)

######@> Defining output folder to the OMs...
if(dir.exists("01_OMs")) {
    print("OK, 01_OMs directory was already created!");
    output.dir <- "01_OMs/"
} else {
    dir.create("01_OMs");
    output.dir <- "01_OMs/"
}

######@> Path to the OMs...
path01 <- "00_SS/"

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
load("05_Basement/tsIndex_ver01.RData")

########################################################################
######@> Importing OMs...

######@>================================================================
######@> WSKJ_EstRec93_Qnt25_h6

######@> Path to the file...
SSdir <- paste0(path01, "WSKJ_EstRec93_Qnt25_h6")

######@> W-SKJ EstRec93 Qnt25_h6 - [Method SS2OM]...
OM1 <- SS2OM(SSdir,
             nsim = 100,
             proyears = 30,
             reps = 100,
             maxF = 3,
             seed = 1,
             interval = 3,
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

######@> Importing data from SS...
OM1.Data <- SS2Data(SSdir,
                    Name = "OM1 Data WSKJ_EstRec93_Qnt25_h6",
                    Common_Name = "Skipjack",
                    Species = "Katsuwonus pelamis",
                    Region = "Western Atlantic Ocean")

#####@> Selecting and looking for the inverse-variance combined
#####@> indices...

####@> Preparing the index combined vector...
index <- c(rep(NA, 29),
    tsIndex$Obs[tsIndex$Fleet == "Averaging Index 03"][1:40])
names(index) <- OM1.Data@Year
cv_index <- c(rep(NA, 29), rep(0.2, 40))
names(cv_index) <- OM1.Data@Year

####@> Populate the index slot in Data object...
OM1.Data@Ind <- array(index, dim = c(1, length(OM1.Data@Year)))
OM1.Data@CV_Ind <- array(cv_index, dim = c(1, length(OM1.Data@Year)))

#####@> Adding OM1.Data in cpars data compartment...
OM1@cpars$Data <- OM1.Data

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
             interval = 3,
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

######@> Importing data from SS...
OM2.Data <- SS2Data(SSdir,
                    Name = "OM2 Data WSKJ_EstRec93_Qnt50_h6",
                    Common_Name = "Skipjack",
                    Species = "Katsuwonus pelamis",
                    Region = "Western Atlantic Ocean")

#####@> Selecting and looking for the inverse-variance combined
#####@> indices...

####@> Preparing the index combined vector...
index <- c(rep(NA, 29),
    tsIndex$Obs[tsIndex$Fleet == "Averaging Index 03"][1:40])
names(index) <- OM2.Data@Year
cv_index <- c(rep(NA, 29), rep(0.2, 40))
names(cv_index) <- OM2.Data@Year

####@> Populate the index slot in Data object...
OM2.Data@Ind <- array(index, dim = c(1, length(OM2.Data@Year)))
OM2.Data@CV_Ind <- array(cv_index, dim = c(1, length(OM2.Data@Year)))

#####@> Adding OM2.Data in cpars data compartment...
OM2@cpars$Data <- OM2.Data

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
             interval = 3,
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

######@> Importing data from SS...
OM3.Data <- SS2Data(SSdir,
                    Name = "OM3 Data WSKJ_EstRec93_Qnt75_h6",
                    Common_Name = "Skipjack",
                    Species = "Katsuwonus pelamis",
                    Region = "Western Atlantic Ocean")

#####@> Selecting and looking for the inverse-variance combined
#####@> indices...

####@> Preparing the index combined vector...
index <- c(rep(NA, 29),
    tsIndex$Obs[tsIndex$Fleet == "Averaging Index 03"][1:40])
names(index) <- OM3.Data@Year
cv_index <- c(rep(NA, 29), rep(0.2, 40))
names(cv_index) <- OM3.Data@Year

####@> Populate the index slot in Data object...
OM3.Data@Ind <- array(index, dim = c(1, length(OM3.Data@Year)))
OM3.Data@CV_Ind <- array(cv_index, dim = c(1, length(OM3.Data@Year)))

#####@> Adding OM3.Data in cpars data compartment...
OM3@cpars$Data <- OM3.Data

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
             interval = 3,
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

######@> Importing data from SS...
OM4.Data <- SS2Data(SSdir,
                    Name = "OM4 Data WSKJ_EstRec93_Qnt25_h7",
                    Common_Name = "Skipjack",
                    Species = "Katsuwonus pelamis",
                    Region = "Western Atlantic Ocean")

#####@> Selecting and looking for the inverse-variance combined
#####@> indices...

####@> Preparing the index combined vector...
index <- c(rep(NA, 29),
    tsIndex$Obs[tsIndex$Fleet == "Averaging Index 03"][1:40])
names(index) <- OM4.Data@Year
cv_index <- c(rep(NA, 29), rep(0.2, 40))
names(cv_index) <- OM4.Data@Year

####@> Populate the index slot in Data object...
OM4.Data@Ind <- array(index, dim = c(1, length(OM4.Data@Year)))
OM4.Data@CV_Ind <- array(cv_index, dim = c(1, length(OM4.Data@Year)))

#####@> Adding OM4.Data in cpars data compartment...
OM4@cpars$Data <- OM4.Data

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
             interval = 3,
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

######@> Importing data from SS...
OM5.Data <- SS2Data(SSdir,
                    Name = "OM5 Data WSKJ_EstRec93_Qnt50_h7",
                    Common_Name = "Skipjack",
                    Species = "Katsuwonus pelamis",
                    Region = "Western Atlantic Ocean")

#####@> Selecting and looking for the inverse-variance combined
#####@> indices...

####@> Preparing the index combined vector...
index <- c(rep(NA, 29),
    tsIndex$Obs[tsIndex$Fleet == "Averaging Index 03"][1:40])
names(index) <- OM5.Data@Year
cv_index <- c(rep(NA, 29), rep(0.2, 40))
names(cv_index) <- OM5.Data@Year

####@> Populate the index slot in Data object...
OM5.Data@Ind <- array(index, dim = c(1, length(OM5.Data@Year)))
OM5.Data@CV_Ind <- array(cv_index, dim = c(1, length(OM5.Data@Year)))

#####@> Adding OM5.Data in cpars data compartment...
OM5@cpars$Data <- OM5.Data

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
             interval = 3,
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


######@> Importing data from SS...
OM6.Data <- SS2Data(SSdir,
                    Name = "OM6 Data WSKJ_EstRec93_Qnt75_h7",
                    Common_Name = "Skipjack",
                    Species = "Katsuwonus pelamis",
                    Region = "Western Atlantic Ocean")

#####@> Selecting and looking for the inverse-variance combined
#####@> indices...

####@> Preparing the index combined vector...
index <- c(rep(NA, 29),
    tsIndex$Obs[tsIndex$Fleet == "Averaging Index 03"][1:40])
names(index) <- OM6.Data@Year
cv_index <- c(rep(NA, 29), rep(0.2, 40))
names(cv_index) <- OM6.Data@Year

####@> Populate the index slot in Data object...
OM6.Data@Ind <- array(index, dim = c(1, length(OM6.Data@Year)))
OM6.Data@CV_Ind <- array(cv_index, dim = c(1, length(OM6.Data@Year)))

#####@> Adding OM6.Data in cpars data compartment...
OM6@cpars$Data <- OM6.Data

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
             interval = 3,
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

######@> Importing data from SS...
OM7.Data <- SS2Data(SSdir,
                    Name = "OM7 Data WSKJ_EstRec93_Qnt25_h8",
                    Common_Name = "Skipjack",
                    Species = "Katsuwonus pelamis",
                    Region = "Western Atlantic Ocean")

#####@> Selecting and looking for the inverse-variance combined
#####@> indices...

####@> Preparing the index combined vector...
index <- c(rep(NA, 29),
    tsIndex$Obs[tsIndex$Fleet == "Averaging Index 03"][1:40])
names(index) <- OM7.Data@Year
cv_index <- c(rep(NA, 29), rep(0.2, 40))
names(cv_index) <- OM7.Data@Year

####@> Populate the index slot in Data object...
OM7.Data@Ind <- array(index, dim = c(1, length(OM7.Data@Year)))
OM7.Data@CV_Ind <- array(cv_index, dim = c(1, length(OM7.Data@Year)))

#####@> Adding OM7.Data in cpars data compartment...
OM7@cpars$Data <- OM7.Data

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
             interval = 3,
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

######@> Importing data from SS...
OM8.Data <- SS2Data(SSdir,
                    Name = "OM8 Data WSKJ_EstRec93_Qnt50_h8",
                    Common_Name = "Skipjack",
                    Species = "Katsuwonus pelamis",
                    Region = "Western Atlantic Ocean")

#####@> Selecting and looking for the inverse-variance combined
#####@> indices...

####@> Preparing the index combined vector...
index <- c(rep(NA, 29),
    tsIndex$Obs[tsIndex$Fleet == "Averaging Index 03"][1:40])
names(index) <- OM8.Data@Year
cv_index <- c(rep(NA, 29), rep(0.2, 40))
names(cv_index) <- OM8.Data@Year

####@> Populate the index slot in Data object...
OM8.Data@Ind <- array(index, dim = c(1, length(OM8.Data@Year)))
OM8.Data@CV_Ind <- array(cv_index, dim = c(1, length(OM8.Data@Year)))

#####@> Adding OM8.Data in cpars data compartment...
OM8@cpars$Data <- OM8.Data

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
             interval = 3,
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

######@> Importing data from SS...
OM9.Data <- SS2Data(SSdir,
                    Name = "OM9 Data WSKJ_EstRec93_Qnt75_h8",
                    Common_Name = "Skipjack",
                    Species = "Katsuwonus pelamis",
                    Region = "Western Atlantic Ocean")

#####@> Selecting and looking for the inverse-variance combined
#####@> indices...

####@> Preparing the index combined vector...
index <- c(rep(NA, 29),
    tsIndex$Obs[tsIndex$Fleet == "Averaging Index 03"][1:40])
names(index) <- OM9.Data@Year
cv_index <- c(rep(NA, 29), rep(0.2, 40))
names(cv_index) <- OM9.Data@Year

####@> Populate the index slot in Data object...
OM9.Data@Ind <- array(index, dim = c(1, length(OM9.Data@Year)))
OM9.Data@CV_Ind <- array(cv_index, dim = c(1, length(OM9.Data@Year)))

#####@> Adding OM9.Data in cpars data compartment...
OM9@cpars$Data <- OM9.Data

#####@> Looking to the OM...
plot_SS2OM(OM9,
           SSdir,
           gender = 1,
           filename = "WSKJ_EstRec93_Qnt75_h8",
           dir = output.dir,
           open_file = TRUE,
           silent = FALSE)

######@>----------------------------------------------------------------
######@> Implementation error 10%

######@> Cloning from the first 9 OMs...
OM10 <- OM1
OM10@Name <- paste0(OM1@Name, "_ImpErr_10%")
OM10@TACFrac <- c(1.1, 1.1)

OM11 <- OM2
OM11@Name <- paste0(OM2@Name, "_ImpErr_10%")
OM11@TACFrac <- c(1.1, 1.1)

OM12 <- OM3
OM12@Name <- paste0(OM3@Name, "_ImpErr_10%")
OM12@TACFrac <- c(1.1, 1.1)

OM13 <- OM4
OM13@Name <- paste0(OM4@Name, "_ImpErr_10%")
OM13@TACFrac <- c(1.1, 1.1)

OM14 <- OM5
OM14@Name <- paste0(OM5@Name, "_ImpErr_10%")
OM14@TACFrac <- c(1.1, 1.1)

OM15 <- OM6
OM15@Name <- paste0(OM6@Name, "_ImpErr_10%")
OM15@TACFrac <- c(1.1, 1.1)

OM16 <- OM7
OM16@Name <- paste0(OM7@Name, "_ImpErr_10%")
OM16@TACFrac <- c(1.1, 1.1)

OM17 <- OM8
OM17@Name <- paste0(OM8@Name, "_ImpErr_10%")
OM17@TACFrac <- c(1.1, 1.1)

OM18 <- OM9
OM18@Name <- paste0(OM9@Name, "_ImpErr_10%")
OM18@TACFrac <- c(1.1, 1.1)

######@>----------------------------------------------------------------
######@> Implementation error 20%

######@> Cloning from the first 9 OMs...
OM19 <- OM1
OM19@Name <- paste0(OM1@Name, "_ImpErr_20%")
OM19@TACFrac <- c(1.2, 1.2)

OM20 <- OM2
OM20@Name <- paste0(OM2@Name, "_ImpErr_20%")
OM20@TACFrac <- c(1.2, 1.2)

OM21 <- OM3
OM21@Name <- paste0(OM3@Name, "_ImpErr_20%")
OM21@TACFrac <- c(1.2, 1.2)

OM22 <- OM4
OM22@Name <- paste0(OM4@Name, "_ImpErr_20%")
OM22@TACFrac <- c(1.2, 1.2)

OM23 <- OM5
OM23@Name <- paste0(OM5@Name, "_ImpErr_20%")
OM23@TACFrac <- c(1.2, 1.2)

OM24 <- OM6
OM24@Name <- paste0(OM6@Name, "_ImpErr_20%")
OM24@TACFrac <- c(1.2, 1.2)

OM25 <- OM7
OM25@Name <- paste0(OM7@Name, "_ImpErr_20%")
OM25@TACFrac <- c(1.2, 1.2)

OM26 <- OM8
OM26@Name <- paste0(OM8@Name, "_ImpErr_20%")
OM26@TACFrac <- c(1.2, 1.2)

OM27 <- OM9
OM27@Name <- paste0(OM9@Name, "_ImpErr_20%")
OM27@TACFrac <- c(1.2, 1.2)

########################################################################
######@> Exporting OMs...

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
OMs[[10]] <- OM10
OMs[[11]] <- OM11
OMs[[12]] <- OM12
OMs[[13]] <- OM13
OMs[[14]] <- OM14
OMs[[15]] <- OM15
OMs[[16]] <- OM16
OMs[[17]] <- OM17
OMs[[18]] <- OM18
OMs[[19]] <- OM19
OMs[[20]] <- OM20
OMs[[21]] <- OM21
OMs[[22]] <- OM22
OMs[[23]] <- OM23
OMs[[24]] <- OM24
OMs[[25]] <- OM25
OMs[[26]] <- OM26
OMs[[27]] <- OM27

######@> Image file...
saveRDS(OMs,
        file = "01_OMs/SS_Operating_models_IVInds_ver01.rds")

######@> OMs Descriptions...
OM.Description <- data.frame(
    OMs = paste0("OM", 1:27),
    Name = rep(c("Reference",
                 "Implementation Catch Error 10%",
                 "Implementation Catch Error 20%"), each = 9),
    Description = sapply(OMs, function(x) x@Name),
    Origem = rep(c("SS3 Uncertainty Grid",
                   "SS3 Uncertainty Grid + Robustness 01",
                   "SS3 Uncertainty Grid + Robustness 02"), each = 9),
    Index = rep("Inverse-Variance Weighting", 27)
)

#####@> Exporting OM Description table...
write.table(OM.Description, "01_OMs/OM_Description_ver00.csv",
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
