########################################################################
## Description: Testing different tunnings in the MPs looking for
## maximizing Yield without crash PGK PM...
##
## Maintainer: Datenkraft (DFKT) - ICCAT W-SKJ MSE
## Author: Rodrigo Sant'Ana
## Created: qui jun 13 14:02:21 2024 (-0300)
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
library(ggmse)

########################################################################
######@> Setup R...

######@> Global environment variables...
## Sys.setenv(Initial_MP_Yr = 2021) ##> NOTE: Not working...
## Sys.getenv("Initial_MP_Yr")

####@> Old method...
Initial_MP_Yr <- 2021

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
######@> Loading data...

######@> Catch for recent years...
Catchdf <- data.frame(
    Year = 2021:2022,
    Catch = c(20048.21, 21377.24))

######@> Loading history objects...
Hists <- sapply(dir("../02_Hists", pattern = "ver02",
    full.names = TRUE), readRDS)

#####@> Setting names for Hists...
names(Hists) <- paste0("OM", sprintf("%03d", 1:27))

#####@> Data object that will be used in the first projection year...
Data <- Hists[[5]]@Data

#####@> Reference case Hists...
HistsRF <- Hists[1:9]

######@> Loading default MPs...
source("../03_script_prepare_MPs_ver00.R")

######@> Loading default PMs...
source("../04_script_prepare_PMs_ver00.R")

########################################################################
######@> Tunning Model-Based MPs...

######@> Changing tunning parameter...

#####@>-----------------------------------------------------------------
#####@> Model based on the SP_01...
modSP01 <- SAMtool::SP(1, Data, prior = list(r = c(0.416, 0.148),
    MSY = c(30000, 0.2)), start = list(dep = 0.98, n = 1))

#####@> SP_01...
SP_01 <- function(x, Data, Data_Lag = 1, Interval = 3, tunepar = 1,
                  mc = c(0.25, 0.1),...) {
    Rec <- new('Rec')
    ## Does TAC need to be updated? (or set a fixed catch if before
    ## Initial_MP_Yr)
    if(SWOMSE::SameTAC(Initial_MP_Yr, Interval, Data)) {
        Rec@TAC <- Data@MPrec[x]
        Rec <- SWOMSE::FixedTAC(Rec, Data) ## use actual catches if they are
                                   ## available
        return(Rec)
    }
    ## Lag Data
    Data <- Lag_Data(Data, Data_Lag)
    Data@Year <- Data@Year[1:(length(Data@Year)-Data_Lag)]
    ## Smooth combined index
    index <- smoothed_index <- Data@Ind[x,]
    smoothed <- stats::smooth(index[!is.na(index)])
    smoothed_index[!is.na(smoothed_index)] <- smoothed
    Data@Ind[x,] <- smoothed_index
    ## apply SP assessment model
    Mod <- SAMtool::SP(x, Data, prior = list(r = c(0.416, 0.148),
        MSY = c(30000, 0.2)), start = list(dep = 0.98, n = 1))
    ## harvest control rule
    ## based on: https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
    Bthresh <- Mod@BMSY
    Blim <- 0.4 * Bthresh
    ## Ftar <- 0.8 * tunepar * Mod@FMSY
    ## Fmin <- 0.1 * tunepar * Mod@FMSY
    Ftar <- tunepar * 0.085534608 ## Mod@FMSY
    Fmin <- 0.1 * Ftar
    Bcurr <- Mod@B[length(Mod@B)]
    if (Bcurr>=Bthresh) {
        Fmort <- Ftar
    } else if (Bcurr>Blim) {
        Fmort <- Ftar * (-0.367 + 1.167 * Bcurr/Bthresh)
    } else {
        Fmort <- Fmin
    }
    TAC <-  Fmort*Bcurr
    Rec@TAC <- TAC
    Rec@TAC <- MaxChange(TAC, Data@MPrec[x], mc)
    Rec
}
class(SP_01) <- 'MP'

#####@> Tunning SP_01...

####@> tunepar = 1.03
SP_01A <- SP_01
formals(SP_01A)$tunepar <- 1.01
class(SP_01A) <- "MP"

####@> tunepar = 1.025
SP_01B <- SP_01
formals(SP_01B)$tunepar <- 1.02
class(SP_01B) <- "MP"

####@> tunepar = 1.02
SP_01C <- SP_01
formals(SP_01C)$tunepar <- 1.03
class(SP_01C) <- "MP"

####@> tunepar = 1.015
SP_01D <- SP_01
formals(SP_01D)$tunepar <- 1.04
class(SP_01D) <- "MP"

######@> Testing tuned MPs...
MSE <- Project(Hists[[6]], MPs = c("SP_01", "SP_01A", "SP_01B",
    "SP_01C", "SP_01D"))

#####@> Observing the results...

####@> Projections...
plot_projection_ts(MSE, "SSB", bbmsy_zones = c(0.4, 0.8, 1))
plot_projection_ts(MSE, "FM", bbmsy_zones = c(0.4, 0.8, 1))
plot_projection_ts(MSE, "Catch", bbmsy_zones = c(0.4, 0.8, 1))

####@> Probability of being in a Green Quadrant...
PGK_short(MSE)
PGK_med(MSE)
PGK_long(MSE)
PGK(MSE)

#####@>-----------------------------------------------------------------
#####@> SP_02...

#####@> MP SP_02...
SP_02 <- function(x, Data, Data_Lag = 1, Interval = 3, tunepar = 1,
                  mc = c(0.2, 0.1), ...) {
    Rec <- new('Rec')
    ## Does TAC need to be updated? (or set a fixed catch if before
    ## Initial_MP_Yr)
    if(SWOMSE::SameTAC(Initial_MP_Yr, Interval, Data)) {
        Rec@TAC <- Data@MPrec[x]
        Rec <- SWOMSE::FixedTAC(Rec, Data) ## use actual catches if they are
                                   ## available
        return(Rec)
    }
    ## Lag Data
    Data <- Lag_Data(Data, Data_Lag)
    Data@Year <- Data@Year[1:(length(Data@Year)-Data_Lag)]
    ## Smooth combined index
    index <- smoothed_index <- Data@Ind[x,]
    smoothed <- stats::smooth(index[!is.na(index)])
    smoothed_index[!is.na(smoothed_index)] <- smoothed
    Data@Ind[x,] <- smoothed_index
    ## apply SP assessment model
    Mod <- SAMtool::SP_SS(x, Data, prior = list(r = c(0.416, 0.148),
        MSY = c(30000, 0.2)), start = list(dep = 0.98, n = 1),
        fix_sigma = FALSE, fix_tau = TRUE)
    ## harvest control rule
    ## based on: https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
    Bthresh <- Mod@BMSY
    Blim <- 0.4 * Bthresh
    ## Ftar <- 0.8 * tunepar * Mod@FMSY
    ## Fmin <- 0.1 * tunepar * Mod@FMSY
    Ftar <- tunepar * 0.085534608
    Fmin <- 0.1 * Ftar
    Bcurr <- Mod@B[length(Mod@B)]
    if (Bcurr>=Bthresh) {
        Fmort <- Ftar
    } else if (Bcurr>Blim) {
        Fmort <- Ftar * (-0.367 + 1.167 * Bcurr/Bthresh)
    } else {
        Fmort <- Fmin
    }
    TAC <-  Fmort*Bcurr
    Rec@TAC <- TAC
    Rec@TAC <- MaxChange(TAC, Data@MPrec[x], mc)
    Rec
}
class(SP_02) <- 'MP'

#####@> Tunning SP_02...

####@> tunepar = 1.03
SP_02A <- SP_02
formals(SP_02A)$tunepar <- 1.06
class(SP_02A) <- "MP"

####@> tunepar = 1.025
SP_02B <- SP_02
formals(SP_02B)$tunepar <- 1.05
class(SP_02B) <- "MP"

####@> tunepar = 1.02
SP_02C <- SP_02
formals(SP_02C)$tunepar <- 1.04
class(SP_02C) <- "MP"

####@> tunepar = 1.015
SP_02D <- SP_02
formals(SP_02D)$tunepar <- 1.03
class(SP_02D) <- "MP"

######@> Testing tuned MPs...
MSE <- Project(Hists[[5]], MPs = c("SP_02", "SP_02A", "SP_02B",
    "SP_02C", "SP_02D"))

#####@> Observing the results...

####@> Projections...
plot_projection_ts(MSE, "SSB", bbmsy_zones = c(0.4, 0.8, 1))
plot_projection_ts(MSE, "FM", bbmsy_zones = c(0.4, 0.8, 1))
plot_projection_ts(MSE, "Catch", bbmsy_zones = c(0.4, 0.8, 1))

####@> Probability of being in a Green Quadrant...
PGK_long(MSE)
PGK_short(MSE)
PGK_med(MSE)

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
