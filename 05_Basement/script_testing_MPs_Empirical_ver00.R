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
######@> Tunning Empirical Index-Based MPs...

######@> Changing tunning parameter...

#####@> Iratio MP...
Iratio_MOD <- function(x, Data, Interval = 3, Data_Lag = 2,
                       tunepar = 1, mc = 0.2, yrs = c(3, 5), ...) {
    Rec <- new("Rec")
    ## Does TAC need to be updated? (or set a fixed catch if before
    ## Initial_MP_Yr)
    if(SWOMSE::SameTAC(Initial_MP_Yr, Interval, Data)) {
        Rec@TAC <- Data@MPrec[x]
        Rec <- SWOMSE::FixedTAC(Rec, Data) ## use actual catches if they are
        ## available
        return(Rec)
    } else {
        ## Lag Data
        Data <- Lag_Data(Data, Data_Lag)
        Data@Year <- Data@Year[1:(length(Data@Year)-Data_Lag)]
        ## Smooth combined index
        index <- smoothed_index <- Data@Ind[x,]
        smoothed <- stats::smooth(index[!is.na(index)])
        smoothed_index[!is.na(smoothed_index)] <- smoothed
        Data@Ind[x,] <- smoothed_index
        ## Historical data...
        yr.ind <- which(Data@Year == max(Data@Year) - yrs[1])
        hist.yrs <- (yr.ind - yrs[2] + 1):yr.ind
        I.hist <- mean(Data@Ind[x, hist.yrs], na.rm = TRUE)
        ## Current data...
        curr.yrs <- (yr.ind + 1):(yr.ind + yrs[1])
        I.curr <- mean(Data@Ind[x, curr.yrs], na.rm = TRUE)
        ## Alpha...
        alpha <- I.curr / I.hist
        ## Catch data...
        Cat <- Data@Cat[x, length(Data@Cat[x, ])]
        Cc <- MSEtool::trlnorm(reps = 1, Cat, Data@CV_Cat[x, 1])
        ## Harvest control rule...
        if(alpha > 1.2) {
            alpha <- 1.2
        } else if(alpha < 0.8) {
            alpha <- 0.8
        } else {
            alpha <- alpha
        }
        ## TAC...
        TAC <- alpha * tunepar * Cc
        TAC <- MSEtool::TACfilter(TAC)
        ## Return TAC...
        Rec@TAC <- TAC
        Rec@TAC <- MaxChange(TAC, Data@MPrec[x], mc)
        return(Rec)
    }
}
class(Iratio_MOD) <- "MP"

#####@> Tunning SP_01...

####@> tunepar = 1.03
Iratio_MODA <- Iratio_MOD
formals(Iratio_MODA)$tunepar <- 1.06
class(Iratio_MODA) <- "MP"

####@> tunepar = 1.025
Iratio_MODB <- Iratio_MOD
formals(Iratio_MODB)$tunepar <- 1.05
class(Iratio_MODB) <- "MP"

####@> tunepar = 1.02
Iratio_MODC <- Iratio_MOD
formals(Iratio_MODC)$tunepar <- 1.04
class(Iratio_MODC) <- "MP"

####@> tunepar = 1.015
Iratio_MODD <- Iratio_MOD
formals(Iratio_MODD)$tunepar <- 1.03
class(Iratio_MODD) <- "MP"

######@> Testing tuned MPs...
MSE <- Project(Hists[[6]], MPs = c("Iratio_MOD", "Iratio_MODA", "Iratio_MODB",
    "Iratio_MODC", "Iratio_MODD"))

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
#####@> CE (Wrap from SWOMSE)...

#####@> MP CE_01...
CE_01 <- function(x, Data, Data_Lag = 1, Interval = 3, tunepar = 1,
               mc = 0.25, yrs = c(5, 3), ...) {
    Rec <- new("Rec")
    ## Does TAC need to be updated? (or set a fixed catch if before
    ## Initial_MP_Yr)
    if (SWOMSE::SameTAC(Initial_MP_Yr, Interval, Data)) {
        Rec@TAC <- Data@MPrec[x]
        Rec <- SWOMSE::FixedTAC(Rec, Data) # use actual catches if they are available
        return(Rec)
    }
    ## Lag Data
    Data <- Lag_Data(Data, Data_Lag)
    ## smooth combined index
    index <- smoothed_index <- Data@Ind[x,]
    smoothed <- stats::smooth(index[!is.na(index)])
    smoothed_index[!is.na(smoothed_index)] <- smoothed
    Data@Ind[x,] <- smoothed_index
    ## Calculate Historical Relative Exploitation Rate
    yr.ind <- which(Data@Year == 2017)
    hist.yrs <- (yr.ind - yrs[1] + 1):yr.ind
    histER <- mean(Data@Cat[x, hist.yrs]) / mean(Data@Ind[x, hist.yrs])
    ## Calculate Current Relative Exploitation Rate
    current_yr <- length(Data@Ind[x,])
    recent_yrs <- (current_yr - yrs[2]+1):current_yr
    curER <- mean(Data@Cat[x, recent_yrs]) / mean(Data@Ind[x,
               recent_yrs])
    ## Control Rule
    histInd <- mean(Data@Ind[x, hist.yrs])
    curInd <- mean(Data@Ind[x, recent_yrs], na.rm = TRUE)
    ind_ratio <- curInd/histInd
    if (ind_ratio >= 0.8) {
        targER <- histER
    } else if (ind_ratio > 0.5) {
        targER <- histER * (-1.4 + 3 * ind_ratio)
    } else {
        targER <- 0.1 * histER
    }
    ## Exploitation Rate Ratio
    ER_ratio <- targER/curER
    TAC <- ER_ratio * tunepar * Data@MPrec[x]
    ## Maximum allowed change in TAC
    Rec@TAC <- SWOMSE::MaxChange(TAC, Data@MPrec[x], mc)
    ## Rec@TAC <- TAC
    Rec
}
class(CE_01) <- "MP"

#####@> Tunning CE_01...

####@> tunepar = 1.03
CE_01A <- CE_01
formals(CE_01A)$tunepar <- 1.03
class(CE_01A) <- "MP"

####@> tunepar = 1.025
CE_01B <- CE_01
formals(CE_01B)$tunepar <- 0.9
class(CE_01B) <- "MP"

####@> tunepar = 1.02
CE_01C <- CE_01
formals(CE_01C)$tunepar <- 0.8
class(CE_01C) <- "MP"

####@> tunepar = 1.015
CE_01D <- CE_01
formals(CE_01D)$tunepar <- 0.7
class(CE_01D) <- "MP"

######@> Testing tuned MPs...
MSE <- Project(Hists[[5]], MPs = c("CE_01", "CE_01A", "CE_01B",
    "CE_01C", "CE_01D"))

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
