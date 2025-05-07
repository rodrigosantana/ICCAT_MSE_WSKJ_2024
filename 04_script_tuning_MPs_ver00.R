########################################################################
## Description: Tuning parameters - Based on Tom Carruthers Github
## Code...
##
## Maintainer: DatenKraft - SCRS / ICCAT (TT Species Group)
## Author: Rodrigo Sant'Ana
## Created: seg set  2 09:49:19 2024 (-0300)
## Version: 0.0.1
##
## URL:
## Doc URL:
##
## Database info:
##
### Commentary:
## https://github.com/Blue-Matter/ClimateTest/
## https://github.com/Blue-Matter/Blue_Shark_MSE/
##
### Code:
########################################################################

########################################################################
######@> Loading R packages...

######@> Package list...
library(openMSE)
library(tidyverse)
library(snowfall)

########################################################################
######@> Setup R...

######@> Global environment variables...
## Sys.setenv(Initial_MP_Yr = 2021) ##> NOTE: Not working...
## Sys.getenv("Initial_MP_Yr")

####@> Old method...
## Initial_MP_Yr <- 2025

######@> Creating catch data.frame for the recent years...
Catchdf <- read("CatchDF.rda")

######@> Creating index data.frame for the recent years...
load("05_Results/tsIndex_ver02.RData")
Cpuedf <- tsIndex %>%
    filter(Year %in% 2021:2022,
           Fleet == "Inverse variance weighted average")

######@> Loading one historical scenario for testing...
Hist05 <- readRDS("03_Hists/OM005_IVInds_ver02.hist")

#####@> Extract Data for testing MPs...
Data <- Hist05@Data

#####@> Loading history objects...
Hists <- sapply(dir("03_Hists",
                    pattern = "[0-9][0-9][0-9]_IVInds_CORRECTED_ver03",
                    full.names = TRUE), readRDS)

####@> Setting names for Hists...
names(Hists) <- paste0("OM", sprintf("%03d", 1:9))

######@> Function to adjust symetric, asymetric or not adjust the TAC...
adjust_TAC <- function(TAC, Last_TAC, mc) {
    delta_TAC <- TAC / Last_TAC
    if(any(is.na(mc))) {
        return(TAC)
    }
    lower_bound <- 1 - ifelse(length(mc) == 2, mc[1], mc)
    upper_bound <- 1 + ifelse(length(mc) == 2, mc[2], mc)
    delta_TAC <- pmax(pmin(delta_TAC, upper_bound), lower_bound)
    return(Last_TAC * delta_TAC)
}

######@> Function to adjust TAC for model based MPs...
adjust_TAC2 <- function(TAC, Last_TAC, mc, B_rel) {
    delta_TAC <- TAC / Last_TAC
    if(any(is.na(mc))) {
        return(TAC)
    }
    lower_bound <- if (length(mc) == 2) 1 - mc[1] else 1 - mc
    upper_bound <- if (length(mc) == 2) 1 + mc[2] else 1 + mc
    if(B_rel > 1) {
        delta_TAC <- pmax(pmin(delta_TAC, upper_bound), lower_bound)
    } else {
        delta_TAC <- pmin(delta_TAC, upper_bound)
    }
    return(Last_TAC * delta_TAC)
}

######@> Defining number of CPUs available for parallel process...
total_cores <- parallel::detectCores()

########################################################################
######@> Management Procedures...

######@> Relative Abundance Index MPs...

#####@> Iratio with asymetrical TAC correction...
IR_01 <- function(x, Data, Interval = 3, Data_Lag = 2,
                  Initial_MP_Yr = 2025, reps =  1,
                  tunepar = 1.306547889462,
                  mc = c(0.25, 0.2),
                       yrs = c(3, 5), ...) {
    Rec <- new("Rec")
    data_year <- max(Data@Year)
    ## Does TAC need to be updated? (or set a fixed catch if before
    ## Initial_MP_Yr)
    if(data_year < max(Catchdf$Year)) {
        Rec <- SWOMSE::FixedTAC(Rec, Data)
        return(Rec)
    } else {
        ## Including new index values for the recent year...
        Data@Year <- seq(min(Data@Year), max(Cpuedf$Year), 1)
        Data@Cat <-
            cbind(Data@Cat,
                  rep(Catchdf$Catch[Catchdf$Year == 2021], 300),
                  rep(Catchdf$Catch[Catchdf$Year == 2022], 300))
        Data@Ind <-
            cbind(Data@Ind,
                  rep(Cpuedf$Obs[Cpuedf$Year == 2021], 300),
                  rep(Cpuedf$Obs[Cpuedf$Year == 2022], 300))
        ## Smooth combined index
        index <- smoothed_index <- Data@Ind[x,]
        smoothed <- stats::smooth(index[!is.na(index)])
        smoothed_index[!is.na(smoothed_index)] <- smoothed
        Data@Ind[x,] <- smoothed_index
        ## Lag Data
        Data <- Lag_Data(Data, Data_Lag)
        ## Data@Year <- Data@Year[1:(length(Data@Year)-Data_Lag)]
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
        ## TAC <- MSEtool::TACfilter(TAC)
        ## Return TAC...
        ## Rec@TAC <- TAC
        Rec@TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
        return(Rec)
    }
}
class(IR_01) <- "MP"

#####@> Iratio with symetrical TAC correction...
IR_02 <- function(x, Data, Interval = 3, Data_Lag = 2,
                  Initial_MP_Yr = 2025, reps =  1,
                  tunepar = 1.52, mc = c(0.2, 0.2),
                  yrs = c(3, 5), ...) {
    Rec <- new("Rec")
    data_year <- max(Data@Year)
    ## Does TAC need to be updated? (or set a fixed catch if before
    ## Initial_MP_Yr)
    if(data_year < max(Catchdf$Year)) {
        Rec <- SWOMSE::FixedTAC(Rec, Data)
        return(Rec)
    } else {
        ## Including new index values for the recent year...
        Data@Year <- seq(min(Data@Year), max(Cpuedf$Year), 1)
        Data@Cat <-
            cbind(Data@Cat,
                  rep(Catchdf$Catch[Catchdf$Year == 2021], 300),
                  rep(Catchdf$Catch[Catchdf$Year == 2022], 300))
        Data@Ind <-
            cbind(Data@Ind,
                  rep(Cpuedf$Obs[Cpuedf$Year == 2021], 300),
                  rep(Cpuedf$Obs[Cpuedf$Year == 2022], 300))
        ## Smooth combined index
        index <- smoothed_index <- Data@Ind[x,]
        smoothed <- stats::smooth(index[!is.na(index)])
        smoothed_index[!is.na(smoothed_index)] <- smoothed
        Data@Ind[x,] <- smoothed_index
        ## Lag Data
        Data <- Lag_Data(Data, Data_Lag)
        ## Data@Year <- Data@Year[1:(length(Data@Year)-Data_Lag)]
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
        ## TAC <- MSEtool::TACfilter(TAC)
        ## Return TAC...
        ## Rec@TAC <- TAC
        Rec@TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
        return(Rec)
    }
}
class(IR_02) <- "MP"

#####@> Iratio without TAC correction...
IR_03 <- function(x, Data, Interval = 3, Data_Lag = 2,
                  Initial_MP_Yr = 2025, reps =  1,
                  tunepar = 1.52, mc = NA,
                  yrs = c(3, 5), ...) {
    Rec <- new("Rec")
    data_year <- max(Data@Year)
    ## Does TAC need to be updated? (or set a fixed catch if before
    ## Initial_MP_Yr)
    if(data_year < max(Catchdf$Year)) {
        Rec <- SWOMSE::FixedTAC(Rec, Data)
        return(Rec)
    } else {
        ## Including new index values for the recent year...
        Data@Year <- seq(min(Data@Year), max(Cpuedf$Year), 1)
        Data@Cat <-
            cbind(Data@Cat,
                  rep(Catchdf$Catch[Catchdf$Year == 2021], 300),
                  rep(Catchdf$Catch[Catchdf$Year == 2022], 300))
        Data@Ind <-
            cbind(Data@Ind,
                  rep(Cpuedf$Obs[Cpuedf$Year == 2021], 300),
                  rep(Cpuedf$Obs[Cpuedf$Year == 2022], 300))
        ## Smooth combined index
        index <- smoothed_index <- Data@Ind[x,]
        smoothed <- stats::smooth(index[!is.na(index)])
        smoothed_index[!is.na(smoothed_index)] <- smoothed
        Data@Ind[x,] <- smoothed_index
        ## Lag Data
        Data <- Lag_Data(Data, Data_Lag)
        ## Data@Year <- Data@Year[1:(length(Data@Year)-Data_Lag)]
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
        ## TAC <- MSEtool::TACfilter(TAC)
        ## Return TAC...
        ## Rec@TAC <- TAC
        Rec@TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
        return(Rec)
    }
}
class(IR_03) <- "MP"

#####@> Constant Exploitation with Control Rule MP with asymetrical TAC
#####@> correction...
CE_01 <- function(x, Data, Data_Lag = 2, Interval = 3, tunepar = 2.1,
                  Initial_MP_Yr = 2025, mc = c(0.25, 0.2),
                  yrs = c(5, 3), ...) {
    Rec <- new("Rec")
    data_year <- max(Data@Year)
    ## Does TAC need to be updated? (or set a fixed catch if before
    ## Initial_MP_Yr)
    if(data_year < max(Catchdf$Year)) {
        Rec <- SWOMSE::FixedTAC(Rec, Data)
        return(Rec)
    } else {
        ## Including new index values for the recent year...
        Data@Year <- seq(min(Data@Year), max(Cpuedf$Year), 1)
        Data@Cat <-
            cbind(Data@Cat,
                  rep(Catchdf$Catch[Catchdf$Year == 2021], 300),
                  rep(Catchdf$Catch[Catchdf$Year == 2022], 300))
        Data@Ind <-
            cbind(Data@Ind,
                  rep(Cpuedf$Obs[Cpuedf$Year == 2021], 300),
                  rep(Cpuedf$Obs[Cpuedf$Year == 2022], 300))
        ## smooth combined index
        index <- smoothed_index <- Data@Ind[x,]
        smoothed <- stats::smooth(index[!is.na(index)])
        smoothed_index[!is.na(smoothed_index)] <- smoothed
        Data@Ind[x,] <- smoothed_index
        ## Lag Data
        Data <- Lag_Data(Data, Data_Lag)
        ## Calculate Historical Relative Exploitation Rate
        yr.ind <- which(Data@Year == 2017)
        hist.yrs <- (yr.ind - yrs[1] + 1):yr.ind
        histER <- mean(Data@Cat[x, hist.yrs]) / mean(Data@Ind[x, hist.yrs])
        ## Calculate Current Relative Exploitation Rate
        current_yr <- length(Data@Ind[x,])
        recent_yrs <- (current_yr - yrs[2]+1):current_yr
        curER <- mean(Data@Cat[x, recent_yrs]) /
            mean(Data@Ind[x, recent_yrs])
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
        ## Rec@TAC <- TAC
        Rec@TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
        return(Rec)
    }
}
class(CE_01) <- "MP"

#####@> Constant Exploitation with Control Rule MP with symetrical TAC
#####@> correction...
CE_02 <- function(x, Data, Data_Lag = 2, Interval = 3, tunepar = 2.1,
                  Initial_MP_Yr = 2025, mc = c(0.2, 0.2),
                  yrs = c(5, 3), ...) {
    Rec <- new("Rec")
    data_year <- max(Data@Year)
    ## Does TAC need to be updated? (or set a fixed catch if before
    ## Initial_MP_Yr)
    if(data_year < max(Catchdf$Year)) {
        Rec <- SWOMSE::FixedTAC(Rec, Data)
        return(Rec)
    } else {
        ## Including new index values for the recent year...
        Data@Year <- seq(min(Data@Year), max(Cpuedf$Year), 1)
        Data@Cat <-
            cbind(Data@Cat,
                  rep(Catchdf$Catch[Catchdf$Year == 2021], 300),
                  rep(Catchdf$Catch[Catchdf$Year == 2022], 300))
        Data@Ind <-
            cbind(Data@Ind,
                  rep(Cpuedf$Obs[Cpuedf$Year == 2021], 300),
                  rep(Cpuedf$Obs[Cpuedf$Year == 2022], 300))
        ## smooth combined index
        index <- smoothed_index <- Data@Ind[x,]
        smoothed <- stats::smooth(index[!is.na(index)])
        smoothed_index[!is.na(smoothed_index)] <- smoothed
        Data@Ind[x,] <- smoothed_index
        ## Lag Data
        Data <- Lag_Data(Data, Data_Lag)
        ## Calculate Historical Relative Exploitation Rate
        yr.ind <- which(Data@Year == 2017)
        hist.yrs <- (yr.ind - yrs[1] + 1):yr.ind
        histER <- mean(Data@Cat[x, hist.yrs]) / mean(Data@Ind[x, hist.yrs])
        ## Calculate Current Relative Exploitation Rate
        current_yr <- length(Data@Ind[x,])
        recent_yrs <- (current_yr - yrs[2]+1):current_yr
        curER <- mean(Data@Cat[x, recent_yrs]) /
            mean(Data@Ind[x, recent_yrs])
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
        ## Rec@TAC <- TAC
        Rec@TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
        return(Rec)
    }
}
class(CE_02) <- "MP"

#####@> Constant Exploitation with Control Rule MP without TAC
#####@> correction...
CE_03 <- function(x, Data, Data_Lag = 2, Interval = 3, tunepar = 2.1,
                  Initial_MP_Yr = 2025, mc = NA,
                  yrs = c(5, 3), ...) {
    Rec <- new("Rec")
    data_year <- max(Data@Year)
    ## Does TAC need to be updated? (or set a fixed catch if before
    ## Initial_MP_Yr)
    if(data_year < max(Catchdf$Year)) {
        Rec <- SWOMSE::FixedTAC(Rec, Data)
        return(Rec)
    } else {
        ## Including new index values for the recent year...
        Data@Year <- seq(min(Data@Year), max(Cpuedf$Year), 1)
        Data@Cat <-
            cbind(Data@Cat,
                  rep(Catchdf$Catch[Catchdf$Year == 2021], 300),
                  rep(Catchdf$Catch[Catchdf$Year == 2022], 300))
        Data@Ind <-
            cbind(Data@Ind,
                  rep(Cpuedf$Obs[Cpuedf$Year == 2021], 300),
                  rep(Cpuedf$Obs[Cpuedf$Year == 2022], 300))
        ## smooth combined index
        index <- smoothed_index <- Data@Ind[x,]
        smoothed <- stats::smooth(index[!is.na(index)])
        smoothed_index[!is.na(smoothed_index)] <- smoothed
        Data@Ind[x,] <- smoothed_index
        ## Lag Data
        Data <- Lag_Data(Data, Data_Lag)
        ## Calculate Historical Relative Exploitation Rate
        yr.ind <- which(Data@Year == 2017)
        hist.yrs <- (yr.ind - yrs[1] + 1):yr.ind
        histER <- mean(Data@Cat[x, hist.yrs]) / mean(Data@Ind[x, hist.yrs])
        ## Calculate Current Relative Exploitation Rate
        current_yr <- length(Data@Ind[x,])
        recent_yrs <- (current_yr - yrs[2]+1):current_yr
        curER <- mean(Data@Cat[x, recent_yrs]) /
            mean(Data@Ind[x, recent_yrs])
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
        ## Rec@TAC <- TAC
        Rec@TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
        return(Rec)
    }
}
class(CE_03) <- "MP"

#####@> Islope1 MP with asymetrical TAC correction...
IS_01 <- function(x, Data, Interval = 3, Data_Lag = 2,
                  Initial_MP_Yr = 2025, reps = 1, lambda = 0.4,
                  tunepar = 1.04, mc = c(0.25, 0.2), yrsmth = 5, ...) {
    Rec <- new("Rec")
    data_year <- max(Data@Year)
    ## Does TAC need to be updated? (or set a fixed catch if before
    ## Initial_MP_Yr)
    if(data_year < max(Catchdf$Year)) {
        Rec <- SWOMSE::FixedTAC(Rec, Data)
        return(Rec)
    } else {
        ## Including new index values for the recent year...
        Data@Year <- seq(min(Data@Year), max(Cpuedf$Year), 1)
        Data@Cat <-
            cbind(Data@Cat,
                  rep(Catchdf$Catch[Catchdf$Year == 2021], 300),
                  rep(Catchdf$Catch[Catchdf$Year == 2022], 300))
        Data@Ind <-
            cbind(Data@Ind,
                  rep(Cpuedf$Obs[Cpuedf$Year == 2021], 300),
                  rep(Cpuedf$Obs[Cpuedf$Year == 2022], 300))
        ## Smooth combined index
        index <- smoothed_index <- Data@Ind[x,]
        smoothed <- stats::smooth(index[!is.na(index)])
        smoothed_index[!is.na(smoothed_index)] <- smoothed
        Data@Ind[x,] <- smoothed_index
        ## Lag Data
        Data <- Lag_Data(Data, Data_Lag)
        ## Data@Year <- Data@Year[1:(length(Data@Year)-Data_Lag)]
        ## Applied Islope1 MP to lagged data
        ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)
        Years <- Data@Year[ind]
        ylast <- (Data@LHYear[1] - Data@Year[1]) + 1
        C_dat <- Data@Cat[x, ind]
        if (is.na(Data@MPrec[x]) || length(Data@Year) == ylast + 1) {
            TACstar <- (1 - xx) *
                MSEtool::trlnorm(reps, mean(C_dat, na.rm = TRUE),
                                 Data@CV_Cat/(yrsmth^0.5))
        } else {
            TACstar <- rep(Data@MPrec[x], reps)
        }
        I_hist <- Data@Ind[x, ind]
        yind <- 1:yrsmth
        slppar <- summary(lm(log(I_hist) ~ yind))$coefficients[2, 1:2]
        if (reps > 1) {
            Islp <- rnorm(reps, slppar[1], slppar[2])
        } else {
            Islp <- slppar[1]
        }
        TAC <- TACstar * (1 + lambda * Islp) * tunepar
        ## Return TAC...
        ## Rec@TAC <- TAC
        Rec@TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
        return(Rec)
    }
}
class(IS_01) <- "MP"

#####@> Islope1 MP with symetrical TAC correction...
IS_02 <- function(x, Data, Interval = 3, Data_Lag = 2,
                  Initial_MP_Yr = 2025, reps = 1, lambda = 0.4,
                  tunepar = 1.04, mc = c(0.2, 0.2), yrsmth = 5, ...) {
    Rec <- new("Rec")
    data_year <- max(Data@Year)
    ## Does TAC need to be updated? (or set a fixed catch if before
    ## Initial_MP_Yr)
    if(data_year < max(Catchdf$Year)) {
        Rec <- SWOMSE::FixedTAC(Rec, Data)
        return(Rec)
    } else {
        ## Including new index values for the recent year...
        Data@Year <- seq(min(Data@Year), max(Cpuedf$Year), 1)
        Data@Cat <-
            cbind(Data@Cat,
                  rep(Catchdf$Catch[Catchdf$Year == 2021], 300),
                  rep(Catchdf$Catch[Catchdf$Year == 2022], 300))
        Data@Ind <-
            cbind(Data@Ind,
                  rep(Cpuedf$Obs[Cpuedf$Year == 2021], 300),
                  rep(Cpuedf$Obs[Cpuedf$Year == 2022], 300))
        ## Smooth combined index
        index <- smoothed_index <- Data@Ind[x,]
        smoothed <- stats::smooth(index[!is.na(index)])
        smoothed_index[!is.na(smoothed_index)] <- smoothed
        Data@Ind[x,] <- smoothed_index
        ## Lag Data
        Data <- Lag_Data(Data, Data_Lag)
        ## Data@Year <- Data@Year[1:(length(Data@Year)-Data_Lag)]
        ## Applied Islope1 MP to lagged data
        ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)
        Years <- Data@Year[ind]
        ylast <- (Data@LHYear[1] - Data@Year[1]) + 1
        C_dat <- Data@Cat[x, ind]
        if (is.na(Data@MPrec[x]) || length(Data@Year) == ylast + 1) {
            TACstar <- (1 - xx) *
                MSEtool::trlnorm(reps, mean(C_dat, na.rm = TRUE),
                                 Data@CV_Cat/(yrsmth^0.5))
        } else {
            TACstar <- rep(Data@MPrec[x], reps)
        }
        I_hist <- Data@Ind[x, ind]
        yind <- 1:yrsmth
        slppar <- summary(lm(log(I_hist) ~ yind))$coefficients[2, 1:2]
        if (reps > 1) {
            Islp <- rnorm(reps, slppar[1], slppar[2])
        } else {
            Islp <- slppar[1]
        }
        TAC <- TACstar * (1 + lambda * Islp) * tunepar
        ## Return TAC...
        ## Rec@TAC <- TAC
        Rec@TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
        return(Rec)
    }
}
class(IS_02) <- "MP"

#####@> Islope1 MP without TAC correction...
IS_03 <- function(x, Data, Interval = 3, Data_Lag = 2,
                  Initial_MP_Yr = 2025, reps = 1, lambda = 0.4,
                  tunepar = 1.04, mc = NA, yrsmth = 5, ...) {
    Rec <- new("Rec")
    data_year <- max(Data@Year)
    ## Does TAC need to be updated? (or set a fixed catch if before
    ## Initial_MP_Yr)
    if(data_year < max(Catchdf$Year)) {
        Rec <- SWOMSE::FixedTAC(Rec, Data)
        return(Rec)
    } else {
        ## Including new index values for the recent year...
        Data@Year <- seq(min(Data@Year), max(Cpuedf$Year), 1)
        Data@Cat <-
            cbind(Data@Cat,
                  rep(Catchdf$Catch[Catchdf$Year == 2021], 300),
                  rep(Catchdf$Catch[Catchdf$Year == 2022], 300))
        Data@Ind <-
            cbind(Data@Ind,
                  rep(Cpuedf$Obs[Cpuedf$Year == 2021], 300),
                  rep(Cpuedf$Obs[Cpuedf$Year == 2022], 300))
        ## Smooth combined index
        index <- smoothed_index <- Data@Ind[x,]
        smoothed <- stats::smooth(index[!is.na(index)])
        smoothed_index[!is.na(smoothed_index)] <- smoothed
        Data@Ind[x,] <- smoothed_index
        ## Lag Data
        Data <- Lag_Data(Data, Data_Lag)
        ## Data@Year <- Data@Year[1:(length(Data@Year)-Data_Lag)]
        ## Applied Islope1 MP to lagged data
        ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)
        Years <- Data@Year[ind]
        ylast <- (Data@LHYear[1] - Data@Year[1]) + 1
        C_dat <- Data@Cat[x, ind]
        if (is.na(Data@MPrec[x]) || length(Data@Year) == ylast + 1) {
            TACstar <- (1 - xx) *
                MSEtool::trlnorm(reps, mean(C_dat, na.rm = TRUE),
                                 Data@CV_Cat/(yrsmth^0.5))
        } else {
            TACstar <- rep(Data@MPrec[x], reps)
        }
        I_hist <- Data@Ind[x, ind]
        yind <- 1:yrsmth
        slppar <- summary(lm(log(I_hist) ~ yind))$coefficients[2, 1:2]
        if (reps > 1) {
            Islp <- rnorm(reps, slppar[1], slppar[2])
        } else {
            Islp <- slppar[1]
        }
        TAC <- TACstar * (1 + lambda * Islp) * tunepar
        ## Return TAC...
        ## Rec@TAC <- TAC
        Rec@TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
        return(Rec)
    }
}
class(IS_03) <- "MP"

#####@> GBslope MP with asymetrical TAC correction...
GB_01 <- function(x, Data, Interval = 3, Data_Lag = 2,
                  Initial_MP_Yr = 2025, reps = 1, lambda = 0.4,
                  tunepar = 1.6, mc = c(0.2, 0.25), yrsmth = 5, ...) {
    Rec <- new("Rec")
    data_year <- max(Data@Year)
    ## Does TAC need to be updated? (or set a fixed catch if before
    ## Initial_MP_Yr)
    if(data_year < max(Catchdf$Year)) {
        Rec <- SWOMSE::FixedTAC(Rec, Data)
        return(Rec)
    } else {
        ## Including new index values for the recent year...
        Data@Year <- seq(min(Data@Year), max(Cpuedf$Year), 1)
        Data@Cat <-
            cbind(Data@Cat,
                  rep(Catchdf$Catch[Catchdf$Year == 2021], 300),
                  rep(Catchdf$Catch[Catchdf$Year == 2022], 300))
        Data@Ind <-
            cbind(Data@Ind,
                  rep(Cpuedf$Obs[Cpuedf$Year == 2021], 300),
                  rep(Cpuedf$Obs[Cpuedf$Year == 2022], 300))
        ## Smooth combined index
        index <- smoothed_index <- Data@Ind[x,]
        smoothed <- stats::smooth(index[!is.na(index)])
        smoothed_index[!is.na(smoothed_index)] <- smoothed
        Data@Ind[x,] <- smoothed_index
        ## Lag Data
        Data <- Lag_Data(Data, Data_Lag)
        ## Data@Year <- Data@Year[1:(length(Data@Year)-Data_Lag)]
        ## Applied Islope1 MP to lagged data
        Catrec <- Data@Cat[x, length(Data@Cat[x, ])]
        ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)
        I_hist <- Data@Ind[x, ind]
        yind <- 1:yrsmth
        slppar <- summary(lm(log(I_hist) ~ yind))$coefficients[2, 1:2]
        if(reps > 1) {
            Islp <- rnorm(reps, slppar[1], slppar[2])
        } else {
            Islp <- slppar[1]
        }
        MuC <- Data@Cat[x, length(Data@Cat[x, ])]
        Cc <- MSEtool::trlnorm(reps, MuC, Data@CV_Cat[x, 1])
        TAC <- Cc * (1 + lambda * Islp) * tunepar
        TAC[TAC > (1.2 * Catrec)] <- 1.2 * Catrec
        TAC[TAC < (0.8 * Catrec)] <- 0.8 * Catrec
        TAC <- MSEtool::TACfilter(TAC)
        ## Return TAC...
        ## Rec@TAC <- TAC
        Rec@TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
        return(Rec)
    }
}
class(GB_01) <- "MP"

#####@> GBslope MP with symetrical TAC correction...
GB_02 <- function(x, Data, Interval = 3, Data_Lag = 2,
                  Initial_MP_Yr = 2025, reps = 1, lambda = 0.4,
                  tunepar = 1.6, mc = c(0.2, 0.2), yrsmth = 5, ...) {
    Rec <- new("Rec")
    data_year <- max(Data@Year)
    ## Does TAC need to be updated? (or set a fixed catch if before
    ## Initial_MP_Yr)
    if(data_year < max(Catchdf$Year)) {
        Rec <- SWOMSE::FixedTAC(Rec, Data)
        return(Rec)
    } else {
        ## Including new index values for the recent year...
        Data@Year <- seq(min(Data@Year), max(Cpuedf$Year), 1)
        Data@Cat <-
            cbind(Data@Cat,
                  rep(Catchdf$Catch[Catchdf$Year == 2021], 300),
                  rep(Catchdf$Catch[Catchdf$Year == 2022], 300))
        Data@Ind <-
            cbind(Data@Ind,
                  rep(Cpuedf$Obs[Cpuedf$Year == 2021], 300),
                  rep(Cpuedf$Obs[Cpuedf$Year == 2022], 300))
        ## Smooth combined index
        index <- smoothed_index <- Data@Ind[x,]
        smoothed <- stats::smooth(index[!is.na(index)])
        smoothed_index[!is.na(smoothed_index)] <- smoothed
        Data@Ind[x,] <- smoothed_index
        ## Lag Data
        Data <- Lag_Data(Data, Data_Lag)
        ## Data@Year <- Data@Year[1:(length(Data@Year)-Data_Lag)]
        ## Applied Islope1 MP to lagged data
        Catrec <- Data@Cat[x, length(Data@Cat[x, ])]
        ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)
        I_hist <- Data@Ind[x, ind]
        yind <- 1:yrsmth
        slppar <- summary(lm(log(I_hist) ~ yind))$coefficients[2, 1:2]
        if(reps > 1) {
            Islp <- rnorm(reps, slppar[1], slppar[2])
        } else {
            Islp <- slppar[1]
        }
        MuC <- Data@Cat[x, length(Data@Cat[x, ])]
        Cc <- MSEtool::trlnorm(reps, MuC, Data@CV_Cat[x, 1])
        TAC <- Cc * (1 + lambda * Islp) * tunepar
        TAC[TAC > (1.2 * Catrec)] <- 1.2 * Catrec
        TAC[TAC < (0.8 * Catrec)] <- 0.8 * Catrec
        TAC <- MSEtool::TACfilter(TAC)
        ## Return TAC...
        ## Rec@TAC <- TAC
        Rec@TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
        return(Rec)
    }
}
class(GB_02) <- "MP"

#####@> GBslope MP without TAC correction...
GB_03 <- function(x, Data, Interval = 3, Data_Lag = 2,
                  Initial_MP_Yr = 2025, reps = 1, lambda = 0.4,
                  tunepar = 1.6, mc = NA, yrsmth = 5, ...) {
    Rec <- new("Rec")
    data_year <- max(Data@Year)
    ## Does TAC need to be updated? (or set a fixed catch if before
    ## Initial_MP_Yr)
    if(data_year < max(Catchdf$Year)) {
        Rec <- SWOMSE::FixedTAC(Rec, Data)
        return(Rec)
    } else {
        ## Including new index values for the recent year...
        Data@Year <- seq(min(Data@Year), max(Cpuedf$Year), 1)
        Data@Cat <-
            cbind(Data@Cat,
                  rep(Catchdf$Catch[Catchdf$Year == 2021], 300),
                  rep(Catchdf$Catch[Catchdf$Year == 2022], 300))
        Data@Ind <-
            cbind(Data@Ind,
                  rep(Cpuedf$Obs[Cpuedf$Year == 2021], 300),
                  rep(Cpuedf$Obs[Cpuedf$Year == 2022], 300))
        ## Smooth combined index
        index <- smoothed_index <- Data@Ind[x,]
        smoothed <- stats::smooth(index[!is.na(index)])
        smoothed_index[!is.na(smoothed_index)] <- smoothed
        Data@Ind[x,] <- smoothed_index
        ## Lag Data
        Data <- Lag_Data(Data, Data_Lag)
        ## Data@Year <- Data@Year[1:(length(Data@Year)-Data_Lag)]
        ## Applied Islope1 MP to lagged data
        Catrec <- Data@Cat[x, length(Data@Cat[x, ])]
        ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)
        I_hist <- Data@Ind[x, ind]
        yind <- 1:yrsmth
        slppar <- summary(lm(log(I_hist) ~ yind))$coefficients[2, 1:2]
        if(reps > 1) {
            Islp <- rnorm(reps, slppar[1], slppar[2])
        } else {
            Islp <- slppar[1]
        }
        MuC <- Data@Cat[x, length(Data@Cat[x, ])]
        Cc <- MSEtool::trlnorm(reps, MuC, Data@CV_Cat[x, 1])
        TAC <- Cc * (1 + lambda * Islp) * tunepar
        TAC[TAC > (1.2 * Catrec)] <- 1.2 * Catrec
        TAC[TAC < (0.8 * Catrec)] <- 0.8 * Catrec
        TAC <- MSEtool::TACfilter(TAC)
        ## Return TAC...
        ## Rec@TAC <- TAC
        Rec@TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
        return(Rec)
    }
}
class(GB_03) <- "MP"

######@> Model-based MPs...

#####@> A Surplus Production model with a Harvest Control Rule that
#####@> linearly reduces F (Schaeffer model - B 80-40) - harvest
#####@> control rule based on: [with TAC fixed first year]
#####@> https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf...
SP_01 <- function(x, Data, Interval = 3, Data_Lag = 2,
                  Initial_MP_Yr = 2025, tunepar = 1.2, mc = NA, ...) {
    Rec <- new('Rec')
    data_year <- max(Data@Year)
    ## Does TAC need to be updated? (or set a fixed catch if before
    ## Initial_MP_Yr)
    if(data_year < max(Catchdf$Year)) {
        Rec <- SWOMSE::FixedTAC(Rec, Data)
        return(Rec)
    } else {
        ## Including new index values for the recent year...
        Data@Year <- seq(min(Data@Year), max(Cpuedf$Year), 1)
        Data@Cat <-
            cbind(Data@Cat,
                  rep(Catchdf$Catch[Catchdf$Year == 2021], 300),
                  rep(Catchdf$Catch[Catchdf$Year == 2022], 300))
        Data@Ind <-
            cbind(Data@Ind,
                  rep(Cpuedf$Obs[Cpuedf$Year == 2021], 300),
                  rep(Cpuedf$Obs[Cpuedf$Year == 2022], 300))
        Data@CV_Ind <-
            cbind(Data@CV_Ind,
                  rep(0.2, 300),
                  rep(0.2, 300))
        ## Smooth combined index
        index <- smoothed_index <- Data@Ind[x,]
        smoothed <- stats::smooth(index[!is.na(index)])
        smoothed_index[!is.na(smoothed_index)] <- smoothed
        Data@Ind[x,] <- smoothed_index
        ## Lag Data
        Data <- Lag_Data(Data, Data_Lag)
        ## Data@Year <- Data@Year[1:(length(Data@Year)-Data_Lag)]
        ## apply SP assessment model
        Mod <- SAMtool::SP(x, Data, prior = list(r = c(0.416, 0.148),
                                                 MSY = c(30000, 0.2)),
                           start = list(dep = 0.98, n = 1))
        ## harvest control rule
        ## based on: https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
        Bthresh <- Mod@BMSY
        Blim <- 0.4 * Bthresh
        ## Ftar <- 0.8 * tunepar * Mod@FMSY
        ## Fmin <- 0.1 * tunepar * Mod@FMSY
        Ftar <- tunepar * 0.10
        Fmin <- 0.1 * Ftar
        Bcurr <- Mod@B[length(Mod@B)]
        if(Bcurr >= Bthresh) {
            Fmort <- Ftar
        } else if(Bcurr > Blim) {
            Fmort <- Ftar * (-0.367 + 1.167 * Bcurr/Bthresh)
        } else {
            Fmort <- Fmin
        }
        TAC <-  Fmort*Bcurr
        ## Return TAC...
        ## Rec@TAC <- TAC
        Rec@TAC <- adjust_TAC2(TAC, Data@MPrec[x], mc,
                               B_rel = Bcurr/Bthresh)
        return(Rec)
    }
}
class(SP_01) <- "MP"

#####@> A Surplus Production Space-State model with a Harvest Control
#####@> Rule that linearly reduces F (Schaeffer model B 80-40) -
#####@> harvest control rule based on:
#####@> https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf...
SP_02 <- function(x, Data, Interval = 3, Data_Lag = 2,
                  Initial_MP_Yr = 2025, tunepar = 1.2, mc = NA, ...) {
    Rec <- new('Rec')
    data_year <- max(Data@Year)
    ## Does TAC need to be updated? (or set a fixed catch if before
    ## Initial_MP_Yr)
    if(data_year < max(Catchdf$Year)) {
        Rec <- SWOMSE::FixedTAC(Rec, Data)
        return(Rec)
    } else {
        ## Including new index values for the recent year...
        Data@Year <- seq(min(Data@Year), max(Cpuedf$Year), 1)
        Data@Cat <-
            cbind(Data@Cat,
                  rep(Catchdf$Catch[Catchdf$Year == 2021], 300),
                  rep(Catchdf$Catch[Catchdf$Year == 2022], 300))
        Data@Ind <-
            cbind(Data@Ind,
                  rep(Cpuedf$Obs[Cpuedf$Year == 2021], 300),
                  rep(Cpuedf$Obs[Cpuedf$Year == 2022], 300))
        Data@CV_Ind <-
            cbind(Data@CV_Ind,
                  rep(0.2, 300),
                  rep(0.2, 300))
        ## Smooth combined index
        index <- smoothed_index <- Data@Ind[x,]
        smoothed <- stats::smooth(index[!is.na(index)])
        smoothed_index[!is.na(smoothed_index)] <- smoothed
        Data@Ind[x,] <- smoothed_index
        ## Lag Data
        Data <- Lag_Data(Data, Data_Lag)
        ## Data@Year <- Data@Year[1:(length(Data@Year)-Data_Lag)]
        ## apply SP assessment model
        Mod <- SAMtool::SP_SS(x, Data,
                              prior = list(r = c(0.416, 0.148),
                                           MSY = c(30000, 0.2)),
                              start = list(dep = 0.98, n = 1),
                              fix_sigma = FALSE, fix_tau = TRUE)
        ## harvest control rule
        ## based on: https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
        Bthresh <- Mod@BMSY
        Blim <- 0.4 * Bthresh
        ## Ftar <- 0.8 * tunepar * Mod@FMSY
        ## Fmin <- 0.1 * tunepar * Mod@FMSY
        Ftar <- tunepar * 0.10
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
        ## Return TAC...
        ## Rec@TAC <- TAC
        Rec@TAC <- adjust_TAC2(TAC, Data@MPrec[x], mc,
                               B_rel = Bcurr/Bthresh)
        return(Rec)
    }
}
class(SP_02) <- "MP"

#####@> A Surplus Production model with a Harvest Control Rule that
#####@> linearly reduces F (Schaeffer model B 100-40) - harvest control
#####@> rule based on:
#####@> https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
#####@> HCR ram Ftar = 1, Fmin = 0.1, Blim = 0.4, BMSY = 1...
SP_03 <- function(x, Data, Interval = 3, Data_Lag = 2,
                  Initial_MP_Yr = 2025, tunepar = 1.2, mc = NA, ...) {
    Rec <- new('Rec')
    data_year <- max(Data@Year)
    ## Does TAC need to be updated? (or set a fixed catch if before
    ## Initial_MP_Yr)
    if(data_year < max(Catchdf$Year)) {
        Rec <- SWOMSE::FixedTAC(Rec, Data)
        return(Rec)
    } else {
        ## Including new index values for the recent year...
        Data@Year <- seq(min(Data@Year), max(Cpuedf$Year), 1)
        Data@Cat <-
            cbind(Data@Cat,
                  rep(Catchdf$Catch[Catchdf$Year == 2021], 300),
                  rep(Catchdf$Catch[Catchdf$Year == 2022], 300))
        Data@Ind <-
            cbind(Data@Ind,
                  rep(Cpuedf$Obs[Cpuedf$Year == 2021], 300),
                  rep(Cpuedf$Obs[Cpuedf$Year == 2022], 300))
        Data@CV_Ind <-
            cbind(Data@CV_Ind,
                  rep(0.2, 300),
                  rep(0.2, 300))
        ## Smooth combined index
        index <- smoothed_index <- Data@Ind[x,]
        smoothed <- stats::smooth(index[!is.na(index)])
        smoothed_index[!is.na(smoothed_index)] <- smoothed
        Data@Ind[x,] <- smoothed_index
        ## Lag Data
        Data <- Lag_Data(Data, Data_Lag)
        ## Data@Year <- Data@Year[1:(length(Data@Year)-Data_Lag)]
        ## apply SP assessment model
        Mod <- SAMtool::SP(x, Data, prior = list(r = c(0.416, 0.148),
                                                 MSY = c(30000, 0.2)),
                           start = list(dep = 0.98, n = 1))
        ## harvest control rule
        ## based on: https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
        Bthresh <- Mod@BMSY
        Blim <- 0.4 * Bthresh
        ## Ftar <- 0.8 * tunepar * Mod@FMSY
        ## Fmin <- 0.1 * tunepar * Mod@FMSY
        Ftar <- tunepar * 0.10
        Fmin <- 0.1 * Ftar
        Bcurr <- Mod@B[length(Mod@B)]
        if(Bcurr >= Bthresh) {
            Fmort <- Ftar
        } else if(Bcurr > Blim) {
            Fmort <- Ftar * (-0.5 + 1.5 * Bcurr/Bthresh)
        } else {
            Fmort <- Fmin
        }
        TAC <-  Fmort*Bcurr
        ## Return TAC...
        ## Rec@TAC <- TAC
        Rec@TAC <- adjust_TAC2(TAC, Data@MPrec[x], mc,
                               B_rel = Bcurr/Bthresh)
        return(Rec)
    }
}
class(SP_03) <- "MP"

#####@> A Surplus Production Space-State model with a Harvest Control
#####@> Rule that linearly reduces F (Schaeffer model) - harvest control
#####@> rule based on:
#####@> https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
#####@> HCR ram Ftar = 1, Fmin = 0.1, Blim = 0.4, BMSY = 1...
SP_04 <- function(x, Data, Interval = 3, Data_Lag = 2,
                  Initial_MP_Yr = 2025, tunepar = 1.2, mc = NA, ...) {
    Rec <- new('Rec')
    data_year <- max(Data@Year)
    ## Does TAC need to be updated? (or set a fixed catch if before
    ## Initial_MP_Yr)
    if(data_year < max(Catchdf$Year)) {
        Rec <- SWOMSE::FixedTAC(Rec, Data)
        return(Rec)
    } else {
        ## Including new index values for the recent year...
        Data@Year <- seq(min(Data@Year), max(Cpuedf$Year), 1)
        Data@Cat <-
            cbind(Data@Cat,
                  rep(Catchdf$Catch[Catchdf$Year == 2021], 300),
                  rep(Catchdf$Catch[Catchdf$Year == 2022], 300))
        Data@Ind <-
            cbind(Data@Ind,
                  rep(Cpuedf$Obs[Cpuedf$Year == 2021], 300),
                  rep(Cpuedf$Obs[Cpuedf$Year == 2022], 300))
        Data@CV_Ind <-
            cbind(Data@CV_Ind,
                  rep(0.2, 300),
                  rep(0.2, 300))
        ## Smooth combined index
        index <- smoothed_index <- Data@Ind[x,]
        smoothed <- stats::smooth(index[!is.na(index)])
        smoothed_index[!is.na(smoothed_index)] <- smoothed
        Data@Ind[x,] <- smoothed_index
        ## Lag Data
        Data <- Lag_Data(Data, Data_Lag)
        ## Data@Year <- Data@Year[1:(length(Data@Year)-Data_Lag)]
        ## apply SP assessment model
        Mod <- SAMtool::SP_SS(x, Data,
                              prior = list(r = c(0.416, 0.148),
                                           MSY = c(30000, 0.2)),
                              start = list(dep = 0.98, n = 1),
                              fix_sigma = FALSE, fix_tau = TRUE)
        ## harvest control rule
        ## based on: https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
        Bthresh <- Mod@BMSY
        Blim <- 0.4 * Bthresh
        ## Ftar <- 0.8 * tunepar * Mod@FMSY
        ## Fmin <- 0.1 * tunepar * Mod@FMSY
        Ftar <- tunepar * 0.10
        Fmin <- 0.1 * Ftar
        Bcurr <- Mod@B[length(Mod@B)]
        if (Bcurr>=Bthresh) {
            Fmort <- Ftar
        } else if (Bcurr>Blim) {
            Fmort <- Ftar * (-0.5 + 1.5 * Bcurr/Bthresh)
        } else {
            Fmort <- Fmin
        }
        TAC <-  Fmort*Bcurr
        ## Return TAC...
        ## Rec@TAC <- TAC
        Rec@TAC <- adjust_TAC2(TAC, Data@MPrec[x], mc,
                               B_rel = Bcurr/Bthresh)
        return(Rec)
    }
}
class(SP_04) <- "MP"

########################################################################
######@> Defining functions for the optimization...

######@> Internal MSE running function for the tune_MP function (Tom
######@> Carruthers, 2024)
results <- list()
int_tune <- function(par, MP_parname, MP, Hist_list, minfunc, parallel) {
    assign("MPtest", get(MP))
    formals(MPtest)[[MP_parname]] <- par
    cat(paste0(MP_parname, " = ", round(par, 6), " \n"))
    class(MPtest) <- "MP"
    if (!parallel) {
        MSE_list <- lapply(Hist_list,
                           function(X) Project(X, MPs = "MPtest"))
    } else {
        sfExport("MPtest", "Hist_list", "minfunc")
        MSE_list <- sfLapply(Hist_list,
                             function(X) Project(X, MPs = "MPtest"))
    }
    PGKw <- minfunc(MSE_list)
    result <- list(par = par, PGKw = PGKw)
    results <<- append(results, list(result))
    return((PGKw - 0.6)^2)
}

######@> A generic function that uses optimize to tune a single MP
######@> parameter to minimize a user-specified function (e.g. squared
######@> distance from a mean yield, PGK = 70%, etc.)...
tune_MP <- function(Hist_list, MP, MP_parname, interval, minfunc,
                    tol = 1E-2, parallel = FALSE) {
    opt <- optimize(int_tune, interval = interval,
                    MP_parname = MP_parname, MP = MP,
                    Hist_list = Hist_list, minfunc = minfunc,
                    tol = tol, parallel = parallel)
    MPout <- get(MP)
    formals(MPout)[MP_parname] <- opt$minimum
    class(MPout) <- "MP"
    ## cat("Resultados de Otimização:\n")
    ## print(results)
    return(list(MPout = MPout, optimization_results = results))
    ## return(MPout)
}

######@> A function that calculates the squared difference between
######@> obtained and target mean catch...
minfunc01 <- function(MSE_list) {
    Bm <- sapply(MSE_list,
                 function(X) { X@SSB[,1, 19:21] /
                                   X@SSB_hist[, X@nyears]})
    Bm <- mean(Bm)
    cat(paste0("Bm = ", round(Bm, 2), "\n"))
    (Bm - 1)^2
    return(Bm)
}

######@> A function that calculates the squared difference between
######@> obtained and target PGK...
minfunc02 <- function(MSE_list) {
    PGKm <- sapply(MSE_list,
                  function(X){ mean(X@SB_SBMSY > 1 & X@F_FMSY < 1)})
    PGKw <- mean(PGKm)
    cat(paste0("PGKw = ", round(PGKw, 6), "\n"))
    (PGKw - 0.7)^2 # PGK = 0.7
    return(PGKw)
}

######@> A function that calculates the squared difference between
######@> obtained and target PGK in the first 6 years...
minfunc03 <- function(MSE_list) {
    PGKm <- sapply(MSE_list, function(X) {
        mean(X@SB_SBMSY[ , , 1:6] > 1 & X@F_FMSY[ , , 1:6] < 1)
    })
    PGKw <- mean(PGKm)
    cat(paste0("PGKw[1-6] = ", round(PGKw, 6), "\n"))
    return(PGKw)
}

## ######@> gpuR example...
## minfunc03_gpu <- function(MSE_list) {
##     SB_SBMSY_gpu <- gpuMatrix(X@SB_SBMSY[,,1:6], nrow(X@SB_SBMSY),
##                               ncol(X@SB_SBMSY))
##     F_FMSY_gpu <- gpuMatrix(X@F_FMSY[,,1:6], nrow(X@F_FMSY), ncol(X@F_FMSY))
##     PGKm_gpu <- colMeans((SB_SBMSY_gpu > 1) & (F_FMSY_gpu < 1))
##     PGKw <- mean(as.numeric(PGKm_gpu[]))  # Converte de volta para CPU
##     cat(paste0("PGKw[1-6] = ", round(PGKw, 6), "\n"))
##     return(PGKw)
## }

########################################################################
######@> Applying tune function...

######@> Projecting MSE...
mse <- Project(Hist05, MPs = "SP_01")

#####@> Applying tune profile for one OM...
teste01 <- tune_MP(Hists, MP = "CE_01", MP_parname = "tunepar",
                   interval = c(1.5, 2), minfunc = minfunc02,
                   tol = 1E-2, parallel = FALSE) ## 1.942351 1.939018

#####@> Parallel not working!!!
## num_cores <- 9
## sfInit(parallel = TRUE, cpus = num_cores)
## library(MSEtool)
## sfLibrary(MSEtool)
## teste02 <- tune_MP(Hist_list = Hists, MP = "CE_01",
##                    MP_parname = "tunepar",
##                    interval = c(0.8, 3), minfunc = minfunc03,
##                    tol = 1E-2, parallel = TRUE)
## sfStop()


######@> Applying tunings profile - For short period (1-6 years)...

#####@> Index-based MPs...
system.time(
    IR_01tune <- tune_MP(Hists, MP = "IR_01", MP_parname = "tunepar",
                         interval = c(0.8, 3), minfunc = minfunc03,
                         tol = 1E-3, parallel = FALSE))

system.time(
    CE_01tune <- tune_MP(Hists, MP = "CE_01", MP_parname = "tunepar",
                         interval = c(0.8, 3), minfunc = minfunc03,
                         tol = 1E-3, parallel = FALSE)
)

system.time(
    IS_01tune <- tune_MP(Hists, MP = "IS_01", MP_parname = "tunepar",
                         interval = c(0.8, 3), minfunc = minfunc03,
                         tol = 1E-3, parallel = FALSE)
)

#####@> Model-based MPs...
SP_01tune <- tune_MP(Hists, MP = "SP_01", MP_parname = "tunepar",
                     interval = c(0.8, 3), minfunc = minfunc03,
                     tol = 1E-3, parallel = FALSE)

SP_02tune <- tune_MP(Hists, MP = "SP_02", MP_parname = "tunepar",
                     interval = c(0.8, 3), minfunc = minfunc03,
                     tol = 1E-3, parallel = FALSE)

SP_03tune <- tune_MP(Hists, MP = "SP_03", MP_parname = "tunepar",
                     interval = c(0.8, 3), minfunc = minfunc03,
                     tol = 1E-3, parallel = FALSE)

SP_04tune <- tune_MP(Hists, MP = "SP_04", MP_parname = "tunepar",
                     interval = c(0.8, 3), minfunc = minfunc03,
                     tol = 1E-3, parallel = FALSE)

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
