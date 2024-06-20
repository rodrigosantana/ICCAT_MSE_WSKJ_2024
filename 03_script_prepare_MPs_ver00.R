########################################################################
## Description: Management Procedures WSKJ MSE...
##
## Maintainer: DatenKraft - SCRS/ICCAT (TT Species Group)
## Author: Rodrigo Sant'Ana
## Created: seg jun 11 17:11:31 2024 (-0300)
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

########################################################################
######@> Setup R...

######@> Global environment variables...
## Sys.setenv(Initial_MP_Yr = 2021) ##> NOTE: Not working...
## Sys.getenv("Initial_MP_Yr")

####@> Old method...
Initial_MP_Yr <- 2021

########################################################################
######@> Management Procedures...

######@> MP Wrapper general function...
MP_wrapper <- function(MP, delay = 2, ...) {
    MP <- substitute(MP)
    MP_call <- as.call(c(MP, x = quote(x), Data = quote(Data),
                         reps = quote(reps), list(...)))
    force(delay)
    MP_body <- bquote({
        Data@Ind <- Data@AddInd[, 1, ]
        Data@CV_Ind <- Data@CV_AddInd[, 1, ]
        if(delay > 0) {
            delay_ind <- seq(ncol(Data@Cat) - delay + 1,
                             ncol(Data@Cat), 1)
            Data@Ind <- Data@Ind[, -delay_ind]
            Data@CV_Ind <- Data@CV_Ind[, -delay_ind]
            Data@Cat <- Data@Cat[, -delay_ind, drop = FALSE]
            Data@CV_Cat <- Data@CV_Cat[, -delay_ind, drop = FALSE]
            Data@Year <- Data@Year[-delay_ind]
        }
        Rec <- .(MP_call)
        return(Rec)
    })
    MP_out <- eval(call("function",
                        as.pairlist(alist(x = 1, Data = , reps = 1)),
                        MP_body))
    class(MP_out) <- "MP"
    return(MP_out)
}

######@> Constant catches MPs...
CC_40kt <- function(x, Data, reps = 1, ...) {
    Rec <- new("Rec")
    ## TAC will be set for data_year+1
    ## ie data_year will be 2020 in projection year 2021
    data_year <- max(Data@Year)
    if(data_year < max(Catchdf$Year)) {
        Rec <- SWOMSE::FixedTAC(Rec, Data)
        return(Rec)
    } else {
        Rec@TAC <- rep(4e4, reps)
        return(Rec)
    }
}
class(CC_40kt) <- "MP"

CC_35kt <- function(x, Data, reps = 1, ...) {
    Rec <- new("Rec")
    ## TAC will be set for data_year+1
    ## ie data_year will be 2020 in projection year 2021
    data_year <- max(Data@Year)
    if(data_year < max(Catchdf$Year)) {
        Rec <- SWOMSE::FixedTAC(Rec, Data)
        return(Rec)
    } else {
        Rec@TAC <- rep(35e3, reps)
        return(Rec)
    }
}
class(CC_35kt) <- "MP"

CC_30kt <- function(x, Data, reps = 1, ...) {
    Rec <- new("Rec")
    ## TAC will be set for data_year+1
    ## ie data_year will be 2020 in projection year 2021
    data_year <- max(Data@Year)
    if(data_year < max(Catchdf$Year)) {
        Rec <- SWOMSE::FixedTAC(Rec, Data)
        return(Rec)
    } else {
        Rec@TAC <- rep(3e4, reps)
        return(Rec)
    }
}
class(CC_30kt) <- "MP"

CC_25kt <- function(x, Data, reps = 1, ...) {
    Rec <- new("Rec")
    ## TAC will be set for data_year+1
    ## ie data_year will be 2020 in projection year 2021
    data_year <- max(Data@Year)
    if(data_year < max(Catchdf$Year)) {
        Rec <- SWOMSE::FixedTAC(Rec, Data)
        return(Rec)
    } else {
        Rec@TAC <- rep(25e3, reps)
        return(Rec)
    }
}
class(CC_25kt) <- "MP"

CC_20kt <- function(x, Data, reps = 1, ...) {
    Rec <- new("Rec")
    ## TAC will be set for data_year+1
    ## ie data_year will be 2020 in projection year 2021
    data_year <- max(Data@Year)
    if(data_year < max(Catchdf$Year)) {
        Rec <- SWOMSE::FixedTAC(Rec, Data)
        return(Rec)
    } else {
        Rec@TAC <- rep(2e4, reps)
        return(Rec)
    }
}
class(CC_20kt) <- "MP"

CC_15kt <- function(x, Data, reps = 1, ...) {
    Rec <- new("Rec")
    ## TAC will be set for data_year+1
    ## ie data_year will be 2020 in projection year 2021
    data_year <- max(Data@Year)
    if(data_year < max(Catchdf$Year)) {
        Rec <- SWOMSE::FixedTAC(Rec, Data)
        return(Rec)
    } else {
        Rec@TAC <- rep(1.5e4, reps)
        return(Rec)
    }
}
class(CC_15kt) <- "MP"

######@> Relative Abundance Index MPs...

#####@> Iratio MP...
Iratio_MOD <- function(x, Data, Interval = 3, Data_Lag = 1, ...) {
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
        ## Applied Iration MP to lagged data
        Rec <- Iratio(x, Data, reps = 1)
        return(Rec)
    }
}
class(Iratio_MOD) <- "MP"

#####@> Islope1 MP...
Islope1_MOD <- function(x, Data, Interval = 3, Data_Lag = 1, ...) {
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
        ## Applied Islope1 MP to lagged data
        Rec <- Islope1(x, Data, reps = 1)
        return(Rec)
    }
}
class(Islope1_MOD) <- "MP"

#####@> GBslope MP...
GBslope_MOD <- function(x, Data, Interval = 3, Data_Lag = 1, ...) {
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
        ## Applied Islope1 MP to lagged data
        Rec <- GB_slope(x, Data, reps = 1)
        return(Rec)
    }
}
class(GBslope_MOD) <- "MP"

#####@> Constant Exploitation with Control Rule MP...
CE <- function(x, Data, Data_Lag = 1, Interval = 3, tunepar = 1,
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
        targER <- histER * ( -1.4 + 3 * ind_ratio)
    } else {
        targER <- 0.1 * histER
    }
    ## Exploitation Rate Ratio
    ER_ratio <- targER/curER
    TAC <- ER_ratio * tunepar * Data@MPrec[x]
    ## Maximum allowed change in TAC
    Rec@TAC <- SWOMSE::MaxChange(TAC, Data@MPrec[x], mc)
    Rec
}
class(CE) <- "MP"

#####@> Stastical catch-at-age models HCRs...
SCA_00 <- make_MP(SCA,
    HCR_ramp,
    OCP_type = "SSB_SSBMSY",
    Ftarget_type = "FMSY",
    LOCP = 0.4,
    TOCP = 1,
    relF_min = 0.1,
    relF_max = 1)

######@> A Statistical Catch-At-Age model with a Harvest Control Rule that
######@> linearly reduces F - harvest control rule
######@> based on:
######@> https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf...
SCA_01 <- function(x, Data, Interval = 3, tunepar = 1.2, ...) {
    Rec <- new('Rec')
    ## apply SP assessment model
    Mod <- SAMtool::SCA(x, Data, SR = "BH", vulneability = "dome",
        CAA_dist = "lognormal")
    ## harvest control rule
    ## based on: https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
    Bthresh <- Mod@BMSY
    Blim <- 0.4 * Bthresh
    Ftar <- tunepar * 0.1
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
    Rec
}
class(SCA_01) <- 'MP'

######@> A Statistical Catch-At-Age model with a Harvest Control Rule that
######@> linearly reduces F - harvest control rule
######@> based on:
######@> https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf...
SCA_02 <- function(x, Data, Interval = 3, tunepar = 1.2, ...) {
    Rec <- new('Rec')
    ## apply SP assessment model
    Mod <- SAMtool::SCA(x, Data, SR = "BH", vulneability = "dome",
        CAA_dist = "lognormal")
    ## harvest control rule
    ## based on: https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
    Bthresh <- Mod@BMSY
    Blim <- 0.4 * Bthresh
    Ftar <- tunepar * 0.1
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
    Rec@TAC <- TAC
    Rec
}
class(SCA_02) <- 'MP'

######@> A Surplus Production model with a Harvest Control Rule that
######@> linearly reduces F (Schaeffer model) - harvest control rule
######@> based on: [with TAC fixed first year]
######@> https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf...
SP_01 <- function(x, Data, Data_Lag = 1, Interval = 3, tunepar = 1.04, ...) {
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
    Rec@TAC <- TAC
    Rec
}
class(SP_01) <- 'MP'

######@> A Surplus Production Space-State model with a Harvest Control
######Rule that linearly reduces F (Schaeffer model) - harvest control
######rule based on:
######@> https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf...
SP_02 <- function(x, Data, Data_Lag = 1, Interval = 3, tunepar = 1.22, ...) {
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
    Rec@TAC <- TAC
    Rec
}
class(SP_02) <- 'MP'


######@> A Surplus Production model with a Harvest Control Rule that
######@> linearly reduces F (Schaeffer model) - harvest control rule
######@> based on:
######@> https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
######@> HCR ram Ftar = 1, Fmin = 0.1, Blim = 0.4, BMSY = 1...
SP_03 <- function(x, Data, Data_Lag = 1, Interval = 3, tunepar = 1.27, ...) {
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
    ## Ftar <- 1 * tunepar * Mod@FMSY
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
    Rec@TAC <- TAC
    Rec
}
class(SP_03) <- 'MP'

######@> A Surplus Production Space-State model with a Harvest Control
######@> Rule that linearly reduces F (Schaeffer model) - harvest control
######@> rule based on:
######@> https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
######@> HCR ram Ftar = 1, Fmin = 0.1, Blim = 0.4, BMSY = 1...
SP_04 <- function(x, Data, Data_Lag = 1, Interval = 3, tunepar = 1.22, ...) {
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
        MSY = c(30000, 0.2)),
        start = list(dep = 0.98, n = 1),
        fix_sigma = FALSE, fix_tau = TRUE)
    ## harvest control rule
    ## based on: https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
    Bthresh <- Mod@BMSY
    Blim <- 0.4 * Bthresh
    ## Ftar <- 1 * tunepar * Mod@FMSY
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
    Rec@TAC <- TAC
    Rec
}
class(SP_04) <- 'MP'

########################################################################
######@> Playing with MPs...

## ######@> Loading example data...
## Hist_005 <- readRDS("02_Hists/OM005_IVInds.hist")

## #####@> Data object that will be used in the first projection year...
## Data <- Hist_005@Data

## ######@> Applying model...
## Mod <- SAMtool::SP_SS(1, Data,
##     prior = list(r = c(0.416, 0.148),
##         MSY = c(30000, 0.2)),
##     start = list(dep = 0.98, n = 1),
##     fix_sigma = FALSE, fix_tau = TRUE)


## Mod2 <- SAMtool::SCA(1, Data, SR = "BH", vulneability = "dome",
##     CAA_dist = "lognormal")

## ######@> Defining Catchdf...
## Catchdf <- data.frame(Year = 2021:2022,
##     Catch = c(20048.21, 21377.24))

## ######@> Testing MPs...

## #####@> CC MP's...

## ####@> 40 kt..
## tempData <- Data
## CC_40kt(1, tempData) # 2021 TAC

## tempData@Year <- 1952:2021
## CC_40kt(1, tempData) # 2022 TAC

## tempData@Year <- 1952:2022
## CC_40kt(1, tempData) # 2023 TAC

## tempData@Year <- 1952:2023
## CC_40kt(1, tempData) # 2024 TAC

## tempData@Year <- 1952:2024
## CC_40kt(1, tempData) # 2025 TAC

## ####@> 35 kt..
## tempData <- Data
## CC_35kt(1, tempData) # 2021 TAC

## tempData@Year <- 1952:2021
## CC_35kt(1, tempData) # 2022 TAC

## tempData@Year <- 1952:2022
## CC_35kt(1, tempData) # 2023 TAC

## tempData@Year <- 1952:2023
## CC_35kt(1, tempData) # 2024 TAC

## tempData@Year <- 1952:2024
## CC_35kt(1, tempData) # 2025 TAC

## ####@> 30 kt..
## tempData <- Data
## CC_30kt(1, tempData) # 2021 TAC

## tempData@Year <- 1952:2021
## CC_30kt(1, tempData) # 2022 TAC

## tempData@Year <- 1952:2022
## CC_30kt(1, tempData) # 2023 TAC

## tempData@Year <- 1952:2023
## CC_30kt(1, tempData) # 2024 TAC

## tempData@Year <- 1952:2024
## CC_30kt(1, tempData) # 2025 TAC

## ####@> 25 kt..
## tempData <- Data
## CC_25kt(1, tempData) # 2021 TAC

## tempData@Year <- 1952:2021
## CC_25kt(1, tempData) # 2022 TAC

## tempData@Year <- 1952:2022
## CC_25kt(1, tempData) # 2023 TAC

## tempData@Year <- 1952:2023
## CC_25kt(1, tempData) # 2024 TAC

## tempData@Year <- 1952:2024
## CC_25kt(1, tempData) # 2025 TAC

## ####@> 20 kt..
## tempData <- Data
## CC_20kt(1, tempData) # 2021 TAC

## tempData@Year <- 1952:2021
## CC_20kt(1, tempData) # 2022 TAC

## tempData@Year <- 1952:2022
## CC_20kt(1, tempData) # 2023 TAC

## tempData@Year <- 1952:2023
## CC_20kt(1, tempData) # 2024 TAC

## tempData@Year <- 1952:2024
## CC_20kt(1, tempData) # 2025 TAC

## ####@> 15 kt..
## tempData <- Data
## CC_15kt(1, tempData) # 2021 TAC

## tempData@Year <- 1952:2021
## CC_15kt(1, tempData) # 2022 TAC

## tempData@Year <- 1952:2022
## CC_15kt(1, tempData) # 2023 TAC

## tempData@Year <- 1952:2023
## CC_15kt(1, tempData) # 2024 TAC

## tempData@Year <- 1952:2024
## CC_15kt(1, tempData) # 2025 TAC

## #####@> Empirical MP's...

## ####@> Iratio..
## tempData <- Data
## Iratio_MOD(1, tempData) # 2021 TAC

## tempData@Year <- 1952:2021
## Iratio_MOD(1, tempData) # 2022 TAC

## tempData@Year <- 1952:2022
## tempData@Ind <- cbind(tempData@Ind,
##     tempData@Ind[ , ncol(tempData@Ind)]) # add another index data point
## tempData@Cat <- cbind(tempData@Cat,
##     tempData@Cat[ , ncol(tempData@Cat)]) # add another catch data point
## Iratio_MOD(1, tempData) # 2023 TAC

## tempData@Year <- 1952:2023
## tempData@Ind <- cbind(tempData@Ind,
##     tempData@Ind[ , ncol(tempData@Ind)]) # add another index data point
## tempData@Cat <- cbind(tempData@Cat,
##     tempData@Cat[ , ncol(tempData@Cat)]) # add another catch data point
## Iratio_MOD(1, tempData) # 2024 TAC

## tempData@Year <- 1952:2024
## tempData@Ind <- cbind(tempData@Ind,
##     tempData@Ind[ , ncol(tempData@Ind)]) # add another index data point
## tempData@Cat <- cbind(tempData@Cat,
##     tempData@Cat[ , ncol(tempData@Cat)]) # add another catch data point
## Iratio_MOD(1, tempData) # 2025 TAC

## ####@> Islope..
## tempData <- Data
## Islope1_MOD(1, tempData) # 2021 TAC

## tempData@Year <- 1952:2021
## Islope1_MOD(1, tempData) # 2022 TAC

## tempData@Year <- 1952:2022
## tempData@Ind <- cbind(tempData@Ind,
##     tempData@Ind[ , ncol(tempData@Ind)]) # add another index data point
## tempData@Cat <- cbind(tempData@Cat,
##     tempData@Cat[ , ncol(tempData@Cat)]) # add another catch data point
## Islope1_MOD(1, tempData) # 2023 TAC

## tempData@Year <- 1952:2023
## tempData@Ind <- cbind(tempData@Ind,
##     tempData@Ind[ , ncol(tempData@Ind)]) # add another index data point
## tempData@Cat <- cbind(tempData@Cat,
##     tempData@Cat[ , ncol(tempData@Cat)]) # add another catch data point
## Islope1_MOD(1, tempData) # 2024 TAC

## tempData@Year <- 1952:2024
## tempData@Ind <- cbind(tempData@Ind,
##     tempData@Ind[ , ncol(tempData@Ind)]) # add another index data point
## tempData@Cat <- cbind(tempData@Cat,
##     tempData@Cat[ , ncol(tempData@Cat)]) # add another catch data point
## Islope1_MOD(1, tempData) # 2025 TAC

## ####@> Islope..
## tempData <- Data
## GBslope_MOD(1, tempData) # 2021 TAC

## tempData@Year <- 1952:2021
## GBslope_MOD(1, tempData) # 2022 TAC

## tempData@Year <- 1952:2022
## tempData@Ind <- cbind(tempData@Ind,
##     tempData@Ind[ , ncol(tempData@Ind)]) # add another index data point
## tempData@Cat <- cbind(tempData@Cat,
##     tempData@Cat[ , ncol(tempData@Cat)]) # add another catch data point
## GBslope_MOD(1, tempData) # 2023 TAC

## tempData@Year <- 1952:2023
## tempData@Ind <- cbind(tempData@Ind,
##     tempData@Ind[ , ncol(tempData@Ind)]) # add another index data point
## tempData@Cat <- cbind(tempData@Cat,
##     tempData@Cat[ , ncol(tempData@Cat)]) # add another catch data point
## GBslope_MOD(1, tempData) # 2024 TAC

## tempData@Year <- 1952:2024
## tempData@Ind <- cbind(tempData@Ind,
##     tempData@Ind[ , ncol(tempData@Ind)]) # add another index data point
## tempData@Cat <- cbind(tempData@Cat,
##     tempData@Cat[ , ncol(tempData@Cat)]) # add another catch data point
## GBslope_MOD(1, tempData) # 2025 TAC

## #####@> Model-based MP's...

## ####@> SP_01..
## tempData <- Data
## SP_01(1, tempData) # 2021 TAC

## tempData@Year <- 1952:2021
## SP_01(1, tempData) # 2022 TAC

## tempData@Year <- 1952:2022
## tempData@Ind <- cbind(tempData@Ind,
##     tempData@Ind[ , ncol(tempData@Ind)]) # add another index data point
## tempData@CV_Ind <- cbind(tempData@CV_Ind,
##     tempData@CV_Ind[ , ncol(tempData@CV_Ind)]) # add another CV data point
## tempData@Cat <- cbind(tempData@Cat,
##     tempData@Cat[ , ncol(tempData@Cat)]) # add another catch data point
## SP_01(1, tempData) # 2023 TAC

## tempData@Year <- 1952:2023
## tempData@Ind <- cbind(tempData@Ind,
##     tempData@Ind[ , ncol(tempData@Ind)]) # add another index data point
## tempData@CV_Ind <- cbind(tempData@CV_Ind,
##     tempData@CV_Ind[ , ncol(tempData@CV_Ind)]) # add another CV data point
## tempData@Cat <- cbind(tempData@Cat,
##     tempData@Cat[ , ncol(tempData@Cat)]) # add another catch data point
## SP_01(1, tempData) # 2024 TAC

## tempData@Year <- 1952:2024
## tempData@Ind <- cbind(tempData@Ind,
##     tempData@Ind[ , ncol(tempData@Ind)]) # add another index data point
## tempData@CV_Ind <- cbind(tempData@CV_Ind,
##     tempData@CV_Ind[ , ncol(tempData@CV_Ind)]) # add another CV data point
## tempData@Cat <- cbind(tempData@Cat,
##     tempData@Cat[ , ncol(tempData@Cat)]) # add another catch data point
## SP_01(1, tempData) # 2025 TAC

## ####@> SP_02..
## tempData <- Data
## SP_02(1, tempData) # 2021 TAC

## tempData@Year <- 1952:2021
## SP_02(1, tempData) # 2022 TAC

## tempData@Year <- 1952:2022
## tempData@Ind <- cbind(tempData@Ind,
##     tempData@Ind[ , ncol(tempData@Ind)]) # add another index data point
## tempData@CV_Ind <- cbind(tempData@CV_Ind,
##     tempData@CV_Ind[ , ncol(tempData@CV_Ind)]) # add another CV data point
## tempData@Cat <- cbind(tempData@Cat,
##     tempData@Cat[ , ncol(tempData@Cat)]) # add another catch data point
## SP_02(1, tempData) # 2023 TAC

## tempData@Year <- 1952:2023
## tempData@Ind <- cbind(tempData@Ind,
##     tempData@Ind[ , ncol(tempData@Ind)]) # add another index data point
## tempData@CV_Ind <- cbind(tempData@CV_Ind,
##     tempData@CV_Ind[ , ncol(tempData@CV_Ind)]) # add another CV data point
## tempData@Cat <- cbind(tempData@Cat,
##     tempData@Cat[ , ncol(tempData@Cat)]) # add another catch data point
## SP_02(1, tempData) # 2024 TAC

## tempData@Year <- 1952:2024
## tempData@Ind <- cbind(tempData@Ind,
##     tempData@Ind[ , ncol(tempData@Ind)]) # add another index data point
## tempData@CV_Ind <- cbind(tempData@CV_Ind,
##     tempData@CV_Ind[ , ncol(tempData@CV_Ind)]) # add another CV data point
## tempData@Cat <- cbind(tempData@Cat,
##     tempData@Cat[ , ncol(tempData@Cat)]) # add another catch data point
## SP_02(1, tempData) # 2025 TAC

## ####@> SP_03..
## tempData <- Data
## SP_03(1, tempData) # 2021 TAC

## tempData@Year <- 1952:2021
## SP_03(1, tempData) # 2022 TAC

## tempData@Year <- 1952:2022
## tempData@Ind <- cbind(tempData@Ind,
##     tempData@Ind[ , ncol(tempData@Ind)]) # add another index data point
## tempData@CV_Ind <- cbind(tempData@CV_Ind,
##     tempData@CV_Ind[ , ncol(tempData@CV_Ind)]) # add another CV data point
## tempData@Cat <- cbind(tempData@Cat,
##     tempData@Cat[ , ncol(tempData@Cat)]) # add another catch data point
## SP_03(1, tempData) # 2023 TAC

## tempData@Year <- 1952:2023
## tempData@Ind <- cbind(tempData@Ind,
##     tempData@Ind[ , ncol(tempData@Ind)]) # add another index data point
## tempData@CV_Ind <- cbind(tempData@CV_Ind,
##     tempData@CV_Ind[ , ncol(tempData@CV_Ind)]) # add another CV data point
## tempData@Cat <- cbind(tempData@Cat,
##     tempData@Cat[ , ncol(tempData@Cat)]) # add another catch data point
## SP_03(1, tempData) # 2024 TAC

## tempData@Year <- 1952:2024
## tempData@Ind <- cbind(tempData@Ind,
##     tempData@Ind[ , ncol(tempData@Ind)]) # add another index data point
## tempData@CV_Ind <- cbind(tempData@CV_Ind,
##     tempData@CV_Ind[ , ncol(tempData@CV_Ind)]) # add another CV data point
## tempData@Cat <- cbind(tempData@Cat,
##     tempData@Cat[ , ncol(tempData@Cat)]) # add another catch data point
## SP_03(1, tempData) # 2025 TAC

## ####@> SP_04..
## tempData <- Data
## SP_04(1, tempData) # 2021 TAC

## tempData@Year <- 1952:2021
## SP_04(1, tempData) # 2022 TAC

## tempData@Year <- 1952:2022
## tempData@Ind <- cbind(tempData@Ind,
##     tempData@Ind[ , ncol(tempData@Ind)]) # add another index data point
## tempData@CV_Ind <- cbind(tempData@CV_Ind,
##     tempData@CV_Ind[ , ncol(tempData@CV_Ind)]) # add another CV data point
## tempData@Cat <- cbind(tempData@Cat,
##     tempData@Cat[ , ncol(tempData@Cat)]) # add another catch data point
## SP_04(1, tempData) # 2023 TAC

## tempData@Year <- 1952:2023
## tempData@Ind <- cbind(tempData@Ind,
##     tempData@Ind[ , ncol(tempData@Ind)]) # add another index data point
## tempData@CV_Ind <- cbind(tempData@CV_Ind,
##     tempData@CV_Ind[ , ncol(tempData@CV_Ind)]) # add another CV data point
## tempData@Cat <- cbind(tempData@Cat,
##     tempData@Cat[ , ncol(tempData@Cat)]) # add another catch data point
## SP_04(1, tempData) # 2024 TAC

## tempData@Year <- 1952:2024
## tempData@Ind <- cbind(tempData@Ind,
##     tempData@Ind[ , ncol(tempData@Ind)]) # add another index data point
## tempData@CV_Ind <- cbind(tempData@CV_Ind,
##     tempData@CV_Ind[ , ncol(tempData@CV_Ind)]) # add another CV data point
## tempData@Cat <- cbind(tempData@Cat,
##     tempData@Cat[ , ncol(tempData@Cat)]) # add another catch data point
## SP_04(1, tempData) # 2025 TAC

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
