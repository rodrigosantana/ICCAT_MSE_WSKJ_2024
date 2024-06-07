########################################################################
## Description: Management Procedures WSKJ MSE...
##
## Maintainer: DatenKraft - SCRS/ICCAT (TT Species Group)
## Author: Rodrigo Sant'Ana
## Created: seg abr 17 17:11:31 2023 (-0300)
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
CC_40kt <- function(x, Data, reps) {
    Rec <- new("Rec")
    Rec@TAC <- rep(4e4, reps)
    return(Rec)
}
class(CC_40kt) <- "MP"

CC_35kt <- function(x, Data, reps) {
    Rec <- new("Rec")
    Rec@TAC <- rep(35e3, reps)
    return(Rec)
}
class(CC_35kt) <- "MP"

CC_30kt <- function(x, Data, reps) {
    Rec <- new("Rec")
    Rec@TAC <- rep(3e4, reps)
    return(Rec)
}
class(CC_30kt) <- "MP"

CC_25kt <- function(x, Data, reps) {
    Rec <- new("Rec")
    Rec@TAC <- rep(25e3, reps)
    return(Rec)
}
class(CC_25kt) <- "MP"

CC_20kt <- function(x, Data, reps) {
    Rec <- new("Rec")
    Rec@TAC <- rep(2e4, reps)
    return(Rec)
}
class(CC_20kt) <- "MP"

CC_15kt <- function(x, Data, reps) {
    Rec <- new("Rec")
    Rec@TAC <- rep(1.5e4, reps)
    return(Rec)
}
class(CC_15kt) <- "MP"

######@> Relative Abundance Index MPs...
Iratio_ <- MP_wrapper(Iratio)
Islope <- MP_wrapper(Islope1, xx = 0, yrsmth = 3, lambda = 0.6)
GBslope <- MP_wrapper(GB_slope, yrsmth = 3)

#####@> Stastical catch-at-age models HCRs...
SCA_100_40_SBMSY <- make_MP(SCA,
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
SCA_01 <- function(x, Data, Interval = 3, tunepar = 1, ...) {
    Rec <- new('Rec')
    ## apply SP assessment model
    Mod <- SAMtool::SCA(x, Data)
    ## harvest control rule
    ## based on: https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
    Bthresh <- Mod@SSBMSY
    Blim <- 0.4 * Bthresh
    Ftar <- 0.8 * tunepar * Mod@FMSY
    Fmin <- 0.1 * tunepar * Mod@FMSY
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
SCA_02 <- function(x, Data, Interval = 3, tunepar = 1, ...) {
    Rec <- new('Rec')
    ## apply SP assessment model
    Mod <- SAMtool::SCA(x, Data)
    ## harvest control rule
    ## based on: https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
    Bthresh <- Mod@BMSY
    Blim <- 0.4 * Bthresh
    Ftar <- 1 * tunepar * Mod@FMSY
    Fmin <- 0.1 * tunepar * Mod@FMSY
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

######@> Assessment models MPs with Specific HCR ramps...
SPSS_100_40_SBMSY <- make_MP(SP_SS,
                             HCR_ramp,
                             OCP_type = "SSB_SSBMSY",
                             Ftarget_type = "FMSY",
                             LOCP = 0.4,
                             TOCP = 1,
                             relF_min = 0.1,
                             relF_max = 1)

######@> Assessment models MPs with Specific HCR ramps...
SP_100_40_SBMSY <- make_MP(SP,
                           HCR_ramp,
                           OCP_type = "SSB_SSBMSY",
                           Ftarget_type = "FMSY",
                           LOCP = 0.4,
                           TOCP = 1,
                           relF_min = 0.1,
                           relF_max = 1)

######@> A Surplus Production model with a Harvest Control Rule that
######@> linearly reduces F (Schaeffer model) - harvest control rule
######@> based on: [with TAC fixed first year]
######@> https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf...
SP_01 <- function(x, Data, Data_Lag = 1, Interval = 3, tunepar = 1, ...) {
    Rec <- new('Rec')
    ## Does TAC need to be updated? (or set a fixed catch if before
    ## Initial_MP_Yr)
    if(SameTAC(Initial_MP_Yr, Interval, Data)) {
        Rec@TAC <- Data@MPrec[x]
        Rec <- FixedTAC(Rec, Data) ## use actual catches if they are
                                   ## available
        return(Rec)
    }
    ## Lag Data
    Data <- Lag_Data(Data, Data_Lag)
    ## apply SP assessment model
    Mod <- SAMtool::SP(x, Data, prior = list(r = c(0.5, 0.3)))
    ## harvest control rule
    ## based on: https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
    Bthresh <- Mod@BMSY
    Blim <- 0.4 * Bthresh
    Ftar <- 0.8 * tunepar * Mod@FMSY
    Fmin <- 0.1 * tunepar * Mod@FMSY
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

######@> A Surplus Production model with a Harvest Control Rule that
######@> linearly reduces F (Schaeffer model) - harvest control rule
######@> based on: [with TAC fixed first year - State-space]
######@> https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf...
SP_02 <- function(x, Data, Data_Lag = 1, Interval = 3, tunepar = 1, ...) {
    Rec <- new('Rec')
    ## Does TAC need to be updated? (or set a fixed catch if before
    ## Initial_MP_Yr)
    if(SameTAC(Initial_MP_Yr, Interval, Data)) {
        Rec@TAC <- Data@MPrec[x]
        Rec <- FixedTAC(Rec, Data) ## use actual catches if they are
                                   ## available
        return(Rec)
    }
    ## Lag Data
    Data <- Lag_Data(Data, Data_Lag)
    ## apply SP assessment model
    Mod <- SAMtool::SP_SS(x, Data, prior = list(r = c(0.5, 0.3)))
    ## harvest control rule
    ## based on: https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
    Bthresh <- Mod@BMSY
    Blim <- 0.4 * Bthresh
    Ftar <- 0.8 * tunepar * Mod@FMSY
    Fmin <- 0.1 * tunepar * Mod@FMSY
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
######@> based on: [without TAC fixed first year]
######@> https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf...
SP_03 <- function(x, Data, Data_Lag = 1, Interval = 3, tunepar = 1, ...) {
    Rec <- new('Rec')
    ## Lag Data
    Data <- Lag_Data(Data, Data_Lag)
    ## apply SP assessment model
    Mod <- SAMtool::SP(x, Data, prior = list(r = c(0.5, 0.3)))
    ## harvest control rule
    ## based on: https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
    Bthresh <- Mod@BMSY
    Blim <- 0.4 * Bthresh
    Ftar <- 0.8 * tunepar * Mod@FMSY
    Fmin <- 0.1 * tunepar * Mod@FMSY
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
class(SP_03) <- 'MP'

######@> A Surplus Production model with a Harvest Control Rule that
######@> linearly reduces F (Schaeffer model) - harvest control rule
######@> based on: [without TAC fixed first year - State-space]
######@> https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf...
SP_04 <- function(x, Data, Data_Lag = 1, Interval = 3, tunepar = 1, ...) {
    Rec <- new('Rec')
    ## Lag Data
    Data <- Lag_Data(Data, Data_Lag)
    ## apply SP assessment model
    Mod <- SAMtool::SP_SS(x, Data, prior = list(r = c(0.5, 0.3)))
    ## harvest control rule
    ## based on: https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
    Bthresh <- Mod@BMSY
    Blim <- 0.4 * Bthresh
    Ftar <- 0.8 * tunepar * Mod@FMSY
    Fmin <- 0.1 * tunepar * Mod@FMSY
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
class(SP_04) <- 'MP'

######@> A Surplus Production model with a Harvest Control Rule that
######@> linearly reduces F (Schaeffer model) - harvest control rule
######@> based on: [without TAC fixed first year and tunepar = 3]
######@> https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf...
SP_05 <- function(x, Data, Data_Lag = 1, Interval = 3, tunepar = 3, ...) {
    Rec <- new('Rec')
    ## Lag Data
    Data <- Lag_Data(Data, Data_Lag)
    ## apply SP assessment model
    Mod <- SAMtool::SP(x, Data, prior = list(r = c(0.5, 0.3)))
    ## harvest control rule
    ## based on: https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
    Bthresh <- Mod@BMSY
    Blim <- 0.4 * Bthresh
    Ftar <- 0.8 * tunepar * Mod@FMSY
    Fmin <- 0.1 * tunepar * Mod@FMSY
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
class(SP_05) <- 'MP'

######@> A Surplus Production model with a Harvest Control Rule that
######@> linearly reduces F (Schaeffer model) - harvest control rule
######@> based on: [without TAC fixed first year and tunepar = 3 - State-space]
######@> https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf...
SP_06 <- function(x, Data, Data_Lag = 1, Interval = 3, tunepar = 3, ...) {
    Rec <- new('Rec')
    ## Lag Data
    Data <- Lag_Data(Data, Data_Lag)
    ## apply SP assessment model
    Mod <- SAMtool::SP_SS(x, Data, prior = list(r = c(0.5, 0.3)))
    ## harvest control rule
    ## based on: https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
    Bthresh <- Mod@BMSY
    Blim <- 0.4 * Bthresh
    Ftar <- 0.8 * tunepar * Mod@FMSY
    Fmin <- 0.1 * tunepar * Mod@FMSY
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
class(SP_06) <- 'MP'

######@> A Surplus Production model with a Harvest Control Rule that
######@> linearly reduces F (Schaeffer model) - harvest control rule
######@> based on:
######@> https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
######@> HCR ram Ftar = 1, Fmin = 0.1, Blim = 0.4, BMSY = 1...
SP_07 <- function(x, Data, Data_Lag = 1, Interval = 3, tunepar = 3, ...) {
    Rec <- new('Rec')
    ## Lag Data
    Data <- Lag_Data(Data, Data_Lag)
    ## apply SP assessment model
    Mod <- SAMtool::SP(x, Data, prior = list(r = c(0.5, 0.3)))
    ## harvest control rule
    ## based on: https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
    Bthresh <- Mod@BMSY
    Blim <- 0.4 * Bthresh
    Ftar <- 1 * tunepar * Mod@FMSY
    Fmin <- 0.1 * tunepar * Mod@FMSY
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
class(SP_07) <- 'MP'

######@> A Surplus Production model with a Harvest Control Rule that
######@> linearly reduces F (Schaeffer model) - harvest control rule
######@> based on:
######@> https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
######@> HCR ram Ftar = 1, Fmin = 0.1, Blim = 0.4, BMSY = 1...
SP_08 <- function(x, Data, Data_Lag = 1, Interval = 3, tunepar = 3, ...) {
    Rec <- new('Rec')
    ## Lag Data
    Data <- Lag_Data(Data, Data_Lag)
    ## apply SP assessment model
    Mod <- SAMtool::SP_SS(x, Data, prior = list(r = c(0.5, 0.3)))
    ## harvest control rule
    ## based on: https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
    Bthresh <- Mod@BMSY
    Blim <- 0.4 * Bthresh
    Ftar <- 1 * tunepar * Mod@FMSY
    Fmin <- 0.1 * tunepar * Mod@FMSY
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
class(SP_08) <- 'MP'

######@> A Surplus Production model with a Harvest Control Rule that
######@> linearly reduces F (Fox model) - harvest control rule
######@> based on:
######@> https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
######@> HCR ram Ftar = 1, Fmin = 0.1, Blim = 0.4, BMSY = 1...
SP_09 <- function(x, Data, Data_Lag = 1, Interval = 3, tunepar = 3, ...) {
    Rec <- new('Rec')
    ## Lag Data
    Data <- Lag_Data(Data, Data_Lag)
    ## apply SP assessment model
    Mod <- SAMtool::SP_Fox(x, Data, prior = list(r = c(0.5, 0.3)))
    ## harvest control rule
    ## based on: https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
    Bthresh <- Mod@BMSY
    Blim <- 0.4 * Bthresh
    Ftar <- 1 * tunepar * Mod@FMSY
    Fmin <- 0.1 * tunepar * Mod@FMSY
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
class(SP_09) <- 'MP'

######@> A Surplus Production model with a Harvest Control Rule that
######@> linearly reduces F (Fox model) - harvest control rule
######@> based on:
######@> https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf...
######@> HCR ram Ftar = 0.8, Fmin = 0.1, Blim = 0.4, BMSY = 1...
SP_10 <- function(x, Data, Data_Lag = 1, Interval = 3, tunepar = 3, ...) {
    Rec <- new('Rec')
    ## Lag Data
    Data <- Lag_Data(Data, Data_Lag)
    ## apply SP assessment model
    Mod <- SAMtool::SP_Fox(x, Data, prior = list(r = c(0.5, 0.3)))
    ## harvest control rule
    ## based on: https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
    Bthresh <- Mod@BMSY
    Blim <- 0.4 * Bthresh
    Ftar <- 0.8 * tunepar * Mod@FMSY
    Fmin <- 0.1 * tunepar * Mod@FMSY
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
class(SP_10) <- 'MP'

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
