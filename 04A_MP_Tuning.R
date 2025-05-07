source("000_load_openMSE_packages.R")

library(openMSE)


# Source pre-tuned MPs
source("03A_Define_MPs.R")

# Load Hist files

HistList <- readRDS("03_Hists/HistList.rda")


# tuning optimization function
optPGK60_4_10 <- function(MSE_list) {
  PGKm <- sapply(MSE_list, function(X) {
    # Years 4 - 10 - index 5:11 becuase first projection year is 2025 before MP is used
    mean(X@SB_SBMSY[ , , 5:11] > 1 & X@F_FMSY[ , , 5:11] < 1)
  })
  PGKw <- mean(PGKm)

  ssq <- (PGKw-0.6)^2
  cat(paste0("*************************\n"))
  cat(paste0("PGKw[4-10] = ", round(PGKw, 4), "\n"))
  cat(paste0("SSQ = ", round(ssq, 5), "\n"))
  cat(paste0("*************************\n"))
  ssq
}
    
# Custom Tuning Function 
DoMPTune <- function(HistList, 
                     MPName, 
                     TuneInterval=c(0.2, 2),
                     TuneFunction=optPGK60_4_10,
                     Data_Lag=1, 
                     ManagementInterval=3, 
                     Initial_MP_Yr=2026,
                     parallel=FALSE) {
  
  tuneMP <- get(MPName)
  formals(tuneMP)$Data_Lag <- Data_Lag
  formals(tuneMP)$Interval <- ManagementInterval
  formals(tuneMP)$Initial_MP_Yr <- Initial_MP_Yr 
  class(tuneMP) <- 'MP'
  
  tunedMP <- MSEtool::tune_MP(HistList, "tuneMP", MP_parname='tunepar', TuneInterval,
                              TuneFunction, 
                              tol=1E-3, 
                              parallel=parallel)
  
  dirName <- paste0('DataLag_', Data_Lag, '_Interval_', ManagementInterval)
  
  if (!dir.exists(file.path('TunedMPs', dirName)))
    dir.create(file.path('TunedMPs', dirName))
  
  filename <- paste0(MPName, '.mp')
  saveRDS(tunedMP, file.path('TunedMPs', dirName, filename))
  
}


# not working for some reason
# setup() # setup parallel processing 
# 
# snowfall::sfExport(list=c('FixedTAC',
#                           'SameTAC',
#                           'adjust_TAC',
#                           'adjust_TAC2')
                   )
# ----- IR MPs ---------
IR_MPs <- c("IR_01", "IR_02", "IR_03")

for (mp in IR_MPs) {
  DoMPTune(HistList, mp)
}


# ----- CE MPs ---------
CE_MPs <- c('CE_01', 'CE_02', 'CE_03')

for (mp in CE_MPs) {
  DoMPTune(HistList, mp)
}

# ----- IS MPs ---------
IS_MPs <- c('IS_01', 'IS_02', 'IS_03')

for (mp in IS_MPs) {
  DoMPTune(HistList, mp)
}

# ----- SP MPs ---------
SP_MPs <- c('SP_01', 'SP_02', 'SP_03', 'SP_04')

for (mp in SP_MPs) {
  DoMPTune(HistList, mp)
}







  