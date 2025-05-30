
# Update to latest Github versions
source("000_load_openMSE_packages.R")

library(openMSE)


# Source pre-tuned MPs
suppressWarnings(rm(list=avail('MP'))) # remove MPs from global enviroment

# Load MPs into global environment
for (fl in list.files("04_MPs", full.names = TRUE))   source(fl)

getMPNames <- function() {
  ls(envir = .GlobalEnv)[vapply(ls(envir = .GlobalEnv),
                                MSEtool:::getclass, logical(1), 
                                classy = 'MP')]
}



# Load Hist files
HistList <- readRDS("03_Hists/HistList.rda")


# tuning optimization function
# optPGK60_4_10 <- function(MSE_list) {
#   PGKm <- sapply(MSE_list, function(X) {
#     # Years 4 - 10 - index 5:11 becuase first projection year is 2025 before MP is used
#     mean(X@SB_SBMSY[ , , 5:11] > 1 & X@F_FMSY[ , , 5:11] < 1)
#   })
#   PGKw <- mean(PGKm)
# 
#   ssq <- (PGKw-0.6)^2
#   cat(paste0("*************************\n"))
#   cat(paste0("PGKw[4-10] = ", round(PGKw, 4), "\n"))
#   cat(paste0("SSQ = ", round(ssq, 5), "\n"))
#   cat(paste0("*************************\n"))
#   ssq
# }

optPGK60_1_30 <- function(MSE_list) {
    PGKm <- sapply(MSE_list, function(X) {
        ## Years 1 - 30 - index 2:31 becuase first projection year is 2025 before MP is used
        mean(X@SB_SBMSY[ , , 2:30] > 1 & X@F_FMSY[ , , 2:30] < 1)
    })
    PGKw <- round(mean(PGKm),3)

    ssq <- (PGKw-0.6)^2
    cat(paste0("*************************\n"))
    cat(paste0("PGKw = ", PGKw, "\n"))
    cat(paste0("SSQ = ", round(ssq, 5), "\n"))
    cat(paste0("*************************\n\n\n"))
    ssq
}

# Custom Tuning Function
DoMPTune <- function(HistList,
                     MPName,
                     TuneInterval=c(0.05, 4),
                     TuneFunction=optPGK60_1_30,
                     Data_Lag=1,
                     ManagementInterval=3,
                     Initial_MP_Yr=2026,
                     tol=1E-2,
                     parallel=FALSE) {

  message('Tuning: ', MPName)
  
  tuneMP <- get(MPName)
  formals(tuneMP)$Data_Lag <- Data_Lag
  formals(tuneMP)$Interval <- ManagementInterval
  formals(tuneMP)$Initial_MP_Yr <- Initial_MP_Yr
  class(tuneMP) <- 'MP'

  assign(MPName, tuneMP, envir=.GlobalEnv)

  tunedMP <- MSEtool::tune_MP(HistList, MP=MPName, MP_parname='tunepar',
                              interval=TuneInterval,
                              minfunc=TuneFunction,
                              tol=tol,
                              parallel=parallel)

  dirName <- paste0('DataLag_', Data_Lag, '_Interval_', ManagementInterval)

  if (!dir.exists('TunedMPs'))
    dir.create('TunedMPs')

  if (!dir.exists(file.path('TunedMPs', dirName)))
    dir.create(file.path('TunedMPs', dirName))

  filename <- paste0(MPName, '1_30.mp')
  saveRDS(tunedMP, file.path('TunedMPs', dirName, filename))

}


# 
AllMPs <- getMPNames()

parallel <- TRUE

if (parallel) {
  setup()
  sfExport(list=list('Catchdf', 
                     'FixedTAC', 
                     'SameTAC', 
                     'adjust_TAC', 
                     'adjust_TAC2'))
}

# ----- CE MPs ---------
CE_MPs <- c('CE1', 'CE2', 'CE3')

for (mp in CE_MPs) {
  DoMPTune(HistList, mp,parallel = parallel)
}

# ----- IR MPs ---------
IR_MPs <- c("IR1", "IR2", 'IR3')
IR_MPs <- c( "IR2", 'IR3')



for (mp in IR_MPs) {
  DoMPTune(HistList, mp, c(0.6, 1.5), parallel = parallel)
}




# ----- IS MPs ---------
# IS_MPs <- c('IS_01', 'IS_02', 'IS_03')
# 
# for (mp in IS_MPs) {
#   DoMPTune(HistList, mp)
# }

# ----- SP MPs ---------
SP_MPs <- c('SP_01', 'SP_02', 'SP_03', 'SP_04')

for (mp in SP_MPs) {
  DoMPTune(HistList, mp)
}







  
