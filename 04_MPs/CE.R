source('03AA_MP_Internal_Functions.R')



# Constant Exploitation  MPs

# arguments:
# - tunepar: tuning value for reference index
# - modifier: base modifier for reference index
# - doSmoother: logical; include smoother?
# - refYear: reference year for historical index
# - mc: max change in TAC between management cycles: c(MaxDecrease, MaxIncrease); NA to ignore
# - yrs: smoothing years: 
#        first value - nyears before refYear (including refYear)
#        second value - nyears before current year (including current year)
# - incHCR: include harvest control rule; logical

ConstantExploitation <- function(x, Data,
                                 Data_Lag = 1, 
                                 Interval = 3,
                                 Initial_MP_Yr = 2026, 
                                 reps =  1, 
                                 tunepar = 1,
                                 modifier = 0.4, 
                                 doSmoother=FALSE,
                                 refYear=2017,
                                 mc = c(0.25, 0.25),
                                 yrs = c(2, 2),
                                 incHCR=FALSE,
                                 ...) {
  Rec <- new("Rec")
  
  # Check if TAC needs to be updated
  if (SameTAC(Initial_MP_Yr, Interval, Data)) {
    Rec@TAC <- Data@MPrec[x]
    Rec <- FixedTAC(Rec, Data) # use actual catches if they are available
    return(Rec)
  }
  
  # Lag Data
  Data <- Lag_Data(Data, Data_Lag)
  Data@Year <- Data@Year[1:(length(Data@Year)-Data_Lag)]
  
  ## smooth combined index
  if (doSmoother) {
    index <- smoothed_index <- Data@Ind[x,]
    smoothed <- stats::smooth(index[!is.na(index)])
    smoothed_index[!is.na(smoothed_index)] <- smoothed
    Data@Ind[x,] <- smoothed_index
  }
 
  ## Calculate Historical Relative Exploitation Rate
  yr.ind <- which(Data@Year == refYear)
  hist.yrs <- (yr.ind - yrs[1] + 1):yr.ind
  
  histER <- mean(Data@Cat[x, hist.yrs]) / mean(Data@Ind[x, hist.yrs]) 
  
  ## Calculate Current Relative Exploitation Rate
  current_yr <- length(Data@Ind[x,])
  recent_yrs <- (current_yr - yrs[2]+1):current_yr
  curER <- mean(Data@Cat[x, recent_yrs]) / mean(Data@Ind[x, recent_yrs]) 
  
  ## Control Rule
  histInd <- mean(Data@Ind[x, hist.yrs]) * modifier * tunepar
  curInd <- mean(Data@Ind[x, recent_yrs], na.rm = TRUE)
  ind_ratio <- curInd/histInd
  
  targER <- histER * ind_ratio
  
  if (incHCR) {
    if (ind_ratio >= 0.8) {
      targER <- histER
    } else if (ind_ratio > 0.5) {
      targER <- histER * (-1.4 + 3 * ind_ratio)
    } else {
      targER <- 0.1 * histER
    }
  }

  ## Exploitation Rate Ratio
  ER_ratio <- targER/curER
  TAC <- ER_ratio  * Data@MPrec[x]
  ## Maximum allowed change in TAC
  ## Rec@TAC <- TAC
  Rec@TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
  Rec
}


# Test different configurations of function arguments
# HistList <- readRDS("03_Hists/HistList.rda")
# 
# Hist <- HistList[[1]]
# 
# MSE <- Project(Hist, MPs=c('CE1', 'CE2', 'CE3'))
# Pplot(MSE)

# Results:
# - doSmoother = TRUE - worse performance
# - incHCR = TRUE worse performance

CE1 <- ConstantExploitation
class(CE1) <- 'MP'

CE2 <- ConstantExploitation
formals(CE2)$mc <- c(0.3,0.3)  # slightly better performance
class(CE2) <- 'MP'



