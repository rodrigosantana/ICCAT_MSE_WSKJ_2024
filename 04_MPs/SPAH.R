source('03AA_MP_Internal_Functions.R')


SurplusProductionStockStatus <- function(x, Data,
                                 Data_Lag = 1, 
                                 Interval = 3,
                                 Initial_MP_Yr = 2026, 
                                 reps =  1, 
                                 tunepar = .7,
                                 mc = c(0.3, 0.3),
                                 useHCR=TRUE,
                                 BMSYTarg=1.3,
                                 BMSYLim=0.6,
                                 delta1=1,
                                 delta2=0.5,
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
  
  # Estimated Stock Status from SP 
  SPrun <- SAMtool::SP(x, Data,
                  prior = list(r = c(0.416, 0.148),
                               MSY = c(30000, 0.2)),
                  start = list(dep = 0.98, n = 1))
  
  # Slope to Target Rule 
  EstB_BMSY <- SPrun@B_BMSY
  EstSS <-  tail(EstB_BMSY,1) 
  
  # ##############################################################################
  # 
  # par(mfrow=c(2,2))
  # B <- c(rowSums(HistList[[om]]@TSdata$Biomass[x,,]), 
  #   MSEList[[om]]@B[x,mm,])
  # B_BMSY <- B/HistList[[om]]@Ref$ByYear$BMSY
  # 
  # 
  # plot(Data@Year, B_BMSY[seq_along(Data@Year)], type='l', ylim=c(0,3))
  # lines(Data@Year, EstB_BMSY[seq_along(Data@Year)], col='red')
  # 
  # plot(Data@Year, Data@Cat[x,], type='l', ylim=c(0,max(Data@Cat[x,])))
  # 
  # plot(Data@Year, Data@Ind[x,], type='l', ylim=c(0,2))
  # 
  # ##############################################################################
  
  # delta <- 1+((EstSS - BMSYTarg)/n* tunepar)
  # delta <- 1+(EstSS - BMSYTarg) * tunepar
  
  prevER <- tail(Data@Cat[x,],1)/tail(EstB_BMSY,2)[1] 
 
  Bratio <- EstSS/BMSYTarg
  delta <- Bratio
  if (useHCR) {
    if (EstSS>=BMSYLim & EstSS<BMSYTarg) {
      a <- (delta1-delta2)/(BMSYTarg-BMSYLim)
      b <- delta2 -a*BMSYLim
      delta <- a*Bratio+b
      # Bratio <- seq(BMSYLim, BMSYTarg, by=0.1)
      # plot(Bratio, a*Bratio+b, type='l', ylim=c(0,1))
    } else if (EstSS<BMSYLim){
      delta <- delta2
    }
  }
  
  newER <- prevER * delta * tunepar
  TAC <- newER*tail(EstB_BMSY,1)[1]
  
  ## Maximum allowed change in TAC
  ## Rec@TAC <- TAC
  if (Bratio>=1) {
    TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
  }
  if(!useHCR) {
    TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
  }
  
  Rec@TAC <- TAC
  Rec
}


# mc - performs better without but may exceed 25%
# useHCR - performs better with
# BMSYTarg
# BMSYLim
# delta1
# delta2 

SPA_Base <- SurplusProductionStockStatus
class(SPA_Base) <- 'MP'

mse <- Project(HistList[[5]], 'SPA_Base')
Pplot(mse)

x=1
Data <- SWOMSE::Trim_Data(mse@PPD[[1]], 2031)


# HistList <- readRDS("03_Hists/HistList.rda")
# 
# 
# 
# 
# i <- 5
# MSE <- Project(HistList[[i]], MPs=c('SPA_Base', 
#                                     'SPA_Base_nomc', 
#                                     'SPA_TargHigh',
#                                     'SPA_TargLow',
#                                     'SPA_TargLimHigh',
#                                     'SPA_TargLimLow'
#                                     ))
# Pplot(MSE)
# 
# SPA_LowTune <- SurplusProductionStockStatus
# formals(SPA_LowTune)$tunepar <- 0.2
# class(SPA_LowTune) <- 'MP'
# 
# MSE <- Project(HistList[[i]], MPs=c('SPA_Base', 
#                                     'SPA_LowTune')
# )
# Pplot(MSE)
# 
# 
# 

