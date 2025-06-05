source('03AA_MP_Internal_Functions.R')


SurplusProductionStockStatus <- function(x, Data,
                                 Data_Lag = 1, 
                                 Interval = 3,
                                 Initial_MP_Yr = 2026, 
                                 reps =  1, 
                                 tunepar = 0.6,
                                 mc = c(0.25, 0.25),
                                 useHCR=TRUE,
                                 BMSYTarg=1.3,
                                 BMSYLim=0.6,
                                 FMSYTarg=0.8,
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
  
  EstF_FMSY <- tail(SPrun@F_FMSY,1)
  EstB_BMSY <- SPrun@B_BMSY
  
  prevER <- tail(Data@Cat[x,],1)/tail(EstB_BMSY,2)[1] 
  
  newER <- prevER * (1+(FMSYTarg-EstF_FMSY)) * tunepar
  
  EstSS <-  tail(EstB_BMSY,1) 
  Bratio <- EstSS/BMSYTarg

  if (useHCR) {
    if (EstSS>=BMSYLim & EstSS<BMSYTarg) {
      a <- (delta1-delta2)/(BMSYTarg-BMSYLim)
      b <- delta2 -a*BMSYLim
      newER <- newER * a*Bratio+b
      # Bratio <- seq(BMSYLim, BMSYTarg, by=0.1)
      # plot(Bratio, a*Bratio+b, type='l', ylim=c(0,1))
    } else if (EstSS<BMSYLim){
      newER <- newER * delta2
    }
  }
  
  TAC <- newER*tail(EstB_BMSY,1)[1]
  
  ## Maximum allowed change in TAC
  ## Rec@TAC <- TAC
  # if (Bratio>=1) {
  TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
  # }
  # if(!useHCR) {
  #   TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
  # }
  
  Rec@TAC <- TAC
  Rec
}


# mc - performs better without but may exceed 25%
# useHCR - performs better with
# BMSYTarg
# BMSYLim
# delta1
# delta2 

SPA_1 <- SurplusProductionStockStatus
class(SPA_1) <- 'MP'

SPA_2 <- SPA_1
formals(SPA_2)$FMSYTarg <- 0.6
class(SPA_2) <- 'MP'

SPA_3 <- SPA_1
formals(SPA_3)$FMSYTarg <- 1
class(SPA_3) <- 'MP'


SPA_4 <- SPA_1
formals(SPA_4)$BMSYTarg <- 1.6
formals(SPA_4)$BMSYLim <- 0.8
class(SPA_4) <- 'MP'

SPA_5 <- SPA_1
formals(SPA_5)$mc <- c(0.4,0.4)
class(SPA_5) <- 'MP'

# SPA_nomc <- SurplusProductionStockStatus
# formals(SPA_nomc)$mc <- NA
# class(SPA_nomc) <- 'MP'
# 
# 
# most pessimist
# mse <- Project(HistList[[5]], c('SPA_Base', 'SPA_nomc'))
# Pplot(mse)
# 
# # most optimistic
# mse2 <- Project(HistList[[7]], c('SPA_Base', 'SPA_nomc'))
# Pplot(mse2)
# 

# which.max(mse@F_FMSY[,1,15])
# 
# x=11
# Data <- SWOMSE::Trim_Data(mse@PPD[[1]], 2040)
# 



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

# 
# # ##############################################################################
# # 
# par(mfrow=c(2,2))
# B <- c(rowSums(HistList[[om]]@TSdata$Biomass[x,,]),
#   mse@B[x,mm,])
# B_BMSY <- B/HistList[[om]]@Ref$ByYear$BMSY
# 
# 
# plot(Data@Year, B_BMSY[seq_along(Data@Year)], type='l', ylim=c(0,3))
# lines(Data@Year, EstB_BMSY[seq_along(Data@Year)], col='red')
# 
# F <- c(HistList[[om]]@TSdata$Find[x,],
#        mse@FM[x,mm,])
# F_FMSY <- F/HistList[[om]]@Ref$ByYear$FMSY
# 
# 
# plot(Data@Year, F_FMSY[seq_along(Data@Year)], type='l', ylim=c(0,1))
# lines(Data@Year, SPrun@F_FMSY[seq_along(Data@Year)], col='red')
# 
# 
# plot(Data@Year, Data@Cat[x,], type='l', ylim=c(0,max(Data@Cat[x,])))
# 
# plot(Data@Year, Data@Ind[x,], type='l', ylim=c(0,2))
# # 
# ##############################################################################

# delta <- 1+((EstSS - BMSYTarg)/n* tunepar)
# delta <- 1+(EstSS - BMSYTarg) * tunepar
