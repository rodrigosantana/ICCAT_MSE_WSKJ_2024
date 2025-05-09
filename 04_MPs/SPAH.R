source('03AA_MP_Internal_Functions.R')
HistList <- readRDS("03_Hists/HistList.rda")

SurplusProductionStockStatus <- function(x, Data,
                                 Data_Lag = 1, 
                                 Interval = 3,
                                 Initial_MP_Yr = 2026, 
                                 reps =  1, 
                                 tunepar = 1,
                                 mc = c(0.25, 0.25),
                                 useHCR=TRUE,
                                 BMSYTarg=1.3,
                                 BMSYLim=0.5,
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
  delta <- 1+(EstSS - BMSYTarg) * tunepar

  # yind <- seq_along(EstSS)
  # summary(lm(EstSS~yind))$coefficients[2, 1:2]
  
  prevER <- tail(Data@Cat[x,],1)/tail(EstB_BMSY,2)[1] 
 
  Bratio <- EstSS/BMSYTarg
  if (useHCR) {
    if (Bratio>=1) {
      delta <- delta
    } else if (EstSS>=BMSYLim) {
      a <- (delta1-delta2)/(BMSYTarg-BMSYLim)
      b <- deltaMin -a*BMSYLim
      delta <- a*Bratio+b
      # Bratio <- seq(BMSYLim, BMSYTarg, by=0.1)
      # plot(Bratio, a*Bratio+b, type='l', ylim=c(0,1))
    } else {
      delta <- delta2
    }
  }
  
  newER <- prevER * delta
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


# mc
# useHCR
# BMSYTarg
# BMSYLim
# delta1
# delta2 


SPA_Base <- SurplusProductionStockStatus
class(SPA_Base) <- 'MP'

SPA_Base_nomc <- SurplusProductionStockStatus
formals(SPA_Base_nomc)$mc <- NA
class(SPA_Base_nomc) <- 'MP'

SPA_NoHCR <- SurplusProductionStockStatus
formals(SPA_NoHCR)$useHCR <- FALSE
class(SPA_NoHCR) <- 'MP'

SPA_NoHCR_nomc <- SurplusProductionStockStatus
formals(SPA_NoHCR_nomc)$useHCR <- FALSE
formals(SPA_NoHCR_nomc)$mc <- NA
class(SPA_NoHCR_nomc) <- 'MP'

i <- 7
MSE <- Project(HistList[[i]], MPs=c('SPA_Base', 
                                    'SPA_Base_nomc', 
                                    'SPA_NoHCR',
                                    'SPA_NoHCR_nomc'
                                    )
               )
Pplot(MSE)






SPAH1 <- SurplusProductionStockStatus
class(SPAH1) <- 'MP'

SPAH2 <- SurplusProductionStockStatus
formals(SPAH2)$BMSYLim <- 1
class(SPAH2) <- 'MP'

SPAH3 <- SurplusProductionStockStatus
formals(SPAH3)$BMSYLim <- 0.8
class(SPAH3) <- 'MP'


# Test different configurations of function arguments


selOMs <- c(1,5,9)

MSEList <- list()
for (i in selOMs) {
  MSEList[[i]] <- Project(HistList[[i]], MPs=c('SPAH1', 'SPAH2', 'SPAH3'))
  Pplot(MSEList[[i]])
}




pyears <- 2025:2054 
cbind(1:30, pyears)

pyears[15]

seq(2026, by = 3, length.out = 30)


Pplot(MSE1)
Pplot(MSE5)
Pplot(MSE9)

MSE5@SB_SBMSY[,1,15] |> which.min()
x <- 56
mm <- 1
om <- 5
Data <- SWOMSE::Trim_Data(MSE5@PPD[[mm]], 2034)
