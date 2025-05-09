source('03AA_MP_Internal_Functions.R')


SurplusProductionStockStatus <- function(x, Data,
                                 Data_Lag = 1, 
                                 Interval = 3,
                                 Initial_MP_Yr = 2026, 
                                 reps =  1, 
                                 tunepar = 0.5,
                                 mc = c(0.25, 0.25),
                                 n=3,
                                 BMSYTarg=1.2,
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
  EstSS <-  tail(EstB_BMSY,2) |> mean()

  delta <- 1+((EstSS - BMSYTarg)/n* tunepar)

  # yind <- seq_along(EstSS)
  # summary(lm(EstSS~yind))$coefficients[2, 1:2]
  
  TAC <- delta  * Data@MPrec[x]
  ## Maximum allowed change in TAC
  ## Rec@TAC <- TAC
  Rec@TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
  Rec
}


SPAH1 <- SurplusProductionStockStatus
class(SPAH1) <- 'MP'

SPAH2 <- SurplusProductionStockStatus
formals(SPAH2)$n <- 4
class(SPAH2) <- 'MP'

SPAH3 <- SurplusProductionStockStatus
formals(SPAH3)$tunepar <- 0.3
class(SPAH3) <- 'MP'


# Test different configurations of function arguments
HistList <- readRDS("03_Hists/HistList.rda")

Hist <- HistList[[1]]

MSE <- Project(Hist, MPs=c('SPAH1', 'SPAH2', 'SPAH3'))
Pplot(MSE)

# Results:
# - doSmoother = TRUE - worse performance
# - incHCR = TRUE worse performance
