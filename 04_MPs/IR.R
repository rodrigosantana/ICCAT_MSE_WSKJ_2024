
source('03AA_MP_Internal_Functions.R')

# Index Ratio MPs 

# mc - maximum change in TAC between management cycles 

IndexRatio <- function(x,
                       Data, 
                       Data_Lag = 1, 
                       Interval = 3,
                       Initial_MP_Yr = 2026, 
                       reps =  1, 
                       tunepar = 0.9, 
                       modifier = 1,
                       mc = c(0.25, 0.25), 
                       yrs = c(3, 3), 
                       refYear=2020,
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
  
  # Set TAC
  ## Smooth combined index
  # index <- smoothed_index <- Data@Ind[x,]
  # smoothed <- stats::smooth(index[!is.na(index)])
  # smoothed_index[!is.na(smoothed_index)] <- smoothed
  # Data@Ind[x,] <- smoothed_index
  
  ## Historical data...
  yr.ind <- which(Data@Year == refYear)
  hist.yrs <- (yr.ind - yrs[1] + 1):yr.ind
  I.hist <- mean(Data@Ind[x, hist.yrs], na.rm = TRUE) * modifier
  
  ## Current data...
  yr.ind <- which.max(Data@Year)
  curr.yrs <- (yr.ind - yrs[2]+ 1):yr.ind 
  I.curr <- mean(Data@Ind[x, curr.yrs], na.rm = TRUE)
  
  # plot(Data@Year, Data@Ind[1,], type='l')
  # lines(Data@Year[hist.yrs],Data@Ind[x, hist.yrs], col='blue', lwd=2)
  # lines(Data@Year[curr.yrs],Data@Ind[x, curr.yrs], col='red', lwd=2)
  #
  ## Alpha...
  alpha <- I.curr / I.hist  * tunepar
  
  ## Catch data...
  Cat <- mean(Data@Cat[x, curr.yrs])
 
  ## Harvest control rule...
  # if(alpha > 1.2) {
  #   alpha <- 1.2
  # } else if(alpha < 0.8) {
  #   alpha <- 0.8
  # } else {
  #   alpha <- alpha
  # }
 
  TAC <- alpha  * Cat
  ## TAC <- MSEtool::TACfilter(TAC)
  ## Return TAC...
  ## Rec@TAC <- TAC
  Rec@TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
  
  Rec
}



# Results:
# - doSmoother = TRUE - worse performance
# - incHCR = TRUE worse performance

IR1 <- IndexRatio
class(IR1) <- 'MP'

IR2 <- IndexRatio
formals(IR2)$mc <- c(0.2,0.2)  
class(IR2) <- 'MP'

IR3 <- IndexRatio
formals(IR3)$mc <- c(0.15,0.15)  
class(IR3) <- 'MP'


# Test different configurations of function arguments
# HistList <- readRDS("03_Hists/HistList.rda")
# 
# Hist <- HistList[[1]]
# 
# MSE <- Project(Hist, MPs=c('IR1', 'IR2', 'IR3'))
# Pplot(MSE)
