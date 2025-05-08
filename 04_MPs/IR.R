
source('03AA_MP_Internal_Functions.R')

# Index Ratio MPs 

# mc - maximum change in TAC between management cycles 

IndexRatio <- function(x, Data, Data_Lag = 1, Interval = 3,
                  Initial_MP_Yr = 2026, reps =  1, 
                  tunepar = 1.143986, mc = c(0.25, 0.2), 
                  yrs = c(3, 5), ...) {
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
  yr.ind <- which(Data@Year == max(Data@Year) - yrs[1])
  hist.yrs <- (yr.ind - yrs[2] + 1):yr.ind
  I.hist <- mean(Data@Ind[x, hist.yrs], na.rm = TRUE) * tunepar
  ## Current data...
  curr.yrs <- (yr.ind + 1):(yr.ind + yrs[1])
  I.curr <- mean(Data@Ind[x, curr.yrs], na.rm = TRUE)
  ## Alpha...
  alpha <- I.curr / I.hist 
  ## Catch data...
  Cat <- Data@Cat[x, length(Data@Cat[x, ])]
  Cc <- MSEtool::trlnorm(reps = 1, Cat, Data@CV_Cat[x, 1])
  
  ## Harvest control rule...
  # if(alpha > 1.2) {
  #   alpha <- 1.2
  # } else if(alpha < 0.8) {
  #   alpha <- 0.8
  # } else {
  #   alpha <- alpha
  # }
  ## TAC...
  TAC <- alpha  * Cc
  ## TAC <- MSEtool::TACfilter(TAC)
  ## Return TAC...
  ## Rec@TAC <- TAC
  Rec@TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
  
  Rec
}
class(IR_01) <- "MP"