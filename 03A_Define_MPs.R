
# Constant Catch -----

CC_40kt <- function(x, Data, Data_Lag = 2, Interval = 3,
                    Initial_MP_Yr = 2025, reps = 1, ...) {
  Rec <- new("Rec")
  ## TAC will be set for data_year+1
  ## ie data_year will be 2020 in projection year 2021
  data_year <- max(Data@Year)
  if(data_year < max(Catchdf$Year)) {
    Rec <- FixedTAC(Rec, Data) ## use actual catches if they are
    ## available
    ## return(Rec)
  } else {
    Rec@TAC <- rep(4e4, reps)
    ## return(Rec)
  }
  return(Rec)
}
class(CC_40kt) <- "MP"

CC_30kt <- function(x, Data, Data_Lag = 2, Interval = 3,
                    Initial_MP_Yr = 2025, reps = 1, ...) {
  Rec <- new("Rec")
  ## TAC will be set for data_year+1
  ## ie data_year will be 2020 in projection year 2021
  data_year <- max(Data@Year)
  if(data_year < max(Catchdf$Year)) {
    Rec <- FixedTAC(Rec, Data)
    return(Rec)
  } else {
    Rec@TAC <- rep(3e4, reps)
    return(Rec)
  }
}
class(CC_30kt) <- "MP"

CC_20kt <- function(x, Data, Data_Lag = 2, Interval = 3,
                    Initial_MP_Yr = 2025, reps = 1, ...) {
  Rec <- new("Rec")
  ## TAC will be set for data_year+1
  ## ie data_year will be 2020 in projection year 2021
  data_year <- max(Data@Year)
  if(data_year < max(Catchdf$Year)) {
    Rec <- FixedTAC(Rec, Data)
    return(Rec)
  } else {
    Rec@TAC <- rep(25e3, reps)
    return(Rec)
  }
}
class(CC_20kt) <- "MP"