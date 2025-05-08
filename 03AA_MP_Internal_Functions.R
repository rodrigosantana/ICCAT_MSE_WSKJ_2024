WSKJ_Data <- readRDS("WSKJ_Data.rda")
# assuming 2025 catch is mean of last 3 years
Assumed2025Catch <- tail(WSKJ_Data@Cat[1,],3) |> mean()

Catchdf <- data.frame(Year=2025, Catch=Assumed2025Catch)


#####@> MP Internal Functions ----
# copied from SWOMSE::FixedTAC
FixedTAC <- function(Rec, Data) {
  df <- Catchdf
  if ((max(Data@Year) + 1) %in% df$Year) {
    ind <- match(max(Data@Year) + 1, df$Year)
    Rec@TAC <- df$Catch[ind]
  }
  Rec
}

SameTAC <- function (Initial_MP_Yr, Interval, Data) {
  if (max(Data@Year) < (Initial_MP_Yr - 1)) 
    return(TRUE)
  Imp_Years <- seq(Initial_MP_Yr, by = Interval, length.out = 30)
  if (!(max(Data@Year) + 1) %in% Imp_Years) 
    return(TRUE)
  FALSE
}


######@> Function to adjust symetric, asymetric or not adjust the TAC...
adjust_TAC <- function(TAC, Last_TAC, mc) {
  delta_TAC <- TAC / Last_TAC
  if(any(is.na(mc))) {
    return(TAC)
  }
  lower_bound <- 1 - ifelse(length(mc) == 2, mc[1], mc)
  upper_bound <- 1 + ifelse(length(mc) == 2, mc[2], mc)
  delta_TAC <- pmax(pmin(delta_TAC, upper_bound), lower_bound)
  return(Last_TAC * delta_TAC)
}

######@> Function to adjust TAC for model based MPs...
adjust_TAC2 <- function(TAC, Last_TAC, mc, B_rel) {
  delta_TAC <- TAC / Last_TAC
  if(any(is.na(mc))) {
    return(TAC)
  }
  lower_bound <- if (length(mc) == 2) 1 - mc[1] else 1 - mc
  upper_bound <- if (length(mc) == 2) 1 + mc[2] else 1 + mc
  if(B_rel > 1) {
    delta_TAC <- pmax(pmin(delta_TAC, upper_bound), lower_bound)
  } else {
    delta_TAC <- pmin(delta_TAC, upper_bound)
  }
  return(Last_TAC * delta_TAC)
}