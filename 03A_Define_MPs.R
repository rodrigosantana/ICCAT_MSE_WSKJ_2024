
#



source('03AA_MP_Internal_Functions.R')

#####@> Constant Catch MPs -----
ConstantCatchGeneric <- function(x, Data, Interval = 3,
                                 Initial_MP_Yr = 2026, TAC=4E4, ...) {
  Rec <- new("Rec")
  
  # Check if TAC needs to be updated
  if (SameTAC(Initial_MP_Yr, Interval, Data)) {
    Rec@TAC <- Data@MPrec[x]
    Rec <- FixedTAC(Rec, Data) # use actual catches if they are available
    return(Rec)
  }

  # Set TAC
  Rec@TAC <-TAC
  
  Rec
}

CC_40kt <- ConstantCatchGeneric
formals(CC_40kt)$TAC <- 4E4

CC_30kt <- ConstantCatchGeneric
formals(CC_30kt)$TAC <- 3E4

CC_20kt <- ConstantCatchGeneric
formals(CC_30kt)$TAC <- 2E4

#####@> Relative Abundance Index MPs ----

#####@> Iratio with asymetrical TAC correction... 
IR_01 <- function(x, Data, Data_Lag = 1, Interval = 3,
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
  index <- smoothed_index <- Data@Ind[x,]
  smoothed <- stats::smooth(index[!is.na(index)])
  smoothed_index[!is.na(smoothed_index)] <- smoothed
  Data@Ind[x,] <- smoothed_index
  
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
  if(alpha > 1.2) {
    alpha <- 1.2
  } else if(alpha < 0.8) {
    alpha <- 0.8
  } else {
    alpha <- alpha
  }
  ## TAC...
  TAC <- alpha  * Cc
  ## TAC <- MSEtool::TACfilter(TAC)
  ## Return TAC...
  ## Rec@TAC <- TAC
  Rec@TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)

  Rec
}
class(IR_01) <- "MP"

#####@> Iratio with symetrical TAC correction... 
IR_02 <- function(x, Data, Data_Lag = 1, Interval = 3,
                  Initial_MP_Yr = 2026, reps =  1,
                  tunepar = 1.143986, mc = c(0.2, 0.2),
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
  
  ## Smooth combined index
  index <- smoothed_index <- Data@Ind[x,]
  smoothed <- stats::smooth(index[!is.na(index)])
  smoothed_index[!is.na(smoothed_index)] <- smoothed
  Data@Ind[x,] <- smoothed_index
  
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
  if(alpha > 1.2) {
    alpha <- 1.2
  } else if(alpha < 0.8) {
    alpha <- 0.8
  } else {
    alpha <- alpha
  }
  ## TAC...
  TAC <- alpha  * Cc
  ## TAC <- MSEtool::TACfilter(TAC)
  ## Return TAC...
  ## Rec@TAC <- TAC
  Rec@TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
  return(Rec)
}
class(IR_02) <- "MP"


#####@> Iratio without TAC correction...
IR_03 <- function(x, Data, Data_Lag = 1, Interval = 3,
                  Initial_MP_Yr = 2026, reps =  1,
                  tunepar = 1.143986, mc = NA,
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
  
  ## Smooth combined index
  index <- smoothed_index <- Data@Ind[x,]
  smoothed <- stats::smooth(index[!is.na(index)])
  smoothed_index[!is.na(smoothed_index)] <- smoothed
  Data@Ind[x,] <- smoothed_index
  
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
  if(alpha > 1.2) {
    alpha <- 1.2
  } else if(alpha < 0.8) {
    alpha <- 0.8
  } else {
    alpha <- alpha
  }
  ## TAC...
  TAC <- alpha * Cc
  ## TAC <- MSEtool::TACfilter(TAC)
  ## Return TAC...
  ## Rec@TAC <- TAC
  Rec@TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
  return(Rec)
  
}
class(IR_03) <- "MP"

# ----- Constant Exploitation -----
#####@> Constant Exploitation with Control Rule MP with asymetrical TAC
#####@> correction...
CE_01 <- function(x, Data, Data_Lag = 1, Interval = 3,
                  Initial_MP_Yr = 2026, reps =  1, 
                  tunepar = 1.792558, mc = c(0.25, 0.2),
                  yrs = c(5, 3), ...) {
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
  index <- smoothed_index <- Data@Ind[x,]
  smoothed <- stats::smooth(index[!is.na(index)])
  smoothed_index[!is.na(smoothed_index)] <- smoothed
  Data@Ind[x,] <- smoothed_index
  
  ## Calculate Historical Relative Exploitation Rate
  yr.ind <- which(Data@Year == 2017)
  hist.yrs <- (yr.ind - yrs[1] + 1):yr.ind
  histER <- mean(Data@Cat[x, hist.yrs]) / mean(Data@Ind[x, hist.yrs])
  ## Calculate Current Relative Exploitation Rate
  current_yr <- length(Data@Ind[x,])
  recent_yrs <- (current_yr - yrs[2]+1):current_yr
  curER <- mean(Data@Cat[x, recent_yrs]) * tunepar /
    mean(Data@Ind[x, recent_yrs]) 
  ## Control Rule
  histInd <- mean(Data@Ind[x, hist.yrs])
  curInd <- mean(Data@Ind[x, recent_yrs], na.rm = TRUE)
  ind_ratio <- curInd/histInd
  if (ind_ratio >= 0.8) {
    targER <- histER
  } else if (ind_ratio > 0.5) {
    targER <- histER * (-1.4 + 3 * ind_ratio)
  } else {
    targER <- 0.1 * histER
  }
  ## Exploitation Rate Ratio
  ER_ratio <- targER/curER
  TAC <- ER_ratio  * Data@MPrec[x]
  ## Maximum allowed change in TAC
  ## Rec@TAC <- TAC
  Rec@TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
  return(Rec)
}
class(CE_01) <- "MP"

#####@> Constant Exploitation with Control Rule MP with symetrical TAC
#####@> correction...
CE_02 <- function(x, Data, Data_Lag = 1, Interval = 3,
                  Initial_MP_Yr = 2026, reps =  1,
                   tunepar = 1.792558, mc = c(0.2, 0.2),
                  yrs = c(5, 3), ...) {
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
  index <- smoothed_index <- Data@Ind[x,]
  smoothed <- stats::smooth(index[!is.na(index)])
  smoothed_index[!is.na(smoothed_index)] <- smoothed
  Data@Ind[x,] <- smoothed_index
  ## Lag Data
  
  ## Calculate Historical Relative Exploitation Rate
  yr.ind <- which(Data@Year == 2017)
  hist.yrs <- (yr.ind - yrs[1] + 1):yr.ind
  histER <- mean(Data@Cat[x, hist.yrs]) / mean(Data@Ind[x, hist.yrs])
  ## Calculate Current Relative Exploitation Rate
  current_yr <- length(Data@Ind[x,])
  recent_yrs <- (current_yr - yrs[2]+1):current_yr
  curER <- mean(Data@Cat[x, recent_yrs]) * tunepar /
    mean(Data@Ind[x, recent_yrs]) 
  ## Control Rule
  histInd <- mean(Data@Ind[x, hist.yrs])
  curInd <- mean(Data@Ind[x, recent_yrs], na.rm = TRUE)
  ind_ratio <- curInd/histInd
  if (ind_ratio >= 0.8) {
    targER <- histER
  } else if (ind_ratio > 0.5) {
    targER <- histER * (-1.4 + 3 * ind_ratio)
  } else {
    targER <- 0.1 * histER
  }
  ## Exploitation Rate Ratio
  ER_ratio <- targER/curER
  TAC <- ER_ratio  * Data@MPrec[x]
  ## Maximum allowed change in TAC
  ## Rec@TAC <- TAC
  Rec@TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
  return(Rec)
  
}
class(CE_02) <- "MP"

#####@> Constant Exploitation with Control Rule MP without TAC
#####@> correction...
CE_03 <- function(x, Data, Data_Lag = 1, Interval = 3,
                  Initial_MP_Yr = 2026, reps =  1,
                  tunepar = 1.792558, mc = NA,
                  yrs = c(5, 3), ...) {
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
  index <- smoothed_index <- Data@Ind[x,]
  smoothed <- stats::smooth(index[!is.na(index)])
  smoothed_index[!is.na(smoothed_index)] <- smoothed
  Data@Ind[x,] <- smoothed_index
  ## Lag Data
  
  ## Calculate Historical Relative Exploitation Rate
  yr.ind <- which(Data@Year == 2017)
  hist.yrs <- (yr.ind - yrs[1] + 1):yr.ind
  histER <- mean(Data@Cat[x, hist.yrs]) / mean(Data@Ind[x, hist.yrs])
  ## Calculate Current Relative Exploitation Rate
  current_yr <- length(Data@Ind[x,])
  recent_yrs <- (current_yr - yrs[2]+1):current_yr
  curER <- mean(Data@Cat[x, recent_yrs]) * tunepar /
    mean(Data@Ind[x, recent_yrs]) 
  ## Control Rule
  histInd <- mean(Data@Ind[x, hist.yrs])
  curInd <- mean(Data@Ind[x, recent_yrs], na.rm = TRUE)
  ind_ratio <- curInd/histInd
  if (ind_ratio >= 0.8) {
    targER <- histER
  } else if (ind_ratio > 0.5) {
    targER <- histER * (-1.4 + 3 * ind_ratio)
  } else {
    targER <- 0.1 * histER
  }
  ## Exploitation Rate Ratio
  ER_ratio <- targER/curER
  TAC <- ER_ratio  * Data@MPrec[x]
  ## Maximum allowed change in TAC
  ## Rec@TAC <- TAC
  Rec@TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
  return(Rec)
  
}
class(CE_03) <- "MP"

# ------ Index Slope MPS -----

#####@> Islope1 MP with asymetrical TAC correction...
IS_01 <- function(x, Data, Data_Lag = 1, Interval = 3,
                  Initial_MP_Yr = 2026, reps =  1, lambda = 0.4,
                  tunepar = 0.9, mc = c(0.25, 0.2), yrsmth = 5, ...) {
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
  
  
  ## Smooth combined index
  index <- smoothed_index <- Data@Ind[x,]
  smoothed <- stats::smooth(index[!is.na(index)])
  smoothed_index[!is.na(smoothed_index)] <- smoothed
  Data@Ind[x,] <- smoothed_index
  
  ## Applied Islope1 MP to lagged data
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)
  Years <- Data@Year[ind]
  ylast <- (Data@LHYear[1] - Data@Year[1]) + 1
  C_dat <- Data@Cat[x, ind]
  if (is.na(Data@MPrec[x]) || length(Data@Year) == ylast + 1) {
    TACstar <- (1 - xx) *
      MSEtool::trlnorm(reps, mean(C_dat, na.rm = TRUE),
                       Data@CV_Cat/(yrsmth^0.5))
  } else {
    TACstar <- rep(Data@MPrec[x], reps)
  }
  I_hist <- Data@Ind[x, ind] * tunepar
  yind <- 1:yrsmth
  slppar <- summary(lm(log(I_hist) ~ yind))$coefficients[2, 1:2]
  if (reps > 1) {
    Islp <- rnorm(reps, slppar[1], slppar[2])
  } else {
    Islp <- slppar[1]
  }
  TAC <- TACstar * (1 + lambda * Islp) 
  ## Return TAC...
  ## Rec@TAC <- TAC
  Rec@TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
  return(Rec)
  
}
class(IS_01) <- "MP"

#####@> Islope1 MP with symetrical TAC correction...
IS_02 <- function(x, Data, Data_Lag = 1, Interval = 3,
                  Initial_MP_Yr = 2026, reps =  1, lambda = 0.4,
                  tunepar = 0.9, mc = c(0.2, 0.2), yrsmth = 5, ...) {
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
  
  ## Smooth combined index
  index <- smoothed_index <- Data@Ind[x,]
  smoothed <- stats::smooth(index[!is.na(index)])
  smoothed_index[!is.na(smoothed_index)] <- smoothed
  Data@Ind[x,] <- smoothed_index
  
  ## Applied Islope1 MP to lagged data
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)
  Years <- Data@Year[ind]
  ylast <- (Data@LHYear[1] - Data@Year[1]) + 1
  C_dat <- Data@Cat[x, ind]
  if (is.na(Data@MPrec[x]) || length(Data@Year) == ylast + 1) {
    TACstar <- (1 - xx) *
      MSEtool::trlnorm(reps, mean(C_dat, na.rm = TRUE),
                       Data@CV_Cat/(yrsmth^0.5))
  } else {
    TACstar <- rep(Data@MPrec[x], reps)
  }
  I_hist <- Data@Ind[x, ind] * tunepar
  yind <- 1:yrsmth
  slppar <- summary(lm(log(I_hist) ~ yind))$coefficients[2, 1:2]
  if (reps > 1) {
    Islp <- rnorm(reps, slppar[1], slppar[2])
  } else {
    Islp <- slppar[1]
  }
  TAC <- TACstar * (1 + lambda * Islp) 
  ## Return TAC...
  ## Rec@TAC <- TAC
  Rec@TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
  return(Rec)
  
}
class(IS_02) <- "MP"

#####@> Islope1 MP without TAC correction...
IS_03 <- function(x, Data, Data_Lag = 1, Interval = 3,
                  Initial_MP_Yr = 2026, reps =  1, lambda = 0.4,
                  tunepar = 0.9, mc = NA, yrsmth = 5, ...) {
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
  
  ## Smooth combined index
  index <- smoothed_index <- Data@Ind[x,]
  smoothed <- stats::smooth(index[!is.na(index)])
  smoothed_index[!is.na(smoothed_index)] <- smoothed
  Data@Ind[x,] <- smoothed_index
  
  ## Applied Islope1 MP to lagged data
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)
  Years <- Data@Year[ind]
  ylast <- (Data@LHYear[1] - Data@Year[1]) + 1
  C_dat <- Data@Cat[x, ind]
  if (is.na(Data@MPrec[x]) || length(Data@Year) == ylast + 1) {
    TACstar <- (1 - xx) *
      MSEtool::trlnorm(reps, mean(C_dat, na.rm = TRUE),
                       Data@CV_Cat/(yrsmth^0.5))
  } else {
    TACstar <- rep(Data@MPrec[x], reps)
  }
  I_hist <- Data@Ind[x, ind]  * tunepar
  yind <- 1:yrsmth
  slppar <- summary(lm(log(I_hist) ~ yind))$coefficients[2, 1:2]
  if (reps > 1) {
    Islp <- rnorm(reps, slppar[1], slppar[2])
  } else {
    Islp <- slppar[1]
  }
  TAC <- TACstar * (1 + lambda * Islp)
  ## Return TAC...
  ## Rec@TAC <- TAC
  Rec@TAC <- adjust_TAC(TAC, Data@MPrec[x], mc)
  return(Rec)
}
class(IS_03) <- "MP"

# ----- SP Models ----

######@> A Surplus Production model with a Harvest Control Rule that
######@> linearly reduces F (Schaeffer model - B 80-40) - harvest
######@> control rule based on: [with TAC fixed first year]
######@> https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf...
SP_01 <- function(x, Data, Data_Lag = 1, Interval = 3,
                  Initial_MP_Yr = 2026, reps =  1, tunepar = 0.8796482, mc = NA, ...) {
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
  
  
  ## Smooth combined index
  index <- smoothed_index <- Data@Ind[x,]
  smoothed <- stats::smooth(index[!is.na(index)])
  smoothed_index[!is.na(smoothed_index)] <- smoothed
  Data@Ind[x,] <- smoothed_index
  
  ## apply SP assessment model
  Mod <- SAMtool::SP(x, Data, prior = list(r = c(0.416, 0.148),
                                           MSY = c(30000, 0.2)),
                     start = list(dep = 0.98, n = 1))
  ## harvest control rule
  ## based on: https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
  Bthresh <- Mod@BMSY
  Blim <- 0.4 * Bthresh
  ## Ftar <- 0.8 * tunepar * Mod@FMSY
  ## Fmin <- 0.1 * tunepar * Mod@FMSY
  Ftar <- tunepar * 0.10
  Fmin <- 0.1 * Ftar
  Bcurr <- Mod@B[length(Mod@B)]
  if(Bcurr >= Bthresh) {
    Fmort <- Ftar
  } else if(Bcurr > Blim) {
    Fmort <- Ftar * (-0.367 + 1.167 * Bcurr/Bthresh)
  } else {
    Fmort <- Fmin
  }
  TAC <-  Fmort*Bcurr
  ## Return TAC...
  ## Rec@TAC <- TAC
  Rec@TAC <- adjust_TAC2(TAC, Data@MPrec[x], mc,
                         B_rel = Bcurr/Bthresh)
  return(Rec)
  
}
class(SP_01) <- "MP"

######@> A Surplus Production Space-State model with a Harvest Control
######@> Rule that linearly reduces F (Schaeffer model B 80-40) -
######@> harvest control rule based on:
######@> https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf...
SP_02 <- function(x, Data, Data_Lag = 1, Interval = 3,
                  Initial_MP_Yr = 2026, reps =  1,
                  tunepar = 0.616231, mc = NA, ...) {
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
  
  ## Smooth combined index
  index <- smoothed_index <- Data@Ind[x,]
  smoothed <- stats::smooth(index[!is.na(index)])
  smoothed_index[!is.na(smoothed_index)] <- smoothed
  Data@Ind[x,] <- smoothed_index
  
  ## apply SP assessment model
  Mod <- SAMtool::SP_SS(x, Data,
                        prior = list(r = c(0.416, 0.148),
                                     MSY = c(30000, 0.2)),
                        start = list(dep = 0.98, n = 1),
                        fix_sigma = FALSE, fix_tau = TRUE)
  ## harvest control rule
  ## based on: https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
  Bthresh <- Mod@BMSY
  Blim <- 0.4 * Bthresh
  ## Ftar <- 0.8 * tunepar * Mod@FMSY
  ## Fmin <- 0.1 * tunepar * Mod@FMSY
  Ftar <- tunepar * 0.10
  Fmin <- 0.1 * Ftar
  Bcurr <- Mod@B[length(Mod@B)]
  if (Bcurr>=Bthresh) {
    Fmort <- Ftar
  } else if (Bcurr>Blim) {
    Fmort <- Ftar * (-0.367 + 1.167 * Bcurr/Bthresh)
  } else {
    Fmort <- Fmin
  }
  TAC <-  Fmort*Bcurr
  ## Return TAC...
  ## Rec@TAC <- TAC
  Rec@TAC <- adjust_TAC2(TAC, Data@MPrec[x], mc,
                         B_rel = Bcurr/Bthresh)
  return(Rec)
 
}
class(SP_02) <- "MP"

######@> A Surplus Production model with a Harvest Control Rule that
######@> linearly reduces F (Schaeffer model B 100-40) - harvest control
######@> rule based on:
######@> https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
######@> HCR ram Ftar = 1, Fmin = 0.1, Blim = 0.4, BMSY = 1...
SP_03 <- function(x, Data, Data_Lag = 1, Interval = 3,
                  Initial_MP_Yr = 2026, reps =  1,
                  tunepar = 0.8796482, mc = NA, ...) {
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
  
  ## Smooth combined index
  index <- smoothed_index <- Data@Ind[x,]
  smoothed <- stats::smooth(index[!is.na(index)])
  smoothed_index[!is.na(smoothed_index)] <- smoothed
  Data@Ind[x,] <- smoothed_index
  
  ## apply SP assessment model
  Mod <- SAMtool::SP(x, Data, prior = list(r = c(0.416, 0.148),
                                           MSY = c(30000, 0.2)),
                     start = list(dep = 0.98, n = 1))
  ## harvest control rule
  ## based on: https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
  Bthresh <- Mod@BMSY
  Blim <- 0.4 * Bthresh
  ## Ftar <- 0.8 * tunepar * Mod@FMSY
  ## Fmin <- 0.1 * tunepar * Mod@FMSY
  Ftar <- tunepar * 0.10
  Fmin <- 0.1 * Ftar
  Bcurr <- Mod@B[length(Mod@B)]
  if(Bcurr >= Bthresh) {
    Fmort <- Ftar
  } else if(Bcurr > Blim) {
    Fmort <- Ftar * (-0.5 + 1.5 * Bcurr/Bthresh)
  } else {
    Fmort <- Fmin
  }
  TAC <-  Fmort*Bcurr
  ## Return TAC...
  ## Rec@TAC <- TAC
  Rec@TAC <- adjust_TAC2(TAC, Data@MPrec[x], mc,
                         B_rel = Bcurr/Bthresh)
  return(Rec)
  
}
class(SP_03) <- "MP"

######@> A Surplus Production Space-State model with a Harvest Control
######@> Rule that linearly reduces F (Schaeffer model) - harvest control
######@> rule based on:
######@> https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
######@> HCR ram Ftar = 1, Fmin = 0.1, Blim = 0.4, BMSY = 1...
SP_04 <- function(x, Data, Data_Lag = 1, Interval = 3,
                  Initial_MP_Yr = 2026, reps =  1, tunepar = 0.616231, mc = NA, ...) {
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
  
  
  ## Smooth combined index
  index <- smoothed_index <- Data@Ind[x,]
  smoothed <- stats::smooth(index[!is.na(index)])
  smoothed_index[!is.na(smoothed_index)] <- smoothed
  Data@Ind[x,] <- smoothed_index

  ## apply SP assessment model
  Mod <- SAMtool::SP_SS(x, Data,
                        prior = list(r = c(0.416, 0.148),
                                     MSY = c(30000, 0.2)),
                        start = list(dep = 0.98, n = 1),
                        fix_sigma = FALSE, fix_tau = TRUE)
  ## harvest control rule
  ## based on: https://www.iccat.int/Documents/Recs/compendiopdf-e/2017-04-e.pdf
  Bthresh <- Mod@BMSY
  Blim <- 0.4 * Bthresh
  ## Ftar <- 0.8 * tunepar * Mod@FMSY
  ## Fmin <- 0.1 * tunepar * Mod@FMSY
  Ftar <- tunepar * 0.10
  Fmin <- 0.1 * Ftar
  Bcurr <- Mod@B[length(Mod@B)]
  if (Bcurr>=Bthresh) {
    Fmort <- Ftar
  } else if (Bcurr>Blim) {
    Fmort <- Ftar * (-0.5 + 1.5 * Bcurr/Bthresh)
  } else {
    Fmort <- Fmin
  }
  TAC <-  Fmort*Bcurr
  ## Return TAC...
  ## Rec@TAC <- TAC
  Rec@TAC <- adjust_TAC2(TAC, Data@MPrec[x], mc,
                         B_rel = Bcurr/Bthresh)
  return(Rec)
}
class(SP_04) <- "MP"
