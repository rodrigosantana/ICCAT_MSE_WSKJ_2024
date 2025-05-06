####@> Function to generate recruitment deviations from truncated log-normal dist
rtnorm <- function(n, mu, sigma, lower, upper) {
  qnorm(
    runif(n,
          pnorm(lower, mu, sigma),
          pnorm(upper, mu, sigma)),
    mu,
    sigma)
}

GenerateRecDevs <- function(OM, truncSD=2, plot=FALSE) {
  SD <- OM@Perr[1]  # sigmaR 
  AC <- OM@AC[1] # lag 1 autocorrelation
  
  mu <- -0.5 * SD^2  * (1 - AC)/sqrt(1 - AC^2)
  lower <- mu-truncSD*SD
  upper <- mu+truncSD*SD
  
  logDevsOrig <- log(OM@cpars$Perr_y)  # recruitment deviations
  dd <- dim(logDevsOrig)
  
  # generate rec devs from truncated dist
  logDevs <- array(rtnorm(OM@nsim*OM@proyears, mu, SD, lower, upper), 
                   dim=c(OM@nsim, OM@proyears))
  
  # add auto-correlation
  lastHistDev <- logDevsOrig[,dd[2]+1 - OM@proyears]
  for (y in 1:ncol(logDevs)) {
    if (y==1) {
      logDevs[,y] <- AC * lastHistDev + logDevs[,y] * (1-AC*AC)^0.5
    } else {
      logDevs[,y] <- AC * logDevs[,y-1] + logDevs[,y] * (1-AC*AC)^0.5
    }
  }
  
  # rescale to SD 
  logDevs <- (scale(logDevs) + mean(logDevs)) * SD
  
  # replace values
  OM@cpars$Perr_y[,(dd[2]+1 - OM@proyears):dd[2]] <- exp(logDevs)
  
  if (plot) {
    par(mfrow=c(1,2))
    histYrs <- (OM@CurrentYr-OM@nyears-1):OM@CurrentYr
    projYrs <- seq(OM@CurrentYr+1, by=1, length.out=OM@proyears)
    years <- c(histYrs, projYears)
    ind <- (dd[2] - length(years) + 1):dd[2]
    matplot(years, t(logDevsOrig[,ind]), type='l', ylim=c(-1.5, 1.5),
            xlab='Year', ylab='log Recruitment Deviations')
    mtext("Normal Dist", 3)
    matplot(years, t(log(OM@cpars$Perr_y[,ind])), type='l', ylim=c(-1.5, 1.5), 
            xlab='Year', ylab='')
    mtext("Truncated Normal Dist", 3)
  }
  OM
}

####@> Function to add WSJK_Data and I_beta = 1 to cpars...
update_cpars <- function(OM) {
  OM@cpars$Data <- WSKJ_Data
  OM@cpars$I_beta <- rep(1, OM@nsim)
  OM <- GenerateRecDevs(OM) # generate rec devs from truncated dist
  OM
}
