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
  
  # rec devs estimated by SS3
  estRecDevs <- log(OM@cpars$Perr_y[1,OM@maxage:(OM@maxage+OM@nyears-1)])
  estRecDevs <- estRecDevs[estRecDevs!=0]
  
  # calculate SD and AC from estimated historical rec devs
  SD <- sd(estRecDevs) # OM@Perr[1] # sigmaR 
  AC <- acf(estRecDevs, plot=FALSE)$acf[2] # OM@AC[1] lag 1 autocorrelation

  mu <- -0.5 * SD^2  * (1 - AC)/sqrt(1 - AC^2)
  lower <- mu-truncSD*SD
  upper <- mu+truncSD*SD
  
  logDevsOrig <- log(OM@cpars$Perr_y)  # recruitment deviations
  dd <- dim(logDevsOrig)
  
  # generate rec devs from truncated dist
  set.seed(OM@seed)
  nsim <- OM@nsim
  nyear <- OM@proyears
  nsamp <- nsim * nyear
  logDevs <- array(rtnorm(nsamp, mu, SD, lower, upper), 
                   dim=c(nsim, nyear))
  

  # add auto-correlation
  lastHistDev <- logDevsOrig[,dd[2]+1 - OM@proyears]
  for (y in 1:ncol(logDevs)) {
    if (y==1) {
      logDevs[,y] <- AC * lastHistDev + logDevs[,y] * (1-AC*AC)^0.5
    } else {
      logDevs[,y] <- AC * logDevs[,y-1] + logDevs[,y] * (1-AC*AC)^0.5
    }
  }
  
  # replace values
  OM@cpars$Perr_y[,(dd[2]+1 - OM@proyears):dd[2]] <- exp(logDevs)
  
  if (plot) {
    df <- data.frame(EstRecDevs=as.vector(estRecDevs))
    df2 <- data.frame(SimRecDevs=as.vector(logDevs))  
    dist <- data.frame(x=c(-1, 1))
    
    mu <- mean(df$EstRecDevs)
    sd <- sd(df$EstRecDevs)
    SDdf <- data.frame(intercept=c(mu - sd,
                                   mu + sd,
                                   mu - 2*sd,
                                   mu + 2*sd),
                       type=c('1 SD', '1 SD', '2 SD', '2 SD')
    )
    
    ggplot(df) + 
      geom_histogram(aes(x = EstRecDevs, y =..density..),
                     breaks = seq(-1, 1, by = 0.075), 
                     colour = "black", 
                     fill = "grey") +
      geom_vline(data=SDdf, aes(xintercept=intercept, color=type), linetype=2) +
      stat_function(fun = dnorm, data=dist, aes(x=x),
                    args = list(mean = mean(df$EstRecDevs), 
                                sd = sd(df$EstRecDevs))) +
      theme_bw()
    
    ggplot(df2) + 
      geom_histogram(aes(x = SimRecDevs, y =..density..),
                     breaks = seq(-1, 1, by = 0.075), 
                     colour = "black", 
                     fill = "grey") +
      stat_function(fun = dnorm, data=dist, aes(x=x),
                    args = list(mean = mean(df2$SimRecDevs), 
                                sd = sd(df2$SimRecDevs))) +
      theme_bw()
    
    # par(mfrow=c(1,2))
    # lim <- range(estRecDevs) |> abs() |> max()
    # hist(estRecDevs, main='Estimated Historical', xlim=c(-lim,lim))
    # hist(logDevs, main='Simulated', xlim=c(-lim,lim))
    # 
    # par(mfrow=c(1,2))
    # histYrs <- (OM@CurrentYr-OM@nyears-1):OM@CurrentYr
    # projYrs <- seq(OM@CurrentYr+1, by=1, length.out=OM@proyears)
    # years <- c(histYrs, projYrs)
    # ind <- (dd[2] - length(years) + 1):dd[2]
    # matplot(years, t(logDevsOrig[,ind]), type='l', ylim=c(-1.5, 1.5),
    #         xlab='Year', ylab='log Recruitment Deviations')
    # mtext("Normal Dist", 3)
    # matplot(years, t(log(OM@cpars$Perr_y[,ind])), type='l', ylim=c(-1.5, 1.5), 
    #         xlab='Year', ylab='')
    # mtext("Truncated Normal Dist", 3)
    
  }
  OM
}

####@> Function to add WSJK_Data and I_beta = 1 to cpars...
update_cpars <- function(OM) {
  OM@cpars$Data <- WSKJ_Data
  OM@cpars$I_beta <- rep(1, OM@nsim)
  OM@cpars$control$ObsCatch <- "Removals"
  OM <- GenerateRecDevs(OM) # generate rec devs from truncated dist
  OM
}
