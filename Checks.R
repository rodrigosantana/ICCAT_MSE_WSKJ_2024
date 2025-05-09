



r = SAMtool::SP(1, WSKJ_Data,
                prior = list(r = c(0.416, 0.148),
                                    MSY = c(30000, 0.2)),
                start = list(dep = 0.98, n = 1))


HistB <- Hist@TSdata$Biomass[1,,] |> rowSums()
HistB_BMSY <- HistB/Hist@Ref$ByYear$BMSY[1,1:73]


plot(HistB_BMSY, type='l', ylim=c(0,3))
lines(r@B_BMSY, col='blue')




MSE <- readRDS('04_MSEs/DataLag_1_Interval_3/IR11_30/001_WSKJ_EstRec93_Qnt25_h6.mse')
# MSE <- readRDS('04_MSEs/DataLag_1_Interval_3/IR11_30/002_WSKJ_EstRec93_Qnt50_h6.mse')
# MSE <- readRDS('04_MSEs/DataLag_1_Interval_3/SP_011_30/002_WSKJ_EstRec93_Qnt50_h6.mse')

par(mfrow=c(3,3))
samps <- sample(1:MSE@nsim, 9)
for (x in samps) {
  DF <- openMSE::get_Biomass(MSE) |> dplyr::filter(Sim==x)
  DF$Value <- DF$Value/MSE@RefPoint$ByYear$BMSY[x,]
  
  r = SAMtool::SP(x, MSE@PPD[[1]],
                  prior = list(r = c(0.416, 0.148),
                               MSY = c(30000, 0.2)),
                  start = list(dep = 0.98, n = 1))
  
  
  plot(DF$Year, DF$Value, type='l', ylim=c(0,4))
  lines(DF$Year[DF$Period=='Historical'], DF$Value[DF$Period=='Historical'], col='green', lwd=2)
  abline(h=1, lty=3)
  lines(as.numeric(names(r@B_BMSY)), r@B_BMSY, col='blue')
}



