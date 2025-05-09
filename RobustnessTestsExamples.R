
# Example code to demonstrate changes to OM for robustness tests

##########################################
#                                        #
# Changes to Historical fishery dynamics #
#                                        #
##########################################

# Changes to historical fishery dynamics :

## Option 1: Modify and run SS3
# 1. change data, parameters, etc in SS3 
# 2. re-run SS3 assessment model
# 3. Import OM with usual process
# OM dynamics will match SS3 output

# Option 2: Modify OM@cpars
# 1. Modify OM@cpars - similar process to Hist below - see ?ValidCpars
# 2. Run `Simulate` with modified OM
# 3. OM fishery dynamics will no longer match SS3 output


################################################################################
#
# NOTE: for Option 2 above and `Hist` object below, only changes to arrays 
# that have a time and age (if appicable) dimension will have an impact on 
# the simulation.
#
# E.g., If you want to modify natural mortalit rate:
#       - Changes to `OM@M` or `OM@cpars$M` will have no effect 
#         (vector length `nsim`)
#       - Change must be made on `OM@cpars$M_ageArray`  
#         ( array  dim = c(nsim, maxage+1, nyears+proyears) )
#        
################################################################################

####################################################
#                                                  #
# Changes to Projection dynamics in `Hist` object  #
#                                                  #
####################################################

HistList <- readRDS("03_Hists/HistList.rda")

om <- 1

Hist <- HistList[[om]]
OM <- Hist@OM

# ----- Recruitment Deviations -----
# OM@cpars$Perr_y
RecDevs <- Hist@SampPars$Stock$Perr_y
dim(RecDevs) # nsim, maxage + nyears+proyears
c(OM@nsim, OM@maxage+OM@nyears+OM@proyears)

ProjInd <- (OM@maxage+OM@nyears+1):dim(RecDevs)[2]
ProjYears <- seq(OM@CurrentYr+1, by=1, length.out=OM@proyears)
length(ProjInd) == OM@proyears

ProjectionRecDevs <- RecDevs[,ProjInd]

prevSD <- sd(log(ProjectionRecDevs))
newSD  <- prevSD * 1.5 # eg increase variabilty by 50%
NewProjectionRecDevs <- exp(rnorm(OM@nsim*OM@proyears, 0, newSD))
NewProjectionRecDevs <- matrix(NewProjectionRecDevs, OM@nsim, OM@proyears)

par(mfrow=c(1,2))
matplot(ProjYears, t(ProjectionRecDevs), type='l', ylab='', ylim=c(0,4))
matplot(ProjYears, t(NewProjectionRecDevs), type='l', ylab='', ylim=c(0,4))

# replace in Hist
Hist@SampPars$Stock$Perr_y[,ProjInd] <- NewProjectionRecDevs
# MSEnewRecDevs <- Project(Hist, ...)


# ----- Growth -----

# Length and Weight-at-Age
LengthAtAge <- Hist@SampPars$Stock$Len_age # OM@nsim, OM@nyears+OM@proyears
WeightAtAge <- Hist@SampPars$Stock$Wt_age # OM@nsim, OM@nyears+OM@proyears
FleetWeightAtAge <- Hist@SampPars$Fleet$Wt_age_C  # OM@nsim, OM@nyears+OM@proyears

# FleetWeightAtAge = empirical weight at age used to calculate Catch Biomas from Catch Numbers

# example - linear increase in K for projection years
Linf <- Hist@SampPars$Stock$Linf[1] # same value for all simulations
K <- Hist@SampPars$Stock$K[1] # same value for all simulations
t0 <- Hist@SampPars$Stock$t0[1] # same value for all simulations

Wa <- Hist@SampPars$Stock$a
Wb <- Hist@SampPars$Stock$b

newK <- seq(K, K*1.5, length.out=OM@proyears)
Ages <- 0:Hist@SampPars$Stock$maxage

newLengthAge <- sapply(seq_along(newK), function(i)
  Linf *  (1-exp(-newK[i] *(Ages-t0)))
  )

par(mfrow=c(1,1))
plot(Ages, LengthAtAge[1,,OM@nyears+1], type='l')
lines(Ages, newLengthAge[,1], col='blue')
lines(Ages, newLengthAge[,OM@proyears], col='green')

# adjust weight-at-age
newWeightAtAge <- Wa*newLengthAge^Wb

par(mfrow=c(1,1))
plot(Ages, WeightAtAge[1,,OM@nyears+1], type='l')
lines(Ages, newWeightAtAge[,1], col='blue')
lines(Ages, newWeightAtAge[,OM@proyears], col='green')



plot(LengthAtAge[1,,OM@nyears], WeightAtAge[1,,OM@nyears], type='l')
lines(LengthAtAge[1,,OM@nyears], FleetWeightAtAge[1,,OM@nyears], col='blue')

mod <- lm(log(FleetWeightAtAge[1,,OM@nyears])~log(LengthAtAge[1,,OM@nyears]))
params <- summary(mod)$coefficients[1:2,1] 

newFleetWeightAtAge <- exp(params[1])*newLengthAge^params[2]


par(mfrow=c(1,1))
plot(Ages, FleetWeightAtAge[1,,OM@nyears+1], type='l')
lines(Ages, newFleetWeightAtAge[,1], col='blue')
lines(Ages, newFleetWeightAtAge[,OM@proyears], col='green')


# Update Arrays - add sim dimension
newLengthAge <- replicate(OM@nsim, newLengthAge) |> aperm(c(3,1,2))
newWeightAtAge <- replicate(OM@nsim, newWeightAtAge) |> aperm(c(3,1,2))
newFleetWeightAtAge <- replicate(OM@nsim, newFleetWeightAtAge) |> aperm(c(3,1,2))

projind <- (OM@nyears+1):(OM@nyears+OM@proyears)
Hist@SampPars$Stock$Len_age[,,projind] <- newLengthAge
Hist@SampPars$Stock$Wt_age[,,projind] <- newWeightAtAge
Hist@SampPars$Fleet$Wt_age_C[,,projind] <- newFleetWeightAtAge

# MSEnewGrowth <- Project(Hist, ...)


