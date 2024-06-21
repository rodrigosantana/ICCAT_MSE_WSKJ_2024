

# ---- 01_script_import_SS3_ver00.R ----

## 1. OM@interval needs to be set to 1 -----

OM1 <- SS2OM(SSdir,
             nsim = 100,
             proyears = 30,
             reps = 100,
             maxF = 3,
             seed = 1,
             # interval = 3,
             interval = 1,
             pstar = 0.5,
             Obs = MSEtool::Precise_Unbiased,
             Imp = MSEtool::Perfect_Imp,
             gender = 1,
             silent = FALSE,
             Name = "OM WSKJ_EstRec93_Qnt25_h6",
             Source = "ICCAT 2022 SKJ assessment",
             Author = "Rodrigo Sant'Ana & Bruno Leite Mourato",
             report = FALSE,
             filename = "OM_WSKJ_EstRec93_Qnt25_h6",
             dir = output.dir,
             open_file = TRUE)

# The reason for this is that we have the code that checks each
# year if the TAC needs to be updated within the MP:
# if (SWOMSE::SameTAC(Initial_MP_Yr, Interval, Data))
# ie the MP code needs to be run every year (OM@interval=1) but
# it will only update the TAC if it's an update year
# messy i know!




## 2. OM@cpars$Data is being updated each OM -----
# This can be a problem when the life-history parameters change between OMs.


# For example:
OM1.Data <- SS2Data(paste0(path01, "WSKJ_EstRec93_Qnt25_h6"),
                    Name = "OM1 Data WSKJ_EstRec93_Qnt25_h6",
                    Common_Name = "Skipjack",
                    Species = "Katsuwonus pelamis",
                    Region = "Western Atlantic Ocean")

OM2.Data <- SS2Data(paste0(path01, "WSKJ_EstRec93_Qnt50_h6"),
                    Name = "OM1 Data WSKJ_EstRec93_Qnt25_h6",
                    Common_Name = "Skipjack",
                    Species = "Katsuwonus pelamis",
                    Region = "Western Atlantic Ocean")


OM1.Data@Mort
OM2.Data@Mort

OM1.Data@vbLinf
OM2.Data@vbLinf

OM1.Data@steep
OM2.Data@steep

# This only matters when you have MPs that use the life-history information in the Data object.
# eg, SCA model that uses M, h, Linf, etc.

# This is probably not an issue now, but to avoid any potential problems, my
# recommendation would be to first create a Data object that
# contains all the actual data/parameters that would be used if the MP was
# applied to it for an TAC recommendation for 2025.
# eg.,

WSJK_Data <- new('Data')

slots_to_copy <- c('Year', 'LHYear', 'Cat', 'CV_Cat')
for (sl in slots_to_copy) {
  slot(WSJK_Data, sl) <- slot(OM1.Data, sl)
}

## ---- NOTE:  ----##
OM1.Data@CAA
##  The CAA data isn't returned from `SS2Data`, I'm not sure why
## You'll need to manually add CAA data if you want to use that for the SCA models
## otherwise it will use simulated CAA data - again, apologies that this isn't clear
# in the openMSE code!

# data manually added:
index <- c(rep(NA, 29),
           tsIndex$Obs[tsIndex$Fleet == "Averaging Index 03"][1:40])
names(index) <- OM1.Data@Year
cv_index <- c(rep(NA, 29), rep(0.2, 40))
names(cv_index) <- OM1.Data@Year

WSJK_Data@Ind <- array(index, dim = c(1, length(index)))
WSJK_Data@CV_Ind <- array(cv_index, dim = c(1, length(index)))

# then loop over OMs and add data:
# eg.
OM1@cpars$Data <- WSJK_Data
OM2@cpars$Data <- WSJK_Data
...

# etc

## 3. Add I_beta to cpars -----
# I found this when I made a plot to compare the index and biomass in the projections
# and found the index was not tracking the biomass AT ALL
#
# this is because the OM estimates the standard deviation, auto-correlation
# AND the hyperstability/depletion parameter internally.
#
# the hyperstability parameter is often poorly estimated, and given that we
# are using the index to set management advice we must assume it tracks the
# abundance.
# We need to set I_beta in cpars
# there was NO WAY you were to know this. This is another thing that needs to
# be improved in MSEtool.

# Need to add I_beta <- 1 to cpars for each OM
# eg.
OM1@cpars$I_beta <- rep(1, OM1@nsim)
OM2@cpars$I_beta <- rep(1, OM2@nsim)
...




# ---- script_testing_MPs_Empirical_ver00.R -----


# OM@interval = 1 fixed the issue where TAC was constant for all years
# if the first TAC set by the MP is 2025:
Initial_MP_Yr <- 2025

# if the MP is setting the first TAC in 2025, we need to hard code the
# catches for 2021:2024

# here i'm assuming 2023 & 2024 catch are mean of 2021:2022
# (I know you said you'll have real data values for these soon, so
# these can be updated once they are available)
Catchdf <- data.frame(
  Year = 2021:2024,
  Catch = c(20048.21, 21377.24, 20712.72, 20712.72))


