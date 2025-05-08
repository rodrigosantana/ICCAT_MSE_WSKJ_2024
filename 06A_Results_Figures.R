library(openMSE)
library(ggplot2)

ManagementOptions <- list.dirs('04_MSEs',
                               recursive = FALSE, 
                               full.names = TRUE)

# Select the DataLag and Management Interval options
MngOption <- ManagementOptions[1] # 'DataLag_1_Interval_3

MPpaths <- list.dirs(MngOption, recursive = FALSE, full.names = TRUE)


MPpaths <- MPpaths[grepl("1_30", MPpaths)] # MPs tuned to 1-30


CalcPGK <- function(DF, Years=2021:2054) {
  yrs <- paste(range(Years), collapse='-')
  fDF <- DF |> dplyr::filter(Year %in% Years) |> 
    dplyr::group_by(MP) |> 
    dplyr::mutate(InGreen=SB_SBMSY>1 & F_FMSY<1) |>
    dplyr::summarise(PGK=mean(InGreen))
  names(fDF)[2] <- paste(names(fDF)[2], yrs)
  fDF
}


# ----- Create Data Frames -----

GetInfo <- function(MPpath) {
  msefiles <- list.files(MPpath, full.names = TRUE)
  if (length(msefiles)<1)
    stop('No files found for: ', basename(MPpath))
  
  OMList <- list()
  for (j in seq_along(msefiles)) {
    mse <- readRDS(msefiles[j])
    t <- strsplit(mse@Name, '_')[[1]][3:4]
    
    ProjectYears <- seq( mse@OM$CurrentYr[1]+1, 
                         by=1, 
                         length.out=mse@proyears)
   
    OMList[[j]] <- data.frame(Sim=1:mse@nsim,
                              MP=basename(MPpath),
                              Growth=t[1],
                              Steepness=t[2],
                              Year=rep(ProjectYears, each=mse@nsim),
                              SB_SBMSY=as.vector(mse@SB_SBMSY),
                              F_FMSY=as.vector(mse@F_FMSY),
                              Removals=as.vector(mse@Removals),
                              TAC=as.vector(mse@TAC),
                              Landings=as.vector(mse@Catch))
  }
  do.call('rbind', OMList)
}


MPList <- lapply(MPpaths, GetInfo)


DF <- do.call('rbind', MPList)

CalcPGK(DF)

DF |> dplyr::filter(Year==2025) |> 
  dplyr::group_by(MP) |>
  dplyr::summarise(mean(SB_SBMSY),
                   min(SB_SBMSY))

# ---- Time Series Plot -----
plotVars <- c('F_FMSY', 'SB_SBMSY', 'TAC', 'Removals', 'Landings')
ribbonCols <- c('darkgray', 'black')

## ---- All OMs Combined ----
DFplot <- DF |> tidyr::pivot_longer(dplyr::all_of(plotVars)) |> 
  dplyr::group_by(Year, MP, name) |>
  dplyr::summarise(Median=median(value),
                   Lower1=quantile(value, 0.025),
                   Upper1=quantile(value, 0.975),
                   Lower2=quantile(value, 0.15),
                   Upper2=quantile(value, 0.85))
                   
DFplot$name <- factor(DFplot$name, levels=plotVars, ordered = TRUE)
DFplot$MP <- gsub('1_30', '', DFplot$MP)

ggplot(DFplot, aes(x=Year)) +
  facet_grid(name~MP, scales='free') +
  expand_limits(y=0) +
  geom_ribbon(aes(ymin = Lower1,
                  ymax = Upper1),
              alpha = 0.5, fill = ribbonCols[1]) + 
  geom_ribbon(aes(ymin = Lower2,
                  ymax = Upper2),
              alpha = 0.5, fill = ribbonCols[2]) +
  geom_line(aes(y=Median)) +
  theme_bw()


## ---- SB/SBMSY By OM ----
DFplotbyOM <- DF |>
  tidyr::pivot_longer(c('Growth', 'Steepness'))  |> 
  dplyr::group_by(Year, MP, name, value) |>
  dplyr::summarise(Median=median(SB_SBMSY),
                   Lower1=quantile(SB_SBMSY, 0.025),
                   Upper1=quantile(SB_SBMSY, 0.975),
                   Lower2=quantile(SB_SBMSY, 0.15),
                   Upper2=quantile(SB_SBMSY, 0.85))

DFplotbyOM$MP <- gsub('1_30', '', DFplotbyOM$MP)

ggplot(DFplotbyOM, aes(x=Year)) +
  facet_grid(value~MP) +
  geom_hline(yintercept = 1, linetype=2) +
  expand_limits(y=0) +
  geom_ribbon(aes(ymin = Lower1,
                  ymax = Upper1),
              alpha = 0.5, fill = ribbonCols[1]) + 
  geom_ribbon(aes(ymin = Lower2,
                  ymax = Upper2),
              alpha = 0.5, fill = ribbonCols[2]) +
  geom_line(aes(y=Median)) +
  theme_bw()



# ---- Kobe Time ----


# ---- Violin ----


CalcPGK(DF)


PGKm <- sapply(MSE_list, function(X) {
  # Years 4 - 10 - index 5:11 becuase first projection year is 2025 before MP is used
  mean(X@SB_SBMSY[ , , 5:11] > 1 & X@F_FMSY[ , , 5:11] < 1)
})
PGKw <- mean(PGKm)


data.frame()


# ----- Time Series Plot -----
