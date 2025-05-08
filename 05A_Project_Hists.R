if(packageVersion("MSEtool") == "3.7.9999") {
  print("MSEtool R package already installed is the correct version")
} else {
  pak::pkg_install("blue-matter/MSEtool")
}

#####@> SAMtool (1.7.0)...
if(packageVersion("SAMtool") == "1.8.1") {
  print("SAMtool R package already installed is the correct version")
} else {
  pak::pkg_install("blue-matter/SAMtool")
}

#####@> openMSE (1.1.1)...
if(packageVersion("openMSE") == "1.0.1") {
  print("openMSE R package already installed is the correct version")
} else {
  pak::pkg_install("blue-matter/openMSE")
}


library(dplyr)
library(openMSE)

source('03AA_MP_Internal_Functions.R')

########################################################################
######@> Setup MSE...

######@> Creating a folder to receive the MSE objects...
if(!dir.exists("04_MSEs"))
  dir.create("04_MSEs")

######@> OM folder...
path01 <- "03_Hists"

######@> MPs folder...
path02 <- "../R-Work/"

########################################################################
######@> Simulating Projections Data...

Hists <- readRDS("03_Hists/HistList.rda")


ManagementOptions <- list.dirs('TunedMPs', recursive = FALSE, full.names = TRUE)

# Select the DataLag and Management Interval options
MngOption <- ManagementOptions[1] # 'DataLag_1_Interval_3'

MSEDir <- file.path('04_MSEs', basename(MngOption))
if (!dir.exists(MSEDir))
  dir.create(MSEDir)


tunedMPs <- list.files(MngOption, full.names = TRUE, recursive = )
MPnames <- gsub('.mp', '', basename(tunedMPs))


################################################################################
#
# Select MPs to Project - some may already be done 
#
################################################################################
ProjectMPs <- MPnames # MPs to project

ProjectMPs <- ProjectMPs[grepl("1_30", ProjectMPs)] # MPs tuned to 1-30

################################################################################


for (i in seq_along(ProjectMPs)) {
  ind <- match(ProjectMPs[i], MPnames)
  mp <- readRDS(tunedMPs[ind])
  assign(ProjectMPs[i], mp, envir=.GlobalEnv)
}


#### Loop Over MPs and then over OMs 
for (mp in ProjectMPs) {
  message('\nProjecting MP: ', mp)
  if (!dir.exists(file.path(MSEDir, mp)))
    dir.create(file.path(MSEDir, mp))
  for (i in seq_along(Hists)) {
    message('OM: ', i, '/', length(Hists))
    Hist <- Hists[[i]]
    nm <- gsub('OM', '', Hist@OM@Name) |> trimws() 
    nm <- paste(sprintf("%03d", i), nm, sep="_")
    nm <- paste0(nm, '.mse')
    MSE <- Project(Hist, MPs = mp, parallel = FALSE, silent=TRUE)
    saveRDS(MSE, file.path(MSEDir, mp, nm))
  }
}

