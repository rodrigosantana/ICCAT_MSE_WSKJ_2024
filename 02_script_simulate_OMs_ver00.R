########################################################################
## Description: Simulating Operating Models based on the new OMs...
##
## Maintainer: Datenkraft - ICCAT/SCRS - TT MSE
## Author: Rodrigo Sant'Ana
## Created: dom ago 25 14:20:53 2024 (-0300)
## Version: 0.0.1
##
## URL:
## Doc URL:
##
## Database info:
##
### Commentary:
##
### Code:
########################################################################

########################################################################
######@> Setup R...

######@> Update openMSE R packages...
######@> Update openMSE R packages...
source("000_load_openMSE_packages.R")

######@> Loading R packages...
library(openMSE)
library(dplyr)
library(reshape2)
library(readr)
library(ggplot2)
library(tidyverse)
library(viridis)
library(janitor)
library(readxl)
# library(ggmse)
library(patchwork)

######@> Instaling MSEextra...
## MSEextra(force = TRUE)
# library(MSEextra)
## library(SWOMSE)

######@> ggplot theme...
extrafont::loadfonts(device = "postscript")
rgb01 <- "black"
rgb02 <- "black"
seta <- grid::arrow(length = grid::unit(0.2, "cm"), type = "open")
seta2 <- grid::arrow(length = grid::unit(0.2, "cm"), type = "open",
                     ends = "both")
my_theme <- function(base_size = 18, base_family = "Helvetica") {
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
        theme(axis.ticks = element_line(colour = rgb01),
              axis.line = element_line(colour = rgb01, linewidth = 0.2),
              axis.text = element_text(colour = rgb02, size = 14),
              axis.title = element_text(size = 18),
              legend.background = element_blank(),
              legend.key = element_blank(),
              panel.background = element_blank(),
              panel.grid = element_line(linetype = "solid",
                                        linewidth = 0.1,
                                        colour = "gray90"),
              plot.background = element_blank(),
              complete = TRUE)
}

######@> Creating a folder to receive the hist objects...
if(!dir.exists("03_Hists"))
    dir.create("03_Hists")

######@> OM folder...
path01 <- "02_OMs/"

######@> MPs folder...
path02 <- "../R-Work/"


########################################################################
######@> Loading data...

######@> Operating models...
OMs <- readRDS(paste0(path01, "SS_Operating_models_IVInds_ver02.rds"),
               refhook = NULL)

#####@> Setting names for OMs...
names(OMs) <- paste0("OM", sprintf("%03d", 1:9))

########################################################################
######@> Simulating Historical Data...

######@> List of objects...
OM_Objects <- names(OMs)

######@> Looping to create historical data...
HistList <- list()
for(i in seq_along(OM_Objects)) {
  OM <- OMs[[i]]
  Hist <- Simulate(OM, parallel = FALSE, silent = FALSE)
  HistList[[i]] <- Hist
  nm <- paste0(OM_Objects[i], "_IVInds_ver02", ".hist")
  saveRDS(Hist, file.path("03_Hists", nm))
}

######@> Examine Recruitment Deviations
sims <- sample(1:100, 3) # random 3 simulations

par(mfrow=c(3,3),  mar=c(1,1,1,1), oma=c(1,1,1,1))
for (i in 1:9) {
  devs <- t(HistList[[i]]@SampPars$Stock$Perr_y[sims,])
  devs <- devs[45:105,]
  matplot(  log(devs), type='l', xlab='', ylab='',
          lwd=1.3, ylim=c(-1, 1))
}





######@> Check openMSE vs SS3 reference points

compareRefs <- function() {
  Histfiles <- list.files("03_Hists")

  ReplistList <- HistList <- list()
  for (i in seq_along(Histfiles)) {
    message("Reading Hist ", i)
    HistList[[i]] <- readRDS(file.path('03_Hists', Histfiles[[i]]))
  }
  
  SSdirs <- list.dirs("01_SS", recursive = FALSE, full.names = TRUE)

  DFList <- list()
  for (i in seq_along(HistList)) {
    Hist <- HistList[[i]]
    omname <- gsub("OM", "", Hist@OM@Name) |> trimws()
    SSind <- which(grepl(omname, basename(SSdirs)))
    replist <- r4ss::SS_output(SSdirs[SSind])
  
    vars <- c("SSB_MSY", "annF_MSY", "Dead_Catch_MSY", "Ret_Catch_MSY")
    df <- replist$derived_quants |> dplyr::filter(Label %in% vars) |> 
      dplyr::select(Label, SS=Value)
    df$Label <- c('SBMSY', 'FMSY', "MSYDead", "MSYRetain")
    
    df$OM <- c(  Hist@Ref$ReferencePoints$SSBMSY[1],
                 Hist@Ref$ReferencePoints$FMSY[1],
                 Hist@Ref$ReferencePoints$MSY[1],
                 Hist@Ref$ReferencePoints$MSY[1])
    df$SS <- round(df$SS,3)
    df$OM <- round(df$OM,3)
    df$ratio <- df$SS/df$OM
    df$OMname <- omname 
    DFList[[i]] <- df
  }
  do.call('rbind', DFList)
}

CompareRefDF <- compareRefs()

# all equal except FMSY consistently lower in SS3
CompareRefDF |> dplyr::group_by(Label) |>
  dplyr::summarise(Mean=mean(ratio))



# Commenting out code below - was used to update Ref Points in Hist Objects
# No longer needed as Ref Points in OM are ok
# Possibly can be deleted but keeping for now

# #####@> Functions to extract some quantities...
# 
# ####@> Reference values...
# get_Quantities <- function(Hists) {
#     temp <- lapply(
#         Hists,
#         function(hist) {
#             cbind.data.frame(
#                 scenario = str_sub(hist@OM@Name, 4, 31),
#                 msy = mean(hist@Ref$ByYear$MSY),
#                 fmsy = mean(hist@Ref$ByYear$FMSY),
#                 ssbmsy = mean(hist@Ref$ByYear$SSBMSY)
#             )})
#     temp <- do.call("rbind.data.frame", temp)
#     return(temp)
# }
# 
# ####@> Trajectory values...
# get_Find <- function(Hist) {
#     value <-
#         as.data.frame(stack(as.data.frame(Hist@TSdata$Find)))$values
#     year <- rep(1952:2020, each = 300)
#     scenario <- Hist@OM@Name
#     sim <- dim(Hist@TSdata$Number)[1]
#     df <- data.frame(
#         scenario = gsub("OM ", "", scenario),
#         year = year,
#         sim = 1:sim,
#         value = value,
#         variable = "Find",
#         period = "Historical",
#         model = "Model 1"
#     )
# }
# 
# get_Quantities02 <- function(Hists, variable = "SSB") {
#     if(variable == "SSB") {
#         temp <- lapply(
#             Hists,
#             function(hist) {
#                 cbind.data.frame(
#                     scenario = str_sub(hist@OM@Name, 4, 31),
#                     temp02 = get_SSB(hist)
#                 )})
#     } else {
#         temp <- lapply(
#             Hists,
#             function(hist) {
#                 cbind.data.frame(
#                     temp02 = get_Find(hist)
#                 )})
#     }
#     temp <- do.call("rbind.data.frame", temp)
#     names(temp) <- c("scenario", "year", "sim", "value", "variable",
#                      "period", "model")
#     return(temp)
# }
# 
# 
# 
# ########################################################################
# ######@> Custom Historical Data...
# 
# ######@> Loading data...
# 
# #####@> Reference quantities from the last stock assessment...
# load("05_Results/Reference_Quantities_ver03.RData")
# 
# #####@> Reference catches from the last stock assessment...
# load("05_Results/Reference_Catches_ver02.RData")
# 
# #####@> Loading history objects...
# Hists <- sapply(dir("03_Hists", pattern = "IVInds_ver02",
#                     full.names = TRUE), readRDS)
# 
# ####@> Setting names for Hists...
# names(Hists) <- paste0("OM", sprintf("%03d", 1:9))

# ######@>----------------------------------------------------------------
# ######@> Comparing values - Reference Quantities...
# 
# ######@> Quantities...
# tab01 <- get_Quantities(Hists)
# 
# #####@> Figures...
# 
# ####@> MSY...
# temp <- tab01 %>%
#     left_join(msy, by = "scenario") %>%
#     select(scenario, "openMSE" = msy, "SS3" = replist1) %>%
#     pivot_longer(names_to = "Method", values_to = "MSY", 2:3)
# 
# p00 <- ggplot(data = temp, aes(x = scenario, y = MSY, fill = Method)) +
#     geom_bar(stat = "identity", position = "dodge", colour = "black",
#              alpha = 0.8) +
#     labs(x = "Uncertainty grid scenario", y = "MSY") +
#     scale_y_continuous(expand = c(0, 0), limits = c(0, 60000)) +
#     scale_fill_grey() +
#     my_theme() +
#     theme(legend.position = "none") +
#     coord_flip()
# p00
# 
# ####@> F_MSY...
# temp <- tab01 %>%
#     left_join(fmsy, by = "scenario") %>%
#     select(scenario, "openMSE" = fmsy, "SS3" = replist1) %>%
#     pivot_longer(names_to = "Method", values_to = "FMSY", 2:3)
# 
# p01 <- ggplot(data = temp, aes(x = scenario, y = FMSY, fill = Method)) +
#     geom_bar(stat = "identity", position = "dodge", colour = "black",
#              alpha = 0.8) +
#     labs(x = "", y = expression(F[MSY])) +
#     scale_y_continuous(expand = c(0, 0), limits = c(0, 1.3)) +
#     scale_fill_grey() +
#     my_theme() +
#     theme(legend.position = "none") +
#     coord_flip()
# p01
# 
# ####@> SSB_MSY...
# temp <- tab01 %>%
#     left_join(ssbmsy, by = "scenario") %>%
#     select(scenario, "openMSE" = ssbmsy, "SS3" = replist1) %>%
#     pivot_longer(names_to = "Method", values_to = "SSBMSY", 2:3)
# 
# p02 <- ggplot(data = temp, aes(x = scenario, y = SSBMSY, fill = Method)) +
#     geom_bar(stat = "identity", position = "dodge", colour = "black",
#              alpha = 0.8) +
#     labs(x = "", y = expression(SSB[MSY])) +
#     scale_y_continuous(expand = c(0, 0), limits = c(0, 120000)) +
#     scale_fill_grey() +
#     my_theme() +
#     coord_flip()
# p02
# 
# #####@> Saving figures...
# plot00 <- p00 | p01 | p02
# ggsave("05_Results/Comp_Ref_Quantities_ver00.png", plot = plot00,
#        device = "png", units = "cm", w = 55, h = 25, dpi = 300,
#        bg = "white")
# 
# #####@> Replacing the values...
# 
# ####@> OM001...
# Hists[[1]]@OM@Name
# Hists[[1]]@Ref$ByYear$MSY[] <-
#     msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt25_h6"]
# Hists[[1]]@Ref$ByYear$FMSY[] <-
#     fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt25_h6"]
# Hists[[1]]@Ref$ByYear$SSBMSY[] <-
#     ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt25_h6"]
# 
# Hists[[1]]@Ref$ReferencePoints$MSY <-
#     msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt25_h6"]
# Hists[[1]]@Ref$ReferencePoints$FMSY <-
#     fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt25_h6"]
# Hists[[1]]@Ref$ReferencePoints$SSBMSY <-
#     ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt25_h6"]
# 
# ####@> OM002...
# Hists[[2]]@OM@Name
# Hists[[2]]@Ref$ByYear$MSY[] <-
#     msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt50_h6"]
# Hists[[2]]@Ref$ByYear$FMSY[] <-
#     fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt50_h6"]
# Hists[[2]]@Ref$ByYear$SSBMSY[] <-
#     ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt50_h6"]
# 
# Hists[[2]]@Ref$ReferencePoints$MSY <-
#     msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt50_h6"]
# Hists[[2]]@Ref$ReferencePoints$FMSY <-
#     fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt50_h6"]
# Hists[[2]]@Ref$ReferencePoints$SSBMSY <-
#     ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt50_h6"]
# 
# ####@> OM003...
# Hists[[3]]@OM@Name
# Hists[[3]]@Ref$ByYear$MSY[] <-
#     msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt75_h6"]
# Hists[[3]]@Ref$ByYear$FMSY[] <-
#     fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt75_h6"]
# Hists[[3]]@Ref$ByYear$SSBMSY[] <-
#     ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt75_h6"]
# 
# Hists[[3]]@Ref$ReferencePoints$MSY <-
#     msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt75_h6"]
# Hists[[3]]@Ref$ReferencePoints$FMSY <-
#     fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt75_h6"]
# Hists[[3]]@Ref$ReferencePoints$SSBMSY <-
#     ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt75_h6"]
# 
# ####@> OM004...
# Hists[[4]]@OM@Name
# Hists[[4]]@Ref$ByYear$MSY[] <-
#     msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt25_h7"]
# Hists[[4]]@Ref$ByYear$FMSY[] <-
#     fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt25_h7"]
# Hists[[4]]@Ref$ByYear$SSBMSY[] <-
#     ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt25_h7"]
# 
# Hists[[4]]@Ref$ReferencePoints$MSY <-
#     msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt25_h7"]
# Hists[[4]]@Ref$ReferencePoints$FMSY <-
#     fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt25_h7"]
# Hists[[4]]@Ref$ReferencePoints$SSBMSY <-
#     ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt25_h7"]
# 
# ####@> OM005...
# Hists[[5]]@OM@Name
# Hists[[5]]@Ref$ByYear$MSY[] <-
#     msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt50_h7"]
# Hists[[5]]@Ref$ByYear$FMSY[] <-
#     fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt50_h7"]
# Hists[[5]]@Ref$ByYear$SSBMSY[] <-
#     ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt50_h7"]
# 
# Hists[[5]]@Ref$ReferencePoints$MSY <-
#     msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt50_h7"]
# Hists[[5]]@Ref$ReferencePoints$FMSY <-
#     fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt50_h7"]
# Hists[[5]]@Ref$ReferencePoints$SSBMSY <-
#     ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt50_h7"]
# 
# ####@> OM006...
# Hists[[6]]@OM@Name
# Hists[[6]]@Ref$ByYear$MSY[] <-
#     msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt75_h7"]
# Hists[[6]]@Ref$ByYear$FMSY[] <-
#     fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt75_h7"]
# Hists[[6]]@Ref$ByYear$SSBMSY[] <-
#     ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt75_h7"]
# 
# Hists[[6]]@Ref$ReferencePoints$MSY <-
#     msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt75_h7"]
# Hists[[6]]@Ref$ReferencePoints$FMSY <-
#     fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt75_h7"]
# Hists[[6]]@Ref$ReferencePoints$SSBMSY <-
#     ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt75_h7"]
# 
# ####@> OM007...
# Hists[[7]]@OM@Name
# Hists[[7]]@Ref$ByYear$MSY[] <-
#     msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt25_h8"]
# Hists[[7]]@Ref$ByYear$FMSY[] <-
#     fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt25_h8"]
# Hists[[7]]@Ref$ByYear$SSBMSY[] <-
#     ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt25_h8"]
# 
# Hists[[7]]@Ref$ReferencePoints$MSY <-
#     msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt25_h8"]
# Hists[[7]]@Ref$ReferencePoints$FMSY <-
#     fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt25_h8"]
# Hists[[7]]@Ref$ReferencePoints$SSBMSY <-
#     ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt25_h8"]
# 
# ####@> OM008...
# Hists[[8]]@OM@Name
# Hists[[8]]@Ref$ByYear$MSY[] <-
#     msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt50_h8"]
# Hists[[8]]@Ref$ByYear$FMSY[] <-
#     fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt50_h8"]
# Hists[[8]]@Ref$ByYear$SSBMSY[] <-
#     ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt50_h8"]
# 
# Hists[[8]]@Ref$ReferencePoints$MSY <-
#     msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt50_h8"]
# Hists[[8]]@Ref$ReferencePoints$FMSY <-
#     fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt50_h8"]
# Hists[[8]]@Ref$ReferencePoints$SSBMSY <-
#     ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt50_h8"]
# 
# ####@> OM009...
# Hists[[9]]@OM@Name
# Hists[[9]]@Ref$ByYear$MSY[] <-
#     msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt75_h8"]
# Hists[[9]]@Ref$ByYear$FMSY[] <-
#     fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt75_h8"]
# Hists[[9]]@Ref$ByYear$SSBMSY[] <-
#     ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt75_h8"]
# 
# Hists[[9]]@Ref$ReferencePoints$MSY <-
#     msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt75_h8"]
# Hists[[9]]@Ref$ReferencePoints$FMSY <-
#     fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt75_h8"]
# Hists[[9]]@Ref$ReferencePoints$SSBMSY <-
#     ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt75_h8"]
# 
# ######@>----------------------------------------------------------------
# ######@> Comparing values - Trajectories...
# 
# ######@> Quantities...
# tab02.SSB <- get_Quantities02(Hists, variable = "SSB") %>%
#     left_join(ssbmsy, by = "scenario") %>%
#     mutate(SSB_SSBMSY = value/replist1)
# tab02.F <- get_Quantities02(Hists, variable = "F")
# tab02.F <- tab02.F %>%
#     left_join(fmsy, by = "scenario") %>%
#     mutate(F_FMSY = value/replist1)
# 
# #####@> Average levels...
# tab02.SSBm <- tab02.SSB %>%
#     group_by(scenario, year) %>%
#     summarise(SSB = mean(value, na.rm = TRUE),
#               SSB_SSBMSY = mean(SSB_SSBMSY, na.rm = TRUE))
# tab02.Fm <- tab02.F %>%
#     group_by(scenario, year) %>%
#     summarise(F = mean(value, na.rm = TRUE),
#               F_FMSY = mean(F_FMSY, na.rm = TRUE))
# 
# #####@> Figures...
# 
# ####@> SSB...
# p03 <- ggplot() +
#     geom_line(data = filter(ssb, Yr %in% 1952:2020),
#               aes(x = Yr, y = replist1, colour = "SS3",
#                   linetype = "SS3"),
#               linewidth = 1) +
#     geom_line(data = tab02.SSBm,
#               aes(x = year, y = SSB, colour = "openMSE",
#                   linetype = "openMSE"),
#               linewidth = 1) +
#     scale_colour_manual(values =
#                             c("SS3" = "red", "openMSE" = "black")) +
#     scale_linetype_manual(values =
#                               c("SS3" = "solid", "openMSE" = "dashed"))+
#     facet_wrap(~scenario) +
#     labs(x = "Year", y = "SSB", colour = "Method", linetype = "Method") +
#     my_theme()
# p03
# 
# ####@> SSB_SSBMSY...
# p04 <- ggplot() +
#     geom_line(data = filter(ssb_ssbmsy, Year %in% 1952:2020),
#               aes(x = Year, y = stock, colour = "SS3",
#                   linetype = "SS3"),
#               linewidth = 1) +
#     geom_line(data = tab02.SSBm,
#               aes(x = year, y = SSB_SSBMSY, colour = "openMSE",
#                   linetype = "openMSE"),
#               linewidth = 1) +
#     scale_colour_manual(values =
#                             c("SS3" = "red", "openMSE" = "black")) +
#     scale_linetype_manual(values =
#                               c("SS3" = "solid", "openMSE" = "dashed"))+
#     facet_wrap(~scenario) +
#     labs(x = "Year", y = expression(SSB/SSB[MSY]),
#          colour = "Method", linetype = "Method") +
#     my_theme()
# p04
# 
# ####@> F_FMSY...
# p05 <- ggplot() +
#     geom_line(data = filter(ffmsy, Yr %in% 1952:2020),
#               aes(x = Yr, y = replist1, colour = "SS3",
#                   linetype = "SS3"),
#               linewidth = 1) +
#     geom_line(data = tab02.Fm,
#               aes(x = year, y = F_FMSY, colour = "openMSE",
#                   linetype = "openMSE"),
#               linewidth = 1) +
#     scale_colour_manual(values =
#                             c("SS3" = "red", "openMSE" = "black")) +
#     scale_linetype_manual(values =
#                               c("SS3" = "solid", "openMSE" = "dashed"))+
#     facet_wrap(~scenario) +
#     labs(x = "Year", y = expression(F/F[MSY]),
#          colour = "Method", linetype = "Method") +
#     my_theme()
# p05
# 
# #####@> Saving figures...
# ggsave("05_Results/Comp_SSB_Quantities_ver00.png", plot = p03,
#        device = "png", units = "cm", w = 35, h = 30, dpi = 300,
#        bg = "white")
# 
# ggsave("05_Results/Comp_SSB_SSBMSY_Quantities_ver00.png", plot = p04,
#        device = "png", units = "cm", w = 35, h = 30, dpi = 300,
#        bg = "white")
# 
# ggsave("05_Results/Comp_F_FMSY_Quantities_ver00.png", plot = p05,
#        device = "png", units = "cm", w = 35, h = 30, dpi = 300,
#        bg = "white")
# 
# #####@> Replacing the values...
# 
# ####@> To equalize the F from openMSE to F from SS3 will be necessary to
# ####@> update Find slot in Hist object to be the product of F_FMSY and
# ####@> FMSY estimated in SS3...
# 
# ###@> Combine values from SS3...
# SS3.F <- ffmsy %>%
#     left_join(select(fmsy, scenario, "FMSY" = replist1),
#               by = "scenario") %>%
#     mutate(F = replist1 * FMSY) %>%
#     filter(Yr %in% 1952:2020)
# 
# ###@> OM001...
# Hists[[1]]@OM@Name
# Hists[[1]]@TSdata$Find <-
#     matrix(
#         rep(SS3.F$F[SS3.F$scenario == "WSKJ_EstRec93_Qnt25_h6"], each = 300),
#         nrow = 300, byrow = FALSE)
# 
# ###@> OM002...
# Hists[[2]]@OM@Name
# Hists[[2]]@TSdata$Find <-
#     matrix(
#         rep(SS3.F$F[SS3.F$scenario == "WSKJ_EstRec93_Qnt50_h6"], each = 300),
#         nrow = 300, byrow = FALSE)
# 
# ###@> OM003...
# Hists[[3]]@OM@Name
# Hists[[3]]@TSdata$Find <-
#     matrix(
#         rep(SS3.F$F[SS3.F$scenario == "WSKJ_EstRec93_Qnt75_h6"], each = 300),
#         nrow = 300, byrow = FALSE)
# 
# ###@> OM004...
# Hists[[4]]@OM@Name
# Hists[[4]]@TSdata$Find <-
#     matrix(
#         rep(SS3.F$F[SS3.F$scenario == "WSKJ_EstRec93_Qnt25_h7"], each = 300),
#         nrow = 300, byrow = FALSE)
# 
# ###@> OM005...
# Hists[[5]]@OM@Name
# Hists[[5]]@TSdata$Find <-
#     matrix(
#         rep(SS3.F$F[SS3.F$scenario == "WSKJ_EstRec93_Qnt50_h7"], each = 300),
#         nrow = 300, byrow = FALSE)
# 
# ###@> OM006...
# Hists[[6]]@OM@Name
# Hists[[6]]@TSdata$Find <-
#     matrix(
#         rep(SS3.F$F[SS3.F$scenario == "WSKJ_EstRec93_Qnt75_h7"], each = 300),
#         nrow = 300, byrow = FALSE)
# 
# ###@> OM007...
# Hists[[7]]@OM@Name
# Hists[[7]]@TSdata$Find <-
#     matrix(
#         rep(SS3.F$F[SS3.F$scenario == "WSKJ_EstRec93_Qnt25_h8"], each = 300),
#         nrow = 300, byrow = FALSE)
# 
# ###@> OM008...
# Hists[[8]]@OM@Name
# Hists[[8]]@TSdata$Find <-
#     matrix(
#         rep(SS3.F$F[SS3.F$scenario == "WSKJ_EstRec93_Qnt50_h8"], each = 300),
#         nrow = 300, byrow = FALSE)
# 
# ###@> OM009...
# Hists[[9]]@OM@Name
# Hists[[9]]@TSdata$Find <-
#     matrix(
#         rep(SS3.F$F[SS3.F$scenario == "WSKJ_EstRec93_Qnt75_h8"], each = 300),
#         nrow = 300, byrow = FALSE)
# 
# ######@> Comparing new trajectories...
# ######@> Quantities...
# tab02.SSB <- get_Quantities02(Hists, variable = "SSB") %>%
#     left_join(ssbmsy, by = "scenario") %>%
#     mutate(SSB_SSBMSY = value/replist1)
# tab02.F <- get_Quantities02(Hists, variable = "F")
# tab02.F <- tab02.F %>%
#     left_join(fmsy, by = "scenario") %>%
#     mutate(F_FMSY = value/replist1)
# 
# #####@> Average levels...
# tab02.SSBm <- tab02.SSB %>%
#     group_by(scenario, year) %>%
#     summarise(SSB = mean(value, na.rm = TRUE),
#               SSB_SSBMSY = mean(SSB_SSBMSY, na.rm = TRUE))
# tab02.Fm <- tab02.F %>%
#     group_by(scenario, year) %>%
#     summarise(F = mean(value, na.rm = TRUE),
#               F_FMSY = mean(F_FMSY, na.rm = TRUE))
# 
# #####@> Figures...
# 
# ####@> SSB...
# p06 <- ggplot() +
#     geom_line(data = filter(ssb, Yr %in% 1952:2020),
#               aes(x = Yr, y = replist1, colour = "SS3",
#                   linetype = "SS3"),
#               linewidth = 1) +
#     geom_line(data = tab02.SSBm,
#               aes(x = year, y = SSB, colour = "openMSE",
#                   linetype = "openMSE"),
#               linewidth = 1) +
#     scale_colour_manual(values =
#                             c("SS3" = "red", "openMSE" = "black")) +
#     scale_linetype_manual(values =
#                               c("SS3" = "solid", "openMSE" = "dashed"))+
#     facet_wrap(~scenario) +
#     labs(x = "Year", y = "SSB", colour = "Method", linetype = "Method") +
#     my_theme()
# p06
# 
# ####@> SSB_SSBMSY...
# p07 <- ggplot() +
#     geom_line(data = filter(ssb_ssbmsy, Year %in% 1952:2020),
#               aes(x = Year, y = stock, colour = "SS3",
#                   linetype = "SS3"),
#               linewidth = 1) +
#     geom_line(data = tab02.SSBm,
#               aes(x = year, y = SSB_SSBMSY, colour = "openMSE",
#                   linetype = "openMSE"),
#               linewidth = 1) +
#     scale_colour_manual(values =
#                             c("SS3" = "red", "openMSE" = "black")) +
#     scale_linetype_manual(values =
#                               c("SS3" = "solid", "openMSE" = "dashed"))+
#     facet_wrap(~scenario) +
#     labs(x = "Year", y = expression(SSB/SSB[MSY]),
#          colour = "Method", linetype = "Method") +
#     my_theme()
# p07
# 
# ####@> F_FMSY...
# p08 <- ggplot() +
#     geom_line(data = filter(ffmsy, Yr %in% 1952:2020),
#               aes(x = Yr, y = replist1, colour = "SS3",
#                   linetype = "SS3"),
#               linewidth = 1) +
#     geom_line(data = tab02.Fm,
#               aes(x = year, y = F_FMSY, colour = "openMSE",
#                   linetype = "openMSE"),
#               linewidth = 1) +
#     scale_colour_manual(values =
#                             c("SS3" = "red", "openMSE" = "black")) +
#     scale_linetype_manual(values =
#                               c("SS3" = "solid", "openMSE" = "dashed"))+
#     facet_wrap(~scenario) +
#     labs(x = "Year", y = expression(F/F[MSY]),
#          colour = "Method", linetype = "Method") +
#     my_theme()
# p08
# 
# #####@> Saving figures...
# ggsave("05_Results/NEW_Comp_SSB_Quantities_ver00.png", plot = p06,
#        device = "png", units = "cm", w = 35, h = 30, dpi = 300,
#        bg = "white")
# 
# ggsave("05_Results/NEW_Comp_SSB_SSBMSY_Quantities_ver00.png", plot = p07,
#        device = "png", units = "cm", w = 35, h = 30, dpi = 300,
#        bg = "white")
# 
# ggsave("05_Results/NEW_Comp_F_FMSY_Quantities_ver00.png", plot = p08,
#        device = "png", units = "cm", w = 35, h = 30, dpi = 300,
#        bg = "white")
# 
# ######@>----------------------------------------------------------------
# ######@> Comparing values - Reference Catches...
# 
# ######@> Quantities...
# tab03 <- get_Landings(Hists)
# 
# #####@> Summary of the Historical simulations...
# temp00 <- tab03 %>%
#     group_by(Model, Year) %>%
#     summarise(m = mean(Value))
# 
# #####@> Merging bases...
# temp01 <- temp00 %>%
#     left_join(catch, by = "Year") %>%
#     mutate(Diff = Obs - m)
# 
# #####@> Looking for differences between OMs...
# 
# ####@> Table...
# SummTab <- tab03 %>%
#     group_by(Model) %>%
#     summarise(m = mean(Value),
#               d = sd(Value))
# 
# ####@> Figure...
# p09 <- ggplot(data = temp01, aes(x = Year, y = m)) +
#     geom_area(fill = "gray", colour = "black") +
#     geom_line(aes(y = Obs), colour = "red") +
#     geom_hline(yintercept = 30000, lty = 2, alpha = 0.4) +
#     facet_wrap(~Model) +
#     scale_y_continuous(expand = c(0, 0), limits = c(0, 50000)) +
#     labs(x = "Year", y = "Catches (mt)") +
#     my_theme()
# p09
# 
# #####@> Saving figures...
# ggsave("05_Results/Comp_Catches_Quantities_ver00.png", plot = p09,
#        device = "png", units = "cm", w = 35, h = 35, dpi = 300,
#        bg = "white")
# 
# ######@> Replacing values...
# for(i in 1:9) {
#     Hists[[i]]@TSdata$Landings[,,1] <- matrix(rep(catch$Half, each = 300),
#                                               nrow = 300, byrow = FALSE)
#     Hists[[i]]@TSdata$Landings[,,2] <- matrix(rep(catch$Half, each = 300),
#                                               nrow = 300, byrow = FALSE)
# }
# 
# ######@> Comparing values again...
# tab04 <- get_Landings(Hists)
# 
# #####@> Summary of the Historical simulations...
# temp00 <- tab04 %>%
#     group_by(Model, Year) %>%
#     summarise(m = mean(Value))
# 
# #####@> Merging bases...
# temp01 <- temp00 %>%
#     left_join(catch, by = "Year") %>%
#     mutate(Diff = Obs - m)
# 
# #####@> Looking for differences between OMs...
# 
# ####@> Table...
# SummTab2 <- tab04 %>%
#     group_by(Model) %>%
#     summarise(m = mean(Value),
#               d = sd(Value))
# 
# ####@> Figure...
# p10 <- ggplot(data = temp01, aes(x = Year, y = m)) +
#     geom_area(fill = "gray", colour = "black") +
#     geom_line(aes(y = Obs), colour = "red") +
#     geom_hline(yintercept = 30000, lty = 2, alpha = 0.4) +
#     facet_wrap(~Model) +
#     scale_y_continuous(expand = c(0, 0), limits = c(0, 50000)) +
#     labs(x = "Year", y = "Catches (mt)") +
#     my_theme()
# p10
# 
# #####@> Saving figures...
# ggsave("05_Results/New_Comp_Catches_Quantities_ver00.png", plot = p10,
#        device = "png", units = "cm", w = 35, h = 35, dpi = 300,
#        bg = "white")
# 
# ########################################################################
# ######@> Exporting corrected simulations...
# 
# ######@> Looping to export / save simulated historical process...
# for(i in 1:9) {
#     nm <- paste0(OM_Objects[i], "_IVInds_CORRECTED_ver03", ".hist")
#     saveRDS(Hists[[i]], file.path("03_Hists", nm))
# }

########################################################################
##
##                  Creative Commons License 4.0
##                       (CC BY-NC-SA 4.0)
##
##  This is a humam-readable summary of (and not a substitute for) the
##  license (https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
##
##  You are free to:
##
##  Share - copy and redistribute the material in any medium or format.
##
##  The licensor cannot revoke these freedoms as long as you follow the
##  license terms.
##
##  Under the following terms:
##
##  Attribution - You must give appropriate credit, provide a link to
##  license, and indicate if changes were made. You may do so in any
##  reasonable manner, but not in any way that suggests the licensor
##  endorses you or your use.
##
##  NonCommercial - You may not use the material for commercial
##  purposes.
##
##  ShareAlike - If you remix, transform, or build upon the material,
##  you must distributive your contributions under the same license
##  as the  original.
##
##  No additional restrictions — You may not apply legal terms or
##  technological measures that legally restrict others from doing
##  anything the license permits.
##
########################################################################
