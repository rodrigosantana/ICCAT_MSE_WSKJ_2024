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

#####@> MSEtool (3.7.2)...
if(packageVersion("MSEtool") == "3.7.2") {
    print("MSEtool R package already installed is the correct version")
} else {
    pak::pkg_install("blue-matter/MSEtool")
}

#####@> SAMtool (1.7.0)...
if(packageVersion("SAMtool") == "1.7.0") {
    print("SAMtool R package already installed is the correct version")
} else {
    pak::pkg_install("blue-matter/SAMtool")
}

#####@> openMSE (1.1.1)...
if(packageVersion("openMSE") == "1.1.1") {
    print("openMSE R package already installed is the correct version")
} else {
    pak::pkg_install("blue-matter/openMSE")
}

#####@> nswo-mse (0.29.0)...
if(packageVersion("SWOMSE") == "0.29.0") {
    print("SWOMSE R package already installed is the correct version")
} else {
    pak::pkg_install("ICCAT/nswo-mse")
}

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
library(ggmse)
library(patchwork)

######@> Instaling MSEextra...
## MSEextra(force = TRUE)
library(MSEextra)
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

#####@> Functions to extract some quantities...

####@> Reference values...
get_Quantities <- function(Hists) {
    temp <- lapply(
        Hists,
        function(hist) {
            cbind.data.frame(
                scenario = str_sub(hist@OM@Name, 4, 31),
                msy = mean(hist@Ref$ByYear$MSY),
                fmsy = mean(hist@Ref$ByYear$FMSY),
                ssbmsy = mean(hist@Ref$ByYear$SSBMSY)
            )})
    temp <- do.call("rbind.data.frame", temp)
    return(temp)
}

####@> Trajectory values...
get_Find <- function(Hist) {
    value <-
        as.data.frame(stack(as.data.frame(Hist@TSdata$Find)))$values
    year <- rep(1952:2020, each = 300)
    scenario <- Hist@OM@Name
    sim <- dim(Hist@TSdata$Number)[1]
    df <- data.frame(
        scenario = gsub("OM ", "", scenario),
        year = year,
        sim = 1:sim,
        value = value,
        variable = "Find",
        period = "Historical",
        model = "Model 1"
    )
}

get_Quantities02 <- function(Hists, variable = "SSB") {
    if(variable == "SSB") {
        temp <- lapply(
            Hists,
            function(hist) {
                cbind.data.frame(
                    scenario = str_sub(hist@OM@Name, 4, 31),
                    temp02 = get_SSB(hist)
                )})
    } else {
        temp <- lapply(
            Hists,
            function(hist) {
                cbind.data.frame(
                    temp02 = get_Find(hist)
                )})
    }
    temp <- do.call("rbind.data.frame", temp)
    names(temp) <- c("scenario", "year", "sim", "value", "variable",
                     "period", "model")
    return(temp)
}

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
for(i in seq_along(OM_Objects)) {
    OM <- OMs[[i]]
    OM@nsim <- 300
    Hist <- Simulate(OM, parallel = FALSE, silent = FALSE)
    nm <- paste0(OM_Objects[i], "_IVInds_ver02", ".hist")
    saveRDS(Hist, file.path("03_Hists", nm))
}

########################################################################
######@> Custom Historical Data...

######@> Loading data...

#####@> Reference quantities from the last stock assessment...
load("05_Results/Reference_Quantities_ver02.RData")

#####@> Loading history objects...
Hists <- sapply(dir("03_Hists", full.names = TRUE), readRDS)

####@> Setting names for Hists...
names(Hists) <- paste0("OM", sprintf("%03d", 1:9))

######@>----------------------------------------------------------------
######@> Comparing values - Reference Quantities...

######@> Quantities...
tab01 <- get_Quantities(Hists)

#####@> Figures...

####@> MSY...
temp <- tab01 %>%
    left_join(msy, by = "scenario") %>%
    select(scenario, "openMSE" = msy, "SS3" = replist1) %>%
    pivot_longer(names_to = "Method", values_to = "MSY", 2:3)

p00 <- ggplot(data = temp, aes(x = scenario, y = MSY, fill = Method)) +
    geom_bar(stat = "identity", position = "dodge", colour = "black",
             alpha = 0.8) +
    labs(x = "Uncertainty grid scenario", y = "MSY") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 60000)) +
    scale_fill_grey() +
    my_theme() +
    theme(legend.position = "none") +
    coord_flip()
p00

####@> F_MSY...
temp <- tab01 %>%
    left_join(fmsy, by = "scenario") %>%
    select(scenario, "openMSE" = fmsy, "SS3" = replist1) %>%
    pivot_longer(names_to = "Method", values_to = "FMSY", 2:3)

p01 <- ggplot(data = temp, aes(x = scenario, y = FMSY, fill = Method)) +
    geom_bar(stat = "identity", position = "dodge", colour = "black",
             alpha = 0.8) +
    labs(x = "", y = expression(F[MSY])) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.3)) +
    scale_fill_grey() +
    my_theme() +
    theme(legend.position = "none") +
    coord_flip()
p01

####@> SSB_MSY...
temp <- tab01 %>%
    left_join(ssbmsy, by = "scenario") %>%
    select(scenario, "openMSE" = ssbmsy, "SS3" = replist1) %>%
    pivot_longer(names_to = "Method", values_to = "SSBMSY", 2:3)

p02 <- ggplot(data = temp, aes(x = scenario, y = SSBMSY, fill = Method)) +
    geom_bar(stat = "identity", position = "dodge", colour = "black",
             alpha = 0.8) +
    labs(x = "", y = expression(SSB[MSY])) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 120000)) +
    scale_fill_grey() +
    my_theme() +
    coord_flip()
p02

#####@> Saving figures...
plot00 <- p00 | p01 | p02
ggsave("05_Results/Comp_Ref_Quantities_ver00.png", plot = plot00,
       device = "png", units = "cm", w = 55, h = 25, dpi = 300,
       bg = "white")

#####@> Replacing the values...

####@> OM001...
Hists[[1]]@OM@Name
Hists[[1]]@Ref$ByYear$MSY[] <-
    msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt25_h6"]
Hists[[1]]@Ref$ByYear$FMSY[] <-
    fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt25_h6"]
Hists[[1]]@Ref$ByYear$SSBMSY[] <-
    ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt25_h6"]

####@> OM002...
Hists[[2]]@OM@Name
Hists[[2]]@Ref$ByYear$MSY[] <-
    msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt50_h6"]
Hists[[2]]@Ref$ByYear$FMSY[] <-
    fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt50_h6"]
Hists[[2]]@Ref$ByYear$SSBMSY[] <-
    ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt50_h6"]

####@> OM003...
Hists[[3]]@OM@Name
Hists[[3]]@Ref$ByYear$MSY[] <-
    msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt75_h6"]
Hists[[3]]@Ref$ByYear$FMSY[] <-
    fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt75_h6"]
Hists[[3]]@Ref$ByYear$SSBMSY[] <-
    ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt75_h6"]

####@> OM004...
Hists[[4]]@OM@Name
Hists[[4]]@Ref$ByYear$MSY[] <-
    msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt25_h7"]
Hists[[4]]@Ref$ByYear$FMSY[] <-
    fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt25_h7"]
Hists[[4]]@Ref$ByYear$SSBMSY[] <-
    ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt25_h7"]

####@> OM005...
Hists[[5]]@OM@Name
Hists[[5]]@Ref$ByYear$MSY[] <-
    msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt50_h7"]
Hists[[5]]@Ref$ByYear$FMSY[] <-
    fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt50_h7"]
Hists[[5]]@Ref$ByYear$SSBMSY[] <-
    ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt50_h7"]

####@> OM006...
Hists[[6]]@OM@Name
Hists[[6]]@Ref$ByYear$MSY[] <-
    msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt75_h7"]
Hists[[6]]@Ref$ByYear$FMSY[] <-
    fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt75_h7"]
Hists[[6]]@Ref$ByYear$SSBMSY[] <-
    ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt75_h7"]

####@> OM007...
Hists[[7]]@OM@Name
Hists[[7]]@Ref$ByYear$MSY[] <-
    msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt25_h8"]
Hists[[7]]@Ref$ByYear$FMSY[] <-
    fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt25_h8"]
Hists[[7]]@Ref$ByYear$SSBMSY[] <-
    ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt25_h8"]

####@> OM008...
Hists[[8]]@OM@Name
Hists[[8]]@Ref$ByYear$MSY[] <-
    msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt50_h8"]
Hists[[8]]@Ref$ByYear$FMSY[] <-
    fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt50_h8"]
Hists[[8]]@Ref$ByYear$SSBMSY[] <-
    ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt50_h8"]

####@> OM009...
Hists[[9]]@OM@Name
Hists[[9]]@Ref$ByYear$MSY[] <-
    msy$replist1[msy$scenario == "WSKJ_EstRec93_Qnt75_h8"]
Hists[[9]]@Ref$ByYear$FMSY[] <-
    fmsy$replist1[fmsy$scenario == "WSKJ_EstRec93_Qnt75_h8"]
Hists[[9]]@Ref$ByYear$SSBMSY[] <-
    ssbmsy$replist1[ssbmsy$scenario == "WSKJ_EstRec93_Qnt75_h8"]

######@>----------------------------------------------------------------
######@> Comparing values - Trajectories - I HAD STOPPED HERE...

######@> Quantities...
tab02.SSB <- get_Quantities02(Hists, variable = "SSB") %>%
    left_join(ssbmsy, by = "scenario") %>%
    mutate(SSB_SSBMSY = value/replist1)
tab02.F <- get_Quantities02(Hists, variable = "F")
tab02.F <- tab02.F %>%
    left_join(fmsy, by = "scenario") %>%
    mutate(F_FMSY = value/replist1)

#####@> Average levels...
tab02.SSBm <- tab02.SSB %>%
    group_by(scenario, year) %>%
    summarise(SSB = mean(value, na.rm = TRUE),
              SSB_SSBMSY = mean(SSB_SSBMSY, na.rm = TRUE))
tab02.Fm <- tab02.F %>%
    group_by(scenario, year) %>%
    summarise(F = mean(value, na.rm = TRUE),
              F_FMSY = mean(F_FMSY, na.rm = TRUE))

#####@> Figures...

####@> SSB...
p03 <- ggplot() +
    geom_line(data = filter(ssb, Yr %in% 1952:2020),
              aes(x = Yr, y = replist1, colour = "SS3",
                  linetype = "SS3"),
              linewidth = 1) +
    geom_line(data = tab02.SSBm,
              aes(x = year, y = SSB, colour = "openMSE",
                  linetype = "openMSE"),
              linewidth = 1) +
    scale_colour_manual(values =
                            c("SS3" = "red", "openMSE" = "black")) +
    scale_linetype_manual(values =
                              c("SS3" = "solid", "openMSE" = "dashed"))+
    facet_wrap(~scenario) +
    labs(x = "Year", y = "SSB", colour = "Method", linetype = "Method") +
    my_theme()
p03

####@> SSB_SSBMSY...
p04 <- ggplot() +
    geom_line(data = filter(ssb_ssbmsy, Year %in% 1952:2020),
              aes(x = Year, y = stock, colour = "SS3",
                  linetype = "SS3"),
              linewidth = 1) +
    geom_line(data = tab02.SSBm,
              aes(x = year, y = SSB_SSBMSY, colour = "openMSE",
                  linetype = "openMSE"),
              linewidth = 1) +
    scale_colour_manual(values =
                            c("SS3" = "red", "openMSE" = "black")) +
    scale_linetype_manual(values =
                              c("SS3" = "solid", "openMSE" = "dashed"))+
    facet_wrap(~scenario) +
    labs(x = "Year", y = expression(SSB/SSB[MSY]),
         colour = "Method", linetype = "Method") +
    my_theme()
p04

####@> F_FMSY...
p05 <- ggplot() +
    geom_line(data = filter(ffmsy, Yr %in% 1952:2020),
              aes(x = Yr, y = replist1, colour = "SS3",
                  linetype = "SS3"),
              linewidth = 1) +
    geom_line(data = tab02.Fm,
              aes(x = year, y = F_FMSY, colour = "openMSE",
                  linetype = "openMSE"),
              linewidth = 1) +
    scale_colour_manual(values =
                            c("SS3" = "red", "openMSE" = "black")) +
    scale_linetype_manual(values =
                              c("SS3" = "solid", "openMSE" = "dashed"))+
    facet_wrap(~scenario) +
    labs(x = "Year", y = expression(F/F[MSY]),
         colour = "Method", linetype = "Method") +
    my_theme()
p05

#####@> Saving figures...
ggsave("05_Results/Comp_SSB_Quantities_ver00.png", plot = p03,
       device = "png", units = "cm", w = 35, h = 30, dpi = 300,
       bg = "white")

ggsave("05_Results/Comp_SSB_SSBMSY_Quantities_ver00.png", plot = p04,
       device = "png", units = "cm", w = 35, h = 30, dpi = 300,
       bg = "white")

ggsave("05_Results/Comp_F_FMSY_Quantities_ver00.png", plot = p05,
       device = "png", units = "cm", w = 35, h = 30, dpi = 300,
       bg = "white")

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
