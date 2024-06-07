########################################################################
## Description: Evaluation abundance indices and biomass trajectory...
##
## Maintainer: DatenKraft - SCRS/ICCAT (TT Species Group)
## Author: Rodrigo Sant'Ana
## Created: ter abr 30 23:54:13 2024 (-0300)
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

######@> Loading R packages...
library(dplyr)
library(reshape2)
library(readr)
library(ggplot2)
library(tidyverse)
library(viridis)
library(janitor)
library(readxl)
library(r4ss)
library(patchwork)
library(cowplot)
library(openMSE)
library(SWOMSE)

######@> Path to the OMs...
path01 <- "00_SS/"

######@> Setting up ggplot theme...

#####@> Importing fonts...
extrafont::loadfonts(device = "postscript")

#####@> ggplot theme...
rgb01 <- "black"
rgb02 <- "black"
seta <- grid::arrow(length = grid::unit(0.2, "cm"), type = "open")
seta2 <- grid::arrow(length = grid::unit(0.2, "cm"), type = "open",
                     ends = "both")
my_theme <- function(base_size = 18, base_family = "Helvetica") {
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
        theme(axis.ticks = element_line(colour = rgb01),
              axis.line = element_line(colour = rgb01, linewidth = 0.2),
              axis.text = element_text(colour = rgb02, size = 18),
              axis.title = element_text(size = 18),
              strip.text = element_text(size = 16,
                                        margin = ggplot2::margin(0.3,
                                                                 0.3,
                                                                 0.3,
                                                                 0.3,
                                                                 "cm"),
                                        face = "bold"),
              legend.background = element_blank(),
              legend.key = element_blank(),
              legend.text = element_text(size = 18),
              legend.title = element_text(size = 20),
              panel.background = element_blank(),
              panel.grid = element_line(linetype = "solid",
                                        linewidth = 0.2,
                                        colour = "gray90"),
              plot.background = element_blank(),
              complete = TRUE)
}

#####@> Testing theme...
df <- data.frame(x = rnorm(10), y = rnorm(10),
                 z = rep(c("A", "B"), each = 5))
ggplot(data = df, aes(x = x, y = y, fill = z)) +
    geom_point(pch = 21, size = 5) +
    facet_grid(~z) +
    my_theme()

########################################################################
######@> Loading Datasets...

######@> Loading SS3 results...
SS01 <- SS_output(paste0(path01, "WSKJ_EstRec93_Qnt25_h6"))
SS02 <- SS_output(paste0(path01, "WSKJ_EstRec93_Qnt50_h6"))
SS03 <- SS_output(paste0(path01, "WSKJ_EstRec93_Qnt75_h6"))
SS04 <- SS_output(paste0(path01, "WSKJ_EstRec93_Qnt25_h7"))
SS05 <- SS_output(paste0(path01, "WSKJ_EstRec93_Qnt50_h7"))
SS06 <- SS_output(paste0(path01, "WSKJ_EstRec93_Qnt75_h7"))
SS07 <- SS_output(paste0(path01, "WSKJ_EstRec93_Qnt25_h8"))
SS08 <- SS_output(paste0(path01, "WSKJ_EstRec93_Qnt50_h8"))
SS09 <- SS_output(paste0(path01, "WSKJ_EstRec93_Qnt75_h8"))

######@> Loading T1NC from ICCAT...

#####@> Download file...
## link <- "https://iccat.int/Data/t1nc_20240131.7z"
## download.file(link, destfile = "05_Basement/T1NC.7z", method = "wget")

#####@> Decompress file...
## zip_file <- "05_Basement/T1NC.7z"
## output_dir <- "05_Basement/"
## command <- sprintf("7z x '%s' -o'%s'", zip_file, output_dir)
## system(command)

#####@> Loading...
t1nc <- read_excel("05_Basement/t1nc-20240131_ALL.xlsx",
                   sheet = "dsT1NC_v1", range = "A4:W96299") %>%
    data.frame

########################################################################
######@> Exploring possibilities...

######@> Extracting trajectories quantities...

#####@> SS01..
timeseries1 <- SS01$Kobe
timeseries1 <- timeseries1[-nrow(timeseries1),]
deriv_quants <- SS01$derived_quants
fmsy <- deriv_quants[deriv_quants$Label == "annF_MSY", 2]
timeseries1$F.Fmsy_rev <- timeseries1$F.Fmsy / fmsy
timeseries1$cenario <- rep("WSKJ_EstRec93_Qnt25_h6",
                           length(timeseries1[, 1]))

#####@> SS02..
timeseries2 <- SS02$Kobe
timeseries2 <- timeseries2[-nrow(timeseries2),]
deriv_quants <- SS02$derived_quants
fmsy <- deriv_quants[deriv_quants$Label == "annF_MSY", 2]
timeseries2$F.Fmsy_rev <- timeseries2$F.Fmsy / fmsy
timeseries2$cenario <- rep("WSKJ_EstRec93_Qnt50_h6",
                          length(timeseries2[, 1]))

#####@> SS03..
timeseries3 <- SS03$Kobe
timeseries3 <- timeseries3[-nrow(timeseries3),]
deriv_quants <- SS03$derived_quants
fmsy <- deriv_quants[deriv_quants$Label == "annF_MSY", 2]
timeseries3$F.Fmsy_rev <- timeseries3$F.Fmsy / fmsy
timeseries3$cenario <- rep("WSKJ_EstRec93_Qnt75_h6",
                          length(timeseries3[, 1]))

#####@> SS04..
timeseries4 <- SS04$Kobe
timeseries4 <- timeseries4[-nrow(timeseries4),]
deriv_quants <- SS01$derived_quants
fmsy <- deriv_quants[deriv_quants$Label == "annF_MSY", 2]
timeseries4$F.Fmsy_rev <- timeseries4$F.Fmsy / fmsy
timeseries4$cenario <- rep("WSKJ_EstRec93_Qnt25_h7",
                          length(timeseries4[, 1]))

#####@> SS05..
timeseries5 <- SS05$Kobe
timeseries5 <- timeseries5[-nrow(timeseries5),]
deriv_quants <- SS05$derived_quants
fmsy <- deriv_quants[deriv_quants$Label == "annF_MSY", 2]
timeseries5$F.Fmsy_rev <- timeseries5$F.Fmsy / fmsy
timeseries5$cenario <- rep("WSKJ_EstRec93_Qnt50_h7",
                          length(timeseries5[, 1]))

#####@> SS06..
timeseries6 <- SS06$Kobe
timeseries6 <- timeseries6[-nrow(timeseries6),]
deriv_quants <- SS06$derived_quants
fmsy <- deriv_quants[deriv_quants$Label == "annF_MSY", 2]
timeseries6$F.Fmsy_rev <- timeseries6$F.Fmsy / fmsy
timeseries6$cenario <- rep("WSKJ_EstRec93_Qnt75_h7",
                          length(timeseries6[, 1]))

#####@> SS07..
timeseries7 <- SS07$Kobe
timeseries7 <- timeseries7[-nrow(timeseries7),]
deriv_quants <- SS07$derived_quants
fmsy <- deriv_quants[deriv_quants$Label == "annF_MSY", 2]
timeseries7$F.Fmsy_rev <- timeseries7$F.Fmsy / fmsy
timeseries7$cenario <- rep("WSKJ_EstRec93_Qnt25_h8",
                          length(timeseries7[, 1]))

#####@> SS08..
timeseries8 <- SS08$Kobe
timeseries8 <- timeseries8[-nrow(timeseries8),]
deriv_quants <- SS08$derived_quants
fmsy <- deriv_quants[deriv_quants$Label == "annF_MSY", 2]
timeseries8$F.Fmsy_rev <- timeseries8$F.Fmsy / fmsy
timeseries8$cenario <- rep("WSKJ_EstRec93_Qnt50_h8",
                          length(timeseries8[, 1]))

#####@> SS09..
timeseries9 <- SS09$Kobe
timeseries9 <- timeseries9[-nrow(timeseries9),]
deriv_quants <- SS09$derived_quants
fmsy <- deriv_quants[deriv_quants$Label == "annF_MSY", 2]
timeseries9$F.Fmsy_rev <- timeseries9$F.Fmsy / fmsy
timeseries9$cenario <- rep("WSKJ_EstRec93_Qnt75_h8",
                          length(timeseries9[, 1]))

#####@> Merging everything...
tsFinal <- rbind(timeseries1, timeseries2, timeseries3, timeseries4,
                 timeseries5, timeseries6, timeseries7, timeseries8,
                 timeseries9) %>%
    mutate(Text = "Stock Assessment")

#####@> Median distribution...
temp <- tsFinal %>%
    group_by(Yr) %>%
    summarise(BBmsy = median(B.Bmsy)) %>%
    mutate(Text = "Stock Assessment")

####@> Visualizing...
p00 <- ggplot(data = tsFinal,
              aes(x = Yr, y = B.Bmsy, colour = cenario)) +
    geom_line(linewidth = 1.2) +
    geom_line(data = temp, aes(x = Yr, y = BBmsy),
              colour = "black", linewidth = 1.5) +
    facet_wrap(~Text, strip.position = "right") +
    scale_colour_grey(start = 0.7, end = 0.7) +
    scale_y_continuous(limits = c(0, 8), expand = c(0, 0)) +
    labs(x = "Year", y = expression(B/B[MSY])) +
    my_theme() +
    theme(legend.position = "none")
p00

######@> Exploring indices used in SS3 models...

#####@> SS01 (all scenarios must be equal)...
tsCPUE <- SS01$cpue %>%
    dplyr::select(Fleet_name, Yr, Obs, SE) %>%
    mutate(Fleet_name = factor(Fleet_name,
        levels = c("BRA_BB_hist", "BB_West", "PS_West",
            "LL_USMX", "HL_RR")))

####@> Visualizing...
p01 <- ggplot(data = tsCPUE,
              aes(x = Yr, y = Obs)) +
    facet_wrap(~Fleet_name, ncol = 3) +
    geom_line(linewidth = 1.2) +
    labs(x = "Year", y = "Scaled abundance index") +
    scale_y_continuous(limits = c(0, 3), expand = c(0, 0)) + 
    my_theme()
p01

ggsave("04_Results/Fig_Indices_SA_ver00.png",
    plot = p01, device = "png", units = "cm", w = 30, h = 25,
    dpi = 300, bg = "white")

#####@> Merging Brazilian BB indices...
tmp00 <- tsCPUE %>%
    filter(Fleet_name %in% c("BRA_BB_hist", "BB_West")) %>%
    mutate(Fleet = "BB_BRA")

#####@> Venezuelan PS indices...
tmp01 <- tsCPUE %>%
    filter(Fleet_name %in% c("PS_West")) %>%
    mutate(Fleet = "PS_VEN")

#####@> United States LL indices...
tmp02 <- tsCPUE %>%
    filter(Fleet_name %in% c("LL_USMX")) %>%
    mutate(Fleet = "LL_USA")

#####@> Averaging Brazilian, Venezuelan and USA...
tmp03 <- tsCPUE %>%
    filter(Fleet_name %in%
           c("BRA_BB_hist", "BB_West", "PS_West", "LL_USMX")) %>%
    group_by(Yr) %>%
    summarise(Obs = mean(Obs, na.rm = TRUE)) %>%
    mutate(Fleet = "Averaging Index 01")

#####@> Weighted average by Brazilian, Venezuelan and USA...
tmp04 <- tsCPUE %>%
    filter(Fleet_name %in%
           c("BRA_BB_hist", "BB_West", "PS_West", "LL_USMX")) %>%
    mutate(w = case_when(
               Fleet_name %in% c("BRA_BB_hist", "BB_West") ~ 2,
               Fleet_name %in% "PS_West" ~ 1,
               Fleet_name %in% "LL_USMX" ~ 1
           )) %>%
    group_by(Yr) %>%
    summarise(Obs = weighted.mean(Obs, w, na.rm = TRUE)) %>%
    mutate(Fleet = "Averaging Index 02")

#####@> Inverse variance weighted mean for the Brazilian, Venezuelan and
######> USA indices...
temp00 <- tsCPUE %>%
    filter(Fleet_name %in%
           c("BRA_BB_hist", "BB_West", "PS_West", "LL_USMX")) %>%
    mutate(New_fleet = case_when(
               Fleet_name %in% c("BRA_BB_hist", "BB_West") ~ "BRA_BB",
               Fleet_name %in% "PS_West" ~ "VEN_PS",
               Fleet_name %in% "LL_USMX" ~ "USMX_LL"
           )) %>%
    select(Yr:New_fleet) %>%
    mutate(w = 1/SE) %>%
    mutate(mu_w = Obs * w)

####@> Dump method 01...
tmp05 <- temp00 %>%
    group_by(Yr) %>%
    summarise(A = sum(mu_w), B = sum(w)) %>%
    mutate(Obs = A/B) %>%
    mutate(Fleet = "Averaging Index 03")

####@> Method 02...
tmp05 <- temp00 %>%
    group_by(Yr) %>%
    summarise(Obs = weighted.mean(Obs, 1/SE, na.rm = TRUE))

####@> Visualizing...
p02 <- ggplot() +
    geom_line(data = tmp00, aes(x = Yr, y = Obs), linewidth = 1.2) +
    geom_line(data = tmp01, aes(x = Yr, y = Obs), linewidth = 1.2) +
    geom_line(data = tmp02, aes(x = Yr, y = Obs), linewidth = 1.2) +
    geom_line(data = tmp03, aes(x = Yr, y = Obs), linewidth = 1.2) +
    geom_line(data = tmp04, aes(x = Yr, y = Obs), linewidth = 1.2) +
    geom_line(data = tmp05, aes(x = Yr, y = Obs), linewidth = 1.2) +
    facet_wrap(~Fleet, strip.position = "right", ncol = 1) +
    scale_y_continuous(limits = c(0, 3), expand = c(0, 0)) +
    scale_x_continuous(limits = c(1952, 2020)) +
    labs(x = "Year", y = "Scaled abundance index") +
    my_theme()
p02

####@> Exporting figure...
plot00 <- p00 | p02
ggsave("04_Results/Fig_Comp_Indices_Trajectory_ver01.png",
       plot = plot00, device = "png", units = "cm", w = 45, h = 35,
       dpi = 300, bg = "white")

######@> Comparing trajectories at the same figure...
p03 <- ggplot() +
    geom_line(data = temp, aes(x = Yr, y = BBmsy),
              colour = "black", linewidth = 1.5) +
    geom_line(data = tmp00, aes(x = Yr, y = Obs), linewidth = 1.2) +
    scale_y_continuous(limits = c(0, 6), expand = c(0, 0)) +
    scale_x_continuous(limits = c(1952, 2020)) +
    labs(x = "Year", y = "Trajectory vs BB_BRA Index") +
    my_theme()
p03

p04 <- ggplot() +
    geom_line(data = temp, aes(x = Yr, y = BBmsy),
              colour = "black", linewidth = 1.5) +
    geom_line(data = tmp01, aes(x = Yr, y = Obs), linewidth = 1.2) +
    scale_y_continuous(limits = c(0, 6), expand = c(0, 0)) +
    scale_x_continuous(limits = c(1952, 2020)) +
    labs(x = "Year", y = "Trajectory vs PS_VEN Index") +
    my_theme()
p04

p05 <- ggplot() +
    geom_line(data = temp, aes(x = Yr, y = BBmsy),
              colour = "black", linewidth = 1.5) +
    geom_line(data = tmp02, aes(x = Yr, y = Obs), linewidth = 1.2) +
    scale_y_continuous(limits = c(0, 6), expand = c(0, 0)) +
    scale_x_continuous(limits = c(1952, 2020)) +
    labs(x = "Year", y = "Trajectory vs LL_USA Index") +
    my_theme()
p05

p06 <- ggplot() +
    geom_line(data = temp, aes(x = Yr, y = BBmsy),
              colour = "black", linewidth = 1.5) +
    geom_line(data = tmp03, aes(x = Yr, y = Obs), linewidth = 1.2) +
    scale_y_continuous(limits = c(0, 6), expand = c(0, 0)) +
    scale_x_continuous(limits = c(1952, 2020)) +
    labs(x = "Year", y = "Trajectory vs Averaging Index 01") +
    my_theme()
p06

p07 <- ggplot() +
    geom_line(data = temp, aes(x = Yr, y = BBmsy),
              colour = "black", linewidth = 1.5) +
    geom_line(data = tmp04, aes(x = Yr, y = Obs), linewidth = 1.2) +
    scale_y_continuous(limits = c(0, 6), expand = c(0, 0)) +
    scale_x_continuous(limits = c(1952, 2020)) +
    labs(x = "Year", y = "Trajectory vs Averaging Index 02") +
    my_theme()
p07

p08 <- ggplot() +
    geom_line(data = temp, aes(x = Yr, y = BBmsy),
        colour = "black", linewidth = 1.5) +
    geom_line(data = tmp05, aes(x = Yr, y = Obs), linewidth = 1.2) +
    scale_y_continuous(limits = c(0, 6), expand = c(0, 0)) +
    scale_x_continuous(limits = c(1952, 2020)) +
    labs(x = "Year", y = "Trajectory vs Averaging Index 03") +
    my_theme()
p08

####@> Exporting figure...
plot01 <- (p03 / p04 / p05) | (p06 / p07 / p08)
ggsave("04_Results/Fig_Comp_Indices_Trajectory_02_ver01.png",
       plot = plot01, device = "png", units = "cm", w = 45, h = 35,
       dpi = 300, bg = "white")

######@> Working with T1NC...

#####@> Filtering W-SKJ data only...
wskj <- t1nc %>%
    filter(Species == "SKJ",
           Stock == "ATW")

#####@> Estimating proportions per FlagName...
tab01 <- wskj %>%
    mutate(flag = ifelse(FlagName == "Brazil", "Brazil", "Other CPCs")) %>%
    group_by(YearC, flag) %>%
    summarise(Catch = sum(Qty_t, na.rm = TRUE)) %>%
    data.frame

####@> Visualizing...
pal <- c("#003f5c", "#2f4b7c", "#665191", "#a05195", "#d45087", "#f95d6a",
         "#ff7c43", "#ffa600")
p08 <- ggplot(data = tab01, aes(x = YearC, y = Catch, fill = flag)) +
    geom_area(colour = "black", position = "stack") +
    labs(x = "Year", y = "Metric tonnes", fill = "") +
    scale_x_continuous(limits = c(1950, 2022),
                       breaks = seq(1950, 2022, 10)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 50000)) +
    scale_fill_manual(values = pal[c(3, 5)]) +
    coord_cartesian(xlim = c(1950, 2020)) +
    my_theme() +
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.position = c(0.15, 0.98))
p08

p09 <- ggplot(data = filter(tab01, YearC >= 2000),
              aes(x = YearC, y = Catch, fill = flag)) +
    geom_area(colour = "black", position = "fill") +
    labs(x = "Year", y = "", fill = "") +
    scale_x_continuous(breaks = seq(2000, 2022, 10),
                       limits = c(2000, 2022),
                       expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), position = "right") +
    scale_fill_manual(values = pal[c(3, 5)]) +
    coord_cartesian(xlim = c(2000, 2022)) +
    my_theme() +
    theme(legend.position = "none",
          plot.margin = margin(t = 10, r = 30, b = 10, l = 10))
p09

####@> Combining figures...
plot02 <- plot_grid(p08, p09, rel_widths = c(1, 0.2))
ggsave("04_Results/Fig_WSKJ_Catches_TS_ver00.png",
       plot = plot02, device = "png", units = "cm", w = 45, h = 25,
       dpi = 300, bg = "white")

########################################################################
######@> Working with different indices inside OMs...

######@> Path to the file...
SSdir <- paste0(path01, "WSKJ_EstRec93_Qnt50_h7")

######@> Path to the output files...
output.dir <- "05_Basement/"

######@> W-SKJ EstRec93 Qnt50_h7 - [Method SS2OM]...
OM5 <- SS2OM(
    SSdir,
    nsim = 100,
    proyears = 30,
    reps = 100,
    maxF = 3,
    seed = 1,
    interval = 3,
    pstar = 0.5,
    Obs = MSEtool::Precise_Unbiased,
    Imp = MSEtool::Perfect_Imp,
    gender = 1,
    silent = FALSE,
    Name = "OM005",
    Source = "ICCAT 2022 SKJ assessment",
    Author = "Rodrigo Sant'Ana & Bruno Leite Mourato",
    report = FALSE,
    filename = "OM_WSKJ_EstRec93_Qnt50_h7",
    dir = output.dir,
    open_file = TRUE)

#####@> Selecting and looking for the combining Brazilian BB indices...

####@> Importing data from SS...
Data <- SS2Data(SSdir,
    Name = "OM005",
    Common_Name = "Skipjack",
    Species = "Katsuwonus pelamis",
    Region = "Western Atlantic Ocean")

####@> Preparing the index combined vector...
index <- vctrs::vec_c(Data@AddInd[1, "BRA_BB_hist", 1:48],
    Data@AddInd[1, "BB_West", 49:69])
cv_index <- vctrs::vec_c(Data@CV_AddInd[1, 5, 1:48],
    Data@CV_AddInd[1, 2, 49:69])
names(cv_index) <- Data@Year

####@> Populate the index slot in Data object...
Data@Ind <- array(index, dim = c(1, length(Data@Year)))
Data@CV_Ind <- array(cv_index, dim = c(1, length(Data@Year)))

####@> Adding Data in cpars data compartment...
OM5.BB_BRA <- OM5
OM5.BB_BRA@cpars$Data <- Data
OM5.BB_BRA@Name <- "OM005 BB BRA"
OM5.BB_BRA@cpars$Data@Name <- "OM005 BB BRA"

####@> Looking to the OM...
plot_SS2OM(OM5.BB_BRA,
           SSdir,
           gender = 1,
           filename = "OM005 BB BRA",
           dir = output.dir,
           open_file = TRUE,
           silent = FALSE)

#####@> Selecting and looking for the PS VEN indices...

####@> Importing data from SS...
Data <- SS2Data(SSdir,
    Name = "OM005",
    Common_Name = "Skipjack",
    Species = "Katsuwonus pelamis",
    Region = "Western Atlantic Ocean")

####@> Preparing the index combined vector...
index <- Data@AddInd[1, "PS_West", ]
cv_index <- Data@CV_AddInd[1, 1, ]
names(cv_index) <- Data@Year

####@> Populate the index slot in Data object...
Data@Ind <- array(index, dim = c(1, length(Data@Year)))
Data@CV_Ind <- array(cv_index, dim = c(1, length(Data@Year)))

####@> Adding Data in cpars data compartment...
OM5.PS_VEN <- OM5
OM5.PS_VEN@cpars$Data <- Data
OM5.PS_VEN@Name <- "OM005 PS VEN"
OM5.PS_VEN@cpars$Data@Name <- "OM005 PS VEN"

####@> Looking to the OM...
plot_SS2OM(OM5.PS_VEN,
    SSdir,
    gender = 1,
    filename = "OM005 PS VEN",
    dir = output.dir,
    open_file = TRUE,
    silent = FALSE)

#####@> Selecting and looking for the LL USA indices...

####@> Importing data from SS...
Data <- SS2Data(SSdir,
    Name = "OM005",
    Common_Name = "Skipjack",
    Species = "Katsuwonus pelamis",
    Region = "Western Atlantic Ocean")

####@> Preparing the index combined vector...
index <- Data@AddInd[1, "LL_USMX", ]
cv_index <- Data@CV_AddInd[1, 3, ]
names(cv_index) <- Data@Year

####@> Populate the index slot in Data object...
Data@Ind <- array(index, dim = c(1, length(Data@Year)))
Data@CV_Ind <- array(cv_index, dim = c(1, length(Data@Year)))

####@> Adding Data in cpars data compartment...
OM5.LL_US <- OM5
OM5.LL_US@cpars$Data <- Data
OM5.LL_US@Name <- "OM005 LL US"
OM5.LL_US@cpars$Data@Name <- "OM005 LL US"

####@> Looking to the OM...
plot_SS2OM(OM5.LL_US,
    SSdir,
    gender = 1,
    filename = "OM005 LL US",
    dir = output.dir,
    open_file = TRUE,
    silent = FALSE)

#####@> Selecting and looking for the Averaging 01 indices...

####@> Importing data from SS...
Data <- SS2Data(SSdir,
    Name = "OM005",
    Common_Name = "Skipjack",
    Species = "Katsuwonus pelamis",
    Region = "Western Atlantic Ocean")

####@> Preparing the index combined vector...
index <- apply(Data@AddInd[1, c(1:3, 5), ], 2, mean, na.rm = TRUE)
cv_index <- apply(Data@CV_AddInd[1, c(1:3, 5), ], 2, mean, na.rm = TRUE)
names(cv_index) <- Data@Year
index[is.na(index)] <- NA
cv_index[is.na(cv_index)] <- NA

####@> Populate the index slot in Data object...
Data@Ind <- array(index, dim = c(1, length(Data@Year)))
Data@CV_Ind <- array(cv_index, dim = c(1, length(Data@Year)))

####@> Adding Data in cpars data compartment...
OM5.AV01 <- OM5
OM5.AV01@cpars$Data <- Data
OM5.AV01@Name <- "OM005 AV01"
OM5.AV01@cpars$Data@Name <- "OM005 AV01"

####@> Looking to the OM...
plot_SS2OM(OM5.AV01,
    SSdir,
    gender = 1,
    filename = "OM005 AV01",
    dir = output.dir,
    open_file = TRUE,
    silent = FALSE)

#####@> Selecting and looking for the Averaging 02 indices...

####@> Importing data from SS...
Data <- SS2Data(SSdir,
    Name = "OM005",
    Common_Name = "Skipjack",
    Species = "Katsuwonus pelamis",
    Region = "Western Atlantic Ocean")

####@> Preparing the index combined vector...
index <- apply(Data@AddInd[1, c(1:3, 5), ], 2,
    function(x) weighted.mean(x, w = c(1, 2, 1, 2), na.rm = TRUE))
cv_index <- apply(Data@CV_AddInd[1, c(1:3, 5), ], 2,
    function(x) weighted.mean(x, w = c(1, 2, 1, 2), na.rm = TRUE))
names(cv_index) <- Data@Year
index[is.na(index)] <- NA
cv_index[is.na(cv_index)] <- NA

####@> Populate the index slot in Data object...
Data@Ind <- array(index, dim = c(1, length(Data@Year)))
Data@CV_Ind <- array(cv_index, dim = c(1, length(Data@Year)))

####@> Adding Data in cpars data compartment...
OM5.AV02 <- OM5
OM5.AV02@cpars$Data <- Data
OM5.AV02@Name <- "OM005 AV02"
OM5.AV02@cpars$Data@Name <- "OM005 AV02"

####@> Looking to the OM...
plot_SS2OM(OM5.AV02,
    SSdir,
    gender = 1,
    filename = "OM005 AV02",
    dir = output.dir,
    open_file = TRUE,
    silent = FALSE)

#####@> Selecting and looking for the Averaging 03 indices...

####@> Importing data from SS...
Data <- SS2Data(SSdir,
    Name = "OM005",
    Common_Name = "Skipjack",
    Species = "Katsuwonus pelamis",
    Region = "Western Atlantic Ocean")

####@> Preparing the index combined vector...
index <- c(rep(NA, 29), as.numeric(tmp05$Obs)[1:40])
names(index) <- Data@Year
cv_index <- c(rep(NA, 29), rep(0.2, 40))
names(cv_index) <- Data@Year
## index[is.na(index)] <- NA
## cv_index[is.na(cv_index)] <- NA

####@> Populate the index slot in Data object...
Data@Ind <- array(index, dim = c(1, length(Data@Year)))
Data@CV_Ind <- array(cv_index, dim = c(1, length(Data@Year)))

####@> Adding Data in cpars data compartment...
OM5.AV03 <- OM5
OM5.AV03@cpars$Data <- Data
OM5.AV03@Name <- "OM005 AV03"
OM5.AV03@cpars$Data@Name <- "OM005 AV03"

####@> Looking to the OM...
plot_SS2OM(OM5.AV03,
    SSdir,
    gender = 1,
    filename = "OM005 AV03",
    dir = output.dir,
    open_file = TRUE,
    silent = FALSE)

######@> Re-evaluating the stock assessment to look which of these
######@> possibilities tracks better the old assessment...

######@> BB BRA...
m0 <- SP_SS(
    x = 1, Data = OM5.BB_BRA@cpars$Data,
    prior = list(r = c(0.416, 0.148), MSY = c(33000, 0.2)),
    start = list(dep = 0.98, n = 1),
    fix_sigma = FALSE, fix_tau = TRUE)

out0 <- data.frame(
    Yr = 1952:2020,
    B.Bmsy = m0@B_BMSY[1:69],
    F.Fmsy = m0@F_FMSY)

p10 <- ggplot() +
    geom_line(data = filter(tsFinal,
        cenario == "WSKJ_EstRec93_Qnt50_h7"),
        aes(x = Yr, y = B.Bmsy), colour = "black", linewidth = 1) +
    geom_line(data = out0,
        aes(x = Yr, y = B.Bmsy), colour = "purple", linewidth = 1) +
    scale_y_continuous(limits = c(0, 6), expand = c(0, 0)) +
    labs(x = "Year", y = expression(B/B[MSY])) +
    my_theme()
p10

p11 <- ggplot() +
    geom_line(data = filter(tsFinal,
        cenario == "WSKJ_EstRec93_Qnt50_h7"),
        aes(x = Yr, y = F.Fmsy), colour = "black", linewidth = 1) +
    geom_line(data = out0,
        aes(x = Yr, y = F.Fmsy), colour = "purple", linewidth = 1) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = "Year", y = expression(F/F[MSY])) +
    my_theme()
p11

######@> PS VEN...
m1 <- SP_SS(
    x = 1, Data = OM5.PS_VEN@cpars$Data,
    prior = list(r = c(0.416, 0.148), MSY = c(33000, 0.2)),
    start = list(dep = 0.98, n = 1),
    fix_sigma = FALSE, fix_tau = TRUE)

out1 <- data.frame(
    Yr = 1952:2020,
    B.Bmsy = m1@B_BMSY[1:69],
    F.Fmsy = m1@F_FMSY)

p12 <- ggplot() +
    geom_line(data = filter(tsFinal,
        cenario == "WSKJ_EstRec93_Qnt50_h7"),
        aes(x = Yr, y = B.Bmsy), colour = "black", linewidth = 1) +
    geom_line(data = out1,
        aes(x = Yr, y = B.Bmsy), colour = "purple", linewidth = 1) +
    scale_y_continuous(limits = c(0, 6), expand = c(0, 0)) +
    labs(x = "Year", y = expression(B/B[MSY])) +
    my_theme()
p12

p13 <- ggplot() +
    geom_line(data = filter(tsFinal,
        cenario == "WSKJ_EstRec93_Qnt50_h7"),
        aes(x = Yr, y = F.Fmsy), colour = "black", linewidth = 1) +
    geom_line(data = out1,
        aes(x = Yr, y = F.Fmsy), colour = "purple", linewidth = 1) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = "Year", y = expression(F/F[MSY])) +
    my_theme()
p13

######@> LL US...
m2 <- SP_SS(
    x = 1, Data = OM5.LL_US@cpars$Data,
    prior = list(r = c(0.416, 0.148), MSY = c(33000, 0.2)),
    start = list(dep = 0.98, n = 1),
    fix_sigma = FALSE, fix_tau = TRUE)

out2 <- data.frame(
    Yr = 1952:2020,
    B.Bmsy = m2@B_BMSY[1:69],
    F.Fmsy = m2@F_FMSY)

p14 <- ggplot() +
    geom_line(data = filter(tsFinal,
        cenario == "WSKJ_EstRec93_Qnt50_h7"),
        aes(x = Yr, y = B.Bmsy), colour = "black", linewidth = 1) +
    geom_line(data = out2,
        aes(x = Yr, y = B.Bmsy), colour = "purple", linewidth = 1) +
    scale_y_continuous(limits = c(0, 6), expand = c(0, 0)) +
    labs(x = "Year", y = expression(B/B[MSY])) +
    my_theme()
p14

p15 <- ggplot() +
    geom_line(data = filter(tsFinal,
        cenario == "WSKJ_EstRec93_Qnt50_h7"),
        aes(x = Yr, y = F.Fmsy), colour = "black", linewidth = 1) +
    geom_line(data = out2,
        aes(x = Yr, y = F.Fmsy), colour = "purple", linewidth = 1) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = "Year", y = expression(F/F[MSY])) +
    my_theme()
p15

######@> AV01...
m3 <- SP_SS(
    x = 1, Data = OM5.AV01@cpars$Data,
    prior = list(r = c(0.416, 0.148), MSY = c(33000, 0.2)),
    start = list(dep = 0.98, n = 1),
    fix_sigma = FALSE, fix_tau = TRUE)

out3 <- data.frame(
    Yr = 1952:2020,
    B.Bmsy = m3@B_BMSY[1:69],
    F.Fmsy = m3@F_FMSY)

p16 <- ggplot() +
    geom_line(data = filter(tsFinal,
        cenario == "WSKJ_EstRec93_Qnt50_h7"),
        aes(x = Yr, y = B.Bmsy), colour = "black", linewidth = 1) +
    geom_line(data = out3,
        aes(x = Yr, y = B.Bmsy), colour = "purple", linewidth = 1) +
    scale_y_continuous(limits = c(0, 6), expand = c(0, 0)) +
    labs(x = "Year", y = expression(B/B[MSY])) +
    my_theme()
p16

p17 <- ggplot() +
    geom_line(data = filter(tsFinal,
        cenario == "WSKJ_EstRec93_Qnt50_h7"),
        aes(x = Yr, y = F.Fmsy), colour = "black", linewidth = 1) +
    geom_line(data = out3,
        aes(x = Yr, y = F.Fmsy), colour = "purple", linewidth = 1) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = "Year", y = expression(F/F[MSY])) +
    my_theme()
p17

######@> AV02...
m4 <- SP_SS(
    x = 1, Data = OM5.AV02@cpars$Data,
    prior = list(r = c(0.416, 0.148), MSY = c(33000, 0.2)),
    start = list(dep = 0.98, n = 1),
    fix_sigma = FALSE, fix_tau = TRUE)

out4 <- data.frame(
    Yr = 1952:2020,
    B.Bmsy = m4@B_BMSY[1:69],
    F.Fmsy = m4@F_FMSY)

p18 <- ggplot() +
    geom_line(data = filter(tsFinal,
        cenario == "WSKJ_EstRec93_Qnt50_h7"),
        aes(x = Yr, y = B.Bmsy), colour = "black", linewidth = 1) +
    geom_line(data = out4,
        aes(x = Yr, y = B.Bmsy), colour = "purple", linewidth = 1) +
    scale_y_continuous(limits = c(0, 6), expand = c(0, 0)) +
    labs(x = "Year", y = expression(B/B[MSY])) +
    my_theme()
p18

p19 <- ggplot() +
    geom_line(data = filter(tsFinal,
        cenario == "WSKJ_EstRec93_Qnt50_h7"),
        aes(x = Yr, y = F.Fmsy), colour = "black", linewidth = 1) +
    geom_line(data = out4,
        aes(x = Yr, y = F.Fmsy), colour = "purple", linewidth = 1) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = "Year", y = expression(F/F[MSY])) +
    my_theme()
p19

######@> AV03...
m5 <- SP_SS(
    x = 1, Data = OM5.AV03@cpars$Data,
    prior = list(r = c(0.416, 0.148), MSY = c(33000, 0.2)),
    start = list(dep = 0.98, n = 1),
    fix_sigma = FALSE, fix_tau = TRUE)

out5 <- data.frame(
    Yr = 1952:2020,
    B.Bmsy = m5@B_BMSY[1:69],
    F.Fmsy = m5@F_FMSY)

p20 <- ggplot() +
    geom_line(data = filter(tsFinal,
        cenario == "WSKJ_EstRec93_Qnt50_h7"),
        aes(x = Yr, y = B.Bmsy), colour = "black", linewidth = 1) +
    geom_line(data = out5,
        aes(x = Yr, y = B.Bmsy), colour = "purple", linewidth = 1) +
    scale_y_continuous(limits = c(0, 6), expand = c(0, 0)) +
    labs(x = "Year", y = expression(B/B[MSY])) +
    my_theme()
p20

p21 <- ggplot() +
    geom_line(data = filter(tsFinal,
        cenario == "WSKJ_EstRec93_Qnt50_h7"),
        aes(x = Yr, y = F.Fmsy), colour = "black", linewidth = 1) +
    geom_line(data = out5,
        aes(x = Yr, y = F.Fmsy), colour = "purple", linewidth = 1) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = "Year", y = expression(F/F[MSY])) +
    my_theme()
p21

######@> Exporting figures...
plot03 <- p10/p11
ggsave("05_Basement/Fig_BB_BRA_ver00.png", plot = plot03,
    device = "png", units = "cm", w = 25, h = 25, dpi = 600,
    bg = "white")

plot04 <- p12/p13
ggsave("05_Basement/Fig_PS_VEN_ver00.png", plot = plot04,
    device = "png", units = "cm", w = 25, h = 25, dpi = 600,
    bg = "white")

plot05 <- p14/p15
ggsave("05_Basement/Fig_LL_US_ver00.png", plot = plot05,
    device = "png", units = "cm", w = 25, h = 25, dpi = 600,
    bg = "white")

plot06 <- p16/p17
ggsave("05_Basement/Fig_AV01_ver00.png", plot = plot06,
    device = "png", units = "cm", w = 25, h = 25, dpi = 600,
    bg = "white")

plot07 <- p18/p19
ggsave("05_Basement/Fig_AV02_ver00.png", plot = plot07,
    device = "png", units = "cm", w = 25, h = 25, dpi = 600,
    bg = "white")

plot08 <- p20/p21
ggsave("05_Basement/Fig_AV03_ver00.png", plot = plot08,
    device = "png", units = "cm", w = 25, h = 25, dpi = 600,
    bg = "white")

########################################################################
######@> Exporting Indices Scenarios...

######@> Concatenating final table...
tsIndex <- combine(tmp00, tmp01, tmp02, tmp03, tmp04,
    select(tmp05, Yr, Obs, Fleet))

######@> Exporting...
save(tsIndex, file = "05_Basement/tsIndex_ver01.RData")

######@> Exporting workspace...
## save.image("05_Basement/Index_Evaluation.RData")

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
