########################################################################
## Description: Preparing updated abundance indices...
##
## Maintainer: DatenKraft - SCRS/ICCAT (TT Species Group)
## Author: Rodrigo Sant'Ana
## Created: dom ago 04 01:36:13 2024 (-0300)
## Version: 0.0.1
##
## URL:
## Doc URL:
##
## Database info: Updated abundance index - YFT SA July 2024
## [SKJEW_CPUE2022_8July2024.xlsx]
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
path01 <- "01_SS/"

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

######@> Função para extrair textos...
extract_obs_se <- function(fleet_name) {
    ## Usando sub() para identificar e extrair 'Obs' ou 'SE' do final da
    ## string
    result <- sub(".*_(Obs|SE)$", "\\1", fleet_name)
    return(result)
}

######@> Função para remover '_Obs' ou '_SE' do final da string...
remove_obs_se <- function(fleet_name) {
    ## Usando sub() para substituir '_Obs' ou '_SE' por uma string vazia
    result <- sub("_(Obs|SE)$", "", fleet_name)
    return(result)
}

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

######@> Loading SS3 results using Ai's method...

#####@> Listing files in SS3 directory...
listF <- list.files(path01)

#####@> Extracting labels...
Labels <- listF

#####@> Number of scenarios in the grid...
Ngrid <- 9

#####@> Looping for extracting values...
GRIDSummary <- NULL
for(i in 1:Ngrid) {
    GRIDSummary[[i]] <- SSsummarize(
        SSgetoutput(keyvec = NULL,
                    dirvec = paste0(path01, listF[i]),
                    getcovar = FALSE,
                    getcomp = FALSE,
                    forecast = FALSE),
        SpawnOutputUnits = "biomass")$quants
}

#####@> Extracting values...

## ssb <- GRIDSummary[[1]][grep("SSB_[0-9]", GRIDSummary[[1]]$Label),]
## ssbmsy <- GRIDSummary[[1]][grep("SSB_MSY", GRIDSummary[[1]]$Label),]
## stock <- ssb[, 1]/ssbmsy[, 1]
## harvest <- GRIDSummary[[1]][grep("")]

####@> SSB...
ssb <- lapply(seq_along(GRIDSummary), function(i) {
    GRIDSummary[[i]] %>%
        filter(grepl("SSB_[0-9]", Label)) %>%
        mutate(scenario = listF[[i]])
})
ssb <- do.call("rbind.data.frame", ssb)

####@> SSB_MSY...
ssbmsy <- lapply(seq_along(GRIDSummary), function(i) {
    GRIDSummary[[i]] %>%
        filter(grepl("SSB_MSY", Label)) %>%
        mutate(scenario = listF[[i]])
})
ssbmsy <- do.call("rbind.data.frame", ssbmsy)

####@> F_FMSY...
ffmsy <- lapply(seq_along(GRIDSummary), function(i) {
    GRIDSummary[[i]] %>%
        filter(grepl("F_[0-9]", Label)) %>%
        mutate(scenario = listF[[i]])
})
ffmsy <- do.call("rbind.data.frame", ffmsy)

####@> FMSY...
fmsy <- lapply(seq_along(GRIDSummary), function(i) {
    GRIDSummary[[i]] %>%
        filter(grepl("annF_MSY", Label)) %>%
        mutate(scenario = listF[[i]])
})
fmsy <- do.call("rbind.data.frame", fmsy)

####@> MSY...
msy <- lapply(seq_along(GRIDSummary), function(i) {
    GRIDSummary[[i]] %>%
        filter(grepl("Dead_Catch_MSY", Label)) %>%
        mutate(scenario = listF[[i]])
})
msy <- do.call("rbind.data.frame", msy)

####@> SSB_SSBMSY...
ssb_ssbmsy <- ssb %>%
    left_join(ssbmsy, by = "scenario") %>%
    mutate(stock = replist1.x/replist1.y) %>%
    select("Year" = Yr.x, scenario, "SSB" = replist1.x,
           "SSB_MSY" = replist1.y, stock)

######@> Loading T1NC from ICCAT...

#####@> Download file...
## link <- "https://iccat.int/Data/t1nc_20240131.7z"
## download.file(link, destfile = "00_Data/T1NC.7z", method = "wget")

#####@> Decompress file...
## zip_file <- "00_Data/T1NC.7z"
## output_dir <- "00_Data/"
## command <- sprintf("7z x '%s' -o'%s'", zip_file, output_dir)
## system(command)

#####@> Loading...
t1nc <- read_excel("00_Data/t1nc-20240131_ALL.xlsx",
                   sheet = "dsT1NC_v1", range = "A4:W96299") %>%
    data.frame

######@> Loading updated indices..
tsCPUEnew <- read_excel("00_Data/SKJEW_CPUE2022_8July2024.xlsx",
                        sheet = 5, range = "A6:M49") %>%
    select(Year,
           "BB_West_Obs" = "Index...2", "BB_West_SE" = "SD",
           "HL_RR_Obs" = "Index...4", "HL_RR_SE" = "SE",
           "US_GOM_Obs" = "Index...6", "US_GOM_SE" = "CV...7",
           "LL_USMX_Obs" = "Index...8", "LL_USMX_SE" = "CV...9",
           "PS_West_Obs" = "Index...10", "PS_West_SE" = "CV...11",
           "BRA_BB_hist_Obs" = "Index...12", "BRA_BB_hist_SE" =
                                                 "CV...13") %>%
    mutate(Text = "Updated Indices 2024") %>%
    data.frame
tsCPUEnew$BRA_BB_hist_Obs <- tsCPUEnew$BRA_BB_hist_Obs /
    mean(tsCPUEnew$BRA_BB_hist_Obs, na.rm = TRUE)
tsCPUEnew$PS_West_Obs <- tsCPUEnew$PS_West_Obs /
    mean(tsCPUEnew$PS_West_Obs, na.rm = TRUE)

######@> Loading determinist results from Ai's folder...

#####@> B/Bmsy...
detBratio <- read_excel("00_Data/KobeMatrices_final_forExeSum_WSKJ.xlsx",
                        sheet = 2, range = "F5:O74")

#####@> F/Fmsy...
detFratio <- read_excel("00_Data/KobeMatrices_final_forExeSum_WSKJ.xlsx",
                        sheet = 2, range = "P5:X74")

######@> Loading MVLN results from Ai's folder...

#####@> B/Bmsy...
mvlBratio <- read_excel("00_Data/KobeMatrices_final_forExeSum_WSKJ.xlsx",
                        sheet = 2, range = "Z5:AC74")

#####@> F/Fmsy...
mvlFratio <- read_excel("00_Data/KobeMatrices_final_forExeSum_WSKJ.xlsx",
                        sheet = 2, range = "AD5:AF74")

########################################################################
######@> Exploring possibilities...

######@> Extracting trajectories quantities...

#####@> SS01..
timeseries1 <- data.frame(Year = 1952:2020,
                          B.Bmsy = SS01$Kobe$B.Bmsy[2:nrow(SS01$Kobe)],
                          F.Fmsy = SS01$Kobe$F.Fmsy[1:nrow(SS01$Kobe)-1])
timeseries1$cenario <- rep("WSKJ_EstRec93_Qnt25_h6",
                           length(timeseries1[, 1]))

#####@> SS02..
timeseries2 <- data.frame(Year = 1952:2020,
                          B.Bmsy = SS02$Kobe$B.Bmsy[2:nrow(SS02$Kobe)],
                          F.Fmsy = SS02$Kobe$F.Fmsy[1:nrow(SS02$Kobe) - 1])
timeseries2$cenario <- rep("WSKJ_EstRec93_Qnt50_h6",
                          length(timeseries2[, 1]))

#####@> SS03..
timeseries3 <- data.frame(Year = 1952:2020,
                          B.Bmsy = SS03$Kobe$B.Bmsy[2:nrow(SS03$Kobe)],
                          F.Fmsy = SS03$Kobe$F.Fmsy[1:nrow(SS03$Kobe) - 1])
timeseries3$cenario <- rep("WSKJ_EstRec93_Qnt75_h6",
                          length(timeseries3[, 1]))

#####@> SS04..
timeseries4 <- data.frame(Year = 1952:2020,
                          B.Bmsy = SS04$Kobe$B.Bmsy[2:nrow(SS04$Kobe)],
                          F.Fmsy = SS04$Kobe$F.Fmsy[1:nrow(SS04$Kobe) - 1])
timeseries4$cenario <- rep("WSKJ_EstRec93_Qnt25_h7",
                          length(timeseries4[, 1]))

#####@> SS05..
timeseries5 <- data.frame(Year = 1952:2020,
                          B.Bmsy = SS05$Kobe$B.Bmsy[2:nrow(SS05$Kobe)],
                          F.Fmsy = SS05$Kobe$F.Fmsy[1:nrow(SS05$Kobe) - 1])
timeseries5$cenario <- rep("WSKJ_EstRec93_Qnt50_h7",
                          length(timeseries5[, 1]))

#####@> SS06..
timeseries6 <- data.frame(Year = 1952:2020,
                          B.Bmsy = SS06$Kobe$B.Bmsy[2:nrow(SS06$Kobe)],
                          F.Fmsy = SS06$Kobe$F.Fmsy[1:nrow(SS06$Kobe) - 1])
timeseries6$cenario <- rep("WSKJ_EstRec93_Qnt75_h7",
                          length(timeseries6[, 1]))

#####@> SS07..
timeseries7 <- data.frame(Year = 1952:2020,
                          B.Bmsy = SS07$Kobe$B.Bmsy[2:nrow(SS07$Kobe)],
                          F.Fmsy = SS07$Kobe$F.Fmsy[1:nrow(SS07$Kobe) - 1])
timeseries7$cenario <- rep("WSKJ_EstRec93_Qnt25_h8",
                          length(timeseries7[, 1]))

#####@> SS08..
timeseries8 <- data.frame(Year = 1952:2020,
                          B.Bmsy = SS08$Kobe$B.Bmsy[2:nrow(SS08$Kobe)],
                          F.Fmsy = SS08$Kobe$F.Fmsy[1:nrow(SS08$Kobe) - 1])
timeseries8$cenario <- rep("WSKJ_EstRec93_Qnt50_h8",
                          length(timeseries8[, 1]))

#####@> SS09..
timeseries9 <- data.frame(Year = 1952:2020,
                          B.Bmsy = SS09$Kobe$B.Bmsy[2:nrow(SS09$Kobe)],
                          F.Fmsy = SS09$Kobe$F.Fmsy[1:nrow(SS09$Kobe) - 1])
timeseries9$cenario <- rep("WSKJ_EstRec93_Qnt75_h8",
                          length(timeseries9[, 1]))

#####@> Merging everything...
tsFinal01 <- rbind(timeseries1, timeseries2, timeseries3, timeseries4,
                   timeseries5, timeseries6, timeseries7, timeseries8,
                   timeseries9) %>%
    mutate(Text = "Method 01",
           cenario2 = gsub("_", "", str_sub(cenario, 15, 22)))

######@> Merging results from Ai's method...
tsFinal02 <- ssb_ssbmsy %>%
    mutate(Year = Year - 1) %>%
    filter(Year %in% 1952:2020)
tmp00 <- ffmsy %>%
    filter(Yr %in% 1952:2020) %>%
    select("Year" = Yr, scenario, "harvest" = replist1)

tsFinal02 <- tsFinal02 %>%
    left_join(tmp00, by = c("scenario", "Year")) %>%
    mutate(Text = "Method 02",
           cenario2 = gsub("_", "", str_sub(scenario, 15, 22)))

########################################################################
######@> Contributions to the total catches...

######@> Estimating the contributions to the total catches...

#####@> Filtering the SKJ data...
skjw <- t1nc %>%
    filter(Species == "SKJ",
           Stock == "ATW")

#####@> Table with proportions by CPC...
tab01 <- skjw %>%
    group_by(YearC) %>%
    summarise(Total = sum(Qty_t, na.rm = TRUE))
temp00 <- skjw %>%
    filter(FleetCode %in% c("BRA", "VEN", "USA")) %>%
    group_by(YearC, FleetCode) %>%
    summarise(Fleet = sum(Qty_t, na.rm = TRUE)) %>%
    pivot_wider(names_from = "FleetCode", values_from = "Fleet") %>%
    mutate_if(is.numeric, ~replace_na(.x, 0))
tab01 <- tab01 %>%
    left_join(temp00, by = "YearC") %>%
    mutate(pBRA = BRA/Total,
           pUSA = USA/Total,
           pVEN = VEN/Total) %>%
    mutate_if(is.numeric, ~replace_na(.x, 0)) %>%
    select(YearC, pBRA:pVEN)

########################################################################
######@> Comparing the deterministic runs from Ai's excel, Ai's SS3
######@> method of importing and I...

######@> For B/Bmsy...

#####@> Converting Ai's structure...
tsDetAi <- detBratio %>%
    pivot_longer(names_to = "cenario2", values_to = "B.Bmsy",
                 2:ncol(.)) %>%
    mutate(Text = "ICCAT Consolidation")

####@> Visualizing...
p00 <- ggplot() +
    geom_line(data = tsFinal01,
              aes(x = Year, y = B.Bmsy, colour = cenario2),
              linewidth = 1.2) +
    geom_line(data = tsFinal02,
              aes(x = Year, y = stock, colour = cenario2),
              linewidth = 1.2) +
    geom_line(data = tsDetAi,
              aes(x = Year, y = B.Bmsy, colour = cenario2),
              linewidth = 1.2) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
    facet_wrap(~Text, strip.position = "right") +
    scale_colour_brewer(palette = 2) +
    scale_y_continuous(limits = c(0, 8), breaks = 0:8, expand = c(0, 0)) +
    labs(x = "Year", y = expression(B/B[MSY]), colour = "") +
    my_theme() +
    theme(legend.position = "bottom")
p00

p01 <- ggplot() +
    geom_line(data = tsFinal01,
              aes(x = Year, y = B.Bmsy, colour = Text),
              linewidth = 1) +
    geom_line(data = tsFinal02,
              aes(x = Year, y = stock, colour = Text),
              linewidth = 1, linetype = "longdash") +
    geom_line(data = tsDetAi,
              aes(x = Year, y = B.Bmsy, colour = Text),
              linewidth = 2, linetype = "dotted") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
    facet_wrap(~cenario2, strip.position = "top") +
    scale_colour_brewer(palette = "Dark2", direction = -1) +
    scale_y_continuous(limits = c(0, 8), breaks = 0:8, expand = c(0, 0)) +
    labs(x = "Year", y = expression(B/B[MSY]), colour = "") +
    my_theme() +
    theme(legend.position = "bottom")
p01

######@> For F/Fmsy...

#####@> Converting Ai's structure...
tsDetAi <- detFratio %>%
    mutate(Year = detBratio$Year) %>%
    pivot_longer(names_to = "cenario2", values_to = "F.Fmsy",
                 1:(ncol(.) - 1)) %>%
    mutate(Text = "ICCAT Consolidation")

####@> Visualizing...
p00a <- ggplot() +
    geom_line(data = tsFinal01,
              aes(x = Year, y = F.Fmsy, colour = cenario2),
              linewidth = 1.2) +
    geom_line(data = tsFinal02,
              aes(x = Year, y = harvest, colour = cenario2),
              linewidth = 1.2) +
    geom_line(data = tsDetAi,
              aes(x = Year, y = F.Fmsy, colour = cenario2),
              linewidth = 1.2) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
    facet_wrap(~Text, strip.position = "right") +
    scale_colour_brewer(palette = 2) +
    scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.5),
                       expand = c(0, 0)) +
    labs(x = "Year", y = expression(F/F[MSY]), colour = "") +
    my_theme() +
    theme(legend.position = "bottom")
p00a

p01b <- ggplot() +
    geom_line(data = tsFinal01,
              aes(x = Year, y = F.Fmsy, colour = Text),
              linewidth = 1) +
    geom_line(data = tsFinal02,
              aes(x = Year, y = harvest, colour = Text),
              linewidth = 1, linetype = "longdash") +
    geom_line(data = tsDetAi,
              aes(x = Year, y = F.Fmsy, colour = Text),
              linewidth = 2, linetype = "dotted") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
    facet_wrap(~cenario2, strip.position = "top") +
    scale_colour_brewer(palette = "Dark2", direction = -1) +
    scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.5),
                       expand = c(0, 0)) +
    labs(x = "Year", y = expression(F/F[MSY]), colour = "") +
    my_theme() +
    theme(legend.position = "bottom")
p01b

######@> Exploring indices used in SS3 models...

#####@> SS01 (all scenarios must be equal)...
tsCPUE <- SS01$cpue %>%
    dplyr::select(Fleet_name, "Year" = Yr, Obs, SE) %>%
    mutate(Fleet_name = factor(Fleet_name,
                               levels = c("BRA_BB_hist", "BB_West",
                                          "PS_West", "LL_USMX",
                                          "HL_RR")),
           Text = "Stock Assessment 2022")

######@> Pivoting tsCPUEnew...
tsCPUEnew <- tsCPUEnew %>%
    pivot_longer(names_to = "Fleet_name", values_to = "Valor", 2:13) %>%
    mutate(Tipo = extract_obs_se(Fleet_name),
           Fleet_name = remove_obs_se(Fleet_name)) %>%
    pivot_wider(names_from = "Tipo", values_from = "Valor")

####@> Visualizing...
p02 <- ggplot() +
    geom_line(data = tsCPUE,
              aes(x = Year, y = Obs), linewidth = 0.5) +
    geom_line(data = tsCPUEnew,
              aes(x = Year, y = Obs), linewidth = 0.3, col = "purple") +
    facet_wrap(~Fleet_name, ncol = 3) +
    labs(x = "Year", y = "Scaled abundance index") +
    ## scale_y_continuous(limits = c(0, 3), expand = c(0, 0)) +
    my_theme()
p02

#####@> Merging Brazilian BB indices...
tmp00 <- tsCPUEnew %>%
    filter(Fleet_name %in% c("BRA_BB_hist", "BB_West")) %>%
    mutate(Fleet = "BB_BRA") %>%
    group_by(Year, Fleet) %>%
    summarise(Obs = mean(Obs, na.rm = TRUE),
              SE = mean(SE, na.rm = TRUE))

#####@> Venezuelan PS indices...
tmp01 <- tsCPUEnew %>%
    filter(Fleet_name %in% c("PS_West")) %>%
    mutate(Fleet = "PS_VEN")

#####@> United States LL indices...
tmp02 <- tsCPUEnew %>%
    filter(Fleet_name %in% c("LL_USMX")) %>%
    mutate(Fleet = "LL_USA")

#####@> Averaging Brazilian, Venezuelan and USA...
tmp03 <- tsCPUEnew %>%
    filter(Fleet_name %in%
           c("BRA_BB_hist", "BB_West", "PS_West", "LL_USMX")) %>%
    group_by(Year) %>%
    summarise(Obs = mean(Obs, na.rm = TRUE)) %>%
    mutate(Fleet = "Simple average")

#####@> Weighted average by Brazilian, Venezuelan and USA... AQUIII
tmp04 <- tsCPUEnew %>%
    filter(Fleet_name %in%
           c("BRA_BB_hist", "BB_West", "PS_West", "LL_USMX")) %>%
    left_join(tab01, by = c("Year" = "YearC")) %>%
    mutate(w = case_when(
               Fleet_name %in% c("BRA_BB_hist", "BB_West") ~ pBRA,
               Fleet_name %in% "PS_West" ~ pVEN,
               Fleet_name %in% "LL_USMX" ~ pUSA
           )) %>%
    group_by(Year) %>%
    summarise(Obs = weighted.mean(Obs, w, na.rm = TRUE)) %>%
    mutate(Fleet = "Average weighted by T1NC")

#####@> Inverse variance weighted mean for the Brazilian, Venezuelan and
######> USA indices...
temp00 <- tsCPUEnew %>%
    filter(Fleet_name %in%
           c("BRA_BB_hist", "BB_West", "PS_West", "LL_USMX")) %>%
    mutate(New_fleet = case_when(
               Fleet_name %in% c("BRA_BB_hist", "BB_West") ~ "BRA_BB",
               Fleet_name %in% "PS_West" ~ "VEN_PS",
               Fleet_name %in% "LL_USMX" ~ "USMX_LL"
           )) %>%
    select(Year:New_fleet) %>%
    mutate(w = 1/SE) %>%
    mutate(mu_w = Obs * w)

####@> Dump method 01...
tmp05 <- temp00 %>%
    group_by(Year) %>%
    summarise(A = sum(mu_w, na.rm = TRUE),
              B = sum(w, na.rm = TRUE)) %>%
    mutate(Obs = A/B) %>%
    mutate(Fleet = "Inverse variance weighted average")

####@> Method 02...
tmp05 <- temp00 %>%
    group_by(Year) %>%
    summarise(Obs = weighted.mean(Obs, 1/SE, na.rm = TRUE)) %>%
    mutate(Fleet = "Inverse variance weighted average")

####@> Visualizing...
p03 <- ggplot() +
    geom_line(data = tmp00, aes(x = Year, y = Obs), linewidth = 0.8) +
    geom_line(data = tmp01, aes(x = Year, y = Obs), linewidth = 0.8) +
    geom_line(data = tmp02, aes(x = Year, y = Obs), linewidth = 0.8) +
    facet_wrap(~Fleet, strip.position = "right", ncol = 1) +
    scale_y_continuous(limits = c(0, 3), expand = c(0, 0)) +
    scale_x_continuous(limits = c(1980, 2024),
                       breaks = seq(1952, 2024, 4)) +
    labs(x = "Year", y = "Scaled abundance index") +
    my_theme()
p03

p04 <- ggplot() +
    geom_line(data = tmp03, aes(x = Year, y = Obs), linewidth = 0.8) +
    geom_line(data = tmp04, aes(x = Year, y = Obs), linewidth = 0.8) +
    geom_line(data = tmp05, aes(x = Year, y = Obs), linewidth = 0.8) +
    facet_wrap(~Fleet, strip.position = "right", ncol = 1) +
    scale_y_continuous(limits = c(0, 3), expand = c(0, 0)) +
    scale_x_continuous(limits = c(1980, 2024),
                       breaks = seq(1952, 2024, 4)) +
    labs(x = "Year", y = "Scaled abundance index") +
    my_theme()
p04

####@> Exporting figure...
plot00 <- p00 / (p03 | p04)
ggsave("05_Results/Fig_Comp_Indices_Trajectory_ver00.png",
       plot = plot00, device = "png", units = "cm", w = 45, h = 55,
       dpi = 300, bg = "white")

ggsave("05_Results/Fig_BBmsy_Trajectories_ver00.png",
       plot = p00, device = "png", units = "cm", w = 35, h = 25,
       dpi = 300, bg = "white")

plot00 <- (p03 | p04)
ggsave("05_Results/Fig_Comp_Indices_Trajectory_02_ver00.png",
       plot = plot00, device = "png", units = "cm", w = 45, h = 55,
       dpi = 300, bg = "white")

######@> Comparing trajectories at the same figure...

#####@> Finding the periods with data overlap both timeseries...
ts00 <- tmp00$Year[!is.na(tmp00$Obs)]
ts01 <- tmp01$Year[!is.na(tmp01$Obs)]
ts02 <- tmp02$Year[!is.na(tmp02$Obs)]
ts03 <- tmp03$Year[!is.na(tmp03$Obs)]
ts04 <- tmp04$Year[!is.na(tmp04$Obs)]
ts05 <- tmp05$Year[!is.na(tmp05$Obs)]

#####@> Preparing deterministic data...
temp <- detBratio %>%
    rowwise() %>%
    mutate(BBmsy = median(c_across(Qnt25h6:Qnt75h8))) %>%
    data.frame

####@> Average overlaped timeseries...
mBBmsy00 <- mean(temp$BBmsy[temp$Year %in% ts00])
mBBmsy01 <- mean(temp$BBmsy[temp$Year %in% ts01])
mBBmsy02 <- mean(temp$BBmsy[temp$Year %in% ts02])
mBBmsy03 <- mean(temp$BBmsy[temp$Year %in% ts03])
mBBmsy04 <- mean(temp$BBmsy[temp$Year %in% ts04])
mBBmsy05 <- mean(temp$BBmsy[temp$Year %in% ts05])

####@> Restimating the Scaled BBmsy...
temp <- temp %>%
    mutate(SBBmsy00 = BBmsy / mBBmsy00,
           SBBmsy01 = BBmsy / mBBmsy01,
           SBBmsy02 = BBmsy / mBBmsy02,
           SBBmsy03 = BBmsy / mBBmsy03,
           SBBmsy04 = BBmsy / mBBmsy04,
           SBBmsy05 = BBmsy / mBBmsy05)

#####@> Preparing index data...
temp2 <- tmp00 %>%
    ungroup() %>%
    mutate(SObs = Obs / mean(tmp00$Obs[tmp00$Year %in% ts00]))
temp3 <- tmp01 %>%
    ungroup() %>%
    mutate(SObs = Obs / mean(tmp01$Obs[tmp01$Year %in% ts01]))
temp4 <- tmp02 %>%
    ungroup() %>%
    mutate(SObs = Obs / mean(tmp02$Obs[tmp02$Year %in% ts02]))
temp5 <- tmp03 %>%
    ungroup() %>%
    mutate(SObs = Obs / mean(tmp03$Obs[tmp03$Year %in% ts03]))
temp6 <- tmp04 %>%
    ungroup() %>%
    mutate(SObs = Obs / mean(tmp04$Obs[tmp04$Year %in% ts04]))
temp7 <- tmp05 %>%
    ungroup() %>%
    mutate(SObs = Obs / mean(tmp05$Obs[tmp05$Year %in% ts05]))

p05 <- ggplot() +
    geom_line(data = temp, aes(x = Year, y = SBBmsy00),
              colour = "black", linewidth = 0.8) +
    geom_line(data = temp2, aes(x = Year, y = SObs),
              linewidth = 1.2, colour = "gray50") +
    scale_y_continuous(limits = c(0, 3), expand = c(0, 0)) +
    scale_x_continuous(limits = c(1950, 2024), ) +
    labs(x = "Year", y = "BB_BRA Index") +
    my_theme()
p05

p06 <- ggplot() +
    geom_line(data = temp, aes(x = Year, y = SBBmsy01),
              colour = "black", linewidth = 0.8) +
    geom_line(data = temp3, aes(x = Year, y = SObs),
              linewidth = 1.2, colour = "gray50") +
    scale_y_continuous(limits = c(0, 3), expand = c(0, 0)) +
    scale_x_continuous(limits = c(1950, 2024)) +
    labs(x = "Year", y = "PS_VEN Index") +
    my_theme()
p06

p07 <- ggplot() +
    geom_line(data = temp, aes(x = Year, y = SBBmsy02),
              colour = "black", linewidth = 0.8) +
    geom_line(data = temp4, aes(x = Year, y = SObs),
              linewidth = 1.2, colour = "gray50") +
    scale_y_continuous(limits = c(0, 3), expand = c(0, 0)) +
    scale_x_continuous(limits = c(1950, 2024)) +
    labs(x = "Year", y = "LL_USA Index") +
    my_theme()
p07

p08 <- ggplot() +
    geom_line(data = temp, aes(x = Year, y = SBBmsy03),
              colour = "black", linewidth = 0.8) +
    geom_line(data = temp5, aes(x = Year, y = SObs),
              linewidth = 1.2, colour = "gray50") +
    scale_y_continuous(limits = c(0, 3), expand = c(0, 0)) +
    scale_x_continuous(limits = c(1950, 2024)) +
    labs(x = "Year", y = "Simple Average Index") +
    my_theme()
p08

p09 <- ggplot() +
    geom_line(data = temp, aes(x = Year, y = SBBmsy04),
              colour = "black", linewidth = 0.8) +
    geom_line(data = temp6, aes(x = Year, y = SObs),
              linewidth = 1.2, colour = "gray50") +
    scale_y_continuous(limits = c(0, 3), expand = c(0, 0)) +
    scale_x_continuous(limits = c(1950, 2024)) +
    labs(x = "Year", y = "Average Weighted by T1NC") +
    my_theme()
p09

p10 <- ggplot() +
    geom_line(data = temp, aes(x = Year, y = SBBmsy05),
        colour = "black", linewidth = 0.8) +
    geom_line(data = temp7, aes(x = Year, y = SObs),
              linewidth = 1.2, colour = "gray50") +
    scale_y_continuous(limits = c(0, 3), expand = c(0, 0)) +
    scale_x_continuous(limits = c(1950, 2024)) +
    labs(x = "Year",
         y = "Inverse Variance Weighted Average") +
    my_theme()
p10

####@> Exporting figure...
plot01 <- (p05 / p06 / p07) | (p08 / p09 / p10)
ggsave("05_Results/Fig_Comp_Scaled_Indices_Trajectory_ver00.png",
       plot = plot01, device = "png", units = "cm", w = 45, h = 35,
       dpi = 300, bg = "white")

######@> Working with T1NC...

#####@> Estimating proportions per FlagName...
tab02 <- skjw %>%
    mutate(flag = ifelse(FlagName == "Brazil", "Brazil", "Other CPCs")) %>%
    group_by(YearC, flag) %>%
    summarise(Catch = sum(Qty_t, na.rm = TRUE)) %>%
    data.frame

####@> Visualizing...
pal <- c("#003f5c", "#2f4b7c", "#665191", "#a05195", "#d45087", "#f95d6a",
         "#ff7c43", "#ffa600")
p11 <- ggplot(data = tab02, aes(x = YearC, y = Catch, fill = flag)) +
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
p11

p12 <- ggplot(data = filter(tab02, YearC >= 2000),
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
p12

####@> Combining figures...
plot02 <- plot_grid(p11, p12, rel_widths = c(1, 0.2))
ggsave("05_Results/Fig_WSKJ_Catches_TS_ver00.png",
       plot = plot02, device = "png", units = "cm", w = 45, h = 25,
       dpi = 300, bg = "white")

########################################################################
######@> Testing OMs imports and compare quantities...

######@> Path to the file (Middle grid scenario)...
SSdir <- paste0(path01, "WSKJ_EstRec93_Qnt50_h7")

######@> Path to the output files...
output.dir <- "06_Basement/"

######@> Importing OM...
OM5 <- SS2OM(SSdir,
             nsim = 100,
             proyears = 30,
             reps = 100,
             maxF = 3,
             seed = 1,
             interval = 1,
             pstar = 0.5,
             Obs = MSEtool::Precise_Unbiased,
             Imp = MSEtool::Perfect_Imp,
             gender = 1,
             silent = FALSE,
             Name = "OM WSKJ_EstRec93_Qnt50_h7",
             Source = "ICCAT 2022 SKJ assessment",
             Author = "Rodrigo Sant'Ana & Bruno Leite Mourato",
             report = FALSE,
             filename = "OM_WSKJ_EstRec93_Qnt50_h7",
             dir = output.dir,
             open_file = TRUE)

#####@> Selecting and looking for the Averaging 01 indices...

####@> Importing data from SS...
Data <- SS2Data(SSdir,
                Name = "OM005",
                Common_Name = "Skipjack",
                Species = "Katsuwonus pelamis",
                Region = "Western Atlantic Ocean")

####@> Preparing the Year data...
Year <- 1952:2020

####@> Preparing the index combined vector...
index <- c(rep(NA, 29), tmp03$Obs[1:40])
names(index) <- 1952:2020
cv_index <- c(rep(NA, 29), rep(0.2, 40))
names(cv_index) <- 1952:2020

####@> Populate the index slot in Data object...
Data@Ind <- array(index, dim = c(1, length(Data@Year)))
Data@CV_Ind <- array(cv_index, dim = c(1, length(Data@Year)))

####@> Adding Data in cpars data compartment...
OM5.AV01 <- OM5
OM5.AV01@cpars$Data <- Data
OM5.AV01@cpars$I_beta <- rep(1, OM5.AV01@nsim)
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

####@> Preparing the Year data...
Year <- 1952:2020

####@> Preparing the index combined vector...
index <- c(rep(NA, 29), tmp04$Obs[1:40])
names(index) <- 1952:2020
cv_index <- c(rep(NA, 29), rep(0.2, 40))
names(cv_index) <- 1952:2020

####@> Populate the index slot in Data object...
Data@Ind <- array(index, dim = c(1, length(Data@Year)))
Data@CV_Ind <- array(cv_index, dim = c(1, length(Data@Year)))

####@> Adding Data in cpars data compartment...
OM5.AV02 <- OM5
OM5.AV02@cpars$Data <- Data
OM5.AV02@cpars$I_beta <- rep(1, OM5.AV02@nsim)
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

####@> Preparing the Year data...
Year <- 1952:2020

####@> Preparing the index combined vector...
index <- c(rep(NA, 29), tmp05$Obs[1:40])
names(index) <- 1952:2020
cv_index <- c(rep(NA, 29), rep(0.2, 40))
names(cv_index) <- 1952:2020

####@> Populate the index slot in Data object...
Data@Ind <- array(index, dim = c(1, length(Data@Year)))
Data@CV_Ind <- array(cv_index, dim = c(1, length(Data@Year)))

####@> Adding Data in cpars data compartment...
OM5.AV03 <- OM5
OM5.AV03@cpars$Data <- Data
OM5.AV03@cpars$I_beta <- rep(1, OM5.AV03@nsim)
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

######@> Simple Average Index...
m0 <- SP_SS(
    x = 1, Data = OM5.AV01@cpars$Data,
    prior = list(r = c(0.416, 0.148), MSY = c(33000, 0.2)),
    start = list(dep = 0.98, n = 1),
    fix_sigma = FALSE, fix_tau = TRUE)

out0 <- data.frame(
    Yr = 1952:2020,
    B.Bmsy = m0@B_BMSY[1:69],
    F.Fmsy = m0@F_FMSY)

p13 <- ggplot() +
    geom_line(data = filter(tsFinal02,
        scenario == "WSKJ_EstRec93_Qnt50_h7"),
        aes(x = Year, y = stock), colour = "black", linewidth = 1) +
    geom_line(data = out0,
        aes(x = Yr, y = B.Bmsy), colour = "purple", linewidth = 1) +
    scale_y_continuous(limits = c(0, 6), expand = c(0, 0)) +
    labs(x = "Year", y = expression(B/B[MSY])) +
    my_theme()
p13

p14 <- ggplot() +
    geom_line(data = filter(tsFinal02,
        scenario == "WSKJ_EstRec93_Qnt50_h7"),
        aes(x = Year, y = harvest), colour = "black", linewidth = 1) +
    geom_line(data = out0,
        aes(x = Yr, y = F.Fmsy), colour = "purple", linewidth = 1) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = "Year", y = expression(F/F[MSY])) +
    my_theme()
p14

######@> Average weighted by T1NC...
m1 <- SP_SS(
    x = 1, Data = OM5.AV02@cpars$Data,
    prior = list(r = c(0.416, 0.148), MSY = c(33000, 0.2)),
    start = list(dep = 0.98, n = 1),
    fix_sigma = FALSE, fix_tau = TRUE)

out1 <- data.frame(
    Yr = 1952:2020,
    B.Bmsy = m1@B_BMSY[1:69],
    F.Fmsy = m1@F_FMSY)

p15 <- ggplot() +
    geom_line(data = filter(tsFinal02,
                            scenario == "WSKJ_EstRec93_Qnt50_h7"),
              aes(x = Year, y = stock), colour = "black", linewidth = 1) +
    geom_line(data = out1,
        aes(x = Yr, y = B.Bmsy), colour = "purple", linewidth = 1) +
    scale_y_continuous(limits = c(0, 6), expand = c(0, 0)) +
    labs(x = "Year", y = expression(B/B[MSY])) +
    my_theme()
p15

p16 <- ggplot() +
    geom_line(data = filter(tsFinal02,
                            scenario == "WSKJ_EstRec93_Qnt50_h7"),
              aes(x = Year, y = harvest), colour = "black", linewidth = 1) +
    geom_line(data = out1,
        aes(x = Yr, y = F.Fmsy), colour = "purple", linewidth = 1) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = "Year", y = expression(F/F[MSY])) +
    my_theme()
p16

######@> Inverse variance weighted average...
m2 <- SP_SS(
    x = 1, Data = OM5.AV03@cpars$Data,
    prior = list(r = c(0.416, 0.148), MSY = c(33000, 0.2)),
    start = list(dep = 0.98, n = 1),
    fix_sigma = FALSE, fix_tau = TRUE)

out2 <- data.frame(
    Yr = 1952:2020,
    B.Bmsy = m2@B_BMSY[1:69],
    F.Fmsy = m2@F_FMSY)

p17 <- ggplot() +
    geom_line(data = filter(tsFinal02,
                            scenario == "WSKJ_EstRec93_Qnt50_h7"),
        aes(x = Year, y = stock), colour = "black", linewidth = 1) +
    geom_line(data = out2,
        aes(x = Yr, y = B.Bmsy), colour = "purple", linewidth = 1) +
    scale_y_continuous(limits = c(0, 6), expand = c(0, 0)) +
    labs(x = "Year", y = expression(B/B[MSY])) +
    my_theme()
p17

p18 <- ggplot() +
    geom_line(data = filter(tsFinal02,
                            scenario == "WSKJ_EstRec93_Qnt50_h7"),
        aes(x = Year, y = harvest), colour = "black", linewidth = 1) +
    geom_line(data = out2,
        aes(x = Yr, y = F.Fmsy), colour = "purple", linewidth = 1) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = "Year", y = expression(F/F[MSY])) +
    my_theme()
p18

######@> Exporting figures...
plot06 <- p13/p14
ggsave("05_Results/Fig_AV01_ver00.png", plot = plot06,
    device = "png", units = "cm", w = 25, h = 25, dpi = 600,
    bg = "white")

plot07 <- p15/p16
ggsave("05_Results/Fig_AV02_ver00.png", plot = plot07,
    device = "png", units = "cm", w = 25, h = 25, dpi = 600,
    bg = "white")

plot08 <- p17/p18
ggsave("05_Results/Fig_AV03_ver00.png", plot = plot08,
    device = "png", units = "cm", w = 25, h = 25, dpi = 600,
    bg = "white")

########################################################################
######@> Exporting Indices Scenarios...

######@> Concatenating final table...
tsIndex <- vctrs::vec_c(
                      dplyr::select(tmp00, Year, Fleet, Obs),
                      dplyr::select(tmp01, Year, Fleet, Obs),
                      dplyr::select(tmp02, Year, Fleet, Obs),
                      dplyr::select(tmp03, Year, Fleet, Obs),
                      dplyr::select(tmp04, Year, Fleet, Obs),
                      dplyr::select(tmp05, Year, Fleet, Obs))

######@> Exporting...
save(tsIndex, file = "05_Results/tsIndex_ver02.RData")

######@> Exporting workspace...
save.image("06_Basement/Index_Evaluation_ver02.RData")

######@> Exporting reference quantities...
save(msy, ssb, ssbmsy, ssb_ssbmsy, fmsy, ffmsy,
     file = "05_Results/Reference_Quantities_ver02.RData")

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
