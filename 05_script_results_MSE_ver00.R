########################################################################
## Description: Compiling Final Results W-SKJ MSE based on the old
## OMs...
##
## Maintainer: Datenkraft - SCRS/ICCAT (TT Species Group)
## Author: Rodrigo Sant'Ana
## Created: dom set 10 21:59:06 2023 (-0300)
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
######@> Loading R packages...

######@> Package list...
library(openMSE)
library(SWOMSE)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggmse)
library(patchwork)
library(viridis)
library(data.table)
library(ggrepel)
library(tidyr)
library(ggridges)

######@> Functions...

#####@> Function to access the quantities for the projections...
get_ts <- function(object,
                   type = c("SSB", "FM", "Catch"),
                   this_year = 0) {
    if (!is(object, "MSE")) {
        stop(
            "`object` must be a MSEtool object of class `mse`",
      "that was created by running `MSEtool::runMSE()`."
      )
    }
    .proj_years <- seq(this_year + 1, this_year + object@proyears)
    years_df <- data.frame(
      year = seq_len(object@proyears), real_year = .proj_years,
      stringsAsFactors = FALSE
    )
    mps <- data.frame(
        mp = seq_along(object@MPs),
        mp_name = object@MPs, stringsAsFactors = FALSE
    )
    ts_data <- list()
    for (i in seq_along(type)) {
        if (type[i] == "C") type[i] <- "Catch" # changed in openMSE
        ts_data[[i]] <- slot(object, type[i]) %>%
            reshape2::melt() %>%
            dplyr::rename(
                       iter = .data$Var1, mp = .data$Var2,
                       value = .data$value, year = .data$Var3
                   ) %>%
            dplyr::left_join(mps, by = "mp") %>%
            dplyr::mutate(Type = type[i])
    }
    ts_data <- dplyr::bind_rows(ts_data)
    ts_data <- dplyr::left_join(ts_data, years_df, by = "year")
    iters <- dim(object@SSB)[1]
    ts_data$iter <- rep(seq_len(iters), length(unique(ts_data$Type))) # parallel messed this up
    ## --------------
    ## historical
    n_hist_years <- dim(object@SSB_hist)[2]
    .hist_years <- seq(this_year - n_hist_years + 1, this_year)
    years_df <- data.frame(
        year = seq_len(n_hist_years), real_year = .hist_years,
        stringsAsFactors = FALSE
    )
    hist_data <- list()
    for (i in seq_along(type)) {
        hist_data[[i]] <- slot(object,
                               ifelse(type[i] == "Catch", "CB_hist",
                                      paste0(type[i], "_hist"))) %>%
            ##apply(c(1, 3), if (type[i] == "FM") max else sum) %>%
            reshape2::melt() %>%
            dplyr::rename(
                       iter = .data$Var1,
                       value = .data$value, year = .data$Var2
                   ) %>%
            dplyr::mutate(Type = type[i])
    }
    hist_data <- dplyr::bind_rows(hist_data)
    hist_data$iter <- rep(seq_len(iters),
                          length(unique(hist_data$Type))) # parallel messed this up
    hist_data <- dplyr::left_join(hist_data, years_df, by = "year")
    hist_data2 <- do.call(
        "rbind",
        replicate(length(mps$mp), hist_data, simplify = FALSE)
    )
    hist_data2[["mp_name"]] <- rep(mps$mp_name, each = nrow(hist_data))
    hist_data2[["mp"]] <- rep(mps$mp, each = nrow(hist_data))
    ts_data <- bind_rows(ts_data, hist_data2)
    ## static ref pts:
    ## dim(object@Misc$MSYRefs$ByYear$SSBMSY)
    all_mps <- tibble(mp = unique(ts_data$mp)[!is.na(unique(ts_data$mp))])
    all_mp_yrs <- expand.grid(mp = all_mps$mp,
                              real_year = sort(unique(hist_data2$real_year)))
    ## object@Misc$MSYRefs$Refs
    ref_ssb_hist <- data.frame(
        ref = object@RefPoint$SSBMSY[, 1, object@nyears],
        iter = seq_len(iters),
        Type = "SSB", stringsAsFactors = FALSE) %>%
        left_join(bind_cols(all_mp_yrs, tibble(Type = rep("SSB", nrow(all_mp_yrs)))), by = "Type")
    ref_f_hist <- data.frame(
        ref = object@RefPoint$FMSY[, 1, object@nyears],
        iter = seq_len(iters), Type = "FM",
        stringsAsFactors = FALSE) %>%
        left_join(bind_cols(all_mp_yrs,
                            tibble(Type = rep("FM", nrow(all_mp_yrs)))), by = "Type")
    all_mp_yrs <- expand.grid(mp = seq_along(mps$mp),
                              real_year = sort(unique(ts_data$real_year)))
    ref_msy <- data.frame(ref = 1, iter = seq_len(iters),
                          Type = "Catch", stringsAsFactors = FALSE) %>%
        left_join(bind_cols(all_mp_yrs,
                            tibble(Type = rep("Catch", nrow(all_mp_yrs)))), by = "Type")
    ## dynamic ref pts for projections:
    x1 <- object@RefPoint$SSBMSY
    x2 <- object@RefPoint$FMSY
    if (length(dim(x1)) < 3) {
        .x1 <- array(NA, dim = c(dim(x1)[1], 1, dim(x1)[2]))
        .x1[,1,] <- x1
        .x2 <- array(NA, dim = c(dim(x2)[1], 1, dim(x2)[2]))
        .x2[,1,] <- x2
    } else {
        .x1 <- x1
        .x2 <- x2
    }
    ref_ssb <- .x1 %>%
        reshape2::melt() %>%
        transmute(iter = .data$Var1,
                  real_year = .data$Var3 + object@OM$CurrentYr[1],
                  mp = .data$Var2,
                  ref = .data$value, Type = "SSB")
    ref_f <- .x2 %>%
        reshape2::melt() %>%
        transmute(iter = .data$Var1,
                  real_year = .data$Var3 +
                      object@OM$CurrentYr[1], mp = .data$Var2,
                  ref = .data$value, Type = "FM")
    refs <- ref_ssb %>% bind_rows(ref_f) %>%
        bind_rows(ref_msy) %>%
        bind_rows(ref_ssb_hist) %>%
        bind_rows(ref_f_hist)
    ts_data <- left_join(ts_data, refs,
                         by = c("iter", "mp", "Type", "real_year")) %>%
        mutate(value = value / ref)
    type[type == "SSB"] <- "B_BMSY"
    type[type == "FM"] <- "F_FMSY"
    type[type == "Catch"] <- "Catch"
    ts_data$Type[ts_data$Type == "FM"] <- "F_FMSY"
    ts_data$Type[ts_data$Type == "SSB"] <- "B_BMSY"
    ts_data$Type[ts_data$Type == "Catch"] <- "Catch"
    ts_data
}

get_ts_quantiles <- function(x, probs = c(0.1, 0.5)) {
    x %>%
        dplyr::group_by(.data$mp_name, .data$real_year, .data$Type) %>%
        dplyr::summarize(
                   median_value = median(.data$value, na.rm = TRUE),
                   u = quantile(.data$value, probs = 1 - probs[2] / 2,
                                na.rm = TRUE),
                   l = quantile(.data$value, probs = probs[2] / 2,
                                na.rm = TRUE),
                   m = quantile(.data$value, probs = 0.50,
                                na.rm = TRUE),
                   uu = quantile(.data$value, probs = 1 - probs[1] / 2,
                                 na.rm = TRUE),
                   ll = quantile(.data$value, probs = probs[1] / 2,
                                 na.rm = TRUE)) %>%
        ungroup()
}

########################################################################
######@> Setup MSE...

######@> Folder destination...
if(dir.exists("04_Results")) {
    print("Folder exists!!")
    folder <- "04_Results/"
} else {
    dir.create("04_Results")
    print("Folder will be created!!")
    folder <- "04_Results/"
}

######@> Path to MSEs...
path00 <- "03_MSEs/"

######@> Path to MPs...
path01 <- "../R-Work/"

######@> Path to PMs...
path02 <- "../R-Work/"

######@> Custom R environment...
options(scipen = 16)
options(digits = 2)

######@> Theme for this project...
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

########################################################################
######@> Loading everything...

######@> Management Procedures...
source(paste0(path01, "02_script_prepare_MPs_ver00.R"))

#####@> Defining MPs...
MPs <- c("NFref", "curE", "CC_15kt", "CC_20kt",
         "CC_25kt", "CC_30kt", "CC_35kt", "CC_40kt",
         "GB_slope", "Iratio", "Islope1",
         "SCA_100_40_SBMSY", "SCA_01", "SCA_02",
         "SP_100_40_SBMSY", "SPSS_100_40_SBMSY",
         "SP_03", "SP_04", "SP_05", "SP_06", "SP_07", "SP_08", "SP_09",
         "SP_10")

######@> Performance Metrics...
source(paste0(path02, "03_script_prepare_PMs_ver00.R"))

#####@> List of PMs...
PMs <- avail('PM')[c(1:23)]

######@> MMSEs model results - No error implementation...
nomes <- paste0("MSE", sprintf("%03d", 1:9))
MSEs <- lapply(1:9, function(i) readRDS(paste0(path00, nomes[i],
                                               "_03.mse")))
names(MSEs) <- paste0("MSE", sprintf("%03d", 1:9))

######@> MMSEs model results - 10% error implementation...
nomes02 <- paste0("MSE", sprintf("%03d", 10:18))
MSEs02 <- lapply(1:9, function(i) readRDS(paste0(path00, nomes02[i],
                                                 "_03.mse")))
names(MSEs02) <- paste0("MSE", sprintf("%03d", 10:18))

######@> MMSEs model results - 20% error implementation...
nomes03 <- paste0("MSE", sprintf("%03d", 19:27))
MSEs03 <- lapply(1:9, function(i) readRDS(paste0(path00, nomes03[i],
                                                 "_03.mse")))
names(MSEs03) <- paste0("MSE", sprintf("%03d", 19:27))

########################################################################
######@> Consolidating Results...

######@>----------------------------------------------------------------
######@> No error implementation - OMs 1-9...

######@> Extracting PMs from all MMSE's objects...

#####@> Listing PMs functions...
PMlist <- list()
for(i in seq_along(PMs)) {
    PMlist[[i]] <- try(get(PMs[i]), silent = TRUE)
}

######@> Extracting values...
PMvalues <- list()
PMvalues2 <- list()
for(i in seq_along(MSEs)) {
    for(j in seq_along(PMs)) {
        pm <- PMlist[[j]](MSEs[[i]])
        om <- paste0("MSE", sprintf("%03d", i))
        mps <- pm@MPs
        val <- data.frame(pm@Prob)
        names(val) <- mps
        val <- val %>%
            mutate(sim = 1:100) %>%
            pivot_longer(names_to = "MP", values_to = "Values", 1:length(MPs))
        nom <- pm@Name
        cap <- pm@Caption
        pm <- PMs[j]
        PMvalues[[j]] <- val %>%
            mutate(OM = om, Name = nom, Caption = cap, PM = pm) %>%
            select(OM, Name, Caption, sim, MP, PM, Values) %>%
            as.data.frame()
    }
    out <- do.call(rbind.data.frame, PMvalues)
    PMvalues2[[i]] <- out
}
PM_output <- do.call(rbind.data.frame, PMvalues2)

#####@> Looking for NAs...
tab01_NA <- PM_output %>%
    filter(is.na(Values)) %>%
    group_by(OM, MP, PM) %>%
    summarise(N = n()) %>%
    as.data.frame()

#####@> Looking to the PMs - Estimating the average for each OMs...
tab01_exp <- PM_output %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SCA_01", "SCA_02",
                     "SP_100_40_SBMSY", "SPSS_100_40_SBMSY", "SP_03",
                     "SP_04", "SP_05", "SP_06", "SP_07", "SP_08",
                     "SP_09", "SP_10")) %>%
    group_by(OM, MP, PM) %>%
    summarise(q50 = mean(Values, na.rm = TRUE)) %>%
    pivot_wider(names_from = "MP", values_from = "q50") %>%
    select(OM:SPSS_100_40_SBMSY, SP_100_40_SBMSY, SP_03:SP_10) %>%
    as.data.frame()

###@> Saving the output table...
tab01_exp_ref <- PM_output %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SCA_01", "SCA_02",
                     "SP_100_40_SBMSY", "SPSS_100_40_SBMSY", "SP_03",
                     "SP_04", "SP_05", "SP_06", "SP_07", "SP_08",
                     "SP_09", "SP_10"))

####@> Exporting table...
write.table(tab01_exp, "04_Results/Table_PM-MP_OMs_001-009_ver00.csv",
            row.names = FALSE, sep = ",", dec = ".")

#####@> Looking to the PMs - Estimating the average between all OMs...
tab01 <- PM_output %>%
    group_by(MP, PM, Name, Caption) %>%
    summarise(q50 = median(Values, na.rm = TRUE)) %>%
    as.data.frame()

####@> Figure PGK ~ POF...
tmp <- tab01 %>%
    filter(PM %in% c("POF", "PGK_short")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SCA_01", "SCA_02",
                     "SP_100_40_SBMSY", "SPSS_100_40_SBMSY", "SP_03",
                     "SP_04", "SP_05", "SP_06", "SP_07", "SP_08",
                     "SP_09", "SP_10")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class = "Status")

p00 <- ggplot(data = tmp, aes(x = PGK_short, y = POF, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 20) +
    labs(x = "Prob. Kobe green quadrante over years 1-3",
         y = "Prob. Overfishing over years 1-30") +
    geom_vline(xintercept = 0.7, linetype = "dashed", colour = "gray") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1)) +
    facet_wrap(~class) +
    my_theme() +
    theme(legend.position = "none")
p00

tmp <- tab01 %>%
    filter(PM %in% c("POF", "PGK_med")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class = "Status")

p01 <- ggplot(data = tmp, aes(x = PGK_med, y = POF, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 50) +
    labs(x = "Prob. Kobe green quadrante over years 4-10",
         y = "Prob. Overfishing over years 1-30") +
    geom_vline(xintercept = 0.7, linetype = "dashed", colour = "gray") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1)) +
    facet_wrap(~class) +
    my_theme() +
    theme(legend.position = "none")
p01

tmp <- tab01 %>%
    filter(PM %in% c("POF", "PGK_long")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class = "Status")

p02 <- ggplot(data = tmp, aes(x = PGK_long, y = POF, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 40) +
    labs(x = "Prob. Kobe green quadrante over years 11-30",
         y = "Prob. Overfishing over years 1-30") +
    geom_vline(xintercept = 0.7, linetype = "dashed", colour = "gray") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1)) +
    facet_wrap(~class) +
    my_theme() +
    theme(legend.position = "none")
p02

tmp <- tab01 %>%
    filter(PM %in% c("POF", "PGK")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class = "Status")

p03 <- ggplot(data = tmp, aes(x = PGK, y = POF, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 40) +
    labs(x = "Prob. Kobe green quadrante over years 1-30",
         y = "Prob. Overfishing over years 1-30") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1)) +
    geom_vline(xintercept = 0.7, linetype = "dashed", colour = "gray") +
    facet_wrap(~class) +
    my_theme() +
    theme(legend.position = "none")
p03

ggsave("Results/Tradeplot_PGK_short-POF_ver01.svg", plot = p00,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Tradeplot_PGK_med-POF_ver01.svg", plot = p01,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Tradeplot_PGK_long-POF_ver01.svg", plot = p02,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Tradeplot_PGK-POF_ver01.svg", plot = p03,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

####@> Figure LRP ~ VarC...
tmp <- tab01 %>%
    filter(PM %in% c("VarCmedium", "LRP_med")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class01 = "Safety", class02 = "Stability")

##@> Correcting SP_05...
tmp$VarCmedium[tmp$VarCmedium > 1] <- 1

p00 <- ggplot(data = tmp, aes(x = VarCmedium, y = LRP_med, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 20) +
    labs(x = "Mean Variation in TAC (%) over years 4-10",
         y = "Prob. SB < 0.4 SBMSY over years 4-10") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1)) +
    geom_hline(yintercept = 0.1, linetype = "dashed", colour = "gray") +
    geom_vline(xintercept = 0.2, linetype = "dashed", colour = "gray") +
    facet_grid(class01~class02) +
    my_theme() +
    theme(legend.position = "none")
p00

tmp <- tab01 %>%
    filter(PM %in% c("VarClong", "LRP_long")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class01 = "Safety", class02 = "Stability")

p01 <- ggplot(data = tmp, aes(x = VarClong, y = LRP_long, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 20) +
    labs(x = "Mean Variation in TAC (%) over years 11-30",
         y = "Prob. SB < 0.4 SBMSY over years 11-30") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1)) +
    geom_hline(yintercept = 0.1, linetype = "dashed", colour = "gray") +
    geom_vline(xintercept = 0.2, linetype = "dashed", colour = "gray") +
    facet_grid(class01~class02) +
    my_theme() +
    theme(legend.position = "none")
p01

tmp <- tab01 %>%
    filter(PM %in% c("VarC", "LRP")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class01 = "Safety", class02 = "Stability")

p02 <- ggplot(data = tmp, aes(x = VarC, y = LRP, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 20) +
    labs(x = "Mean Variation in TAC (%) over years 1-30",
         y = "Prob. SB < 0.4 SBMSY over years 1-30") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1)) +
    geom_hline(yintercept = 0.1, linetype = "dashed", colour = "gray") +
    geom_vline(xintercept = 0.2, linetype = "dashed", colour = "gray") +
    facet_grid(class01~class02) +
    my_theme() +
    theme(legend.position = "none")
p02

ggsave("Results/Tradeplot_VarCmed-LRP_med_ver01.svg", plot = p00,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Tradeplot_VarClong-LRP_long_ver01.svg", plot = p01,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Tradeplot_VarC-LRP_ver01.svg", plot = p02,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

####@> Figure PGK ~ AvC...
tmp <- tab01 %>%
    filter(PM %in% c("AvC_short", "PGK_short")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class = "Status", class2 = "Yield")

p00 <- ggplot(data = tmp, aes(x = PGK_short, y = AvC_short, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 20) +
    labs(x = "Prob. Kobe green quadrante over years 1-3",
         y = "Median Catch (t) over years 1-3") +
    geom_vline(xintercept = 0.7, linetype = "dashed", colour = "gray") +
    scale_y_continuous(limits = c(0, 60000)) +
    scale_x_continuous(limits = c(0, 1)) +
    facet_grid(class2~class) +
    my_theme() +
    theme(legend.position = "none")
p00

tmp <- tab01 %>%
    filter(PM %in% c("AvC_med", "PGK_med")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class = "Status", class2 = "Yield")

p01 <- ggplot(data = tmp, aes(x = PGK_med, y = AvC_med, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 50) +
    labs(x = "Prob. Kobe green quadrante over years 4-10",
         y = "Median Catch (t) over years 4-10") +
    geom_vline(xintercept = 0.7, linetype = "dashed", colour = "gray") +
    scale_y_continuous(limits = c(0, 40000)) +
    scale_x_continuous(limits = c(0, 1)) +
    facet_grid(class2~class) +
    my_theme() +
    theme(legend.position = "none")
p01

tmp <- tab01 %>%
    filter(PM %in% c("AvC_long", "PGK_long")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class = "Status", class2 = "Yield")

p02 <- ggplot(data = tmp, aes(x = PGK_long, y = AvC_long, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 40) +
    labs(x = "Prob. Kobe green quadrante over years 11-30",
         y = "Median Catch (t) over years 11-30") +
    geom_vline(xintercept = 0.7, linetype = "dashed", colour = "gray") +
    scale_y_continuous(limits = c(0, 40000)) +
    scale_x_continuous(limits = c(0, 1)) +
    facet_wrap(~class) +
    my_theme() +
    theme(legend.position = "none")
p02

ggsave("Results/Tradeplot_PGK_short-AvC_short_ver01.svg", plot = p00,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Tradeplot_PGK_med-AvC_med_ver01.svg", plot = p01,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Tradeplot_PGK_long-AvC_long_ver01.svg", plot = p02,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

######@> Consolidating general probabilities for a heatmap table...
ProbTab <- tab01 %>%
    filter(PM %in% c("LRP", "LRP_long", "LRP_med", "LRP_short", "PGK",
                     "PGK_long", "PGK_med", "PGK_short", "PNOF", "POF",
                     "VarC","VarClong", "VarCmedium", "nLRP", "nLRP_long",
                     "nLRP_med", "nLRP_short", "rAvC_long", "rAvC_med",
                     "rAvC_short", "AvC_short", "AvC_med", "AvC_long")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    as.data.frame()

##@> Correcting SP_05...
ProbTab$VarCmedium[ProbTab$VarCmedium > 1] <- 1

#####@> Exporting ProbTab - Combined all 9 MOMs...
write.table(ProbTab, file = "Results/Table_PMs_FULL_ver01.csv",
            row.names = FALSE, sep = ",", dec = ".")

######@> Heatmaps for ProbTab...

#####@> Figure for Yield...
tmp <- ProbTab %>%
    select(MP:AvC_short) %>%
    pivot_longer(names_to = "PM", values_to = "value", 2:4) %>%
    mutate(class = "Yield",
           breaks = cut(value,
                        breaks = quantile(value, seq(0, 1, 0.2)),
                        include.lowest = TRUE, right = FALSE),
           PM = factor(PM, levels = c("AvC_short",
                                      "AvC_med",
                                      "AvC_long")),
           MP = factor(MP, levels = c("CC_20kt", "CC_30kt", "CC_40kt",
                                      "GB_slope", "Iratio", "Islope1",
                                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                                      "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                                      "SP_04", "SP_05", "SP_06")),
           class2 = ifelse(value >= 20000, "Up", "Down")) ## Up average
## 4 years
p00 <- ggplot(data = tmp, aes(x = PM, y = MP, fill = class2)) +
    geom_tile(colour = "black", alpha = 0) +
    geom_text(aes(label = round(value, 0)), size = 5,
              fontface = "bold") +
    scale_fill_brewer(palette = "Greys", direction = 1) +
    facet_wrap(~class) +
    labs(x = "Performance Metrics", y = "Management Procedures") +
    my_theme() +
    theme(legend.position = "none")
p00

#####@> Figure for Status...
tmp <- ProbTab %>%
    select(MP, PGK:PGK_short) %>%
    pivot_longer(names_to = "PM", values_to = "value", 2:5) %>%
    mutate(class = "Status",
           breaks = cut(value,
                        breaks = quantile(value, seq(0, 1, 0.2)),
                        include.lowest = TRUE, right = FALSE),
           PM = factor(PM, levels = c("PGK_short",
                                      "PGK_med",
                                      "PGK_long",
                                      "PGK")),
           MP = factor(MP, levels = c("CC_20kt", "CC_30kt", "CC_40kt",
                                      "GB_slope", "Iratio", "Islope1",
                                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                                      "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                                      "SP_04", "SP_05", "SP_06")),
           class2 = ifelse(value >= 0.7, "Up", "Down"))
p01 <- ggplot(data = tmp, aes(x = PM, y = MP, fill = class2)) +
    geom_tile(colour = "black", alpha = 0.8) +
    geom_text(aes(label = round(value, 3)), size = 5,
              fontface = "bold") +
    scale_fill_brewer(palette = "Greys", direction = 1) +
    facet_wrap(~class) +
    labs(x = "Performance Metrics", y = "Management Procedures") +
    my_theme() +
    theme(legend.position = "none")
p01

tmp <- ProbTab %>%
    select(MP, PNOF:POF) %>%
    pivot_longer(names_to = "PM", values_to = "value", 2:3) %>%
    mutate(class = "Status",
           breaks = cut(value,
                        breaks = quantile(value, seq(0, 1, 0.2)),
                        include.lowest = TRUE, right = FALSE),
           PM = factor(PM, levels = c("PNOF",
                                      "POF")),
           MP = factor(MP, levels = c("CC_20kt", "CC_30kt", "CC_40kt",
                                      "GB_slope", "Iratio", "Islope1",
                                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                                      "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                                      "SP_04", "SP_05", "SP_06")),
           class2 = ifelse(value >= 0.7, "Up", "Down"))
p01b <- ggplot(data = tmp, aes(x = PM, y = MP)) +
    geom_tile(colour = "black", fill = "white", alpha = 0.1) +
    geom_text(aes(label = round(value, 3)), size = 5,
              fontface = "bold") +
    ## scale_fill_brewer(palette = "Greys", direction = 1) +
    facet_wrap(~class) +
    labs(x = "Performance Metrics", y = "Management Procedures") +
    my_theme() +
    theme(legend.position = "none")
p01b

#####@> Figure for Safety...
tmp <- ProbTab %>%
    select(MP, LRP:LRP_short) %>%
    pivot_longer(names_to = "PM", values_to = "value", 2:5) %>%
    mutate(class = "Safety",
           breaks = cut(value,
                        breaks = quantile(value, seq(0, 1, 0.2)),
                        include.lowest = TRUE, right = FALSE),
           PM = factor(PM, levels = c("LRP_short",
                                      "LRP_med",
                                      "LRP_long",
                                      "LRP")),
           MP = factor(MP, levels = c("CC_20kt", "CC_30kt", "CC_40kt",
                                      "GB_slope", "Iratio", "Islope1",
                                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                                      "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                                      "SP_04", "SP_05", "SP_06")),
           class2 = ifelse(value <= 0.1, "Up", "Down"))
p02 <- ggplot(data = tmp, aes(x = PM, y = MP, fill = class2)) +
    geom_tile(colour = "black", alpha = 0.8) +
    geom_text(aes(label = round(value, 3)), size = 5,
              fontface = "bold") +
    scale_fill_brewer(palette = "Greys", direction = 1) +
    facet_wrap(~class) +
    labs(x = "Performance Metrics", y = "Management Procedures") +
    my_theme() +
    theme(legend.position = "none")
p02

#####@> Figure for Stability...
tmp <- ProbTab %>%
    select(MP, VarC:VarCmedium) %>%
    pivot_longer(names_to = "PM", values_to = "value", 2:4) %>%
    mutate(class = "Stability",
           breaks = cut(value,
                        breaks = quantile(value, seq(0, 1, 0.2)),
                        include.lowest = TRUE, right = FALSE),
           PM = factor(PM, levels = c("VarCmedium",
                                      "VarClong",
                                      "VarC")),
           MP = factor(MP, levels = c("CC_20kt", "CC_30kt", "CC_40kt",
                                      "GB_slope", "Iratio", "Islope1",
                                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                                      "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                                      "SP_04", "SP_05", "SP_06")),
           class2 = ifelse(value <= 0.2, "Up", "Down"))

p03 <- ggplot(data = tmp, aes(x = PM, y = MP, fill = class2)) +
    geom_tile(colour = "black", alpha = 0.8) +
    geom_text(aes(label = round(value, 3)), size = 5,
              fontface = "bold") +
    scale_fill_brewer(palette = "Greys", direction = 1) +
    facet_wrap(~class) +
    labs(x = "Performance Metrics", y = "Management Procedures") +
    my_theme() +
    theme(legend.position = "none")
p03

ggsave("Results/Table_PM-MP_Yield_ver03.svg", plot = p00,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 30, h = 25)

ggsave("Results/Table_PM-MP_Status_ver02.svg", plot = p01,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 30, h = 25)

ggsave("Results/Table_PM-MP_Status_B_ver02.svg", plot = p01b,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 30, h = 25)

ggsave("Results/Table_PM-MP_Safety_ver02.svg", plot = p02,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 30, h = 25)

ggsave("Results/Table_PM-MP_Stability_ver02.svg", plot = p03,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 30, h = 25)

######@> Extracting References Points from MMSE objects...

#####@> Extracting SB_SBMSY values...
tmp01 <- list()
for(i in seq_along(MSEs)) {
    om <- names(MSEs)[i]
    output <- MSEs[[i]]@SB_SBMSY %>%
        structure(
            dimnames = list(
                Sim = 1:MSEs[[i]]@nsim,
                MPs = MSEs[[i]]@MPs,
                Year = 2020 + 1:MSEs[[i]]@proyears
            )
        ) %>%
        reshape2::melt() %>%
        mutate(OM = om)
    tmp01[[i]] <- output
}
SB_SBMSY_output <- do.call(rbind.data.frame, tmp01)

####@> Looking to the RP - Estimating the average between all MOMs...
tab02 <- SB_SBMSY_output %>%
    filter(MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                      "GB_slope", "Iratio", "Islope1",
                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                      "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                      "SP_04", "SP_05", "SP_06")) %>%
    group_by(MPs, Year) %>%
    summarise(SB_SBMSY = mean(value, na.rm = TRUE),
              SB_SBMSY_Lwr = quantile(value, probs = 0.025,
                                      na.rm = TRUE),
              SB_SBMSY_Upr = quantile(value, probs = 0.975,
                                      na.rm = TRUE),
              SB_SBMSY_Lwr2 = quantile(value, probs = 0.1,
                                      na.rm = TRUE),
              SB_SBMSY_Upr2 = quantile(value, probs = 0.9,
                                      na.rm = TRUE)) %>%
    as.data.frame()

###@> Figure Reference Points...
p01 <- ggplot(data =
                  filter(tab02,
                         MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr,
                    ymax = SB_SBMSY_Upr),
                alpha = 0.5, fill = "#95B634") +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr2,
                    ymax = SB_SBMSY_Upr2),
                alpha = 0.2, fill = "#006750") +
    geom_line(aes(x = Year, y = SB_SBMSY), linewidth = 1,
              colour = "#006750") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 1) +
    labs(x = "Year", y = expression(SB/SB[MSY])) +
    my_theme()
p01

p02 <- ggplot(data =
                  filter(tab02,
                         MPs %in% c("Iratio", "Islope1", "GB_slope"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr,
                    ymax = SB_SBMSY_Upr),
                alpha = 0.5, fill = "#86BEDA") +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr2,
                    ymax = SB_SBMSY_Upr2),
                alpha = 0.2, fill = "#005690") +
    geom_line(aes(x = Year, y = SB_SBMSY), linewidth = 1,
              colour = "#005690") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 1) +
    labs(x = "Year", y = expression(SB/SB[MSY])) +
    my_theme()
p02

p03 <- ggplot(data =
                  filter(tab02,
                         MPs %in% c("SCA_100_40_SBMSY",
                                    "SP_100_40_SBMSY",
                                    "SPSS_100_40_SBMSY", "SP_01",
                                    "SP_02", "SP_03",
                                    "SP_04", "SP_05", "SP_06"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr,
                    ymax = SB_SBMSY_Upr),
                alpha = 0.5, fill = "#533D8B") +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr2,
                    ymax = SB_SBMSY_Upr2),
                alpha = 0.2, fill = "#5C1374") +
    geom_line(aes(x = Year, y = SB_SBMSY), linewidth = 1,
              colour = "#5C1374") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 2) +
    labs(x = "Year", y = expression(SB/SB[MSY])) +
    my_theme()
p03

ggsave("Results/Trajectory_SB_SBMSY_CC_ver01.svg", plot = p01,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Trajectory_SB_SBMSY_Index_ver01.svg", plot = p02,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Trajectory_SB_SBMSY_Model_ver01.svg", plot = p03,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

#####@> Extracting F_FMSY values...
tmp02 <- list()
for(i in seq_along(MSEs)) {
    om <- names(MSEs)[i]
    output <- MSEs[[i]]@F_FMSY %>%
        structure(
            dimnames = list(
                Sim = 1:MSEs[[i]]@nsim,
                MPs = MSEs[[i]]@MPs,
                Year = 2020 + 1:MSEs[[i]]@proyears
            )
        ) %>%
        reshape2::melt() %>%
        mutate(OM = om) %>%
        as.data.frame()
    tmp02[[i]] <- output
}
F_FMSY_output <- do.call(rbind.data.frame, tmp02)

####@> Looking to the RP - Estimating the average between all MOMs...
tab03 <- F_FMSY_output %>%
    filter(MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                      "GB_slope", "Iratio", "Islope1",
                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                      "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                      "SP_04", "SP_05", "SP_06")) %>%
    group_by(MPs, Year) %>%
    summarise(F_FMSY = mean(value, na.rm = TRUE),
              F_FMSY_Lwr = quantile(value, probs = 0.025,
                                    na.rm = TRUE),
              F_FMSY_Upr = quantile(value, probs = 0.975,
                                    na.rm = TRUE),
              F_FMSY_Lwr2 = quantile(value, probs = 0.1,
                                       na.rm = TRUE),
              F_FMSY_Upr2 = quantile(value, probs = 0.9,
                                       na.rm = TRUE)) %>%
    as.data.frame()

###@> Figure Reference Points...
p01 <- ggplot(data =
                  filter(tab03,
                         MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr,
                    ymax = F_FMSY_Upr),
                alpha = 0.5, fill = "#95B634") +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr2,
                    ymax = F_FMSY_Upr2),
                alpha = 0.2, fill = "#006750") +
    geom_line(aes(x = Year, y = F_FMSY), linewidth = 1,
              colour = "#006750") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 1) +
    labs(x = "Year", y = expression(F/F[MSY])) +
    my_theme()
p01

p02 <- ggplot(data =
                  filter(tab03,
                         MPs %in% c("Iratio", "Islope1", "GB_slope"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr,
                    ymax = F_FMSY_Upr),
                alpha = 0.5, fill = "#86BEDA") +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr2,
                    ymax = F_FMSY_Upr2),
                alpha = 0.2, fill = "#005690") +
    geom_line(aes(x = Year, y = F_FMSY), linewidth = 1,
              colour = "#005690") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 1) +
    labs(x = "Year", y = expression(F/F[MSY])) +
    my_theme()
p02

p03 <- ggplot(data =
                  filter(tab03,
                         MPs %in% c("SCA_100_40_SBMSY",
                                    "SP_100_40_SBMSY",
                                    "SPSS_100_40_SBMSY", "SP_01",
                                    "SP_02", "SP_03",
                                    "SP_04", "SP_05", "SP_06"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr,
                    ymax = F_FMSY_Upr),
                alpha = 0.5, fill = "#533D8B") +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr2,
                    ymax = F_FMSY_Upr2),
                alpha = 0.2, fill = "#5C1374") +
    geom_line(aes(x = Year, y = F_FMSY), linewidth = 1,
              colour = "#5C1374") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 2) +
    labs(x = "Year", y = expression(F/F[MSY])) +
    my_theme()
p03

ggsave("Results/Trajectory_F_FMSY_CC_ver01.svg", plot = p01,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Trajectory_F_FMSY_Index_ver01.svg", plot = p02,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Trajectory_F_FMSY_Model_ver01.svg", plot = p03,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

#####@> Extracting Catches values...
tmp03 <- list()
for(i in seq_along(MSEs)) {
    om <- names(MSEs)[i]
    output <- MSEs[[i]]@Catch %>%
        structure(
            dimnames = list(
                Sim = 1:MSEs[[i]]@nsim,
                MPs = MSEs[[i]]@MPs,
                Year = 2020 + 1:MSEs[[i]]@proyears
            )
        ) %>%
        reshape2::melt() %>%
        mutate(MOM = om) %>%
        as.data.frame()
    tmp03[[i]] <- output
}
Catch_output <- do.call(rbind.data.frame, tmp03)

####@> Looking to the RP - Estimating the average between all MOMs...
tab04 <- Catch_output %>%
    filter(MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                      "GB_slope", "Iratio", "Islope1",
                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                      "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                      "SP_04", "SP_05", "SP_06")) %>%
    group_by(MPs, Year) %>%
    summarise(Catch = mean(value, na.rm = TRUE),
              Catch_Lwr = quantile(value, probs = 0.025,
                                   na.rm = TRUE),
              Catch_Upr = quantile(value, probs = 0.975,
                                   na.rm = TRUE),
              Catch_Lwr2 = quantile(value, probs = 0.1,
                                   na.rm = TRUE),
              Catch_Upr2 = quantile(value, probs = 0.9,
                                   na.rm = TRUE)) %>%
    as.data.frame()

###@> Figure Reference Points...
p01 <- ggplot(data =
                  filter(tab04,
                         MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr,
                    ymax = Catch_Upr),
                alpha = 0.5, fill = "#95B634") +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr2,
                    ymax = Catch_Upr2),
                alpha = 0.2, fill = "#006750") +
    geom_line(aes(x = Year, y = Catch), linewidth = 1,
              colour = "#006750") +
    ## geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 1) +
    labs(x = "Year", y = "Catch (t)") +
    my_theme()
p01

p02 <- ggplot(data =
                  filter(tab04,
                         MPs %in% c("Iratio", "Islope1", "GB_slope"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr,
                    ymax = Catch_Upr),
                alpha = 0.5, fill = "#86BEDA") +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr2,
                    ymax = Catch_Upr2),
                alpha = 0.2, fill = "#005690") +
    geom_line(aes(x = Year, y = Catch), linewidth = 1,
              colour = "#005690") +
    ## geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 1) +
    labs(x = "Year", y = "Catch (t)") +
    my_theme()
p02

p03 <- ggplot(data =
                  filter(tab04,
                         MPs %in% c("SCA_100_40_SBMSY",
                                    "SP_100_40_SBMSY",
                                    "SPSS_100_40_SBMSY", "SP_01",
                                    "SP_02", "SP_03",
                                    "SP_04", "SP_05", "SP_06"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr,
                    ymax = Catch_Upr),
                alpha = 0.5, fill = "#533D8B") +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr2,
                    ymax = Catch_Upr2),
                alpha = 0.2, fill = "#5C1374") +
    geom_line(aes(x = Year, y = Catch), linewidth = 1,
              colour = "#5C1374") +
    ## geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 2) +
    labs(x = "Year", y = "Catch (t)") +
    my_theme()
p03

ggsave("Results/Trajectory_Catch_CC_ver01.svg", plot = p01,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Trajectory_Catch_Index_ver01.svg", plot = p02,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Trajectory_Catch_Model_ver01.svg", plot = p03,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

#####@> Extracting TAC values...
tmp04 <- list()
for(i in seq_along(MSEs)) {
    om <- names(MSEs)[i]
    output <- MSEs[[i]]@TAC %>%
        structure(
            dimnames = list(
                Sim = 1:MSEs[[i]]@nsim,
                MPs = MSEs[[i]]@MPs,
                Year = 2020 + 1:MSEs[[i]]@proyears
            )
        ) %>%
        reshape2::melt() %>%
        mutate(MOM = om) %>%
        as.data.frame()
    tmp04[[i]] <- output
}
TAC_output <- do.call(rbind.data.frame, tmp04)

#####@> Looking for NAs...
tab05_NA <- TAC_output %>%
    filter(is.na(value)) %>%
    group_by(Sim, MOM, MPs) %>%
    summarise(N = n()) %>%
    as.data.frame()

####@> Looking to the RP - Estimating the average between all MOMs...
tab05 <- TAC_output %>%
    filter(MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                      "GB_slope", "Iratio", "Islope1",
                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                      "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                      "SP_04", "SP_05", "SP_06")) %>%
    group_by(MPs, Year) %>%
    summarise(TAC = median(value, na.rm = TRUE),
              TAC_Lwr = quantile(value, probs = 0.025,
                                 na.rm = TRUE),
              TAC_Upr = quantile(value, probs = 0.975,
                                 na.rm = TRUE),
              TAC_Lwr2 = quantile(value, probs = 0.1,
                                 na.rm = TRUE),
              TAC_Upr2 = quantile(value, probs = 0.9,
                                 na.rm = TRUE)) %>%
    as.data.frame()

###@> Figure Reference Points...
p01 <- ggplot(data =
                  filter(tab05,
                         MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = TAC_Lwr,
                    ymax = TAC_Upr),
                alpha = 0.5, fill = "#95B634") +
    geom_ribbon(aes(x = Year,
                    ymin = TAC_Lwr2,
                    ymax = TAC_Upr2),
                alpha = 0.2, fill = "#006750") +
    geom_line(aes(x = Year, y = TAC), linewidth = 1,
              colour = "#006750") +
    ## geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 1) +
    labs(x = "Year", y = "TAC (t)") +
    my_theme()
p01

p02 <- ggplot(data =
                  filter(tab05,
                         MPs %in% c("Iratio", "Islope1", "GB_slope"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = TAC_Lwr,
                    ymax = TAC_Upr),
                alpha = 0.5, fill = "#86BEDA") +
    geom_ribbon(aes(x = Year,
                    ymin = TAC_Lwr2,
                    ymax = TAC_Upr2),
                alpha = 0.2, fill = "#005690") +
    geom_line(aes(x = Year, y = TAC), linewidth = 1,
              colour = "#005690") +
    ## geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 1) +
    labs(x = "Year", y = "TAC (t)") +
    my_theme()
p02

p03 <- ggplot(data =
                  filter(tab05,
                         MPs %in% c("SCA_100_40_SBMSY",
                                    "SP_100_40_SBMSY",
                                    "SPSS_100_40_SBMSY", "SP_01",
                                    "SP_02", "SP_03",
                                    "SP_04", "SP_05", "SP_06"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = TAC_Lwr,
                    ymax = TAC_Upr),
                alpha = 0.5, fill = "#533D8B") +
    geom_ribbon(aes(x = Year,
                    ymin = TAC_Lwr2,
                    ymax = TAC_Upr2),
                alpha = 0.2, fill = "#5C1374") +
    geom_line(aes(x = Year, y = TAC), linewidth = 1,
              colour = "#5C1374") +
    ## geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 2) +
    labs(x = "Year", y = "TAC (t)") +
    my_theme()
p03

ggsave("Results/Trajectory_TAC_CC_ver01.svg", plot = p01,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Trajectory_TAC_Index_ver01.svg", plot = p02,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Trajectory_TAC_Model_ver01.svg", plot = p03,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

######@> Creating the Kobe plots (by Year and general)...

#####@> Merging SB/SBMSY and F/FMSY...
df00 <- SB_SBMSY_output %>%
    left_join(F_FMSY_output, by = c("Sim", "MPs", "Year", "OM")) %>%
    filter(MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                      "GB_slope", "Iratio", "Islope1",
                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                      "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                      "SP_04", "SP_05", "SP_06")) %>%
    rename("SB_SBMSY" = value.x, "F_FMSY" = value.y)

#####@> Estimating Kobe values for the last year...
tmp <- df00 %>%
    filter(Year == 2053) %>%
    group_by(Sim, MPs) %>%
    summarise(SB_SBMSY = mean(SB_SBMSY, na.rm = TRUE),
              F_FMSY = mean(F_FMSY, na.rm = TRUE)) %>%
    mutate(OF = ifelse(F_FMSY <= 1, FALSE, TRUE),
           OFD = ifelse(SB_SBMSY <= 1, TRUE, FALSE))

#####@> Percentage of cases...
valdf <- tmp %>%
    group_by(MPs) %>%
    summarise(BL = sum(OF == FALSE & OFD == TRUE)/100 * 100,
              BR = sum(OF == FALSE & OFD == FALSE)/100 * 100,
              TL = sum(OF == TRUE & OFD == TRUE)/100 * 100,
              TR = sum(OF == TRUE & OFD == FALSE)/100 * 100)
valdf <- valdf %>% tidyr::pivot_longer(., cols = 2:5)
valdf$x <- -Inf
valdf$y <- -Inf
valdf$y[valdf$name == "TL"] <- Inf
valdf$y[valdf$name == "TR"] <- Inf
valdf$x[valdf$name == "BR"] <- Inf
valdf$x[valdf$name == "TR"] <- Inf
valdf$value <- round(valdf$value, 2)
valdf$value <- paste0(valdf$value, "%")
valdf$hjustvar <- -2
valdf$vjustvar <- -2
valdf$hjustvar[valdf$name == "TL"] <- -1
valdf$hjustvar[valdf$name == "TR"] <- 2
valdf$hjustvar[valdf$name == "BL"] <- -1
valdf$hjustvar[valdf$name == "BR"] <- 2
valdf$vjustvar[valdf$name == "TL"] <- 2
valdf$vjustvar[valdf$name == "TR"] <- 2
valdf$vjustvar[valdf$name == "BL"] <- -2
valdf$vjustvar[valdf$name == "BR"] <- -2

#####@> Colors to figure...
kobe_df <- bind_rows(data.frame(x = c(0, 0, 1, 1),
                                y = c(0, 1, 1, 0),
                                fill = "bl"),
                     data.frame(x = c(1, 1, 2, 2),
                                y = c(0, 1, 1, 0),
                                fill = "br"),
                     data.frame(x = c(0, 0, 1, 1),
                                y = c(1, 2, 2, 1),
                                fill = "tl"),
                     data.frame(x = c(1, 1, 2, 2),
                                y = c(1, 2, 2, 1),
                                fill = "tr"))
kobe_df$alpha <- 0.3

#####@> Figure dots...
p00 <- ggplot() +
    geom_polygon(data = kobe_df,
                 aes(x = x, y = y, fill = fill, alpha = alpha)) +
    scale_fill_manual(values = c("#F8DC7A", "#67C18B", "#D8775D",
                                 "#FDBD56")) +
    geom_point(data = tmp,
               aes(x = SB_SBMSY, y = F_FMSY), alpha = 0.5) +
    geom_text(data = valdf, fontface = "bold", size = 4,
              aes(x = x, y = y, label = value,
                  hjust = hjustvar, vjust = vjustvar),
              colour = "black") +
    expand_limits(x = c(0, 2), y = c(0, 2)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 2)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 2)) +
    geom_hline(yintercept = 1, color = "darkgray", linetype = 2) +
    geom_vline(xintercept = 1, color = "darkgray", linetype = 2) +
    labs(x = expression(SB/SB[MSY]), y = expression(F/F[MSY])) +
    facet_wrap(~MPs) +
    my_theme() +
    theme(legend.position = "none")
p00

ggsave("Results/Kobe_plot_Last_Year_2053_ver01.svg", plot = p00,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

######@> Estimating Kobe values per year...
tmp <- df00 %>%
    group_by(Sim, Year, MPs) %>%
    summarise(SB_SBMSY = mean(SB_SBMSY, na.rm = TRUE),
              F_FMSY = mean(F_FMSY, na.rm = TRUE)) %>%
    mutate(OF = ifelse(F_FMSY <= 1, FALSE, TRUE),
           OFD = ifelse(SB_SBMSY <= 1, TRUE, FALSE))

#####@> Proportions by year...
valdf <- tmp %>%
    group_by(Year, MPs) %>%
    summarise(Yellow = sum(OF == FALSE & OFD == TRUE)/100 * 100,
              Green = sum(OF == FALSE & OFD == FALSE)/100 * 100,
              Red = sum(OF == TRUE & OFD == TRUE)/100 * 100,
              Orange = sum(OF == TRUE & OFD == FALSE)/100 * 100) %>%
    pivot_longer(names_to = "Cond", values_to = "Perc", 3:6) %>%
    mutate(Cond = factor(Cond, levels = c("Green", "Yellow", "Orange",
                                          "Red")))

p00 <- ggplot() +
    geom_bar(data = filter(valdf, Year %in% 2021:2053),
             aes(x = Year, y = Perc, fill = Cond),
             stat = "identity", colour = "black",
             width = 1) +
    facet_wrap(~MPs) +
    scale_fill_manual(values = c("#67C18B", "#F8DC7A", "#FDBD56",
                                 "#D8775D")) +
    labs(x = "Year", y = "%") +
    my_theme() +
    theme(legend.position = "none")
p00

p01 <- ggplot() +
    geom_area(data = filter(valdf, Year %in% 2021:2053),
              aes(x = Year, y = Perc, fill = Cond),
              stat = "identity", colour = "black") +
    facet_wrap(~MPs) +
    scale_fill_manual(values = c("#67C18B", "#F8DC7A", "#FDBD56",
                                 "#D8775D")) +
    labs(x = "Year", y = "%") +
    my_theme() +
    theme(legend.position = "none")
p01

ggsave("Results/Kobe_plot_by_year_Bar_ver01.svg", plot = p00,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 30)

ggsave("Results/Kobe_plot_by_year_Area_ver01.svg", plot = p01,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 30)

######@>----------------------------------------------------------------
######@> 10% error implementation - OMs 10-18...

######@> Extracting PMs from all MMSE's objects...

#####@> Listing PMs functions...
PMlist <- list()
for(i in seq_along(PMs)) {
    PMlist[[i]] <- try(get(PMs[i]), silent = TRUE)
}

######@> Extracting values...
PMvalues <- list()
PMvalues2 <- list()
for(i in seq_along(MSEs02)) {
    for(j in seq_along(PMs)) {
        pm <- PMlist[[j]](MSEs02[[i]])
        om <- names(MSEs02[i])
        mps <- pm@MPs
        val <- data.frame(pm@Prob)
        names(val) <- mps
        val <- val %>%
            mutate(sim = 1:100) %>%
            pivot_longer(names_to = "MP", values_to = "Values", 1:20)
        nom <- pm@Name
        cap <- pm@Caption
        pm <- PMs[j]
        PMvalues[[j]] <- val %>%
            mutate(OM = om, Name = nom, Caption = cap, PM = pm) %>%
            select(OM, Name, Caption, sim, MP, PM, Values) %>%
            as.data.frame()
    }
    out <- do.call(rbind.data.frame, PMvalues)
    PMvalues2[[i]] <- out
}
PM_output <- do.call(rbind.data.frame, PMvalues2)

#####@> Looking for NAs...
tab01_NA <- PM_output %>%
    filter(is.na(Values)) %>%
    group_by(OM, MP, PM) %>%
    summarise(N = n()) %>%
    as.data.frame()

#####@> Looking to the PMs - Estimating the average for each OMs...
tab01_exp <- PM_output %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    group_by(OM, MP, PM) %>%
    summarise(q50 = mean(Values, na.rm = TRUE)) %>%
    pivot_wider(names_from = "MP", values_from = "q50") %>%
    select(OM:SPSS_100_40_SBMSY, SP_100_40_SBMSY, SP_01:SP_06) %>%
    as.data.frame()

###@> Saving the output table...
tab01_exp_10 <- PM_output %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06"))

####@> Exporting table...
write.table(tab01_exp, "Results/Table_PM-MP_OMs_10-18_ver00.csv",
            row.names = FALSE, sep = ",", dec = ".")

#####@> Looking to the PMs - Estimating the average between all OMs...
tab01 <- PM_output %>%
    group_by(MP, PM, Name, Caption) %>%
    summarise(q50 = mean(Values, na.rm = TRUE)) %>%
    as.data.frame()

####@> Figure PGK ~ POF...
tmp <- tab01 %>%
    filter(PM %in% c("POF", "PGK_short")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class = "Status")

p00 <- ggplot(data = tmp, aes(x = PGK_short, y = POF, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 20) +
    labs(x = "Prob. Kobe green quadrante over years 1-3",
         y = "Prob. Overfishing over years 1-30") +
    geom_vline(xintercept = 0.7, linetype = "dashed", colour = "gray") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1)) +
    facet_wrap(~class) +
    my_theme() +
    theme(legend.position = "none")
p00

tmp <- tab01 %>%
    filter(PM %in% c("POF", "PGK_med")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class = "Status")

p01 <- ggplot(data = tmp, aes(x = PGK_med, y = POF, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 50) +
    labs(x = "Prob. Kobe green quadrante over years 4-10",
         y = "Prob. Overfishing over years 1-30") +
    geom_vline(xintercept = 0.7, linetype = "dashed", colour = "gray") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1)) +
    facet_wrap(~class) +
    my_theme() +
    theme(legend.position = "none")
p01

tmp <- tab01 %>%
    filter(PM %in% c("POF", "PGK_long")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class = "Status")

p02 <- ggplot(data = tmp, aes(x = PGK_long, y = POF, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 40) +
    labs(x = "Prob. Kobe green quadrante over years 11-30",
         y = "Prob. Overfishing over years 1-30") +
    geom_vline(xintercept = 0.7, linetype = "dashed", colour = "gray") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1)) +
    facet_wrap(~class) +
    my_theme() +
    theme(legend.position = "none")
p02

tmp <- tab01 %>%
    filter(PM %in% c("POF", "PGK")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class = "Status")

p03 <- ggplot(data = tmp, aes(x = PGK, y = POF, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 40) +
    labs(x = "Prob. Kobe green quadrante over years 1-30",
         y = "Prob. Overfishing over years 1-30") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1)) +
    geom_vline(xintercept = 0.7, linetype = "dashed", colour = "gray") +
    facet_wrap(~class) +
    my_theme() +
    theme(legend.position = "none")
p03

ggsave("Results/Tradeplot_PGK_short-POF_ver01.svg", plot = p00,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Tradeplot_PGK_med-POF_ver01.svg", plot = p01,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Tradeplot_PGK_long-POF_ver01.svg", plot = p02,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Tradeplot_PGK-POF_ver01.svg", plot = p03,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

####@> Figure LRP ~ VarC...
tmp <- tab01 %>%
    filter(PM %in% c("VarCmedium", "LRP_med")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class01 = "Safety", class02 = "Stability")

##@> Correcting SP_05...
tmp$VarCmedium[tmp$VarCmedium > 1] <- 1

p00 <- ggplot(data = tmp, aes(x = VarCmedium, y = LRP_med, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 20) +
    labs(x = "Mean Variation in TAC (%) over years 4-10",
         y = "Prob. SB < 0.4 SBMSY over years 4-10") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1)) +
    geom_hline(yintercept = 0.1, linetype = "dashed", colour = "gray") +
    geom_vline(xintercept = 0.2, linetype = "dashed", colour = "gray") +
    facet_grid(class01~class02) +
    my_theme() +
    theme(legend.position = "none")
p00

tmp <- tab01 %>%
    filter(PM %in% c("VarClong", "LRP_long")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class01 = "Safety", class02 = "Stability")

p01 <- ggplot(data = tmp, aes(x = VarClong, y = LRP_long, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 20) +
    labs(x = "Mean Variation in TAC (%) over years 11-30",
         y = "Prob. SB < 0.4 SBMSY over years 11-30") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1)) +
    geom_hline(yintercept = 0.1, linetype = "dashed", colour = "gray") +
    geom_vline(xintercept = 0.2, linetype = "dashed", colour = "gray") +
    facet_grid(class01~class02) +
    my_theme() +
    theme(legend.position = "none")
p01

tmp <- tab01 %>%
    filter(PM %in% c("VarC", "LRP")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class01 = "Safety", class02 = "Stability")

p02 <- ggplot(data = tmp, aes(x = VarC, y = LRP, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 20) +
    labs(x = "Mean Variation in TAC (%) over years 1-30",
         y = "Prob. SB < 0.4 SBMSY over years 1-30") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1)) +
    geom_hline(yintercept = 0.1, linetype = "dashed", colour = "gray") +
    geom_vline(xintercept = 0.2, linetype = "dashed", colour = "gray") +
    facet_grid(class01~class02) +
    my_theme() +
    theme(legend.position = "none")
p02

ggsave("Results/Tradeplot_VarCmed-LRP_med_ver01.svg", plot = p00,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Tradeplot_VarClong-LRP_long_ver01.svg", plot = p01,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Tradeplot_VarC-LRP_ver01.svg", plot = p02,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

####@> Figure PGK ~ AvC...
tmp <- tab01 %>%
    filter(PM %in% c("AvC_short", "PGK_short")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class = "Status", class2 = "Yield")

p00 <- ggplot(data = tmp, aes(x = PGK_short, y = AvC_short, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 20) +
    labs(x = "Prob. Kobe green quadrante over years 1-3",
         y = "Median Catch (t) over years 1-3") +
    geom_vline(xintercept = 0.7, linetype = "dashed", colour = "gray") +
    scale_y_continuous(limits = c(0, 60000)) +
    scale_x_continuous(limits = c(0, 1)) +
    facet_grid(class2~class) +
    my_theme() +
    theme(legend.position = "none")
p00

tmp <- tab01 %>%
    filter(PM %in% c("AvC_med", "PGK_med")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class = "Status", class2 = "Yield")

p01 <- ggplot(data = tmp, aes(x = PGK_med, y = AvC_med, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 50) +
    labs(x = "Prob. Kobe green quadrante over years 4-10",
         y = "Median Catch (t) over years 4-10") +
    geom_vline(xintercept = 0.7, linetype = "dashed", colour = "gray") +
    scale_y_continuous(limits = c(0, 40000)) +
    scale_x_continuous(limits = c(0, 1)) +
    facet_grid(class2~class) +
    my_theme() +
    theme(legend.position = "none")
p01

tmp <- tab01 %>%
    filter(PM %in% c("AvC_long", "PGK_long")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class = "Status", class2 = "Yield")

p02 <- ggplot(data = tmp, aes(x = PGK_long, y = AvC_long, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 40) +
    labs(x = "Prob. Kobe green quadrante over years 11-30",
         y = "Median Catch (t) over years 11-30") +
    geom_vline(xintercept = 0.7, linetype = "dashed", colour = "gray") +
    scale_y_continuous(limits = c(0, 40000)) +
    scale_x_continuous(limits = c(0, 1)) +
    facet_wrap(~class) +
    my_theme() +
    theme(legend.position = "none")
p02

ggsave("Results/Tradeplot_PGK_short-AvC_short_ver01.svg", plot = p00,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Tradeplot_PGK_med-AvC_med_ver01.svg", plot = p01,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Tradeplot_PGK_long-AvC_long_ver01.svg", plot = p02,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

######@> Consolidating general probabilities for a heatmap table...
ProbTab <- tab01 %>%
    filter(PM %in% c("LRP", "LRP_long", "LRP_med", "LRP_short", "PGK",
                     "PGK_long", "PGK_med", "PGK_short", "PNOF", "POF",
                     "VarC","VarClong", "VarCmedium", "nLRP", "nLRP_long",
                     "nLRP_med", "nLRP_short", "rAvC_long", "rAvC_med",
                     "rAvC_short", "AvC_short", "AvC_med", "AvC_long")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    as.data.frame()

##@> Correcting SP_05...
ProbTab$VarCmedium[ProbTab$VarCmedium > 1] <- 1

#####@> Exporting ProbTab - Combined all 9 MOMs...
write.table(ProbTab, file = "Results/Table_PMs_FULL_ver01.csv",
            row.names = FALSE, sep = ",", dec = ".")

######@> Heatmaps for ProbTab...

#####@> Figure for Yield...
tmp <- ProbTab %>%
    select(MP:AvC_short) %>%
    pivot_longer(names_to = "PM", values_to = "value", 2:4) %>%
    mutate(class = "Yield",
           breaks = cut(value,
                        breaks = quantile(value, seq(0, 1, 0.2)),
                        include.lowest = TRUE, right = FALSE),
           PM = factor(PM, levels = c("AvC_short",
                                      "AvC_med",
                                      "AvC_long")),
           MP = factor(MP, levels = c("CC_20kt", "CC_30kt", "CC_40kt",
                                      "GB_slope", "Iratio", "Islope1",
                                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                                      "SPSS_100_40_SBMSY", "SP_01",
                                      "SP_02", "SP_03",
                                      "SP_04", "SP_05", "SP_06")),
           class2 = ifelse(value >= 20000, "Up", "Down"))
p00 <- ggplot(data = tmp, aes(x = PM, y = MP, fill = class2)) +
    geom_tile(colour = "black", alpha = 0) +
    geom_text(aes(label = round(value, 0)), size = 5,
              fontface = "bold") +
    scale_fill_brewer(palette = "Greys", direction = 1) +
    facet_wrap(~class) +
    labs(x = "Performance Metrics", y = "Management Procedures") +
    my_theme() +
    theme(legend.position = "none")
p00

#####@> Figure for Status...
tmp <- ProbTab %>%
    select(MP, PGK:PGK_short) %>%
    pivot_longer(names_to = "PM", values_to = "value", 2:5) %>%
    mutate(class = "Status",
           breaks = cut(value,
                        breaks = quantile(value, seq(0, 1, 0.2)),
                        include.lowest = TRUE, right = FALSE),
           PM = factor(PM, levels = c("PGK_short",
                                      "PGK_med",
                                      "PGK_long",
                                      "PGK")),
           MP = factor(MP, levels = c("CC_20kt", "CC_30kt", "CC_40kt",
                                      "GB_slope", "Iratio", "Islope1",
                                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                                      "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                                      "SP_04", "SP_05", "SP_06")),
           class2 = ifelse(value >= 0.7, "Up", "Down"))
p01 <- ggplot(data = tmp, aes(x = PM, y = MP, fill = class2)) +
    geom_tile(colour = "black", alpha = 0.8) +
    geom_text(aes(label = round(value, 3)), size = 5,
              fontface = "bold") +
    scale_fill_brewer(palette = "Greys", direction = 1) +
    facet_wrap(~class) +
    labs(x = "Performance Metrics", y = "Management Procedures") +
    my_theme() +
    theme(legend.position = "none")
p01

tmp <- ProbTab %>%
    select(MP, PNOF:POF) %>%
    pivot_longer(names_to = "PM", values_to = "value", 2:3) %>%
    mutate(class = "Status",
           breaks = cut(value,
                        breaks = quantile(value, seq(0, 1, 0.2)),
                        include.lowest = TRUE, right = FALSE),
           PM = factor(PM, levels = c("PNOF",
                                      "POF")),
           MP = factor(MP, levels = c("CC_20kt", "CC_30kt", "CC_40kt",
                                      "GB_slope", "Iratio", "Islope1",
                                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                                      "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                                      "SP_04", "SP_05", "SP_06")),
           class2 = ifelse(value >= 0.7, "Up", "Down"))
p01b <- ggplot(data = tmp, aes(x = PM, y = MP)) +
    geom_tile(colour = "black", fill = "white", alpha = 0.1) +
    geom_text(aes(label = round(value, 3)), size = 5,
              fontface = "bold") +
    ## scale_fill_brewer(palette = "Greys", direction = 1) +
    facet_wrap(~class) +
    labs(x = "Performance Metrics", y = "Management Procedures") +
    my_theme() +
    theme(legend.position = "none")
p01b

#####@> Figure for Safety...
tmp <- ProbTab %>%
    select(MP, LRP:LRP_short, nLRP:nLRP_short) %>%
    pivot_longer(names_to = "PM", values_to = "value", 2:5) %>%
    mutate(class = "Safety",
           breaks = cut(value,
                        breaks = quantile(value, seq(0, 1, 0.2)),
                        include.lowest = TRUE, right = FALSE),
           PM = factor(PM, levels = c("LRP_short",
                                      "LRP_med",
                                      "LRP_long",
                                      "LRP")),
           MP = factor(MP, levels = c("CC_20kt", "CC_30kt", "CC_40kt",
                                      "GB_slope", "Iratio", "Islope1",
                                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                                      "SPSS_100_40_SBMSY", "SP_01",
                                      "SP_02", "SP_03",
                                      "SP_04", "SP_05", "SP_06")),
           class2 = ifelse(value <= 0.1, "Up", "Down"))
p02 <- ggplot(data = tmp, aes(x = PM, y = MP, fill = class2)) +
    geom_tile(colour = "black", alpha = 0.8) +
    geom_text(aes(label = round(value, 3)), size = 5,
              fontface = "bold") +
    scale_fill_brewer(palette = "Greys", direction = 1) +
    facet_wrap(~class) +
    labs(x = "Performance Metrics", y = "Management Procedures") +
    my_theme() +
    theme(legend.position = "none")
p02

#####@> Figure for Stability...
tmp <- ProbTab %>%
    select(MP, VarC:VarCmedium) %>%
    pivot_longer(names_to = "PM", values_to = "value", 2:4) %>%
    mutate(class = "Stability",
           breaks = cut(value,
                        breaks = quantile(value, seq(0, 1, 0.2)),
                        include.lowest = TRUE, right = FALSE),
           PM = factor(PM, levels = c("VarCmedium",
                                      "VarClong",
                                      "VarC")),
           MP = factor(MP, levels = c("CC_20kt", "CC_30kt", "CC_40kt",
                                      "GB_slope", "Iratio", "Islope1",
                                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                                      "SPSS_100_40_SBMSY", "SP_01",
                                      "SP_02", "SP_03",
                                      "SP_04", "SP_05", "SP_06")),
           class2 = ifelse(value <= 0.2, "Up", "Down"))
p03 <- ggplot(data = tmp, aes(x = PM, y = MP, fill = class2)) +
    geom_tile(colour = "black", alpha = 0.8) +
    geom_text(aes(label = round(value, 3)), size = 5,
              fontface = "bold") +
    scale_fill_brewer(palette = "Greys", direction = 1) +
    facet_wrap(~class) +
    labs(x = "Performance Metrics", y = "Management Procedures") +
    my_theme() +
    theme(legend.position = "none")
p03

ggsave("Results/Table_PM-MP_Yield_10ImpErro_ver03.png", plot = p00,
       device = "png", dpi = 300, bg = "white", unit = "cm",
       w = 30, h = 25)

ggsave("Results/Table_PM-MP_Status_10ImpErro_ver02.png", plot = p01,
       device = "png", dpi = 300, bg = "white", unit = "cm",
       w = 30, h = 25)

ggsave("Results/Table_PM-MP_Status_B_10ImpErro_ver02.png", plot = p01b,
       device = "png", dpi = 300, bg = "white", unit = "cm",
       w = 30, h = 25)

ggsave("Results/Table_PM-MP_Safety_10ImpErro_ver02.png", plot = p02,
       device = "png", dpi = 300, bg = "white", unit = "cm",
       w = 30, h = 25)

ggsave("Results/Table_PM-MP_Stability_10ImpErro_ver02.png", plot = p03,
       device = "png", dpi = 300, bg = "white", unit = "cm",
       w = 30, h = 25)

######@> Extracting References Points from MMSE objects...

#####@> Extracting SB_SBMSY values...
tmp01 <- list()
for(i in seq_along(MSEs02)) {
    om <- names(MSEs02)[i]
    output <- MSEs02[[i]]@SB_SBMSY %>%
        structure(
            dimnames = list(
                Sim = 1:MSEs02[[i]]@nsim,
                MPs = MSEs02[[i]]@MPs,
                Year = 2020 + 1:MSEs02[[i]]@proyears
            )
        ) %>%
        reshape2::melt() %>%
        mutate(OM = om)
    tmp01[[i]] <- output
}
SB_SBMSY_output <- do.call(rbind.data.frame, tmp01)

####@> Looking to the RP - Estimating the average between all MOMs...
tab02 <- SB_SBMSY_output %>%
    filter(MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                      "GB_slope", "Iratio", "Islope1",
                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                      "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                      "SP_04", "SP_05", "SP_06")) %>%
    group_by(MPs, Year) %>%
    summarise(SB_SBMSY = mean(value, na.rm = TRUE),
              SB_SBMSY_Lwr = quantile(value, probs = 0.025,
                                      na.rm = TRUE),
              SB_SBMSY_Upr = quantile(value, probs = 0.975,
                                      na.rm = TRUE),
              SB_SBMSY_Lwr2 = quantile(value, probs = 0.1,
                                      na.rm = TRUE),
              SB_SBMSY_Upr2 = quantile(value, probs = 0.9,
                                      na.rm = TRUE)) %>%
    as.data.frame()

###@> Figure Reference Points...
p01 <- ggplot(data =
                  filter(tab02,
                         MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr,
                    ymax = SB_SBMSY_Upr),
                alpha = 0.5, fill = "#95B634") +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr2,
                    ymax = SB_SBMSY_Upr2),
                alpha = 0.2, fill = "#006750") +
    geom_line(aes(x = Year, y = SB_SBMSY), linewidth = 1,
              colour = "#006750") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 1) +
    labs(x = "Year", y = expression(SB/SB[MSY])) +
    my_theme()
p01

p02 <- ggplot(data =
                  filter(tab02,
                         MPs %in% c("Iratio", "Islope1", "GB_slope"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr,
                    ymax = SB_SBMSY_Upr),
                alpha = 0.5, fill = "#86BEDA") +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr2,
                    ymax = SB_SBMSY_Upr2),
                alpha = 0.2, fill = "#005690") +
    geom_line(aes(x = Year, y = SB_SBMSY), linewidth = 1,
              colour = "#005690") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 1) +
    labs(x = "Year", y = expression(SB/SB[MSY])) +
    my_theme()
p02

p03 <- ggplot(data =
                  filter(tab02,
                         MPs %in% c("SCA_100_40_SBMSY",
                                    "SP_100_40_SBMSY",
                                    "SPSS_100_40_SBMSY", "SP_01",
                                    "SP_02", "SP_03",
                                    "SP_04", "SP_05", "SP_06"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr,
                    ymax = SB_SBMSY_Upr),
                alpha = 0.5, fill = "#533D8B") +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr2,
                    ymax = SB_SBMSY_Upr2),
                alpha = 0.2, fill = "#5C1374") +
    geom_line(aes(x = Year, y = SB_SBMSY), linewidth = 1,
              colour = "#5C1374") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 2) +
    labs(x = "Year", y = expression(SB/SB[MSY])) +
    my_theme()
p03

ggsave("Results/Trajectory_SB_SBMSY_CC_ver01.svg", plot = p01,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Trajectory_SB_SBMSY_Index_ver01.svg", plot = p02,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Trajectory_SB_SBMSY_Model_ver01.svg", plot = p03,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

#####@> Extracting F_FMSY values...
tmp02 <- list()
for(i in seq_along(MSEs02)) {
    om <- names(MSEs02)[i]
    output <- MSEs02[[i]]@F_FMSY %>%
        structure(
            dimnames = list(
                Sim = 1:MSEs02[[i]]@nsim,
                MPs = MSEs02[[i]]@MPs,
                Year = 2020 + 1:MSEs02[[i]]@proyears
            )
        ) %>%
        reshape2::melt() %>%
        mutate(OM = om) %>%
        as.data.frame()
    tmp02[[i]] <- output
}
F_FMSY_output <- do.call(rbind.data.frame, tmp02)

####@> Looking to the RP - Estimating the average between all MOMs...
tab03 <- F_FMSY_output %>%
    filter(MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                      "GB_slope", "Iratio", "Islope1",
                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                      "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                      "SP_04", "SP_05", "SP_06")) %>%
    group_by(MPs, Year) %>%
    summarise(F_FMSY = mean(value, na.rm = TRUE),
              F_FMSY_Lwr = quantile(value, probs = 0.025,
                                    na.rm = TRUE),
              F_FMSY_Upr = quantile(value, probs = 0.975,
                                    na.rm = TRUE),
              F_FMSY_Lwr2 = quantile(value, probs = 0.1,
                                       na.rm = TRUE),
              F_FMSY_Upr2 = quantile(value, probs = 0.9,
                                       na.rm = TRUE)) %>%
    as.data.frame()

###@> Figure Reference Points...
p01 <- ggplot(data =
                  filter(tab03,
                         MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr,
                    ymax = F_FMSY_Upr),
                alpha = 0.5, fill = "#95B634") +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr2,
                    ymax = F_FMSY_Upr2),
                alpha = 0.2, fill = "#006750") +
    geom_line(aes(x = Year, y = F_FMSY), linewidth = 1,
              colour = "#006750") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 1) +
    labs(x = "Year", y = expression(F/F[MSY])) +
    my_theme()
p01

p02 <- ggplot(data =
                  filter(tab03,
                         MPs %in% c("Iratio", "Islope1", "GB_slope"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr,
                    ymax = F_FMSY_Upr),
                alpha = 0.5, fill = "#86BEDA") +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr2,
                    ymax = F_FMSY_Upr2),
                alpha = 0.2, fill = "#005690") +
    geom_line(aes(x = Year, y = F_FMSY), linewidth = 1,
              colour = "#005690") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 1) +
    labs(x = "Year", y = expression(F/F[MSY])) +
    my_theme()
p02

p03 <- ggplot(data =
                  filter(tab03,
                         MPs %in% c("SCA_100_40_SBMSY",
                                    "SP_100_40_SBMSY",
                                    "SPSS_100_40_SBMSY", "SP_01",
                                    "SP_02", "SP_03",
                                    "SP_04", "SP_05", "SP_06"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr,
                    ymax = F_FMSY_Upr),
                alpha = 0.5, fill = "#533D8B") +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr2,
                    ymax = F_FMSY_Upr2),
                alpha = 0.2, fill = "#5C1374") +
    geom_line(aes(x = Year, y = F_FMSY), linewidth = 1,
              colour = "#5C1374") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 2) +
    labs(x = "Year", y = expression(F/F[MSY])) +
    my_theme()
p03

ggsave("Results/Trajectory_F_FMSY_CC_ver01.svg", plot = p01,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Trajectory_F_FMSY_Index_ver01.svg", plot = p02,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Trajectory_F_FMSY_Model_ver01.svg", plot = p03,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

#####@> Extracting Catches values...
tmp03 <- list()
for(i in seq_along(MSEs02)) {
    om <- names(MSEs02)[i]
    output <- MSEs02[[i]]@Catch %>%
        structure(
            dimnames = list(
                Sim = 1:MSEs02[[i]]@nsim,
                MPs = MSEs02[[i]]@MPs,
                Year = 2020 + 1:MSEs02[[i]]@proyears
            )
        ) %>%
        reshape2::melt() %>%
        mutate(MOM = om) %>%
        as.data.frame()
    tmp03[[i]] <- output
}
Catch_output <- do.call(rbind.data.frame, tmp03)

####@> Looking to the RP - Estimating the average between all MOMs...
tab04 <- Catch_output %>%
    filter(MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                      "GB_slope", "Iratio", "Islope1",
                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                      "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                      "SP_04", "SP_05", "SP_06")) %>%
    group_by(MPs, Year) %>%
    summarise(Catch = mean(value, na.rm = TRUE),
              Catch_Lwr = quantile(value, probs = 0.025,
                                   na.rm = TRUE),
              Catch_Upr = quantile(value, probs = 0.975,
                                   na.rm = TRUE),
              Catch_Lwr2 = quantile(value, probs = 0.1,
                                   na.rm = TRUE),
              Catch_Upr2 = quantile(value, probs = 0.9,
                                   na.rm = TRUE)) %>%
    as.data.frame()

###@> Figure Reference Points...
p01 <- ggplot(data =
                  filter(tab04,
                         MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr,
                    ymax = Catch_Upr),
                alpha = 0.5, fill = "#95B634") +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr2,
                    ymax = Catch_Upr2),
                alpha = 0.2, fill = "#006750") +
    geom_line(aes(x = Year, y = Catch), linewidth = 1,
              colour = "#006750") +
    ## geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 1) +
    labs(x = "Year", y = "Catch (t)") +
    my_theme()
p01

p02 <- ggplot(data =
                  filter(tab04,
                         MPs %in% c("Iratio", "Islope1", "GB_slope"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr,
                    ymax = Catch_Upr),
                alpha = 0.5, fill = "#86BEDA") +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr2,
                    ymax = Catch_Upr2),
                alpha = 0.2, fill = "#005690") +
    geom_line(aes(x = Year, y = Catch), linewidth = 1,
              colour = "#005690") +
    ## geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 1) +
    labs(x = "Year", y = "Catch (t)") +
    my_theme()
p02

p03 <- ggplot(data =
                  filter(tab04,
                         MPs %in% c("SCA_100_40_SBMSY",
                                    "SP_100_40_SBMSY",
                                    "SPSS_100_40_SBMSY", "SP_01",
                                    "SP_02", "SP_03",
                                    "SP_04", "SP_05", "SP_06"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr,
                    ymax = Catch_Upr),
                alpha = 0.5, fill = "#533D8B") +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr2,
                    ymax = Catch_Upr2),
                alpha = 0.2, fill = "#5C1374") +
    geom_line(aes(x = Year, y = Catch), linewidth = 1,
              colour = "#5C1374") +
    ## geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 2) +
    labs(x = "Year", y = "Catch (t)") +
    my_theme()
p03

ggsave("Results/Trajectory_Catch_CC_ver01.svg", plot = p01,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Trajectory_Catch_Index_ver01.svg", plot = p02,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Trajectory_Catch_Model_ver01.svg", plot = p03,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

#####@> Extracting TAC values...
tmp04 <- list()
for(i in seq_along(MSEs02)) {
    om <- names(MSEs02)[i]
    output <- MSEs02[[i]]@TAC %>%
        structure(
            dimnames = list(
                Sim = 1:MSEs02[[i]]@nsim,
                MPs = MSEs02[[i]]@MPs,
                Year = 2020 + 1:MSEs02[[i]]@proyears
            )
        ) %>%
        reshape2::melt() %>%
        mutate(MOM = om) %>%
        as.data.frame()
    tmp04[[i]] <- output
}
TAC_output <- do.call(rbind.data.frame, tmp04)

####@> Looking to the RP - Estimating the average between all MOMs...
tab05 <- TAC_output %>%
    filter(MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                      "GB_slope", "Iratio", "Islope1",
                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                      "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                      "SP_04", "SP_05", "SP_06")) %>%
    group_by(MPs, Year) %>%
    summarise(TAC = median(value, na.rm = TRUE),
              TAC_Lwr = quantile(value, probs = 0.025,
                                 na.rm = TRUE),
              TAC_Upr = quantile(value, probs = 0.975,
                                 na.rm = TRUE),
              TAC_Lwr2 = quantile(value, probs = 0.1,
                                 na.rm = TRUE),
              TAC_Upr2 = quantile(value, probs = 0.9,
                                 na.rm = TRUE)) %>%
    as.data.frame()

###@> Figure Reference Points...
p01 <- ggplot(data =
                  filter(tab05,
                         MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = TAC_Lwr,
                    ymax = TAC_Upr),
                alpha = 0.5, fill = "#95B634") +
    geom_ribbon(aes(x = Year,
                    ymin = TAC_Lwr2,
                    ymax = TAC_Upr2),
                alpha = 0.2, fill = "#006750") +
    geom_line(aes(x = Year, y = TAC), linewidth = 1,
              colour = "#006750") +
    ## geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 1) +
    labs(x = "Year", y = "TAC (t)") +
    my_theme()
p01

p02 <- ggplot(data =
                  filter(tab05,
                         MPs %in% c("Iratio", "Islope1", "GB_slope"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = TAC_Lwr,
                    ymax = TAC_Upr),
                alpha = 0.5, fill = "#86BEDA") +
    geom_ribbon(aes(x = Year,
                    ymin = TAC_Lwr2,
                    ymax = TAC_Upr2),
                alpha = 0.2, fill = "#005690") +
    geom_line(aes(x = Year, y = TAC), linewidth = 1,
              colour = "#005690") +
    ## geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 1) +
    labs(x = "Year", y = "TAC (t)") +
    my_theme()
p02

p03 <- ggplot(data =
                  filter(tab05,
                         MPs %in% c("SCA_100_40_SBMSY",
                                    "SP_100_40_SBMSY",
                                    "SPSS_100_40_SBMSY", "SP_01",
                                    "SP_02", "SP_03",
                                    "SP_04", "SP_05", "SP_06"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = TAC_Lwr,
                    ymax = TAC_Upr),
                alpha = 0.5, fill = "#533D8B") +
    geom_ribbon(aes(x = Year,
                    ymin = TAC_Lwr2,
                    ymax = TAC_Upr2),
                alpha = 0.2, fill = "#5C1374") +
    geom_line(aes(x = Year, y = TAC), linewidth = 1,
              colour = "#5C1374") +
    ## geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 2) +
    labs(x = "Year", y = "TAC (t)") +
    my_theme()
p03

ggsave("Results/Trajectory_TAC_CC_ver01.svg", plot = p01,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Trajectory_TAC_Index_ver01.svg", plot = p02,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Trajectory_TAC_Model_ver01.svg", plot = p03,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

######@> Creating the Kobe plots (by Year and general)...

#####@> Merging SB/SBMSY and F/FMSY...
df00 <- SB_SBMSY_output %>%
    left_join(F_FMSY_output, by = c("Sim", "MPs", "Year", "OM")) %>%
    filter(MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                      "GB_slope", "Iratio", "Islope1",
                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                      "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                      "SP_04", "SP_05", "SP_06")) %>%
    rename("SB_SBMSY" = value.x, "F_FMSY" = value.y)

#####@> Estimating Kobe values for the last year...
tmp <- df00 %>%
    filter(Year == 2053) %>%
    group_by(Sim, MPs) %>%
    summarise(SB_SBMSY = mean(SB_SBMSY, na.rm = TRUE),
              F_FMSY = mean(F_FMSY, na.rm = TRUE)) %>%
    mutate(OF = ifelse(F_FMSY <= 1, FALSE, TRUE),
           OFD = ifelse(SB_SBMSY <= 1, TRUE, FALSE))

#####@> Percentage of cases...
valdf <- tmp %>%
    group_by(MPs) %>%
    summarise(BL = sum(OF == FALSE & OFD == TRUE)/100 * 100,
              BR = sum(OF == FALSE & OFD == FALSE)/100 * 100,
              TL = sum(OF == TRUE & OFD == TRUE)/100 * 100,
              TR = sum(OF == TRUE & OFD == FALSE)/100 * 100)
valdf <- valdf %>% tidyr::pivot_longer(., cols = 2:5)
valdf$x <- -Inf
valdf$y <- -Inf
valdf$y[valdf$name == "TL"] <- Inf
valdf$y[valdf$name == "TR"] <- Inf
valdf$x[valdf$name == "BR"] <- Inf
valdf$x[valdf$name == "TR"] <- Inf
valdf$value <- round(valdf$value, 2)
valdf$value <- paste0(valdf$value, "%")
valdf$hjustvar <- -2
valdf$vjustvar <- -2
valdf$hjustvar[valdf$name == "TL"] <- -1
valdf$hjustvar[valdf$name == "TR"] <- 2
valdf$hjustvar[valdf$name == "BL"] <- -1
valdf$hjustvar[valdf$name == "BR"] <- 2
valdf$vjustvar[valdf$name == "TL"] <- 2
valdf$vjustvar[valdf$name == "TR"] <- 2
valdf$vjustvar[valdf$name == "BL"] <- -2
valdf$vjustvar[valdf$name == "BR"] <- -2

#####@> Colors to figure...
kobe_df <- bind_rows(data.frame(x = c(0, 0, 1, 1),
                                y = c(0, 1, 1, 0),
                                fill = "bl"),
                     data.frame(x = c(1, 1, 2, 2),
                                y = c(0, 1, 1, 0),
                                fill = "br"),
                     data.frame(x = c(0, 0, 1, 1),
                                y = c(1, 2, 2, 1),
                                fill = "tl"),
                     data.frame(x = c(1, 1, 2, 2),
                                y = c(1, 2, 2, 1),
                                fill = "tr"))
kobe_df$alpha <- 0.3

#####@> Figure dots...
p00 <- ggplot() +
    geom_polygon(data = kobe_df,
                 aes(x = x, y = y, fill = fill, alpha = alpha)) +
    scale_fill_manual(values = c("#F8DC7A", "#67C18B", "#D8775D",
                                 "#FDBD56")) +
    geom_point(data = tmp,
               aes(x = SB_SBMSY, y = F_FMSY), alpha = 0.5) +
    geom_text(data = valdf, fontface = "bold", size = 4,
              aes(x = x, y = y, label = value,
                  hjust = hjustvar, vjust = vjustvar),
              colour = "black") +
    expand_limits(x = c(0, 2), y = c(0, 2)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 2)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 2)) +
    geom_hline(yintercept = 1, color = "darkgray", linetype = 2) +
    geom_vline(xintercept = 1, color = "darkgray", linetype = 2) +
    labs(x = expression(SB/SB[MSY]), y = expression(F/F[MSY])) +
    facet_wrap(~MPs) +
    my_theme() +
    theme(legend.position = "none")
p00

ggsave("Results/Kobe_plot_Last_Year_2053_ver01.svg", plot = p00,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

######@> Estimating Kobe values per year...
tmp <- df00 %>%
    group_by(Sim, Year, MPs) %>%
    summarise(SB_SBMSY = mean(SB_SBMSY, na.rm = TRUE),
              F_FMSY = mean(F_FMSY, na.rm = TRUE)) %>%
    mutate(OF = ifelse(F_FMSY <= 1, FALSE, TRUE),
           OFD = ifelse(SB_SBMSY <= 1, TRUE, FALSE))

#####@> Proportions by year...
valdf <- tmp %>%
    group_by(Year, MPs) %>%
    summarise(Yellow = sum(OF == FALSE & OFD == TRUE)/100 * 100,
              Green = sum(OF == FALSE & OFD == FALSE)/100 * 100,
              Red = sum(OF == TRUE & OFD == TRUE)/100 * 100,
              Orange = sum(OF == TRUE & OFD == FALSE)/100 * 100) %>%
    pivot_longer(names_to = "Cond", values_to = "Perc", 3:6) %>%
    mutate(Cond = factor(Cond, levels = c("Green", "Yellow", "Orange",
                                          "Red")))

p00 <- ggplot() +
    geom_bar(data = filter(valdf, Year %in% 2021:2053),
             aes(x = Year, y = Perc, fill = Cond),
             stat = "identity", colour = "black",
             width = 1) +
    facet_wrap(~MPs) +
    scale_fill_manual(values = c("#67C18B", "#F8DC7A", "#FDBD56",
                                 "#D8775D")) +
    labs(x = "Year", y = "%") +
    my_theme() +
    theme(legend.position = "none")
p00

p01 <- ggplot() +
    geom_area(data = filter(valdf, Year %in% 2021:2053),
              aes(x = Year, y = Perc, fill = Cond),
              stat = "identity", colour = "black") +
    facet_wrap(~MPs) +
    scale_fill_manual(values = c("#67C18B", "#F8DC7A", "#FDBD56",
                                 "#D8775D")) +
    labs(x = "Year", y = "%") +
    my_theme() +
    theme(legend.position = "none")
p01

ggsave("Results/Kobe_plot_by_year_Bar_ver01.svg", plot = p00,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 30)

ggsave("Results/Kobe_plot_by_year_Area_ver01.svg", plot = p01,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 30)

######@>----------------------------------------------------------------
######@> 20% error implementation - OMs 19-27...

######@> Extracting PMs from all MMSE's objects...

#####@> Listing PMs functions...
PMlist <- list()
for(i in seq_along(PMs)) {
    PMlist[[i]] <- try(get(PMs[i]), silent = TRUE)
}

######@> Extracting values...
PMvalues <- list()
PMvalues2 <- list()
for(i in seq_along(MSEs03)) {
    for(j in seq_along(PMs)) {
        pm <- PMlist[[j]](MSEs03[[i]])
        om <- names(MSEs03[i])
        mps <- pm@MPs
        val <- data.frame(pm@Prob)
        names(val) <- mps
        val <- val %>%
            mutate(sim = 1:100) %>%
            pivot_longer(names_to = "MP", values_to = "Values", 1:20)
        nom <- pm@Name
        cap <- pm@Caption
        pm <- PMs[j]
        PMvalues[[j]] <- val %>%
            mutate(OM = om, Name = nom, Caption = cap, PM = pm) %>%
            select(OM, Name, Caption, sim, MP, PM, Values) %>%
            as.data.frame()
    }
    out <- do.call(rbind.data.frame, PMvalues)
    PMvalues2[[i]] <- out
}
PM_output <- do.call(rbind.data.frame, PMvalues2)

#####@> Looking for NAs...
tab01_NA <- PM_output %>%
    filter(is.na(Values)) %>%
    group_by(OM, MP, PM) %>%
    summarise(N = n()) %>%
    as.data.frame()

#####@> Looking to the PMs - Estimating the average for each OMs...
tab01_exp <- PM_output %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    group_by(OM, MP, PM) %>%
    summarise(q50 = mean(Values, na.rm = TRUE)) %>%
    pivot_wider(names_from = "MP", values_from = "q50") %>%
    select(OM:SPSS_100_40_SBMSY, SP_100_40_SBMSY, SP_01:SP_06) %>%
    as.data.frame()

###@> Saving the output table...
tab01_exp_20 <- PM_output %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06"))

####@> Exporting table...
write.table(tab01_exp, "Results/Table_PM-MP_OMs_19-27_ver00.csv",
            row.names = FALSE, sep = ",", dec = ".")

#####@> Looking to the PMs - Estimating the average between all OMs...
tab01 <- PM_output %>%
    group_by(MP, PM, Name, Caption) %>%
    summarise(q50 = mean(Values, na.rm = TRUE)) %>%
    as.data.frame()

####@> Figure PGK ~ POF...
tmp <- tab01 %>%
    filter(PM %in% c("POF", "PGK_short")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class = "Status")

p00 <- ggplot(data = tmp, aes(x = PGK_short, y = POF, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 20) +
    labs(x = "Prob. Kobe green quadrante over years 1-3",
         y = "Prob. Overfishing over years 1-30") +
    geom_vline(xintercept = 0.7, linetype = "dashed", colour = "gray") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1)) +
    facet_wrap(~class) +
    my_theme() +
    theme(legend.position = "none")
p00

tmp <- tab01 %>%
    filter(PM %in% c("POF", "PGK_med")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class = "Status")

p01 <- ggplot(data = tmp, aes(x = PGK_med, y = POF, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 50) +
    labs(x = "Prob. Kobe green quadrante over years 4-10",
         y = "Prob. Overfishing over years 1-30") +
    geom_vline(xintercept = 0.7, linetype = "dashed", colour = "gray") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1)) +
    facet_wrap(~class) +
    my_theme() +
    theme(legend.position = "none")
p01

tmp <- tab01 %>%
    filter(PM %in% c("POF", "PGK_long")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class = "Status")

p02 <- ggplot(data = tmp, aes(x = PGK_long, y = POF, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 40) +
    labs(x = "Prob. Kobe green quadrante over years 11-30",
         y = "Prob. Overfishing over years 1-30") +
    geom_vline(xintercept = 0.7, linetype = "dashed", colour = "gray") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1)) +
    facet_wrap(~class) +
    my_theme() +
    theme(legend.position = "none")
p02

tmp <- tab01 %>%
    filter(PM %in% c("POF", "PGK")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class = "Status")

p03 <- ggplot(data = tmp, aes(x = PGK, y = POF, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 40) +
    labs(x = "Prob. Kobe green quadrante over years 1-30",
         y = "Prob. Overfishing over years 1-30") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1)) +
    geom_vline(xintercept = 0.7, linetype = "dashed", colour = "gray") +
    facet_wrap(~class) +
    my_theme() +
    theme(legend.position = "none")
p03

ggsave("Results/Tradeplot_PGK_short-POF_ver01.svg", plot = p00,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Tradeplot_PGK_med-POF_ver01.svg", plot = p01,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Tradeplot_PGK_long-POF_ver01.svg", plot = p02,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Tradeplot_PGK-POF_ver01.svg", plot = p03,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

####@> Figure LRP ~ VarC...
tmp <- tab01 %>%
    filter(PM %in% c("VarCmedium", "LRP_med")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class01 = "Safety", class02 = "Stability")

##@> Correcting SP_05...
tmp$VarCmedium[tmp$VarCmedium > 1] <- 1

p00 <- ggplot(data = tmp, aes(x = VarCmedium, y = LRP_med, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 20) +
    labs(x = "Mean Variation in TAC (%) over years 4-10",
         y = "Prob. SB < 0.4 SBMSY over years 4-10") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1)) +
    geom_hline(yintercept = 0.1, linetype = "dashed", colour = "gray") +
    geom_vline(xintercept = 0.2, linetype = "dashed", colour = "gray") +
    facet_grid(class01~class02) +
    my_theme() +
    theme(legend.position = "none")
p00

tmp <- tab01 %>%
    filter(PM %in% c("VarClong", "LRP_long")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class01 = "Safety", class02 = "Stability")

p01 <- ggplot(data = tmp, aes(x = VarClong, y = LRP_long, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 20) +
    labs(x = "Mean Variation in TAC (%) over years 11-30",
         y = "Prob. SB < 0.4 SBMSY over years 11-30") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1)) +
    geom_hline(yintercept = 0.1, linetype = "dashed", colour = "gray") +
    geom_vline(xintercept = 0.2, linetype = "dashed", colour = "gray") +
    facet_grid(class01~class02) +
    my_theme() +
    theme(legend.position = "none")
p01

tmp <- tab01 %>%
    filter(PM %in% c("VarC", "LRP")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class01 = "Safety", class02 = "Stability")

p02 <- ggplot(data = tmp, aes(x = VarC, y = LRP, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 20) +
    labs(x = "Mean Variation in TAC (%) over years 1-30",
         y = "Prob. SB < 0.4 SBMSY over years 1-30") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1)) +
    geom_hline(yintercept = 0.1, linetype = "dashed", colour = "gray") +
    geom_vline(xintercept = 0.2, linetype = "dashed", colour = "gray") +
    facet_grid(class01~class02) +
    my_theme() +
    theme(legend.position = "none")
p02

ggsave("Results/Tradeplot_VarCmed-LRP_med_ver01.svg", plot = p00,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Tradeplot_VarClong-LRP_long_ver01.svg", plot = p01,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Tradeplot_VarC-LRP_ver01.svg", plot = p02,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

####@> Figure PGK ~ AvC...
tmp <- tab01 %>%
    filter(PM %in% c("AvC_short", "PGK_short")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class = "Status", class2 = "Yield")

p00 <- ggplot(data = tmp, aes(x = PGK_short, y = AvC_short, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 20) +
    labs(x = "Prob. Kobe green quadrante over years 1-3",
         y = "Median Catch (t) over years 1-3") +
    geom_vline(xintercept = 0.7, linetype = "dashed", colour = "gray") +
    scale_y_continuous(limits = c(0, 60000)) +
    scale_x_continuous(limits = c(0, 1)) +
    facet_grid(class2~class) +
    my_theme() +
    theme(legend.position = "none")
p00

tmp <- tab01 %>%
    filter(PM %in% c("AvC_med", "PGK_med")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class = "Status", class2 = "Yield")

p01 <- ggplot(data = tmp, aes(x = PGK_med, y = AvC_med, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 50) +
    labs(x = "Prob. Kobe green quadrante over years 4-10",
         y = "Median Catch (t) over years 4-10") +
    geom_vline(xintercept = 0.7, linetype = "dashed", colour = "gray") +
    scale_y_continuous(limits = c(0, 40000)) +
    scale_x_continuous(limits = c(0, 1)) +
    facet_grid(class2~class) +
    my_theme() +
    theme(legend.position = "none")
p01

tmp <- tab01 %>%
    filter(PM %in% c("AvC_long", "PGK_long")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    mutate(class = "Status", class2 = "Yield")

p02 <- ggplot(data = tmp, aes(x = PGK_long, y = AvC_long, fill = MP)) +
    geom_point(size = 4, colour = "black", pch = 21) +
    geom_label_repel(aes(label = MP), force = 40) +
    labs(x = "Prob. Kobe green quadrante over years 11-30",
         y = "Median Catch (t) over years 11-30") +
    geom_vline(xintercept = 0.7, linetype = "dashed", colour = "gray") +
    scale_y_continuous(limits = c(0, 40000)) +
    scale_x_continuous(limits = c(0, 1)) +
    facet_wrap(~class) +
    my_theme() +
    theme(legend.position = "none")
p02

ggsave("Results/Tradeplot_PGK_short-AvC_short_ver01.svg", plot = p00,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Tradeplot_PGK_med-AvC_med_ver01.svg", plot = p01,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Tradeplot_PGK_long-AvC_long_ver01.svg", plot = p02,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

######@> Consolidating general probabilities for a heatmap table...
ProbTab <- tab01 %>%
    filter(PM %in% c("LRP", "LRP_long", "LRP_med", "LRP_short", "PGK",
                     "PGK_long", "PGK_med", "PGK_short", "PNOF", "POF",
                     "VarC","VarClong", "VarCmedium", "nLRP", "nLRP_long",
                     "nLRP_med", "nLRP_short", "rAvC_long", "rAvC_med",
                     "rAvC_short", "AvC_short", "AvC_med", "AvC_long")) %>%
    filter(MP %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                     "GB_slope", "Iratio", "Islope1",
                     "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                     "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                     "SP_04", "SP_05", "SP_06")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    as.data.frame()

##@> Correcting SP_05...
ProbTab$VarCmedium[ProbTab$VarCmedium > 1] <- 1

#####@> Exporting ProbTab - Combined all 9 MOMs...
write.table(ProbTab, file = "Results/Table_PMs_FULL_ver01.csv",
            row.names = FALSE, sep = ",", dec = ".")

######@> Heatmaps for ProbTab...

#####@> Figure for Yield...
tmp <- ProbTab %>%
    select(MP:AvC_short) %>%
    pivot_longer(names_to = "PM", values_to = "value", 2:4) %>%
    mutate(class = "Yield",
           breaks = cut(value,
                        breaks = quantile(value, seq(0, 1, 0.2)),
                        include.lowest = TRUE, right = FALSE),
           PM = factor(PM, levels = c("AvC_short",
                                      "AvC_med",
                                      "AvC_long")),
           MP = factor(MP, levels = c("CC_20kt", "CC_30kt", "CC_40kt",
                                      "GB_slope", "Iratio", "Islope1",
                                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                                      "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                                      "SP_04", "SP_05", "SP_06")),
           class2 = ifelse(value >= 20000, "Up", "Down")) ## Up average
## 4 years
p00 <- ggplot(data = tmp, aes(x = PM, y = MP, fill = class2)) +
    geom_tile(colour = "black", alpha = 0) +
    geom_text(aes(label = round(value, 0)), size = 5,
              fontface = "bold") +
    scale_fill_brewer(palette = "Greys", direction = 1) +
    facet_wrap(~class) +
    labs(x = "Performance Metrics", y = "Management Procedures") +
    my_theme() +
    theme(legend.position = "none")
p00

#####@> Figure for Status...
tmp <- ProbTab %>%
    select(MP, PGK:PGK_short) %>%
    pivot_longer(names_to = "PM", values_to = "value", 2:5) %>%
    mutate(class = "Status",
           breaks = cut(value,
                        breaks = quantile(value, seq(0, 1, 0.2)),
                        include.lowest = TRUE, right = FALSE),
           PM = factor(PM, levels = c("PGK_short",
                                      "PGK_med",
                                      "PGK_long",
                                      "PGK")),
           MP = factor(MP, levels = c("CC_20kt", "CC_30kt", "CC_40kt",
                                      "GB_slope", "Iratio", "Islope1",
                                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                                      "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                                      "SP_04", "SP_05", "SP_06")),
           class2 = ifelse(value >= 0.7, "Up", "Down"))
p01 <- ggplot(data = tmp, aes(x = PM, y = MP, fill = class2)) +
    geom_tile(colour = "black", alpha = 0.8) +
    geom_text(aes(label = round(value, 3)), size = 5,
              fontface = "bold") +
    scale_fill_brewer(palette = "Greys", direction = 1) +
    facet_wrap(~class) +
    labs(x = "Performance Metrics", y = "Management Procedures") +
    my_theme() +
    theme(legend.position = "none")
p01

tmp <- ProbTab %>%
    select(MP, PNOF:POF) %>%
    pivot_longer(names_to = "PM", values_to = "value", 2:3) %>%
    mutate(class = "Status",
           breaks = cut(value,
                        breaks = quantile(value, seq(0, 1, 0.2)),
                        include.lowest = TRUE, right = FALSE),
           PM = factor(PM, levels = c("PNOF",
                                      "POF")),
           MP = factor(MP, levels = c("CC_20kt", "CC_30kt", "CC_40kt",
                                      "GB_slope", "Iratio", "Islope1",
                                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                                      "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                                      "SP_04", "SP_05", "SP_06")),
           class2 = ifelse(value >= 0.7, "Up", "Down"))
p01b <- ggplot(data = tmp, aes(x = PM, y = MP)) +
    geom_tile(colour = "black", fill = "white", alpha = 0.1) +
    geom_text(aes(label = round(value, 3)), size = 5,
              fontface = "bold") +
    ## scale_fill_brewer(palette = "Greys", direction = 1) +
    facet_wrap(~class) +
    labs(x = "Performance Metrics", y = "Management Procedures") +
    my_theme() +
    theme(legend.position = "none")
p01b

#####@> Figure for Safety...
tmp <- ProbTab %>%
    select(MP, LRP:LRP_short) %>%
    pivot_longer(names_to = "PM", values_to = "value", 2:5) %>%
    mutate(class = "Safety",
           breaks = cut(value,
                        breaks = quantile(value, seq(0, 1, 0.2)),
                        include.lowest = TRUE, right = FALSE),
           PM = factor(PM, levels = c("LRP_short",
                                      "LRP_med",
                                      "LRP_long",
                                      "LRP")),
           MP = factor(MP, levels = c("CC_20kt", "CC_30kt", "CC_40kt",
                                      "GB_slope", "Iratio", "Islope1",
                                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                                      "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                                      "SP_04", "SP_05", "SP_06")),
           class2 = ifelse(value <= 0.1, "Up", "Down"))
p02 <- ggplot(data = tmp, aes(x = PM, y = MP, fill = class2)) +
    geom_tile(colour = "black", alpha = 0.8) +
    geom_text(aes(label = round(value, 3)), size = 5,
              fontface = "bold") +
    scale_fill_brewer(palette = "Greys", direction = 1) +
    facet_wrap(~class) +
    labs(x = "Performance Metrics", y = "Management Procedures") +
    my_theme() +
    theme(legend.position = "none")
p02

#####@> Figure for Stability...
tmp <- ProbTab %>%
    select(MP, VarC:VarCmedium) %>%
    pivot_longer(names_to = "PM", values_to = "value", 2:4) %>%
    mutate(class = "Stability",
           breaks = cut(value,
                        breaks = quantile(value, seq(0, 1, 0.2)),
                        include.lowest = TRUE, right = FALSE),
           PM = factor(PM, levels = c("VarCmedium",
                                      "VarClong",
                                      "VarC")),
           MP = factor(MP, levels = c("CC_20kt", "CC_30kt", "CC_40kt",
                                      "GB_slope", "Iratio", "Islope1",
                                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                                      "SPSS_100_40_SBMSY", "SP_01",
                                      "SP_02", "SP_03",
                                      "SP_04", "SP_05", "SP_06")),
           class2 = ifelse(value <= 0.2, "Up", "Down"))
p03 <- ggplot(data = tmp, aes(x = PM, y = MP, fill = class2)) +
    geom_tile(colour = "black", alpha = 0.8) +
    geom_text(aes(label = round(value, 3)), size = 5,
              fontface = "bold") +
    scale_fill_brewer(palette = "Greys", direction = 1) +
    facet_wrap(~class) +
    labs(x = "Performance Metrics", y = "Management Procedures") +
    my_theme() +
    theme(legend.position = "none")
p03

ggsave("Results/Table_PM-MP_Yield_20ImpError_ver03.png", plot = p00,
       device = "png", dpi = 300, bg = "white", unit = "cm",
       w = 30, h = 25)

ggsave("Results/Table_PM-MP_Status_20ImpError_ver02.png", plot = p01,
       device = "png", dpi = 300, bg = "white", unit = "cm",
       w = 30, h = 25)

ggsave("Results/Table_PM-MP_Status_B_20ImpError_ver02.png", plot = p01b,
       device = "png", dpi = 300, bg = "white", unit = "cm",
       w = 30, h = 25)

ggsave("Results/Table_PM-MP_Safety_20ImpError_ver02.png", plot = p02,
       device = "png", dpi = 300, bg = "white", unit = "cm",
       w = 30, h = 25)

ggsave("Results/Table_PM-MP_Stability_20ImpError_ver02.png", plot = p03,
       device = "png", dpi = 300, bg = "white", unit = "cm",
       w = 30, h = 25)

######@> Extracting References Points from MMSE objects...

#####@> Extracting SB_SBMSY values...
tmp01 <- list()
for(i in seq_along(MSEs03)) {
    om <- names(MSEs03)[i]
    output <- MSEs03[[i]]@SB_SBMSY %>%
        structure(
            dimnames = list(
                Sim = 1:MSEs03[[i]]@nsim,
                MPs = MSEs03[[i]]@MPs,
                Year = 2020 + 1:MSEs03[[i]]@proyears
            )
        ) %>%
        reshape2::melt() %>%
        mutate(OM = om)
    tmp01[[i]] <- output
}
SB_SBMSY_output <- do.call(rbind.data.frame, tmp01)

####@> Looking to the RP - Estimating the average between all MOMs...
tab02 <- SB_SBMSY_output %>%
    filter(MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                      "GB_slope", "Iratio", "Islope1",
                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                      "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                      "SP_04", "SP_05", "SP_06")) %>%
    group_by(MPs, Year) %>%
    summarise(SB_SBMSY = mean(value, na.rm = TRUE),
              SB_SBMSY_Lwr = quantile(value, probs = 0.025,
                                      na.rm = TRUE),
              SB_SBMSY_Upr = quantile(value, probs = 0.975,
                                      na.rm = TRUE),
              SB_SBMSY_Lwr2 = quantile(value, probs = 0.1,
                                      na.rm = TRUE),
              SB_SBMSY_Upr2 = quantile(value, probs = 0.9,
                                      na.rm = TRUE)) %>%
    as.data.frame()

###@> Figure Reference Points...
p01 <- ggplot(data =
                  filter(tab02,
                         MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr,
                    ymax = SB_SBMSY_Upr),
                alpha = 0.5, fill = "#95B634") +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr2,
                    ymax = SB_SBMSY_Upr2),
                alpha = 0.2, fill = "#006750") +
    geom_line(aes(x = Year, y = SB_SBMSY), linewidth = 1,
              colour = "#006750") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 1) +
    labs(x = "Year", y = expression(SB/SB[MSY])) +
    my_theme()
p01

p02 <- ggplot(data =
                  filter(tab02,
                         MPs %in% c("Iratio", "Islope1", "GB_slope"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr,
                    ymax = SB_SBMSY_Upr),
                alpha = 0.5, fill = "#86BEDA") +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr2,
                    ymax = SB_SBMSY_Upr2),
                alpha = 0.2, fill = "#005690") +
    geom_line(aes(x = Year, y = SB_SBMSY), linewidth = 1,
              colour = "#005690") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 1) +
    labs(x = "Year", y = expression(SB/SB[MSY])) +
    my_theme()
p02

p03 <- ggplot(data =
                  filter(tab02,
                         MPs %in% c("SCA_100_40_SBMSY",
                                    "SP_100_40_SBMSY",
                                    "SPSS_100_40_SBMSY", "SP_01",
                                    "SP_02", "SP_03",
                                    "SP_04", "SP_05", "SP_06"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr,
                    ymax = SB_SBMSY_Upr),
                alpha = 0.5, fill = "#533D8B") +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr2,
                    ymax = SB_SBMSY_Upr2),
                alpha = 0.2, fill = "#5C1374") +
    geom_line(aes(x = Year, y = SB_SBMSY), linewidth = 1,
              colour = "#5C1374") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 2) +
    labs(x = "Year", y = expression(SB/SB[MSY])) +
    my_theme()
p03

ggsave("Results/Trajectory_SB_SBMSY_CC_ver01.svg", plot = p01,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Trajectory_SB_SBMSY_Index_ver01.svg", plot = p02,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Trajectory_SB_SBMSY_Model_ver01.svg", plot = p03,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

#####@> Extracting F_FMSY values...
tmp02 <- list()
for(i in seq_along(MSEs03)) {
    om <- names(MSEs03)[i]
    output <- MSEs03[[i]]@F_FMSY %>%
        structure(
            dimnames = list(
                Sim = 1:MSEs03[[i]]@nsim,
                MPs = MSEs03[[i]]@MPs,
                Year = 2020 + 1:MSEs03[[i]]@proyears
            )
        ) %>%
        reshape2::melt() %>%
        mutate(OM = om) %>%
        as.data.frame()
    tmp02[[i]] <- output
}
F_FMSY_output <- do.call(rbind.data.frame, tmp02)

####@> Looking to the RP - Estimating the average between all MOMs...
tab03 <- F_FMSY_output %>%
    filter(MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                      "GB_slope", "Iratio", "Islope1",
                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                      "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                      "SP_04", "SP_05", "SP_06")) %>%
    group_by(MPs, Year) %>%
    summarise(F_FMSY = mean(value, na.rm = TRUE),
              F_FMSY_Lwr = quantile(value, probs = 0.025,
                                    na.rm = TRUE),
              F_FMSY_Upr = quantile(value, probs = 0.975,
                                    na.rm = TRUE),
              F_FMSY_Lwr2 = quantile(value, probs = 0.1,
                                       na.rm = TRUE),
              F_FMSY_Upr2 = quantile(value, probs = 0.9,
                                       na.rm = TRUE)) %>%
    as.data.frame()

###@> Figure Reference Points...
p01 <- ggplot(data =
                  filter(tab03,
                         MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr,
                    ymax = F_FMSY_Upr),
                alpha = 0.5, fill = "#95B634") +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr2,
                    ymax = F_FMSY_Upr2),
                alpha = 0.2, fill = "#006750") +
    geom_line(aes(x = Year, y = F_FMSY), linewidth = 1,
              colour = "#006750") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 1) +
    labs(x = "Year", y = expression(F/F[MSY])) +
    my_theme()
p01

p02 <- ggplot(data =
                  filter(tab03,
                         MPs %in% c("Iratio", "Islope1", "GB_slope"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr,
                    ymax = F_FMSY_Upr),
                alpha = 0.5, fill = "#86BEDA") +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr2,
                    ymax = F_FMSY_Upr2),
                alpha = 0.2, fill = "#005690") +
    geom_line(aes(x = Year, y = F_FMSY), linewidth = 1,
              colour = "#005690") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 1) +
    labs(x = "Year", y = expression(F/F[MSY])) +
    my_theme()
p02

p03 <- ggplot(data =
                  filter(tab03,
                         MPs %in% c("SCA_100_40_SBMSY",
                                    "SP_100_40_SBMSY",
                                    "SPSS_100_40_SBMSY", "SP_01",
                                    "SP_02", "SP_03",
                                    "SP_04", "SP_05", "SP_06"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr,
                    ymax = F_FMSY_Upr),
                alpha = 0.5, fill = "#533D8B") +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr2,
                    ymax = F_FMSY_Upr2),
                alpha = 0.2, fill = "#5C1374") +
    geom_line(aes(x = Year, y = F_FMSY), linewidth = 1,
              colour = "#5C1374") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 2) +
    labs(x = "Year", y = expression(F/F[MSY])) +
    my_theme()
p03

ggsave("Results/Trajectory_F_FMSY_CC_ver01.svg", plot = p01,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Trajectory_F_FMSY_Index_ver01.svg", plot = p02,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Trajectory_F_FMSY_Model_ver01.svg", plot = p03,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

#####@> Extracting Catches values...
tmp03 <- list()
for(i in seq_along(MSEs03)) {
    om <- names(MSEs03)[i]
    output <- MSEs03[[i]]@Catch %>%
        structure(
            dimnames = list(
                Sim = 1:MSEs03[[i]]@nsim,
                MPs = MSEs03[[i]]@MPs,
                Year = 2020 + 1:MSEs03[[i]]@proyears
            )
        ) %>%
        reshape2::melt() %>%
        mutate(MOM = om) %>%
        as.data.frame()
    tmp03[[i]] <- output
}
Catch_output <- do.call(rbind.data.frame, tmp03)

####@> Looking to the RP - Estimating the average between all MOMs...
tab04 <- Catch_output %>%
    filter(MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                      "GB_slope", "Iratio", "Islope1",
                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                      "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                      "SP_04", "SP_05", "SP_06")) %>%
    group_by(MPs, Year) %>%
    summarise(Catch = mean(value, na.rm = TRUE),
              Catch_Lwr = quantile(value, probs = 0.025,
                                   na.rm = TRUE),
              Catch_Upr = quantile(value, probs = 0.975,
                                   na.rm = TRUE),
              Catch_Lwr2 = quantile(value, probs = 0.1,
                                   na.rm = TRUE),
              Catch_Upr2 = quantile(value, probs = 0.9,
                                   na.rm = TRUE)) %>%
    as.data.frame()

###@> Figure Reference Points...
p01 <- ggplot(data =
                  filter(tab04,
                         MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr,
                    ymax = Catch_Upr),
                alpha = 0.5, fill = "#95B634") +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr2,
                    ymax = Catch_Upr2),
                alpha = 0.2, fill = "#006750") +
    geom_line(aes(x = Year, y = Catch), linewidth = 1,
              colour = "#006750") +
    ## geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 1) +
    labs(x = "Year", y = "Catch (t)") +
    my_theme()
p01

p02 <- ggplot(data =
                  filter(tab04,
                         MPs %in% c("Iratio", "Islope1", "GB_slope"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr,
                    ymax = Catch_Upr),
                alpha = 0.5, fill = "#86BEDA") +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr2,
                    ymax = Catch_Upr2),
                alpha = 0.2, fill = "#005690") +
    geom_line(aes(x = Year, y = Catch), linewidth = 1,
              colour = "#005690") +
    ## geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 1) +
    labs(x = "Year", y = "Catch (t)") +
    my_theme()
p02

p03 <- ggplot(data =
                  filter(tab04,
                         MPs %in% c("SCA_100_40_SBMSY",
                                    "SP_100_40_SBMSY",
                                    "SPSS_100_40_SBMSY", "SP_01",
                                    "SP_02", "SP_03",
                                    "SP_04", "SP_05", "SP_06"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr,
                    ymax = Catch_Upr),
                alpha = 0.5, fill = "#533D8B") +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr2,
                    ymax = Catch_Upr2),
                alpha = 0.2, fill = "#5C1374") +
    geom_line(aes(x = Year, y = Catch), linewidth = 1,
              colour = "#5C1374") +
    ## geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 2) +
    labs(x = "Year", y = "Catch (t)") +
    my_theme()
p03

ggsave("Results/Trajectory_Catch_CC_ver01.svg", plot = p01,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Trajectory_Catch_Index_ver01.svg", plot = p02,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Trajectory_Catch_Model_ver01.svg", plot = p03,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

#####@> Extracting TAC values...
tmp04 <- list()
for(i in seq_along(MSEs03)) {
    om <- names(MSEs03)[i]
    output <- MSEs03[[i]]@TAC %>%
        structure(
            dimnames = list(
                Sim = 1:MSEs03[[i]]@nsim,
                MPs = MSEs03[[i]]@MPs,
                Year = 2020 + 1:MSEs03[[i]]@proyears
            )
        ) %>%
        reshape2::melt() %>%
        mutate(MOM = om) %>%
        as.data.frame()
    tmp04[[i]] <- output
}
TAC_output <- do.call(rbind.data.frame, tmp04)

####@> Looking to the RP - Estimating the average between all MOMs...
tab05 <- TAC_output %>%
    filter(MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                      "GB_slope", "Iratio", "Islope1",
                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                      "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                      "SP_04", "SP_05", "SP_06")) %>%
    group_by(MPs, Year) %>%
    summarise(TAC = median(value, na.rm = TRUE),
              TAC_Lwr = quantile(value, probs = 0.025,
                                 na.rm = TRUE),
              TAC_Upr = quantile(value, probs = 0.975,
                                 na.rm = TRUE),
              TAC_Lwr2 = quantile(value, probs = 0.1,
                                 na.rm = TRUE),
              TAC_Upr2 = quantile(value, probs = 0.9,
                                 na.rm = TRUE)) %>%
    as.data.frame()

###@> Figure Reference Points...
p01 <- ggplot(data =
                  filter(tab05,
                         MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = TAC_Lwr,
                    ymax = TAC_Upr),
                alpha = 0.5, fill = "#95B634") +
    geom_ribbon(aes(x = Year,
                    ymin = TAC_Lwr2,
                    ymax = TAC_Upr2),
                alpha = 0.2, fill = "#006750") +
    geom_line(aes(x = Year, y = TAC), linewidth = 1,
              colour = "#006750") +
    ## geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 1) +
    labs(x = "Year", y = "TAC (t)") +
    my_theme()
p01

p02 <- ggplot(data =
                  filter(tab05,
                         MPs %in% c("Iratio", "Islope1", "GB_slope"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = TAC_Lwr,
                    ymax = TAC_Upr),
                alpha = 0.5, fill = "#86BEDA") +
    geom_ribbon(aes(x = Year,
                    ymin = TAC_Lwr2,
                    ymax = TAC_Upr2),
                alpha = 0.2, fill = "#005690") +
    geom_line(aes(x = Year, y = TAC), linewidth = 1,
              colour = "#005690") +
    ## geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 1) +
    labs(x = "Year", y = "TAC (t)") +
    my_theme()
p02

p03 <- ggplot(data =
                  filter(tab05,
                         MPs %in% c("SCA_100_40_SBMSY",
                                    "SP_100_40_SBMSY",
                                    "SPSS_100_40_SBMSY", "SP_01",
                                    "SP_02", "SP_03",
                                    "SP_04", "SP_05", "SP_06"),
                         Year %in% 2021:2053)) +
    geom_ribbon(aes(x = Year,
                    ymin = TAC_Lwr,
                    ymax = TAC_Upr),
                alpha = 0.5, fill = "#533D8B") +
    geom_ribbon(aes(x = Year,
                    ymin = TAC_Lwr2,
                    ymax = TAC_Upr2),
                alpha = 0.2, fill = "#5C1374") +
    geom_line(aes(x = Year, y = TAC), linewidth = 1,
              colour = "#5C1374") +
    ## geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_wrap(~MPs, ncol = 2) +
    labs(x = "Year", y = "TAC (t)") +
    my_theme()
p03

ggsave("Results/Trajectory_TAC_CC_ver01.svg", plot = p01,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Trajectory_TAC_Index_ver01.svg", plot = p02,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

ggsave("Results/Trajectory_TAC_Model_ver01.svg", plot = p03,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

######@> Creating the Kobe plots (by Year and general)...

#####@> Merging SB/SBMSY and F/FMSY...
df00 <- SB_SBMSY_output %>%
    left_join(F_FMSY_output, by = c("Sim", "MPs", "Year", "OM")) %>%
    filter(MPs %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                      "GB_slope", "Iratio", "Islope1",
                      "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                      "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                      "SP_04", "SP_05", "SP_06")) %>%
    rename("SB_SBMSY" = value.x, "F_FMSY" = value.y)

#####@> Estimating Kobe values for the last year...
tmp <- df00 %>%
    filter(Year == 2053) %>%
    group_by(Sim, MPs) %>%
    summarise(SB_SBMSY = mean(SB_SBMSY, na.rm = TRUE),
              F_FMSY = mean(F_FMSY, na.rm = TRUE)) %>%
    mutate(OF = ifelse(F_FMSY <= 1, FALSE, TRUE),
           OFD = ifelse(SB_SBMSY <= 1, TRUE, FALSE))

#####@> Percentage of cases...
valdf <- tmp %>%
    group_by(MPs) %>%
    summarise(BL = sum(OF == FALSE & OFD == TRUE)/100 * 100,
              BR = sum(OF == FALSE & OFD == FALSE)/100 * 100,
              TL = sum(OF == TRUE & OFD == TRUE)/100 * 100,
              TR = sum(OF == TRUE & OFD == FALSE)/100 * 100)
valdf <- valdf %>% tidyr::pivot_longer(., cols = 2:5)
valdf$x <- -Inf
valdf$y <- -Inf
valdf$y[valdf$name == "TL"] <- Inf
valdf$y[valdf$name == "TR"] <- Inf
valdf$x[valdf$name == "BR"] <- Inf
valdf$x[valdf$name == "TR"] <- Inf
valdf$value <- round(valdf$value, 2)
valdf$value <- paste0(valdf$value, "%")
valdf$hjustvar <- -2
valdf$vjustvar <- -2
valdf$hjustvar[valdf$name == "TL"] <- -1
valdf$hjustvar[valdf$name == "TR"] <- 2
valdf$hjustvar[valdf$name == "BL"] <- -1
valdf$hjustvar[valdf$name == "BR"] <- 2
valdf$vjustvar[valdf$name == "TL"] <- 2
valdf$vjustvar[valdf$name == "TR"] <- 2
valdf$vjustvar[valdf$name == "BL"] <- -2
valdf$vjustvar[valdf$name == "BR"] <- -2

#####@> Colors to figure...
kobe_df <- bind_rows(data.frame(x = c(0, 0, 1, 1),
                                y = c(0, 1, 1, 0),
                                fill = "bl"),
                     data.frame(x = c(1, 1, 2, 2),
                                y = c(0, 1, 1, 0),
                                fill = "br"),
                     data.frame(x = c(0, 0, 1, 1),
                                y = c(1, 2, 2, 1),
                                fill = "tl"),
                     data.frame(x = c(1, 1, 2, 2),
                                y = c(1, 2, 2, 1),
                                fill = "tr"))
kobe_df$alpha <- 0.3

#####@> Figure dots...
p00 <- ggplot() +
    geom_polygon(data = kobe_df,
                 aes(x = x, y = y, fill = fill, alpha = alpha)) +
    scale_fill_manual(values = c("#F8DC7A", "#67C18B", "#D8775D",
                                 "#FDBD56")) +
    geom_point(data = tmp,
               aes(x = SB_SBMSY, y = F_FMSY), alpha = 0.5) +
    geom_text(data = valdf, fontface = "bold", size = 4,
              aes(x = x, y = y, label = value,
                  hjust = hjustvar, vjust = vjustvar),
              colour = "black") +
    expand_limits(x = c(0, 2), y = c(0, 2)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 2)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 2)) +
    geom_hline(yintercept = 1, color = "darkgray", linetype = 2) +
    geom_vline(xintercept = 1, color = "darkgray", linetype = 2) +
    labs(x = expression(SB/SB[MSY]), y = expression(F/F[MSY])) +
    facet_wrap(~MPs) +
    my_theme() +
    theme(legend.position = "none")
p00

ggsave("Results/Kobe_plot_Last_Year_2053_ver01.svg", plot = p00,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

######@> Estimating Kobe values per year...
tmp <- df00 %>%
    group_by(Sim, Year, MPs) %>%
    summarise(SB_SBMSY = mean(SB_SBMSY, na.rm = TRUE),
              F_FMSY = mean(F_FMSY, na.rm = TRUE)) %>%
    mutate(OF = ifelse(F_FMSY <= 1, FALSE, TRUE),
           OFD = ifelse(SB_SBMSY <= 1, TRUE, FALSE))

#####@> Proportions by year...
valdf <- tmp %>%
    group_by(Year, MPs) %>%
    summarise(Yellow = sum(OF == FALSE & OFD == TRUE)/100 * 100,
              Green = sum(OF == FALSE & OFD == FALSE)/100 * 100,
              Red = sum(OF == TRUE & OFD == TRUE)/100 * 100,
              Orange = sum(OF == TRUE & OFD == FALSE)/100 * 100) %>%
    pivot_longer(names_to = "Cond", values_to = "Perc", 3:6) %>%
    mutate(Cond = factor(Cond, levels = c("Green", "Yellow", "Orange",
                                          "Red")))

p00 <- ggplot() +
    geom_bar(data = filter(valdf, Year %in% 2021:2053),
             aes(x = Year, y = Perc, fill = Cond),
             stat = "identity", colour = "black",
             width = 1) +
    facet_wrap(~MPs) +
    scale_fill_manual(values = c("#67C18B", "#F8DC7A", "#FDBD56",
                                 "#D8775D")) +
    labs(x = "Year", y = "%") +
    my_theme() +
    theme(legend.position = "none")
p00

p01 <- ggplot() +
    geom_area(data = filter(valdf, Year %in% 2021:2053),
              aes(x = Year, y = Perc, fill = Cond),
              stat = "identity", colour = "black") +
    facet_wrap(~MPs) +
    scale_fill_manual(values = c("#67C18B", "#F8DC7A", "#FDBD56",
                                 "#D8775D")) +
    labs(x = "Year", y = "%") +
    my_theme() +
    theme(legend.position = "none")
p01

ggsave("Results/Kobe_plot_by_year_Bar_ver01.svg", plot = p00,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 30)

ggsave("Results/Kobe_plot_by_year_Area_ver01.svg", plot = p01,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 30)

########################################################################
######@> Projections...

######@>----------------------------------------------------------------
######@> Projections OMs 1 - 9...

#####@> Extraction the projections...
projDB <- data.frame(OM = NULL,
                     iter = NULL,
                     mp = NULL,
                     year = NULL,
                     value = NULL,
                     mp_name = NULL,
                     Type = NULL,
                     id_year = NULL,
                     real_year = NULL,
                     ref = NULL)
for(i in 1:9) {
    print(paste0("Compiling data from OM ", i))
    out <- get_ts(MSEs[[i]], this_year = 2020)
    projDB <- rbind(projDB, data.frame(OM = i,
                                       iter = out$iter,
                                       mp = out$mp,
                                       year = out$year,
                                       value = out$value,
                                       mp_name = out$mp_name,
                                       Type = out$Type,
                                       id_year = out$real_year,
                                       real_year = out$real_year,
                                       ref = out$ref))
}

#####@> Estimating the summary quantities...
summProj <- get_ts_quantiles(projDB)

#####@> Segregating base by Type and Year...
histBBMSY <- filter(summProj, Type == "B_BMSY",
                    real_year %in% 1952:2020)
projBBMSY <- filter(summProj, Type == "B_BMSY",
                    real_year %in% 2020:2053,
                    mp_name %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                                   "GB_slope", "Iratio", "Islope1",
                                   "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                                   "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                                   "SP_04", "SP_05", "SP_06"))
histFFMSY <- filter(summProj, Type == "F_FMSY",
                    real_year %in% 1952:2020)
projFFMSY <- filter(summProj, Type == "F_FMSY",
                    real_year %in% 2020:2053,
                    mp_name %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                                   "GB_slope", "Iratio", "Islope1",
                                   "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                                   "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                                   "SP_04", "SP_05", "SP_06"))
histCatch <- filter(summProj, Type == "Catch",
                    real_year %in% 1952:2020)
projCatch <- filter(summProj, Type == "Catch",
                    real_year %in% 2020:2053,
                    mp_name %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                                   "GB_slope", "Iratio", "Islope1",
                                   "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                                   "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                                   "SP_04", "SP_05", "SP_06"))

p00 <- ggplot() +
    geom_ribbon(data = histBBMSY,
                aes(x = real_year, ymin = l, ymax = u),
                alpha = 0.9, fill = "lightgray") +
    geom_line(data = histBBMSY, aes(x = real_year, y = ll),
              linetype = "dashed", colour = "lightgray") +
    geom_line(data = histBBMSY, aes(x = real_year, y = uu),
              linetype = "dashed", colour = "lightgray") +
    geom_line(data = histBBMSY, aes(x = real_year, y = m),
              linewidth = 1.2, colour = "gray50") +
    geom_ribbon(data = projBBMSY,
                aes(x = real_year, ymin = l, ymax = u,
                    fill = mp_name),
                alpha = 0.1) +
    geom_line(data = projBBMSY,
              aes(x = real_year, y = m, colour = mp_name),
              linewidth = 1.2) +
    geom_text_repel(data = filter(projBBMSY, real_year == 2053),
                    aes(x = real_year, y = m, colour = mp_name,
                        label = mp_name), nudge_x = 4, force = 5) +
    geom_hline(yintercept = 1, colour = "green") +
    geom_hline(yintercept = 0.4, colour = "red") +
    geom_vline(xintercept = 2020, colour = "lightgray",
               linetype = "dashed") +
    scale_y_continuous(breaks = seq(0, 4, 0.4), limits = c(0, 4),
                       expand = c(0, 0)) +
    scale_x_continuous(breaks = seq(1952, 2060, 4), limits = c(1952, 2060)) +
    labs(x = "Year",
         y = "Spawning Stock Biomass relative to MSY levels",
         colour = "") +
    theme_pbs(base_size = 18) +
    theme(legend.position = "none")
p00

p01 <- ggplot() +
    geom_ribbon(data = histFFMSY,
                aes(x = real_year, ymin = l, ymax = u),
                alpha = 0.9, fill = "lightgray") +
    geom_line(data = histFFMSY, aes(x = real_year, y = ll),
              linetype = "dashed", colour = "lightgray") +
    geom_line(data = histFFMSY, aes(x = real_year, y = uu),
              linetype = "dashed", colour = "lightgray") +
    geom_line(data = histFFMSY, aes(x = real_year, y = m),
              linewidth = 1.2, colour = "gray50") +
    geom_ribbon(data = projFFMSY,
                aes(x = real_year, ymin = l, ymax = u,
                    fill = mp_name),
                alpha = 0.1) +
    geom_line(data = projFFMSY,
              aes(x = real_year, y = m, colour = mp_name),
              linewidth = 1.2) +
    geom_text_repel(data = filter(projFFMSY, real_year == 2050),
                    aes(x = real_year, y = m, colour = mp_name,
                        label = mp_name), nudge_x = 4) +
    geom_hline(yintercept = 1, colour = "red") +
    ## geom_hline(yintercept = 0.4, colour = "red") +
    geom_vline(xintercept = 2020, colour = "lightgray",
               linetype = "dashed") +
    scale_y_continuous(breaks = seq(0, 7, 1), limits = c(0, 8),
                       expand = c(0, 0)) +
    scale_x_continuous(breaks = seq(1952, 2060, 4)) +
    labs(x = "Year",
         y = "F/FMSY",
         colour = "") +
    theme_pbs(base_size = 18) +
    theme(legend.position = "none")
p01

p02 <- ggplot() +
    geom_ribbon(data = histCatch,
                aes(x = real_year, ymin = l, ymax = u),
                alpha = 0.9, fill = "lightgray") +
    geom_line(data = histCatch, aes(x = real_year, y = ll),
              linetype = "dashed", colour = "lightgray") +
    geom_line(data = histCatch, aes(x = real_year, y = uu),
              linetype = "dashed", colour = "lightgray") +
    geom_line(data = histCatch, aes(x = real_year, y = m),
              linewidth = 1.2, colour = "gray50") +
    geom_ribbon(data = projCatch,
                aes(x = real_year, ymin = l, ymax = u,
                    fill = mp_name),
                alpha = 0.1) +
    geom_line(data = projCatch,
              aes(x = real_year, y = m, colour = mp_name),
              linewidth = 1.2) +
    geom_text_repel(data = filter(projCatch, real_year == 2050),
                    aes(x = real_year, y = m, colour = mp_name,
                        label = mp_name), nudge_x = 5, nudge_y = 2) +
    ## geom_hline(yintercept = 1, colour = "green") +
    ## geom_hline(yintercept = 0.4, colour = "red") +
    geom_vline(xintercept = 2020, colour = "lightgray",
               linetype = "dashed") +
    scale_y_continuous(breaks = seq(0, 80000, 10000), limits = c(0, 80000),
                       expand = c(0, 0)) +
    scale_x_continuous(breaks = seq(1952, 2060, 4)) +
    labs(x = "Year",
         y = "Catches (t)",
         colour = "") +
    theme_pbs(base_size = 18) +
    theme(legend.position = "none")
p02

######@> Catch...
catch <- readxl::read_excel("../WSKJ_Input_DB_ver00.xlsx",
                            sheet = "CATCH")

## "CC_20kt", "CC_30kt", "CC_40kt",
## "GB_slope", "Iratio", "Islope1",
## "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
## "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
## "SP_04", "SP_05", "SP_06"

CMPIslope <- filter(projCatch, mp_name == "Islope1")
CMPIratio <- filter(projCatch, mp_name == "Iratio")
CMPGB_slope <- filter(projCatch, mp_name == "GB_slope")
CMPCC_20kt <- filter(projCatch, mp_name == "CC_20kt")
CMPCC_30kt <- filter(projCatch, mp_name == "CC_30kt")
CMPCC_40kt <- filter(projCatch, mp_name == "CC_40kt")
CMPSCA <- filter(projCatch, mp_name == "SCA_100_40_SBMSY")
CMPSP <- filter(projCatch, mp_name == "SP_100_40_SBMSY")
CMPSPSS <- filter(projCatch, mp_name == "SPSS_100_40_SBMSY")
CMPSP_01 <- filter(projCatch, mp_name == "SP_01")
CMPSP_02 <- filter(projCatch, mp_name == "SP_02")
CMPSP_03 <- filter(projCatch, mp_name == "SP_03")
CMPSP_04 <- filter(projCatch, mp_name == "SP_04")
CMPSP_05 <- filter(projCatch, mp_name == "SP_05")
CMPSP_06 <- filter(projCatch, mp_name == "SP_06")

#####@> Correction...
CMPIslope$u <- CMPIslope$u + (catch$Catch[catch$Year==2020] -
                              CMPIslope$u[1])
CMPIslope$uu <- CMPIslope$uu + (catch$Catch[catch$Year==2020] -
                                CMPIslope$uu[1])
CMPIslope$m <- CMPIslope$m + (catch$Catch[catch$Year==2020] -
                              CMPIslope$m[1])
CMPIslope$l <- CMPIslope$l + (catch$Catch[catch$Year==2020] -
                              CMPIslope$l[1])
CMPIslope$ll <- CMPIslope$ll + (catch$Catch[catch$Year==2020] -
                                CMPIslope$ll[1])

CMPIratio$u <- CMPIratio$u + (catch$Catch[catch$Year==2020] -
                              CMPIratio$u[1])
CMPIratio$uu <- CMPIratio$uu + (catch$Catch[catch$Year==2020] -
                                CMPIratio$uu[1])
CMPIratio$m <- CMPIratio$m + (catch$Catch[catch$Year==2020] -
                              CMPIratio$m[1])
CMPIratio$l <- CMPIratio$l + (catch$Catch[catch$Year==2020] -
                              CMPIratio$l[1])
CMPIratio$ll <- CMPIratio$ll + (catch$Catch[catch$Year==2020] -
                                CMPIratio$ll[1])

CMPGB_slope$u <- CMPGB_slope$u + (catch$Catch[catch$Year==2020] -
                              CMPGB_slope$u[1])
CMPGB_slope$uu <- CMPGB_slope$uu + (catch$Catch[catch$Year==2020] -
                                CMPGB_slope$uu[1])
CMPGB_slope$m <- CMPGB_slope$m + (catch$Catch[catch$Year==2020] -
                              CMPGB_slope$m[1])
CMPGB_slope$l <- CMPGB_slope$l + (catch$Catch[catch$Year==2020] -
                              CMPGB_slope$l[1])
CMPGB_slope$ll <- CMPGB_slope$ll + (catch$Catch[catch$Year==2020] -
                                CMPGB_slope$ll[1])

CMPCC_20kt$u[1] <- CMPCC_20kt$u[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_20kt$u[1])
CMPCC_20kt$uu[1] <- CMPCC_20kt$uu[1] + (catch$Catch[catch$Year==2020] -
                                        CMPCC_20kt$uu[1])
CMPCC_20kt$m[1] <- CMPCC_20kt$m[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_20kt$m[1])
CMPCC_20kt$l[1] <- CMPCC_20kt$l[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_20kt$l[1])
CMPCC_20kt$ll[1] <- CMPCC_20kt$ll[1] + (catch$Catch[catch$Year==2020] -
                                        CMPCC_20kt$ll[1])

CMPCC_30kt$u[1] <- CMPCC_30kt$u[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_30kt$u[1])
CMPCC_30kt$uu[1] <- CMPCC_30kt$uu[1] + (catch$Catch[catch$Year==2020] -
                                        CMPCC_30kt$uu[1])
CMPCC_30kt$m[1] <- CMPCC_30kt$m[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_30kt$m[1])
CMPCC_30kt$l[1] <- CMPCC_30kt$l[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_30kt$l[1])
CMPCC_30kt$ll[1] <- CMPCC_30kt$ll[1] + (catch$Catch[catch$Year==2020] -
                                        CMPCC_30kt$ll[1])

CMPCC_40kt$u[1] <- CMPCC_40kt$u[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_40kt$u[1])
CMPCC_40kt$uu[1] <- CMPCC_40kt$uu[1] + (catch$Catch[catch$Year==2020] -
                                        CMPCC_40kt$uu[1])
CMPCC_40kt$m[1] <- CMPCC_40kt$m[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_40kt$m[1])
CMPCC_40kt$l[1] <- CMPCC_40kt$l[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_40kt$l[1])
CMPCC_40kt$ll[1] <- CMPCC_40kt$ll[1] + (catch$Catch[catch$Year==2020] -
                                        CMPCC_40kt$ll[1])

CMPSCA$u[1] <- CMPSCA$u[1] + (catch$Catch[catch$Year==2020] -
                              CMPSCA$u[1])
CMPSCA$uu[1] <- CMPSCA$uu[1] + (catch$Catch[catch$Year==2020] -
                                CMPSCA$uu[1])
CMPSCA$m[1] <- CMPSCA$m[1] + (catch$Catch[catch$Year==2020] -
                              CMPSCA$m[1])
CMPSCA$l[1] <- CMPSCA$l[1] + (catch$Catch[catch$Year==2020] -
                              CMPSCA$l[1])
CMPSCA$ll[1] <- CMPSCA$ll[1] + (catch$Catch[catch$Year==2020] -
                                CMPSCA$ll[1])

CMPSP$u[1] <- CMPSP$u[1] + (catch$Catch[catch$Year==2020] -
                            CMPSP$u[1])
CMPSP$uu[1] <- CMPSP$uu[1] + (catch$Catch[catch$Year==2020] -
                              CMPSP$uu[1])
CMPSP$m[1] <- CMPSP$m[1] + (catch$Catch[catch$Year==2020] -
                            CMPSP$m[1])
CMPSP$l[1] <- CMPSP$l[1] + (catch$Catch[catch$Year==2020] -
                            CMPSP$l[1])
CMPSP$ll[1] <- CMPSP$ll[1] + (catch$Catch[catch$Year==2020] -
                              CMPSP$ll[1])

CMPSPSS$u[1] <- CMPSPSS$u[1] + (catch$Catch[catch$Year==2020] -
                            CMPSPSS$u[1])
CMPSPSS$uu[1] <- CMPSPSS$uu[1] + (catch$Catch[catch$Year==2020] -
                              CMPSPSS$uu[1])
CMPSPSS$m[1] <- CMPSPSS$m[1] + (catch$Catch[catch$Year==2020] -
                            CMPSPSS$m[1])
CMPSPSS$l[1] <- CMPSPSS$l[1] + (catch$Catch[catch$Year==2020] -
                            CMPSPSS$l[1])
CMPSPSS$ll[1] <- CMPSPSS$ll[1] + (catch$Catch[catch$Year==2020] -
                              CMPSPSS$ll[1])

CMPSP_01$u[1] <- CMPSP_01$u[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_01$u[1])
CMPSP_01$uu[1] <- CMPSP_01$uu[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_01$uu[1])
CMPSP_01$m[1] <- CMPSP_01$m[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_01$m[1])
CMPSP_01$l[1] <- CMPSP_01$l[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_01$l[1])
CMPSP_01$ll[1] <- CMPSP_01$ll[1] + (catch$Catch[catch$Year==2020] -
                                    CMPSP_01$ll[1])


CMPSP_02$u[1] <- CMPSP_02$u[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_02$u[1])
CMPSP_02$uu[1] <- CMPSP_02$uu[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_02$uu[1])
CMPSP_02$m[1] <- CMPSP_02$m[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_02$m[1])
CMPSP_02$l[1] <- CMPSP_02$l[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_02$l[1])
CMPSP_02$ll[1] <- CMPSP_02$ll[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_02$ll[1])

CMPSP_03$u[1] <- CMPSP_03$u[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_03$u[1])
CMPSP_03$uu[1] <- CMPSP_03$uu[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_03$uu[1])
CMPSP_03$m[1] <- CMPSP_03$m[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_03$m[1])
CMPSP_03$l[1] <- CMPSP_03$l[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_03$l[1])
CMPSP_03$ll[1] <- CMPSP_03$ll[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_03$ll[1])

CMPSP_04$u[1] <- CMPSP_04$u[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_04$u[1])
CMPSP_04$uu[1] <- CMPSP_04$uu[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_04$uu[1])
CMPSP_04$m[1] <- CMPSP_04$m[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_04$m[1])
CMPSP_04$l[1] <- CMPSP_04$l[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_04$l[1])
CMPSP_04$ll[1] <- CMPSP_04$ll[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_04$ll[1])

CMPSP_05$u[1] <- CMPSP_05$u[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_05$u[1])
CMPSP_05$uu[1] <- CMPSP_05$uu[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_05$uu[1])
CMPSP_05$m[1] <- CMPSP_05$m[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_05$m[1])
CMPSP_05$l[1] <- CMPSP_05$l[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_05$l[1])
CMPSP_05$ll[1] <- CMPSP_05$ll[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_05$ll[1])

CMPSP_06$u[1] <- CMPSP_06$u[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_06$u[1])
CMPSP_06$uu[1] <- CMPSP_06$uu[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_06$uu[1])
CMPSP_06$m[1] <- CMPSP_06$m[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_06$m[1])
CMPSP_06$l[1] <- CMPSP_06$l[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_06$l[1])
CMPSP_06$ll[1] <- CMPSP_06$ll[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_06$ll[1])

#####@> Merging...
projCatch02 <- gtools::smartbind(data.frame(CMPIslope),
                                 data.frame(CMPIratio),
                                 data.frame(CMPGB_slope),
                                 data.frame(CMPCC_20kt),
                                 data.frame(CMPCC_30kt),
                                 data.frame(CMPCC_40kt),
                                 data.frame(CMPSCA),
                                 data.frame(CMPSP),
                                 data.frame(CMPSPSS),
                                 data.frame(CMPSP_01),
                                 data.frame(CMPSP_02),
                                 data.frame(CMPSP_03),
                                 data.frame(CMPSP_04),
                                 data.frame(CMPSP_05),
                                 data.frame(CMPSP_06))

p02 <- ggplot() +
    geom_line(data = catch, aes(x = Year, y = Catch),
              linewidth = 1.2, colour = "gray50") +
    geom_ribbon(data = projCatch02,
                aes(x = real_year, ymin = l, ymax = u,
                    fill = mp_name),
                alpha = 0.1) +
    geom_line(data = projCatch02,
              aes(x = real_year, y = m, colour = mp_name),
              linewidth = 1.2) +
    geom_text_repel(data = filter(projCatch02, real_year == 2053),
                    aes(x = real_year, y = m, colour = mp_name,
                        label = mp_name), nudge_x = 6, nudge_y = 500,
                    force = 5) +
    ## geom_hline(yintercept = 1, colour = "green") +
    ## geom_hline(yintercept = 0.4, colour = "red") +
    geom_vline(xintercept = 2020, colour = "lightgray",
               linetype = "dashed") +
    scale_y_continuous(breaks = seq(0, 80000, 10000), limits = c(0, 80000),
                       expand = c(0, 0)) +
    scale_x_continuous(breaks = seq(1952, 2060, 4)) +
    ## scale_fill_viridis(discrete = TRUE, option = 10) +
    ## scale_colour_viridis(discrete = TRUE, option = 10) +
    labs(x = "Year",
         y = "Catches (t)",
         colour = "") +
    theme_pbs(base_size = 18) +
    theme(legend.position = "none")
p02


#####@> Exporting figures...
ggsave("Results/No_Error_Imp/TS_Catches_Rescaled_ver00.svg",
       plot = p02, device = "svg", units = "cm", w = 45, h = 30,
       bg = "white", dpi = 600)

ggsave("Results/No_Error_Imp/TS_BBMSY_ver00.svg",
       plot = p00, device = "svg", units = "cm", w = 45, h = 30,
       bg = "white", dpi = 600)

######@>----------------------------------------------------------------
######@> Projections OMs 10 - 18...

#####@> Extraction the projections...
projDB <- data.frame(OM = NULL,
                     iter = NULL,
                     mp = NULL,
                     year = NULL,
                     value = NULL,
                     mp_name = NULL,
                     Type = NULL,
                     id_year = NULL,
                     real_year = NULL,
                     ref = NULL)
for(i in 1:9) {
    print(paste0("Compiling data from OM ", i))
    out <- get_ts(MSEs02[[i]], this_year = 2020)
    projDB <- rbind(projDB, data.frame(OM = i + 9,
                                       iter = out$iter,
                                       mp = out$mp,
                                       year = out$year,
                                       value = out$value,
                                       mp_name = out$mp_name,
                                       Type = out$Type,
                                       id_year = out$real_year,
                                       real_year = out$real_year,
                                       ref = out$ref))
}

#####@> Estimating the summary quantities...
summProj <- get_ts_quantiles(projDB)

#####@> Segregating base by Type and Year...
histBBMSY <- filter(summProj, Type == "B_BMSY",
                    real_year %in% 1952:2020)
projBBMSY <- filter(summProj, Type == "B_BMSY",
                    real_year %in% 2020:2053,
                    mp_name %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                                   "GB_slope", "Iratio", "Islope1",
                                   "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                                   "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                                   "SP_04", "SP_05", "SP_06"))
histFFMSY <- filter(summProj, Type == "F_FMSY",
                    real_year %in% 1952:2020)
projFFMSY <- filter(summProj, Type == "F_FMSY",
                    real_year %in% 2020:2053,
                    mp_name %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                                   "GB_slope", "Iratio", "Islope1",
                                   "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                                   "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                                   "SP_04", "SP_05", "SP_06"))
histCatch <- filter(summProj, Type == "Catch",
                    real_year %in% 1952:2020)
projCatch <- filter(summProj, Type == "Catch",
                    real_year %in% 2020:2053,
                    mp_name %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                                   "GB_slope", "Iratio", "Islope1",
                                   "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                                   "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                                   "SP_04", "SP_05", "SP_06"))

p00 <- ggplot() +
    geom_ribbon(data = histBBMSY,
                aes(x = real_year, ymin = l, ymax = u),
                alpha = 0.9, fill = "lightgray") +
    geom_line(data = histBBMSY, aes(x = real_year, y = ll),
              linetype = "dashed", colour = "lightgray") +
    geom_line(data = histBBMSY, aes(x = real_year, y = uu),
              linetype = "dashed", colour = "lightgray") +
    geom_line(data = histBBMSY, aes(x = real_year, y = m),
              linewidth = 1.2, colour = "gray50") +
    geom_ribbon(data = projBBMSY,
                aes(x = real_year, ymin = l, ymax = u,
                    fill = mp_name),
                alpha = 0.1) +
    geom_line(data = projBBMSY,
              aes(x = real_year, y = m, colour = mp_name),
              linewidth = 1.2) +
    geom_text_repel(data = filter(projBBMSY, real_year == 2053),
                    aes(x = real_year, y = m, colour = mp_name,
                        label = mp_name), nudge_x = 4, force = 5) +
    geom_hline(yintercept = 1, colour = "green") +
    geom_hline(yintercept = 0.4, colour = "red") +
    geom_vline(xintercept = 2020, colour = "lightgray",
               linetype = "dashed") +
    scale_y_continuous(breaks = seq(0, 4, 0.4), limits = c(0, 4),
                       expand = c(0, 0)) +
    scale_x_continuous(breaks = seq(1952, 2060, 4), limits = c(1952, 2060)) +
    labs(x = "Year",
         y = "Spawning Stock Biomass relative to MSY levels",
         colour = "") +
    theme_pbs(base_size = 18) +
    theme(legend.position = "none")
p00

p01 <- ggplot() +
    geom_ribbon(data = histFFMSY,
                aes(x = real_year, ymin = l, ymax = u),
                alpha = 0.9, fill = "lightgray") +
    geom_line(data = histFFMSY, aes(x = real_year, y = ll),
              linetype = "dashed", colour = "lightgray") +
    geom_line(data = histFFMSY, aes(x = real_year, y = uu),
              linetype = "dashed", colour = "lightgray") +
    geom_line(data = histFFMSY, aes(x = real_year, y = m),
              linewidth = 1.2, colour = "gray50") +
    geom_ribbon(data = projFFMSY,
                aes(x = real_year, ymin = l, ymax = u,
                    fill = mp_name),
                alpha = 0.1) +
    geom_line(data = projFFMSY,
              aes(x = real_year, y = m, colour = mp_name),
              linewidth = 1.2) +
    geom_text_repel(data = filter(projFFMSY, real_year == 2050),
                    aes(x = real_year, y = m, colour = mp_name,
                        label = mp_name), nudge_x = 4) +
    geom_hline(yintercept = 1, colour = "red") +
    ## geom_hline(yintercept = 0.4, colour = "red") +
    geom_vline(xintercept = 2020, colour = "lightgray",
               linetype = "dashed") +
    scale_y_continuous(breaks = seq(0, 7, 1), limits = c(0, 8),
                       expand = c(0, 0)) +
    scale_x_continuous(breaks = seq(1952, 2060, 4)) +
    labs(x = "Year",
         y = "F/FMSY",
         colour = "") +
    theme_pbs(base_size = 18) +
    theme(legend.position = "none")
p01

p02 <- ggplot() +
    geom_ribbon(data = histCatch,
                aes(x = real_year, ymin = l, ymax = u),
                alpha = 0.9, fill = "lightgray") +
    geom_line(data = histCatch, aes(x = real_year, y = ll),
              linetype = "dashed", colour = "lightgray") +
    geom_line(data = histCatch, aes(x = real_year, y = uu),
              linetype = "dashed", colour = "lightgray") +
    geom_line(data = histCatch, aes(x = real_year, y = m),
              linewidth = 1.2, colour = "gray50") +
    geom_ribbon(data = projCatch,
                aes(x = real_year, ymin = l, ymax = u,
                    fill = mp_name),
                alpha = 0.1) +
    geom_line(data = projCatch,
              aes(x = real_year, y = m, colour = mp_name),
              linewidth = 1.2) +
    geom_text_repel(data = filter(projCatch, real_year == 2050),
                    aes(x = real_year, y = m, colour = mp_name,
                        label = mp_name), nudge_x = 5, nudge_y = 2) +
    ## geom_hline(yintercept = 1, colour = "green") +
    ## geom_hline(yintercept = 0.4, colour = "red") +
    geom_vline(xintercept = 2020, colour = "lightgray",
               linetype = "dashed") +
    scale_y_continuous(breaks = seq(0, 80000, 10000), limits = c(0, 80000),
                       expand = c(0, 0)) +
    scale_x_continuous(breaks = seq(1952, 2060, 4)) +
    labs(x = "Year",
         y = "Catches (t)",
         colour = "") +
    theme_pbs(base_size = 18) +
    theme(legend.position = "none")
p02

######@> Catch...
catch <- readxl::read_excel("../WSKJ_Input_DB_ver00.xlsx",
                            sheet = "CATCH")

## "CC_20kt", "CC_30kt", "CC_40kt",
## "GB_slope", "Iratio", "Islope1",
## "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
## "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
## "SP_04", "SP_05", "SP_06"

CMPIslope <- filter(projCatch, mp_name == "Islope1")
CMPIratio <- filter(projCatch, mp_name == "Iratio")
CMPGB_slope <- filter(projCatch, mp_name == "GB_slope")
CMPCC_20kt <- filter(projCatch, mp_name == "CC_20kt")
CMPCC_30kt <- filter(projCatch, mp_name == "CC_30kt")
CMPCC_40kt <- filter(projCatch, mp_name == "CC_40kt")
CMPSCA <- filter(projCatch, mp_name == "SCA_100_40_SBMSY")
CMPSP <- filter(projCatch, mp_name == "SP_100_40_SBMSY")
CMPSPSS <- filter(projCatch, mp_name == "SPSS_100_40_SBMSY")
CMPSP_01 <- filter(projCatch, mp_name == "SP_01")
CMPSP_02 <- filter(projCatch, mp_name == "SP_02")
CMPSP_03 <- filter(projCatch, mp_name == "SP_03")
CMPSP_04 <- filter(projCatch, mp_name == "SP_04")
CMPSP_05 <- filter(projCatch, mp_name == "SP_05")
CMPSP_06 <- filter(projCatch, mp_name == "SP_06")

#####@> Correction...
CMPIslope$u <- CMPIslope$u + (catch$Catch[catch$Year==2020] -
                              CMPIslope$u[1])
CMPIslope$uu <- CMPIslope$uu + (catch$Catch[catch$Year==2020] -
                                CMPIslope$uu[1])
CMPIslope$m <- CMPIslope$m + (catch$Catch[catch$Year==2020] -
                              CMPIslope$m[1])
CMPIslope$l <- CMPIslope$l + (catch$Catch[catch$Year==2020] -
                              CMPIslope$l[1])
CMPIslope$ll <- CMPIslope$ll + (catch$Catch[catch$Year==2020] -
                                CMPIslope$ll[1])

CMPIratio$u <- CMPIratio$u + (catch$Catch[catch$Year==2020] -
                              CMPIratio$u[1])
CMPIratio$uu <- CMPIratio$uu + (catch$Catch[catch$Year==2020] -
                                CMPIratio$uu[1])
CMPIratio$m <- CMPIratio$m + (catch$Catch[catch$Year==2020] -
                              CMPIratio$m[1])
CMPIratio$l <- CMPIratio$l + (catch$Catch[catch$Year==2020] -
                              CMPIratio$l[1])
CMPIratio$ll <- CMPIratio$ll + (catch$Catch[catch$Year==2020] -
                                CMPIratio$ll[1])

CMPGB_slope$u <- CMPGB_slope$u + (catch$Catch[catch$Year==2020] -
                              CMPGB_slope$u[1])
CMPGB_slope$uu <- CMPGB_slope$uu + (catch$Catch[catch$Year==2020] -
                                CMPGB_slope$uu[1])
CMPGB_slope$m <- CMPGB_slope$m + (catch$Catch[catch$Year==2020] -
                              CMPGB_slope$m[1])
CMPGB_slope$l <- CMPGB_slope$l + (catch$Catch[catch$Year==2020] -
                              CMPGB_slope$l[1])
CMPGB_slope$ll <- CMPGB_slope$ll + (catch$Catch[catch$Year==2020] -
                                CMPGB_slope$ll[1])

CMPCC_20kt$u[1] <- CMPCC_20kt$u[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_20kt$u[1])
CMPCC_20kt$uu[1] <- CMPCC_20kt$uu[1] + (catch$Catch[catch$Year==2020] -
                                        CMPCC_20kt$uu[1])
CMPCC_20kt$m[1] <- CMPCC_20kt$m[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_20kt$m[1])
CMPCC_20kt$l[1] <- CMPCC_20kt$l[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_20kt$l[1])
CMPCC_20kt$ll[1] <- CMPCC_20kt$ll[1] + (catch$Catch[catch$Year==2020] -
                                        CMPCC_20kt$ll[1])

CMPCC_30kt$u[1] <- CMPCC_30kt$u[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_30kt$u[1])
CMPCC_30kt$uu[1] <- CMPCC_30kt$uu[1] + (catch$Catch[catch$Year==2020] -
                                        CMPCC_30kt$uu[1])
CMPCC_30kt$m[1] <- CMPCC_30kt$m[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_30kt$m[1])
CMPCC_30kt$l[1] <- CMPCC_30kt$l[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_30kt$l[1])
CMPCC_30kt$ll[1] <- CMPCC_30kt$ll[1] + (catch$Catch[catch$Year==2020] -
                                        CMPCC_30kt$ll[1])

CMPCC_40kt$u[1] <- CMPCC_40kt$u[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_40kt$u[1])
CMPCC_40kt$uu[1] <- CMPCC_40kt$uu[1] + (catch$Catch[catch$Year==2020] -
                                        CMPCC_40kt$uu[1])
CMPCC_40kt$m[1] <- CMPCC_40kt$m[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_40kt$m[1])
CMPCC_40kt$l[1] <- CMPCC_40kt$l[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_40kt$l[1])
CMPCC_40kt$ll[1] <- CMPCC_40kt$ll[1] + (catch$Catch[catch$Year==2020] -
                                        CMPCC_40kt$ll[1])

CMPSCA$u[1] <- CMPSCA$u[1] + (catch$Catch[catch$Year==2020] -
                              CMPSCA$u[1])
CMPSCA$uu[1] <- CMPSCA$uu[1] + (catch$Catch[catch$Year==2020] -
                                CMPSCA$uu[1])
CMPSCA$m[1] <- CMPSCA$m[1] + (catch$Catch[catch$Year==2020] -
                              CMPSCA$m[1])
CMPSCA$l[1] <- CMPSCA$l[1] + (catch$Catch[catch$Year==2020] -
                              CMPSCA$l[1])
CMPSCA$ll[1] <- CMPSCA$ll[1] + (catch$Catch[catch$Year==2020] -
                                CMPSCA$ll[1])

CMPSP$u[1] <- CMPSP$u[1] + (catch$Catch[catch$Year==2020] -
                            CMPSP$u[1])
CMPSP$uu[1] <- CMPSP$uu[1] + (catch$Catch[catch$Year==2020] -
                              CMPSP$uu[1])
CMPSP$m[1] <- CMPSP$m[1] + (catch$Catch[catch$Year==2020] -
                            CMPSP$m[1])
CMPSP$l[1] <- CMPSP$l[1] + (catch$Catch[catch$Year==2020] -
                            CMPSP$l[1])
CMPSP$ll[1] <- CMPSP$ll[1] + (catch$Catch[catch$Year==2020] -
                              CMPSP$ll[1])

CMPSPSS$u[1] <- CMPSPSS$u[1] + (catch$Catch[catch$Year==2020] -
                            CMPSPSS$u[1])
CMPSPSS$uu[1] <- CMPSPSS$uu[1] + (catch$Catch[catch$Year==2020] -
                              CMPSPSS$uu[1])
CMPSPSS$m[1] <- CMPSPSS$m[1] + (catch$Catch[catch$Year==2020] -
                            CMPSPSS$m[1])
CMPSPSS$l[1] <- CMPSPSS$l[1] + (catch$Catch[catch$Year==2020] -
                            CMPSPSS$l[1])
CMPSPSS$ll[1] <- CMPSPSS$ll[1] + (catch$Catch[catch$Year==2020] -
                              CMPSPSS$ll[1])

CMPSP_01$u[1] <- CMPSP_01$u[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_01$u[1])
CMPSP_01$uu[1] <- CMPSP_01$uu[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_01$uu[1])
CMPSP_01$m[1] <- CMPSP_01$m[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_01$m[1])
CMPSP_01$l[1] <- CMPSP_01$l[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_01$l[1])
CMPSP_01$ll[1] <- CMPSP_01$ll[1] + (catch$Catch[catch$Year==2020] -
                                    CMPSP_01$ll[1])


CMPSP_02$u[1] <- CMPSP_02$u[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_02$u[1])
CMPSP_02$uu[1] <- CMPSP_02$uu[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_02$uu[1])
CMPSP_02$m[1] <- CMPSP_02$m[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_02$m[1])
CMPSP_02$l[1] <- CMPSP_02$l[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_02$l[1])
CMPSP_02$ll[1] <- CMPSP_02$ll[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_02$ll[1])

CMPSP_03$u[1] <- CMPSP_03$u[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_03$u[1])
CMPSP_03$uu[1] <- CMPSP_03$uu[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_03$uu[1])
CMPSP_03$m[1] <- CMPSP_03$m[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_03$m[1])
CMPSP_03$l[1] <- CMPSP_03$l[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_03$l[1])
CMPSP_03$ll[1] <- CMPSP_03$ll[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_03$ll[1])

CMPSP_04$u[1] <- CMPSP_04$u[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_04$u[1])
CMPSP_04$uu[1] <- CMPSP_04$uu[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_04$uu[1])
CMPSP_04$m[1] <- CMPSP_04$m[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_04$m[1])
CMPSP_04$l[1] <- CMPSP_04$l[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_04$l[1])
CMPSP_04$ll[1] <- CMPSP_04$ll[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_04$ll[1])

CMPSP_05$u[1] <- CMPSP_05$u[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_05$u[1])
CMPSP_05$uu[1] <- CMPSP_05$uu[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_05$uu[1])
CMPSP_05$m[1] <- CMPSP_05$m[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_05$m[1])
CMPSP_05$l[1] <- CMPSP_05$l[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_05$l[1])
CMPSP_05$ll[1] <- CMPSP_05$ll[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_05$ll[1])

CMPSP_06$u[1] <- CMPSP_06$u[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_06$u[1])
CMPSP_06$uu[1] <- CMPSP_06$uu[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_06$uu[1])
CMPSP_06$m[1] <- CMPSP_06$m[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_06$m[1])
CMPSP_06$l[1] <- CMPSP_06$l[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_06$l[1])
CMPSP_06$ll[1] <- CMPSP_06$ll[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_06$ll[1])

#####@> Merging...
projCatch02 <- gtools::smartbind(data.frame(CMPIslope),
                                 data.frame(CMPIratio),
                                 data.frame(CMPGB_slope),
                                 data.frame(CMPCC_20kt),
                                 data.frame(CMPCC_30kt),
                                 data.frame(CMPCC_40kt),
                                 data.frame(CMPSCA),
                                 data.frame(CMPSP),
                                 data.frame(CMPSPSS),
                                 data.frame(CMPSP_01),
                                 data.frame(CMPSP_02),
                                 data.frame(CMPSP_03),
                                 data.frame(CMPSP_04),
                                 data.frame(CMPSP_05),
                                 data.frame(CMPSP_06))

p02 <- ggplot() +
    geom_line(data = catch, aes(x = Year, y = Catch),
              linewidth = 1.2, colour = "gray50") +
    geom_ribbon(data = projCatch02,
                aes(x = real_year, ymin = l, ymax = u,
                    fill = mp_name),
                alpha = 0.1) +
    geom_line(data = projCatch02,
              aes(x = real_year, y = m, colour = mp_name),
              linewidth = 1.2) +
    geom_text_repel(data = filter(projCatch02, real_year == 2053),
                    aes(x = real_year, y = m, colour = mp_name,
                        label = mp_name), nudge_x = 6, nudge_y = 500,
                    force = 5) +
    ## geom_hline(yintercept = 1, colour = "green") +
    ## geom_hline(yintercept = 0.4, colour = "red") +
    geom_vline(xintercept = 2020, colour = "lightgray",
               linetype = "dashed") +
    scale_y_continuous(breaks = seq(0, 80000, 10000), limits = c(0, 80000),
                       expand = c(0, 0)) +
    scale_x_continuous(breaks = seq(1952, 2060, 4)) +
    ## scale_fill_viridis(discrete = TRUE, option = 10) +
    ## scale_colour_viridis(discrete = TRUE, option = 10) +
    labs(x = "Year",
         y = "Catches (t)",
         colour = "") +
    theme_pbs(base_size = 18) +
    theme(legend.position = "none")
p02


#####@> Exporting figures...
ggsave("Results/10_Error_Imp/TS_Catches_Rescaled_10ImpError_ver00.png",
       plot = p02, device = "png", units = "cm", w = 45, h = 30,
       bg = "white", dpi = 600)

ggsave("Results/10_Error_Imp/TS_BBMSY_10ImpError_ver00.png",
       plot = p00, device = "png", units = "cm", w = 45, h = 30,
       bg = "white", dpi = 600)

######@>----------------------------------------------------------------
######@> Projections OMs 19 - 27...

#####@> Extraction the projections...
projDB <- data.frame(OM = NULL,
                     iter = NULL,
                     mp = NULL,
                     year = NULL,
                     value = NULL,
                     mp_name = NULL,
                     Type = NULL,
                     id_year = NULL,
                     real_year = NULL,
                     ref = NULL)
for(i in 1:9) {
    print(paste0("Compiling data from OM ", i))
    out <- get_ts(MSEs03[[i]], this_year = 2020)
    projDB <- rbind(projDB, data.frame(OM = i + 18,
                                       iter = out$iter,
                                       mp = out$mp,
                                       year = out$year,
                                       value = out$value,
                                       mp_name = out$mp_name,
                                       Type = out$Type,
                                       id_year = out$real_year,
                                       real_year = out$real_year,
                                       ref = out$ref))
}

#####@> Estimating the summary quantities...
summProj <- get_ts_quantiles(projDB)

#####@> Segregating base by Type and Year...
histBBMSY <- filter(summProj, Type == "B_BMSY",
                    real_year %in% 1952:2020)
projBBMSY <- filter(summProj, Type == "B_BMSY",
                    real_year %in% 2020:2053,
                    mp_name %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                                   "GB_slope", "Iratio", "Islope1",
                                   "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                                   "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                                   "SP_04", "SP_05", "SP_06"))
histFFMSY <- filter(summProj, Type == "F_FMSY",
                    real_year %in% 1952:2020)
projFFMSY <- filter(summProj, Type == "F_FMSY",
                    real_year %in% 2020:2053,
                    mp_name %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                                   "GB_slope", "Iratio", "Islope1",
                                   "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                                   "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                                   "SP_04", "SP_05", "SP_06"))
histCatch <- filter(summProj, Type == "Catch",
                    real_year %in% 1952:2020)
projCatch <- filter(summProj, Type == "Catch",
                    real_year %in% 2020:2053,
                    mp_name %in% c("CC_20kt", "CC_30kt", "CC_40kt",
                                   "GB_slope", "Iratio", "Islope1",
                                   "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
                                   "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
                                   "SP_04", "SP_05", "SP_06"))

p00 <- ggplot() +
    geom_ribbon(data = histBBMSY,
                aes(x = real_year, ymin = l, ymax = u),
                alpha = 0.9, fill = "lightgray") +
    geom_line(data = histBBMSY, aes(x = real_year, y = ll),
              linetype = "dashed", colour = "lightgray") +
    geom_line(data = histBBMSY, aes(x = real_year, y = uu),
              linetype = "dashed", colour = "lightgray") +
    geom_line(data = histBBMSY, aes(x = real_year, y = m),
              linewidth = 1.2, colour = "gray50") +
    geom_ribbon(data = projBBMSY,
                aes(x = real_year, ymin = l, ymax = u,
                    fill = mp_name),
                alpha = 0.1) +
    geom_line(data = projBBMSY,
              aes(x = real_year, y = m, colour = mp_name),
              linewidth = 1.2) +
    geom_text_repel(data = filter(projBBMSY, real_year == 2053),
                    aes(x = real_year, y = m, colour = mp_name,
                        label = mp_name), nudge_x = 4, force = 5) +
    geom_hline(yintercept = 1, colour = "green") +
    geom_hline(yintercept = 0.4, colour = "red") +
    geom_vline(xintercept = 2020, colour = "lightgray",
               linetype = "dashed") +
    scale_y_continuous(breaks = seq(0, 4, 0.4), limits = c(0, 4),
                       expand = c(0, 0)) +
    scale_x_continuous(breaks = seq(1952, 2060, 4), limits = c(1952, 2060)) +
    labs(x = "Year",
         y = "Spawning Stock Biomass relative to MSY levels",
         colour = "") +
    theme_pbs(base_size = 18) +
    theme(legend.position = "none")
p00

p01 <- ggplot() +
    geom_ribbon(data = histFFMSY,
                aes(x = real_year, ymin = l, ymax = u),
                alpha = 0.9, fill = "lightgray") +
    geom_line(data = histFFMSY, aes(x = real_year, y = ll),
              linetype = "dashed", colour = "lightgray") +
    geom_line(data = histFFMSY, aes(x = real_year, y = uu),
              linetype = "dashed", colour = "lightgray") +
    geom_line(data = histFFMSY, aes(x = real_year, y = m),
              linewidth = 1.2, colour = "gray50") +
    geom_ribbon(data = projFFMSY,
                aes(x = real_year, ymin = l, ymax = u,
                    fill = mp_name),
                alpha = 0.1) +
    geom_line(data = projFFMSY,
              aes(x = real_year, y = m, colour = mp_name),
              linewidth = 1.2) +
    geom_text_repel(data = filter(projFFMSY, real_year == 2050),
                    aes(x = real_year, y = m, colour = mp_name,
                        label = mp_name), nudge_x = 4) +
    geom_hline(yintercept = 1, colour = "red") +
    ## geom_hline(yintercept = 0.4, colour = "red") +
    geom_vline(xintercept = 2020, colour = "lightgray",
               linetype = "dashed") +
    scale_y_continuous(breaks = seq(0, 7, 1), limits = c(0, 8),
                       expand = c(0, 0)) +
    scale_x_continuous(breaks = seq(1952, 2060, 4)) +
    labs(x = "Year",
         y = "F/FMSY",
         colour = "") +
    theme_pbs(base_size = 18) +
    theme(legend.position = "none")
p01

p02 <- ggplot() +
    geom_ribbon(data = histCatch,
                aes(x = real_year, ymin = l, ymax = u),
                alpha = 0.9, fill = "lightgray") +
    geom_line(data = histCatch, aes(x = real_year, y = ll),
              linetype = "dashed", colour = "lightgray") +
    geom_line(data = histCatch, aes(x = real_year, y = uu),
              linetype = "dashed", colour = "lightgray") +
    geom_line(data = histCatch, aes(x = real_year, y = m),
              linewidth = 1.2, colour = "gray50") +
    geom_ribbon(data = projCatch,
                aes(x = real_year, ymin = l, ymax = u,
                    fill = mp_name),
                alpha = 0.1) +
    geom_line(data = projCatch,
              aes(x = real_year, y = m, colour = mp_name),
              linewidth = 1.2) +
    geom_text_repel(data = filter(projCatch, real_year == 2050),
                    aes(x = real_year, y = m, colour = mp_name,
                        label = mp_name), nudge_x = 5, nudge_y = 2) +
    ## geom_hline(yintercept = 1, colour = "green") +
    ## geom_hline(yintercept = 0.4, colour = "red") +
    geom_vline(xintercept = 2020, colour = "lightgray",
               linetype = "dashed") +
    scale_y_continuous(breaks = seq(0, 80000, 10000), limits = c(0, 80000),
                       expand = c(0, 0)) +
    scale_x_continuous(breaks = seq(1952, 2060, 4)) +
    labs(x = "Year",
         y = "Catches (t)",
         colour = "") +
    theme_pbs(base_size = 18) +
    theme(legend.position = "none")
p02

######@> Catch...
catch <- readxl::read_excel("../WSKJ_Input_DB_ver00.xlsx",
                            sheet = "CATCH")

## "CC_20kt", "CC_30kt", "CC_40kt",
## "GB_slope", "Iratio", "Islope1",
## "SCA_100_40_SBMSY", "SP_100_40_SBMSY",
## "SPSS_100_40_SBMSY", "SP_01", "SP_02", "SP_03",
## "SP_04", "SP_05", "SP_06"

CMPIslope <- filter(projCatch, mp_name == "Islope1")
CMPIratio <- filter(projCatch, mp_name == "Iratio")
CMPGB_slope <- filter(projCatch, mp_name == "GB_slope")
CMPCC_20kt <- filter(projCatch, mp_name == "CC_20kt")
CMPCC_30kt <- filter(projCatch, mp_name == "CC_30kt")
CMPCC_40kt <- filter(projCatch, mp_name == "CC_40kt")
CMPSCA <- filter(projCatch, mp_name == "SCA_100_40_SBMSY")
CMPSP <- filter(projCatch, mp_name == "SP_100_40_SBMSY")
CMPSPSS <- filter(projCatch, mp_name == "SPSS_100_40_SBMSY")
CMPSP_01 <- filter(projCatch, mp_name == "SP_01")
CMPSP_02 <- filter(projCatch, mp_name == "SP_02")
CMPSP_03 <- filter(projCatch, mp_name == "SP_03")
CMPSP_04 <- filter(projCatch, mp_name == "SP_04")
CMPSP_05 <- filter(projCatch, mp_name == "SP_05")
CMPSP_06 <- filter(projCatch, mp_name == "SP_06")

#####@> Correction...
CMPIslope$u <- CMPIslope$u + (catch$Catch[catch$Year==2020] -
                              CMPIslope$u[1])
CMPIslope$uu <- CMPIslope$uu + (catch$Catch[catch$Year==2020] -
                                CMPIslope$uu[1])
CMPIslope$m <- CMPIslope$m + (catch$Catch[catch$Year==2020] -
                              CMPIslope$m[1])
CMPIslope$l <- CMPIslope$l + (catch$Catch[catch$Year==2020] -
                              CMPIslope$l[1])
CMPIslope$ll <- CMPIslope$ll + (catch$Catch[catch$Year==2020] -
                                CMPIslope$ll[1])

CMPIratio$u <- CMPIratio$u + (catch$Catch[catch$Year==2020] -
                              CMPIratio$u[1])
CMPIratio$uu <- CMPIratio$uu + (catch$Catch[catch$Year==2020] -
                                CMPIratio$uu[1])
CMPIratio$m <- CMPIratio$m + (catch$Catch[catch$Year==2020] -
                              CMPIratio$m[1])
CMPIratio$l <- CMPIratio$l + (catch$Catch[catch$Year==2020] -
                              CMPIratio$l[1])
CMPIratio$ll <- CMPIratio$ll + (catch$Catch[catch$Year==2020] -
                                CMPIratio$ll[1])

CMPGB_slope$u <- CMPGB_slope$u + (catch$Catch[catch$Year==2020] -
                              CMPGB_slope$u[1])
CMPGB_slope$uu <- CMPGB_slope$uu + (catch$Catch[catch$Year==2020] -
                                CMPGB_slope$uu[1])
CMPGB_slope$m <- CMPGB_slope$m + (catch$Catch[catch$Year==2020] -
                              CMPGB_slope$m[1])
CMPGB_slope$l <- CMPGB_slope$l + (catch$Catch[catch$Year==2020] -
                              CMPGB_slope$l[1])
CMPGB_slope$ll <- CMPGB_slope$ll + (catch$Catch[catch$Year==2020] -
                                CMPGB_slope$ll[1])

CMPCC_20kt$u[1] <- CMPCC_20kt$u[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_20kt$u[1])
CMPCC_20kt$uu[1] <- CMPCC_20kt$uu[1] + (catch$Catch[catch$Year==2020] -
                                        CMPCC_20kt$uu[1])
CMPCC_20kt$m[1] <- CMPCC_20kt$m[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_20kt$m[1])
CMPCC_20kt$l[1] <- CMPCC_20kt$l[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_20kt$l[1])
CMPCC_20kt$ll[1] <- CMPCC_20kt$ll[1] + (catch$Catch[catch$Year==2020] -
                                        CMPCC_20kt$ll[1])

CMPCC_30kt$u[1] <- CMPCC_30kt$u[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_30kt$u[1])
CMPCC_30kt$uu[1] <- CMPCC_30kt$uu[1] + (catch$Catch[catch$Year==2020] -
                                        CMPCC_30kt$uu[1])
CMPCC_30kt$m[1] <- CMPCC_30kt$m[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_30kt$m[1])
CMPCC_30kt$l[1] <- CMPCC_30kt$l[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_30kt$l[1])
CMPCC_30kt$ll[1] <- CMPCC_30kt$ll[1] + (catch$Catch[catch$Year==2020] -
                                        CMPCC_30kt$ll[1])

CMPCC_40kt$u[1] <- CMPCC_40kt$u[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_40kt$u[1])
CMPCC_40kt$uu[1] <- CMPCC_40kt$uu[1] + (catch$Catch[catch$Year==2020] -
                                        CMPCC_40kt$uu[1])
CMPCC_40kt$m[1] <- CMPCC_40kt$m[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_40kt$m[1])
CMPCC_40kt$l[1] <- CMPCC_40kt$l[1] + (catch$Catch[catch$Year==2020] -
                                      CMPCC_40kt$l[1])
CMPCC_40kt$ll[1] <- CMPCC_40kt$ll[1] + (catch$Catch[catch$Year==2020] -
                                        CMPCC_40kt$ll[1])

CMPSCA$u[1] <- CMPSCA$u[1] + (catch$Catch[catch$Year==2020] -
                              CMPSCA$u[1])
CMPSCA$uu[1] <- CMPSCA$uu[1] + (catch$Catch[catch$Year==2020] -
                                CMPSCA$uu[1])
CMPSCA$m[1] <- CMPSCA$m[1] + (catch$Catch[catch$Year==2020] -
                              CMPSCA$m[1])
CMPSCA$l[1] <- CMPSCA$l[1] + (catch$Catch[catch$Year==2020] -
                              CMPSCA$l[1])
CMPSCA$ll[1] <- CMPSCA$ll[1] + (catch$Catch[catch$Year==2020] -
                                CMPSCA$ll[1])

CMPSP$u[1] <- CMPSP$u[1] + (catch$Catch[catch$Year==2020] -
                            CMPSP$u[1])
CMPSP$uu[1] <- CMPSP$uu[1] + (catch$Catch[catch$Year==2020] -
                              CMPSP$uu[1])
CMPSP$m[1] <- CMPSP$m[1] + (catch$Catch[catch$Year==2020] -
                            CMPSP$m[1])
CMPSP$l[1] <- CMPSP$l[1] + (catch$Catch[catch$Year==2020] -
                            CMPSP$l[1])
CMPSP$ll[1] <- CMPSP$ll[1] + (catch$Catch[catch$Year==2020] -
                              CMPSP$ll[1])

CMPSPSS$u[1] <- CMPSPSS$u[1] + (catch$Catch[catch$Year==2020] -
                            CMPSPSS$u[1])
CMPSPSS$uu[1] <- CMPSPSS$uu[1] + (catch$Catch[catch$Year==2020] -
                              CMPSPSS$uu[1])
CMPSPSS$m[1] <- CMPSPSS$m[1] + (catch$Catch[catch$Year==2020] -
                            CMPSPSS$m[1])
CMPSPSS$l[1] <- CMPSPSS$l[1] + (catch$Catch[catch$Year==2020] -
                            CMPSPSS$l[1])
CMPSPSS$ll[1] <- CMPSPSS$ll[1] + (catch$Catch[catch$Year==2020] -
                              CMPSPSS$ll[1])

CMPSP_01$u[1] <- CMPSP_01$u[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_01$u[1])
CMPSP_01$uu[1] <- CMPSP_01$uu[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_01$uu[1])
CMPSP_01$m[1] <- CMPSP_01$m[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_01$m[1])
CMPSP_01$l[1] <- CMPSP_01$l[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_01$l[1])
CMPSP_01$ll[1] <- CMPSP_01$ll[1] + (catch$Catch[catch$Year==2020] -
                                    CMPSP_01$ll[1])


CMPSP_02$u[1] <- CMPSP_02$u[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_02$u[1])
CMPSP_02$uu[1] <- CMPSP_02$uu[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_02$uu[1])
CMPSP_02$m[1] <- CMPSP_02$m[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_02$m[1])
CMPSP_02$l[1] <- CMPSP_02$l[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_02$l[1])
CMPSP_02$ll[1] <- CMPSP_02$ll[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_02$ll[1])

CMPSP_03$u[1] <- CMPSP_03$u[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_03$u[1])
CMPSP_03$uu[1] <- CMPSP_03$uu[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_03$uu[1])
CMPSP_03$m[1] <- CMPSP_03$m[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_03$m[1])
CMPSP_03$l[1] <- CMPSP_03$l[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_03$l[1])
CMPSP_03$ll[1] <- CMPSP_03$ll[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_03$ll[1])

CMPSP_04$u[1] <- CMPSP_04$u[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_04$u[1])
CMPSP_04$uu[1] <- CMPSP_04$uu[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_04$uu[1])
CMPSP_04$m[1] <- CMPSP_04$m[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_04$m[1])
CMPSP_04$l[1] <- CMPSP_04$l[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_04$l[1])
CMPSP_04$ll[1] <- CMPSP_04$ll[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_04$ll[1])

CMPSP_05$u[1] <- CMPSP_05$u[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_05$u[1])
CMPSP_05$uu[1] <- CMPSP_05$uu[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_05$uu[1])
CMPSP_05$m[1] <- CMPSP_05$m[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_05$m[1])
CMPSP_05$l[1] <- CMPSP_05$l[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_05$l[1])
CMPSP_05$ll[1] <- CMPSP_05$ll[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_05$ll[1])

CMPSP_06$u[1] <- CMPSP_06$u[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_06$u[1])
CMPSP_06$uu[1] <- CMPSP_06$uu[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_06$uu[1])
CMPSP_06$m[1] <- CMPSP_06$m[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_06$m[1])
CMPSP_06$l[1] <- CMPSP_06$l[1] + (catch$Catch[catch$Year==2020] -
                                CMPSP_06$l[1])
CMPSP_06$ll[1] <- CMPSP_06$ll[1] + (catch$Catch[catch$Year==2020] -
                                  CMPSP_06$ll[1])

#####@> Merging...
projCatch02 <- gtools::smartbind(data.frame(CMPIslope),
                                 data.frame(CMPIratio),
                                 data.frame(CMPGB_slope),
                                 data.frame(CMPCC_20kt),
                                 data.frame(CMPCC_30kt),
                                 data.frame(CMPCC_40kt),
                                 data.frame(CMPSCA),
                                 data.frame(CMPSP),
                                 data.frame(CMPSPSS),
                                 data.frame(CMPSP_01),
                                 data.frame(CMPSP_02),
                                 data.frame(CMPSP_03),
                                 data.frame(CMPSP_04),
                                 data.frame(CMPSP_05),
                                 data.frame(CMPSP_06))

p02 <- ggplot() +
    geom_line(data = catch, aes(x = Year, y = Catch),
              linewidth = 1.2, colour = "gray50") +
    geom_ribbon(data = projCatch02,
                aes(x = real_year, ymin = l, ymax = u,
                    fill = mp_name),
                alpha = 0.1) +
    geom_line(data = projCatch02,
              aes(x = real_year, y = m, colour = mp_name),
              linewidth = 1.2) +
    geom_text_repel(data = filter(projCatch02, real_year == 2053),
                    aes(x = real_year, y = m, colour = mp_name,
                        label = mp_name), nudge_x = 6, nudge_y = 500,
                    force = 5) +
    ## geom_hline(yintercept = 1, colour = "green") +
    ## geom_hline(yintercept = 0.4, colour = "red") +
    geom_vline(xintercept = 2020, colour = "lightgray",
               linetype = "dashed") +
    scale_y_continuous(breaks = seq(0, 80000, 10000), limits = c(0, 80000),
                       expand = c(0, 0)) +
    scale_x_continuous(breaks = seq(1952, 2060, 4)) +
    ## scale_fill_viridis(discrete = TRUE, option = 10) +
    ## scale_colour_viridis(discrete = TRUE, option = 10) +
    labs(x = "Year",
         y = "Catches (t)",
         colour = "") +
    theme_pbs(base_size = 18) +
    theme(legend.position = "none")
p02

#####@> Exporting figures...
ggsave("Results/20_Error_Imp/TS_Catches_Rescaled_20ImpError_ver00.png",
       plot = p02, device = "png", units = "cm", w = 45, h = 30,
       bg = "white", dpi = 600)

ggsave("Results/20_Error_Imp/TS_BBMSY_20ImpError_ver00.png",
       plot = p00, device = "png", units = "cm", w = 45, h = 30,
       bg = "white", dpi = 600)

######@>----------------------------------------------------------------
######@> Preparing comparisons between Scenarios (Perfect, 10%, 20%)...

######@> Including references for each scenario...
tab01_exp_ref$scenario <- "Perfect TAC Implementation"
tab01_exp_10$scenario <- "10% TAC Error Implementation"
tab01_exp_20$scenario <- "20% TAC Error Implementation"

######@> Combining datasets...
output <- gtools::smartbind(tab01_exp_ref, tab01_exp_10, tab01_exp_20)

######@> Boxplot...
p00 <- ggplot(data = filter(output, PM %in% c("AvC_long", "AvC_med",
                                              "AvC_short")),
              aes(x = PM, y = Values, fill = scenario)) +
    geom_boxplot(outlier.shape = 21, outlier.fill = "gray10",
                 outlier.colour = "gray10", outlier.alpha = 0.05,
                 outlier.size = 0.4, varwidth = FALSE) +
    facet_wrap(~MP, scales = "free", ncol = 4) +
    labs(x = "Performance Metrics", y = "Mediam Catch (t)", fill = "") +
    my_theme() +
    theme(legend.position = c(0.89, 0.08))
p00

p01 <- ggplot(data = filter(output, PM %in% c("LRP_long", "LRP_med",
                                              "LRP_short", "LRP")),
              aes(x = PM, y = Values, fill = scenario)) +
    geom_boxplot(outlier.shape = 21, outlier.fill = "gray10",
                 outlier.colour = "gray10", outlier.alpha = 0.05,
                 outlier.size = 0.4, varwidth = FALSE) +
    facet_wrap(~MP, scales = "free", ncol = 4) +
    labs(x = "Performance Metrics", y = "Prob. SB < 0.4SBMSY", fill =
                                                                   "") +
    scale_y_continuous(limits = c(0, 1)) +
    my_theme() +
    theme(legend.position = c(0.89, 0.08))
p01

p02 <- ggplot(data = filter(output, PM %in% c("PGK_long", "PGK_med",
                                              "PGK_short", "PGK")),
              aes(x = PM, y = Values, fill = scenario)) +
    geom_boxplot(outlier.shape = 21, outlier.fill = "gray10",
                 outlier.colour = "gray10", outlier.alpha = 0.05,
                 outlier.size = 0.4, varwidth = FALSE) +
    facet_wrap(~MP, scales = "free", ncol = 4) +
    labs(x = "Performance Metrics", y = "Prob. Kobe Green Quadrant", fill = "") +
    my_theme() +
    scale_y_continuous(limits = c(0, 1)) +
    theme(legend.position = c(0.89, 0.08))
p02

p03 <- ggplot(data = filter(output, PM %in% c("VarClong", "VarCmedium",
                                              "VarC")),
              aes(x = PM, y = Values, fill = scenario)) +
    geom_boxplot(outlier.shape = 21, outlier.fill = "gray10",
                 outlier.colour = "gray10", outlier.alpha = 0.05,
                 outlier.size = 0.4, varwidth = FALSE) +
    facet_wrap(~MP, scales = "free", ncol = 4) +
    labs(x = "Performance Metrics", y = "Mean Variation in TAC (%)",
         fill = "") +
    scale_y_continuous(limits = c(0, 1)) +
    my_theme() +
    theme(legend.position = c(0.89, 0.08))
p03

ggsave("Boxplot_AvC_ver00.svg", plot = p00, device = "svg",
       units = "cm", width = 48, height = 33, dpi = 300,
       bg = "white")

ggsave("Boxplot_LRP_ver00.svg", plot = p01, device = "svg",
       units = "cm", width = 48, height = 33, dpi = 300,
       bg = "white")

ggsave("Boxplot_PGK_ver00.svg", plot = p02, device = "svg",
       units = "cm", width = 48, height = 33, dpi = 300,
       bg = "white")

ggsave("Boxplot_VarC_ver00.svg", plot = p03, device = "svg",
       units = "cm", width = 48, height = 33, dpi = 300,
       bg = "white")

######@> Ridges...
p00 <- ggplot(data = filter(output, PM %in% c("AvC_long", "AvC_med",
                                              "AvC_short")),
              aes(x = Values, y = PM, fill = scenario, colour = scenario)) +
    stat_density_ridges(quantile_lines = TRUE, alpha = 0.5) +
    facet_wrap(~MP, scales = "free_x", ncol = 5) +
    labs(y = "Performance Metrics", x = "Mediam Catch (t)", fill = "",
         colour = "") +
    scale_fill_hue(direction = -1) +
    scale_colour_hue(direction = -1) +
    my_theme() +
    theme(legend.position = "bottom")
p00

p01 <- ggplot(data = filter(output, PM %in% c("LRP_long", "LRP_med",
                                              "LRP_short", "LRP")),
              aes(x = Values, y = PM, fill = scenario, colour = scenario)) +
    stat_density_ridges(quantile_lines = TRUE, alpha = 0.5) +
    facet_wrap(~MP, scales = "free_x", ncol = 5) +
    labs(y = "Performance Metrics", x = "Prob. SB < 0.4SBMSY",
         fill = "", colour = "") +
    scale_x_continuous(limits = c(0, 1)) +
    scale_fill_hue(direction = -1) +
    scale_colour_hue(direction = -1) +
    my_theme() +
    theme(legend.position = "bottom")
p01

p02 <- ggplot(data = filter(output, PM %in% c("PGK_long", "PGK_med",
                                              "PGK_short", "PGK")),
              aes(x = Values, y = PM, fill = scenario, colour = scenario)) +
    stat_density_ridges(quantile_lines = TRUE, alpha = 0.5) +
    facet_wrap(~MP, scales = "free_x", ncol = 5) +
    labs(y = "Performance Metrics", x = "Prob. Kobe Green Quadrant",
         fill = "", colour = "") +
    scale_fill_hue(direction = -1) +
    scale_colour_hue(direction = -1) +
    my_theme() +
    scale_x_continuous(limits = c(0, 1)) +
    theme(legend.position = "none")
p02

p03 <- ggplot(data = filter(output, PM %in% c("VarClong", "VarCmedium",
                                              "VarC")),
              aes(x = Values, y = PM, fill = scenario, colour = scenario)) +
    stat_density_ridges(quantile_lines = TRUE, alpha = 0.5) +
    facet_wrap(~MP, scales = "free_x", ncol = 5) +
    labs(y = "Performance Metrics", x = "Mean Variation in TAC (%)",
         fill = "", colour = "") +
    scale_x_continuous(limits = c(0, 1)) +
    scale_fill_hue(direction = -1) +
    scale_colour_hue(direction = -1) +
    my_theme() +
    theme(legend.position = "none")
p03

ggsave("Ridges_AvC_ver00.svg", plot = p00, device = "svg",
       units = "cm", width = 48, height = 33, dpi = 300,
       bg = "white")

ggsave("Ridges_LRP_ver00.svg", plot = p01, device = "svg",
       units = "cm", width = 48, height = 33, dpi = 300,
       bg = "white")

ggsave("Ridges_PGK_ver00.svg", plot = p02, device = "svg",
       units = "cm", width = 48, height = 33, dpi = 300,
       bg = "white")

ggsave("Ridges_VarC_ver00.svg", plot = p03, device = "svg",
       units = "cm", width = 48, height = 33, dpi = 300,
       bg = "white")


######@> Spider plot...

########################################################################
######@> HCR ramp plot...

######@> HCR Ramp 100-40...
Brel <- seq(0, 2, length.out = 200)
Frel <- HCRlin(Brel, 0.4, 1, 0.1, 1)

df <- data.frame(Brel, Frel)

#####@> Colors to figure...
kobe_df <- bind_rows(data.frame(x = c(0, 0, 1, 1),
                                y = c(0, 1, 1, 0),
                                fill = "bl"),
                     data.frame(x = c(1, 1, 2, 2),
                                y = c(0, 1, 1, 0),
                                fill = "br"),
                     data.frame(x = c(0, 0, 1, 1),
                                y = c(1, 2, 2, 1),
                                fill = "tl"),
                     data.frame(x = c(1, 1, 2, 2),
                                y = c(1, 2, 2, 1),
                                fill = "tr"))
kobe_df$alpha <- 0.3

#####@> Figure dots...
p00 <- ggplot() +
    geom_polygon(data = kobe_df,
                 aes(x = x, y = y, fill = fill, alpha = alpha)) +
    scale_fill_manual(values = c("#F8DC7A", "#67C18B", "#D8775D",
                                 "#FDBD56")) +
    geom_line(data = df, aes(x = Brel, y = Frel), linewidth = 1) +
    expand_limits(x = c(0, 2), y = c(0, 2)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 2)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 2)) +
    geom_hline(yintercept = 1, color = "darkgray", linetype = 2) +
    geom_vline(xintercept = 1, color = "darkgray", linetype = 2) +
    labs(x = "", y = "") +
    ## labs(x = expression(SB/SB[MSY]), y = expression(F/F[MSY])) +
    my_theme() +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank())
p00

ggsave("HCR_Ramp.svg", plot = p00,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)

######@> HCR Ramp 80-40...
Brel <- seq(0, 2, length.out = 200)
Frel <- HCRlin(Brel, 0.4, 1, 0.1, 0.8)

df <- data.frame(Brel, Frel)

#####@> Colors to figure...
kobe_df <- bind_rows(data.frame(x = c(0, 0, 1, 1),
                                y = c(0, 1, 1, 0),
                                fill = "bl"),
                     data.frame(x = c(1, 1, 2, 2),
                                y = c(0, 1, 1, 0),
                                fill = "br"),
                     data.frame(x = c(0, 0, 1, 1),
                                y = c(1, 2, 2, 1),
                                fill = "tl"),
                     data.frame(x = c(1, 1, 2, 2),
                                y = c(1, 2, 2, 1),
                                fill = "tr"))
kobe_df$alpha <- 0.3

#####@> Figure dots...
p00 <- ggplot() +
    geom_polygon(data = kobe_df,
                 aes(x = x, y = y, fill = fill, alpha = alpha)) +
    scale_fill_manual(values = c("#F8DC7A", "#67C18B", "#D8775D",
                                 "#FDBD56")) +
    geom_line(data = df, aes(x = Brel, y = Frel), linewidth = 1) +
    expand_limits(x = c(0, 2), y = c(0, 2)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 2)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 2)) +
    geom_hline(yintercept = 1, color = "darkgray", linetype = 2) +
    geom_vline(xintercept = 1, color = "darkgray", linetype = 2) +
    labs(x = "", y = "") +
    ## labs(x = expression(SB/SB[MSY]), y = expression(F/F[MSY])) +
    my_theme() +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank())
p00

ggsave("HCR_Ramp_80-10.svg", plot = p00,
       device = "svg", dpi = 300, bg = "white", unit = "cm",
       w = 35, h = 35)



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
