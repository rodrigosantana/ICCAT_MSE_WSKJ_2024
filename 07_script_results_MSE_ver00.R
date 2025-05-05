########################################################################
## Description: Compiling Final Results W-SKJ MSE based on the new
## OMs...
##
## Maintainer: Datenkraft - SCRS/ICCAT (TT Species Group)
## Author: Rodrigo Sant'Ana
## Created: seg set 2 21:59:06 2024 (-0300)
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
library(unicol)
library(unikn)
library(patchwork)
library(data.table)
library(ggrepel)
library(tidyr)
library(ggridges)
library(gt)

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

#####@> Colorpal set to this work...
my_pal <- unikn::usecol(caltech_1, n = 16)

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
    ts_data$iter <- rep(seq_len(iters),
                        length(unique(ts_data$Type))) # parallel messed this up
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
        left_join(bind_cols(all_mp_yrs,
                            tibble(Type = rep("SSB",
                                              nrow(all_mp_yrs)))),
                  by = "Type")
    ref_f_hist <- data.frame(
        ref = object@RefPoint$FMSY[, 1, object@nyears],
        iter = seq_len(iters), Type = "FM",
        stringsAsFactors = FALSE) %>%
        left_join(bind_cols(all_mp_yrs,
                            tibble(Type = rep("FM", nrow(all_mp_yrs)))),
                  by = "Type")
    all_mp_yrs <- expand.grid(mp = seq_along(mps$mp),
                              real_year = sort(unique(ts_data$real_year)))
    ref_msy <- data.frame(ref = 1, iter = seq_len(iters),
                          Type = "Catch", stringsAsFactors = FALSE) %>%
        left_join(bind_cols(all_mp_yrs,
                            tibble(Type = rep("Catch",
                                              nrow(all_mp_yrs)))),
                  by = "Type")
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

max_value_by_3years <- function(data, year_column, value_column) {
    max_by_period <- data %>%
        mutate(period = floor((!!sym(year_column) -
                               min(!!sym(year_column))) / 3)) %>%
        group_by(period) %>%
        summarise(max_value = max(!!sym(value_column), na.rm = TRUE)) %>%
        ungroup()
    max_values_repeated <- rep(max_by_period$max_value, each = 3)
    full_years <- seq(min(data[[year_column]]), max(data[[year_column]]))
    result <- data.frame(year = full_years,
                         max_value = max_values_repeated[1:length(full_years)])
    return(result)
}

########################################################################
######@> Setup MSE...

######@> Folder destination...
if(dir.exists("05_Results")) {
    print("Folder exists!!")
    folder <- "05_Results/"
} else {
    dir.create("05_Results")
    print("Folder will be created!!")
    folder <- "05_Results/"
}

######@> Path to MSEs...
path00 <- "04_MSEs/"

######@> Path to MPs...
path01 <- "../R-Work/"

######@> Path to PMs...
path02 <- "../R-Work/"

######@> Custom R environment...
options(scipen = 16)
options(digits = 5)

########################################################################
######@> Loading everything...

######@> Management Procedures...
source(paste0(path01, "03_script_prepare_MPs_ver00.R"))

#####@> Defining MPs...
MPs <- c("FMSYref", "FMSYref75", "FMSYref110", ## Reference FMSY
         "CE_01", "CE_02", "CE_03", ## Index-based: Exploitation rate
         "IR_01", "IR_02", "IR_03", ## Index-based: I_ratio
         "IS_01", "IS_02", "IS_03", ## Index-based: I_slope
         "SP_01", "SP_02", "SP_03", "SP_04") ## Model-based: SP

######@> Performance Metrics...
source(paste0(path02, "06_script_prepare_PMs_ver00.R"))

#####@> List of PMs...
PMs <- avail("PM")[c(1:23)]

######@> Loading tuning process data...
load("05_Results/Tune_PGK01-06.RData")
load("05_Results/Tune_PGK_FULL.RData")

########################################################################
######@> Loading MSEs...

######@> MSEs FULL...
MSEs <- sapply(dir("04_MSEs",
                   pattern = "[0-9][0-9][0-9]FULL_IVInds_CORRECTED_00",
                   full.names = TRUE), readRDS)
names(MSEs) <- paste0("MSE", sprintf("%03d", 1:9))

######@> MSEs SPSS...
MSEsSPSS <- sapply(dir("04_MSEs",
                       pattern = "_IVInds_CORRECTED_Model_SPSS_01",
                       full.names = TRUE), readRDS)
names(MSEsSPSS) <- paste0("MSE", sprintf("%03d", 1:9))

########################################################################
######@> Consolidating Results...

######@> Preparing tune data...

#####@> Results from PGK 1-6 years...
results01 <- list("IR" = resultsIR, "CE" = resultsCE, "IS" = resultsIS,
                  "SP01" = resultsSP01, "SP03" = resultsSP03)

#####@> Results from PGK 30 years...
results02 <- list("IR" = resultsIR02, "CE" = resultsCE02,
                  "IS" = resultsIS02,
                  "SP01" = resultsSP0102, "SP03" = resultsSP0302)

#####@> Function to extract data from tune profiles...
extract_profile <- function(list_profiles, ...) {
    temp00 <- lapply(names(list_profiles), function(name) {
        df <- do.call("rbind.data.frame",
                      list_profiles[[name]]$optimization_results)
        df$profile_name <- name
        df$order <- seq_len(nrow(df))
        return(df)
    })
    temp00 <- do.call("rbind.data.frame", temp00)
    return(temp00)
}

#####@> Profiles extractions...
profilePGK1.6 <- extract_profile(results01)
profilePGK30 <- extract_profile(results02)

#####@> Tuning profile optimization...
pAA <- ggplot(data = filter(profilePGK1.6, profile_name != "IS"),
              aes(x = par, y = PGKw)) +
    geom_line() +
    geom_point(aes(fill = ifelse(abs(PGKw - 0.7) < 0.01,
                                 "green", "white")),
               pch = 21, size = 4) +
    geom_hline(yintercept = 0.7, lty = 2, alpha = 0.5) +
    facet_wrap(~profile_name, scales = "free_x") +
    labs(x = "Tune parameter",
         y = expression(PGK[1-6*~years])) +
    scale_x_continuous(labels =
                           scales::number_format(accuracy = 0.001)) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_manual(values = c("green", "white")) +
    my_theme() +
    theme(legend.position = "none")
pAA

pAB <- ggplot(data = filter(profilePGK30, profile_name != "IS"),
              aes(x = par, y = PGKw)) +
    geom_line() +
    geom_point(aes(fill = ifelse(abs(PGKw - 0.7) < 0.01,
                                 "green", "white")),
               pch = 21, size = 4) +
    geom_hline(yintercept = 0.7, lty = 2, alpha = 0.5) +
    facet_wrap(~profile_name, scales = "free_x") +
    labs(x = "Tune parameter",
         y = expression(PGK[1-30*~years])) +
    scale_x_continuous(labels =
                           scales::number_format(accuracy = 0.001)) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_manual(values = c("green", "white")) +
    my_theme() +
    theme(legend.position = "none")
pAB

####@> Exporting graphics...
ggsave("05_Results/SCRS_Doc/Fig00_Profile_PGK_1_6_ver00.png",
       plot = pAA, device = "png", units = "cm", w = 25, h = 25,
       dpi = 600, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig01_Profile_PGK_1_30_ver00.png",
       plot = pAB, device = "png", units = "cm", w = 25, h = 25,
       dpi = 600, bg = "white")

######@>----------------------------------------------------------------
######@> No error implementation - OMs 1-9...

######@> Extracting PMs from all MMSE's objects...

#####@> Listing PMs functions - Models default - Necessary exclude SP_02
#####@> and SP_04 MPs and include these MPs fitted separately...

#####@> Subsetting MPs...
MPsSub <- MPs[c(1:9, 13, 15)]
MSEsSub <- lapply(MSEs, function(x) Sub(x, MPs = MPsSub))

#####@> Prepare PMs...
PMlist <- list()
for(i in seq_along(PMs)) {
    PMlist[[i]] <- try(get(PMs[i]), silent = TRUE)
}

######@> Extracting values - old version...
## PMvalues <- list()
## PMvalues2 <- list()
## for(i in seq_along(MSEs)) {
##     for(j in seq_along(PMs)) {
##         pm <- PMlist[[j]](MSEs[[i]])
##         om <- paste0("MSE", sprintf("%03d", i))
##         mps <- pm@MPs
##         val <- data.frame(pm@Prob)
##         names(val) <- mps
##         val <- val %>%
##             mutate(sim = 1:300) %>%
##             pivot_longer(names_to = "MP", values_to = "Values", 1:length(mps))
##         nom <- pm@Name
##         cap <- pm@Caption
##         pm <- PMs[j]
##         PMvalues[[j]] <- val %>%
##             mutate(OM = om, Name = nom, Caption = cap, PM = pm) %>%
##             select(OM, Name, Caption, sim, MP, PM, Values) %>%
##             as.data.frame()
##     }
##     out <- do.call(rbind.data.frame, PMvalues)
##     PMvalues2[[i]] <- out
## }
## PM_output <- do.call(rbind.data.frame, PMvalues2)

######@> Extracting values using optimized process...

#####@> MSEs Sub...
temp00 <- map2(MSEsSub, seq_along(MSEsSub), function(mse, i) {
    om <- paste0("MSE", sprintf("%03d", i))
    map2(PMlist, PMs, function(pm_fun, pm_name) {
        pm <- pm_fun(mse)
        mps <- pm@MPs
        val <- data.frame(pm@Prob)
        names(val) <- mps
        val %>%
            mutate(sim = 1:300) %>%
            pivot_longer(names_to = "MP", values_to = "Values",
                         1:length(mps)) %>%
            mutate(OM = om,
                   Name = pm@Name,
                   Caption = pm@Caption,
                   PM = pm_name) %>%
            select(OM, Name, Caption, sim, MP, PM, Values) %>%
            as.data.frame()
    }) %>%
        bind_rows()
}) %>%
    bind_rows()

#####@> MSEs SPSS...
temp01 <- map2(MSEsSPSS, seq_along(MSEsSPSS), function(mse, i) {
    om <- paste0("MSE", sprintf("%03d", i))
    map2(PMlist, PMs, function(pm_fun, pm_name) {
        pm <- pm_fun(mse)
        mps <- pm@MPs
        val <- data.frame(pm@Prob)
        names(val) <- mps
        val %>%
            mutate(sim = 1:300) %>%
            pivot_longer(names_to = "MP", values_to = "Values",
                         1:length(mps)) %>%
            mutate(OM = om,
                   Name = pm@Name,
                   Caption = pm@Caption,
                   PM = pm_name) %>%
            select(OM, Name, Caption, sim, MP, PM, Values) %>%
            as.data.frame()
    }) %>%
        bind_rows()
}) %>%
    bind_rows()

#####@> Combining results...
PM_output <- bind_rows(temp00, temp01)

######@> Comparing outputs...
## comparacao <- all.equal(PM_output, PM_output2)
## if (isTRUE(comparacao)) {
##     message("Os outputs são iguais!")
## } else {
##     message("Os outputs são diferentes: ", comparacao)
## }

#####@> Looking for NAs...
tab01_NA <- PM_output %>%
    filter(is.na(Values)) %>%
    group_by(OM, MP, PM) %>%
    summarise(N = n()) %>%
    as.data.frame()

######@> Table with probabilities...
tab01 <- PM_output %>%
    group_by(MP, PM) %>%
    summarise(m = mean(Values, na.rm = TRUE),
              d = sd(Values, na.rm = TRUE),
              q2.5 = quantile(Values, probs = 0.025, na.rm = TRUE),
              q50 = quantile(Values, probs = 0.5, na.rm = TRUE),
              q97.5 = quantile(Values, probs = 0.975, na.rm = TRUE)) %>%
    mutate(Class = case_when(
               grepl("FMSY", MP) ~ "Reference",
               grepl("IR", MP) ~ "Index-based: Iratio",
               grepl("CE", MP) ~ "Index-based: Exploitation rate",
               grepl("SP_01|SP_03", MP) ~
                   "Model-based: Surplus Production",
               grepl("SP_02|SP_04", MP) ~
                   "Model-based: State-Space Surplus Production"))

######@> Figure...

#####@> PGK...
temp <- PM_output %>%
    filter(PM %in% c("PGK", "PGK_long", "PGK_med", "PGK_short")) %>%
    mutate(PM = factor(PM, levels = c("PGK_short", "PGK_med",
                                      "PGK_long", "PGK")),
           MP = factor(MP, levels = c("FMSYref", "FMSYref75",
                                      "FMSYref110", "IR_01", "IR_02",
                                      "IR_03", "CE_01", "CE_02",
                                      "CE_03", "SP_01", "SP_02",
                                      "SP_03", "SP_04"))) %>%
    mutate(Class = case_when(
               grepl("FMSY", MP) ~ "Reference",
               grepl("IR", MP) ~ "Index-based: Iratio",
               grepl("CE", MP) ~ "Index-based: Exploitation rate",
               grepl("SP_01|SP_03", MP) ~
                   "Model-based: Surplus Production",
               grepl("SP_02|SP_04", MP) ~
                   "Model-based: State-Space Surplus Production"))

p00 <- ggplot(data = temp,
              aes(x = PM, y = Values, fill = Class,
                  colour = Class)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5,
                scale = "width") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.2,
                 colour = "darkred") +
    facet_wrap(~MP, ncol = 3) +
    geom_hline(yintercept = 0.7, linetype = "dashed", colour = "gray10") +
    scale_fill_manual(values = my_pal[c(3, 5, 8, 13, 15)]) +
    scale_colour_manual(values = my_pal[c(3, 5, 8, 13, 15)]) +
    my_theme() +
    labs(x = "Performance Metric",
         y = expression(PGK[Status])) +
    theme(legend.position = "none")
p00

p00a <- ggplot(data = filter(temp,
                             Class == "Reference"),
              aes(x = PM, y = Values, fill = Class,
                  colour = Class)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5,
                scale = "width") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.2,
                 colour = "darkred") +
    facet_grid(Class~MP) +
    geom_hline(yintercept = 0.7, linetype = "dashed",
               colour = "gray10") +
    scale_fill_manual(values = my_pal[3]) +
    scale_colour_manual(values = my_pal[3]) +
    my_theme() +
    labs(x = "Performance Metric",
         y = expression(PGK[Status])) +
    theme(legend.position = "none")
p00a

p00b <- ggplot(data = filter(temp,
                             Class == "Index-based: Iratio"),
               aes(x = PM, y = Values, fill = Class,
                   colour = Class)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5,
                scale = "width") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.2,
                 colour = "darkred") +
    facet_grid(Class~MP) +
    geom_hline(yintercept = 0.7, linetype = "dashed",
               colour = "gray10") +
    scale_fill_manual(values = my_pal[5]) +
    scale_colour_manual(values = my_pal[5]) +
    my_theme() +
    labs(x = "Performance Metric",
         y = expression(PGK[Status])) +
    theme(legend.position = "none")
p00b

p00c <- ggplot(data = filter(temp,
                             Class == "Index-based: Exploitation rate"),
               aes(x = PM, y = Values, fill = Class,
                   colour = Class)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5,
                scale = "width") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.2,
                 colour = "darkred") +
    facet_grid(Class~MP) +
    geom_hline(yintercept = 0.7, linetype = "dashed",
               colour = "gray10") +
    scale_fill_manual(values = my_pal[8]) +
    scale_colour_manual(values = my_pal[8]) +
    my_theme() +
    labs(x = "Performance Metric",
         y = expression(PGK[Status])) +
    theme(legend.position = "none")
p00c

p00d <- ggplot(data = filter(temp,
                             Class == "Model-based: Surplus Production"),
               aes(x = PM, y = Values, fill = Class,
                   colour = Class)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5,
                scale = "width") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.2,
                 colour = "darkred") +
    facet_grid(Class~MP) +
    geom_hline(yintercept = 0.7, linetype = "dashed",
               colour = "gray10") +
    scale_fill_manual(values = my_pal[13]) +
    scale_colour_manual(values = my_pal[13]) +
    my_theme() +
    labs(x = "Performance Metric",
         y = expression(PGK[Status])) +
    theme(legend.position = "none")
p00d

p00e <-
    ggplot(data =
               filter(temp,
                      Class == "Model-based: State-Space Surplus Production"),
               aes(x = PM, y = Values, fill = Class,
                   colour = Class)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5,
                scale = "width") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.2,
                 colour = "darkred") +
    facet_grid(Class~MP) +
    geom_hline(yintercept = 0.7, linetype = "dashed",
               colour = "gray10") +
    scale_fill_manual(values = my_pal[15]) +
    scale_colour_manual(values = my_pal[15]) +
    my_theme() +
    labs(x = "Performance Metric",
         y = expression(PGK[Status])) +
    theme(legend.position = "none")
p00e

#####@> Exporting figures...
ggsave("05_Results/SCRS_Doc/Fig02_Violin_PGK_RefMPs_Reference_ver00.png",
       plot = p00a, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig02_Violin_PGK_IRMPs_Reference_ver00.png",
       plot = p00b, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig02_Violin_PGK_CEMPs_Reference_ver00.png",
       plot = p00c, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig02_Violin_PGK_SPMPs_Reference_ver00.png",
       plot = p00d, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig02_Violin_PGK_SPSSMPs_Reference_ver00.png",
       plot = p00e, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

#####@> VarC...
temp <- PM_output %>%
    filter(PM %in% c("VarC", "VarClong", "VarCmedium")) %>%
    mutate(PM = factor(PM, levels = c("VarCmedium", "VarClong",
                                      "VarC")),
           MP = factor(MP, levels = c("FMSYref", "FMSYref75",
                                      "FMSYref110", "IR_01", "IR_02",
                                      "IR_03", "CE_01", "CE_02",
                                      "CE_03", "SP_01", "SP_02",
                                      "SP_03", "SP_04"))) %>%
    mutate(Class = case_when(
               grepl("FMSY", MP) ~ "Reference",
               grepl("IR", MP) ~ "Index-based: Iratio",
               grepl("CE", MP) ~ "Index-based: Exploitation rate",
               grepl("SP_01|SP_03", MP) ~
                   "Model-based: Surplus Production",
               grepl("SP_02|SP_04", MP) ~
                   "Model-based: State-Space Surplus Production"))

p01 <- ggplot(data = temp,
              aes(x = PM, y = Values, fill = Class,
                  colour = Class)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5,
                scale = "width") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.2,
                 colour = "darkred") +
    facet_wrap(~MP, ncol = 3) +
    geom_hline(yintercept = 0.25, linetype = "dashed", colour = "gray10") +
    scale_fill_manual(values = my_pal[c(3, 5, 8, 13, 15)]) +
    scale_colour_manual(values = my_pal[c(3, 5, 8, 13, 15)]) +
    coord_cartesian(ylim = c(0, 2)) +
    my_theme() +
    labs(x = "Performance Metric",
         y = expression(VarC[Stability])) +
    theme(legend.position = "none")
p01

## p01a <- ggplot(data = filter(temp,
##                             !MP %in% c("IS_01", "IS_02", "IS_03"),
##                             Class == "Reference"),
##               aes(x = PM, y = Values, fill = Class,
##                   colour = Class)) +
##     geom_violin(draw_quantiles = c(0.1, 0.9), alpha = 0.5,
##                 scale = "width") +
##     stat_summary(fun = mean, geom = "crossbar", width = 0.2,
##                  colour = "darkred") +
##     facet_grid(Class~MP) +
##     geom_hline(yintercept = 0.25, linetype = "dashed",
##                colour = "gray10") +
##     coord_cartesian(ylim = c(0, 2)) +
##     scale_fill_manual(values = my_pal[3]) +
##     scale_colour_manual(values = my_pal[3]) +
##     my_theme() +
##     labs(x = "Performance Metric",
##          y = expression(VarC[Stability])) +
##     theme(legend.position = "none")
## p01a

p01b <- ggplot(data = filter(temp,
                             Class == "Index-based: Iratio"),
               aes(x = PM, y = Values, fill = Class,
                   colour = Class)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5,
                scale = "width") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.2,
                 colour = "darkred") +
    facet_grid(Class~MP) +
    geom_hline(yintercept = 0.25, linetype = "dashed",
               colour = "gray10") +
    coord_cartesian(ylim = c(0, 2)) +
    scale_fill_manual(values = my_pal[5]) +
    scale_colour_manual(values = my_pal[5]) +
    my_theme() +
    labs(x = "Performance Metric",
         y = expression(VarC[Stability])) +
    theme(legend.position = "none")
p01b

p01c <- ggplot(data = filter(temp,
                             Class == "Index-based: Exploitation rate"),
               aes(x = PM, y = Values, fill = Class,
                   colour = Class)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5,
                scale = "width") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.2,
                 colour = "darkred") +
    facet_grid(Class~MP) +
    geom_hline(yintercept = 0.25, linetype = "dashed",
               colour = "gray10") +
    coord_cartesian(ylim = c(0, 2)) +
    scale_fill_manual(values = my_pal[8]) +
    scale_colour_manual(values = my_pal[8]) +
    my_theme() +
    labs(x = "Performance Metric",
         y = expression(VarC[Stability])) +
    theme(legend.position = "none")
p01c

p01d <- ggplot(data = filter(temp,
                             Class == "Model-based: Surplus Production"),
               aes(x = PM, y = Values, fill = Class,
                   colour = Class)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5,
                scale = "width") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.2,
                 colour = "darkred") +
    facet_grid(Class~MP) +
    geom_hline(yintercept = 0.25, linetype = "dashed",
               colour = "gray10") +
    coord_cartesian(ylim = c(0, 2)) +
    scale_fill_manual(values = my_pal[13]) +
    scale_colour_manual(values = my_pal[13]) +
    my_theme() +
    labs(x = "Performance Metric",
         y = expression(VarC[Stability])) +
    theme(legend.position = "none")
p01d

p01e <-
    ggplot(data =
               filter(temp,
                      Class == "Model-based: State-Space Surplus Production"),
               aes(x = PM, y = Values, fill = Class,
                   colour = Class)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5,
                scale = "width") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.2,
                 colour = "darkred") +
    facet_grid(Class~MP) +
    geom_hline(yintercept = 0.25, linetype = "dashed",
               colour = "gray10") +
    coord_cartesian(ylim = c(0, 2)) +
    scale_fill_manual(values = my_pal[15]) +
    scale_colour_manual(values = my_pal[15]) +
    my_theme() +
    labs(x = "Performance Metric",
         y = expression(VarC[Stability])) +
    theme(legend.position = "none")
p01e

#####@> Exporting figures...
## ggsave("05_Results/SCRS_Doc/Fig02_Violin_VarC_RefMPs_Reference_ver00.png",
##        plot = p01a, device = "png", units = "cm", w = 40, h = 25,
##        dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig02_Violin_VarC_IRMPs_Reference_ver00.png",
       plot = p01b, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig02_Violin_VarC_CEMPs_Reference_ver00.png",
       plot = p01c, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig02_Violin_VarC_SPMPs_Reference_ver00.png",
       plot = p01d, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig02_Violin_VarC_SPSSMPs_Reference_ver00.png",
       plot = p01e, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

#####@> LRP...
temp <- PM_output %>%
    filter(PM %in% c("LRP", "LRP_short", "LRP_med", "LRP_long")) %>%
    mutate(PM = factor(PM, levels = c("LRP_short", "LRP_med",
                                      "LRP_long", "LRP")),
           MP = factor(MP, levels = c("FMSYref", "FMSYref75",
                                      "FMSYref110", "IR_01", "IR_02",
                                      "IR_03", "CE_01", "CE_02",
                                      "CE_03", "SP_01", "SP_02",
                                      "SP_03", "SP_04"))) %>%
    mutate(Class = case_when(
               grepl("FMSY", MP) ~ "Reference",
               grepl("IR", MP) ~ "Index-based: Iratio",
               grepl("CE", MP) ~ "Index-based: Exploitation rate",
               grepl("SP_01|SP_03", MP) ~
                   "Model-based: Surplus Production",
               grepl("SP_02|SP_04", MP) ~
                   "Model-based: State-Space Surplus Production"))

p02 <- ggplot(data = temp,
              aes(x = PM, y = Values, fill = Class,
                  colour = Class)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5,
                scale = "width") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.2,
                 colour = "darkred") +
    facet_wrap(~MP, ncol = 3) +
    geom_hline(yintercept = 0.1, linetype = "dashed",
               colour = "gray10") +
    coord_cartesian(ylim = c(0, 1)) +
    scale_fill_manual(values = my_pal[c(3, 5, 8, 13, 15)]) +
    scale_colour_manual(values = my_pal[c(3, 5, 8, 13, 15)]) +
    my_theme() +
    labs(x = "Performance Metric",
         y = expression(LRP[Safety])) +
    theme(legend.position = "none")
p02

p02a <- ggplot(data = filter(temp,
                            Class == "Reference"),
              aes(x = PM, y = Values, fill = Class,
                  colour = Class)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5,
                scale = "width") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.2,
                 colour = "darkred") +
    facet_grid(Class~MP) +
    geom_hline(yintercept = 0.1, linetype = "dashed",
               colour = "gray10") +
    coord_cartesian(ylim = c(0, 1)) +
    scale_fill_manual(values = my_pal[3]) +
    scale_colour_manual(values = my_pal[3]) +
    my_theme() +
    labs(x = "Performance Metric",
         y = expression(LRP[Safety])) +
    theme(legend.position = "none")
p02a

p02b <- ggplot(data = filter(temp,
                             Class == "Index-based: Iratio"),
               aes(x = PM, y = Values, fill = Class,
                   colour = Class)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5,
                scale = "width") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.2,
                 colour = "darkred") +
    facet_grid(Class~MP) +
    geom_hline(yintercept = 0.1, linetype = "dashed",
               colour = "gray10") +
    coord_cartesian(ylim = c(0, 1)) +
    scale_fill_manual(values = my_pal[5]) +
    scale_colour_manual(values = my_pal[5]) +
    my_theme() +
    labs(x = "Performance Metric",
         y = expression(LRP[Safety])) +
    theme(legend.position = "none")
p02b

p02c <- ggplot(data = filter(temp,
                             Class == "Index-based: Exploitation rate"),
               aes(x = PM, y = Values, fill = Class,
                   colour = Class)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5,
                scale = "width") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.2,
                 colour = "darkred") +
    facet_grid(Class~MP) +
    geom_hline(yintercept = 0.1, linetype = "dashed",
               colour = "gray10") +
    coord_cartesian(ylim = c(0, 1)) +
    scale_fill_manual(values = my_pal[8]) +
    scale_colour_manual(values = my_pal[8]) +
    my_theme() +
    labs(x = "Performance Metric",
         y = expression(LRP[Safety])) +
    theme(legend.position = "none")
p02c

p02d <- ggplot(data = filter(temp,
                             Class == "Model-based: Surplus Production"),
               aes(x = PM, y = Values, fill = Class,
                   colour = Class)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5,
                scale = "width") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.2,
                 colour = "darkred") +
    facet_grid(Class~MP) +
    geom_hline(yintercept = 0.1, linetype = "dashed",
               colour = "gray10") +
    coord_cartesian(ylim = c(0, 1)) +
    scale_fill_manual(values = my_pal[13]) +
    scale_colour_manual(values = my_pal[13]) +
    my_theme() +
    labs(x = "Performance Metric",
         y = expression(LRP[Safety])) +
    theme(legend.position = "none")
p02d

p02e <-
    ggplot(data =
               filter(temp,
                      Class == "Model-based: State-Space Surplus Production"),
               aes(x = PM, y = Values, fill = Class,
                   colour = Class)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5,
                scale = "width") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.2,
                 colour = "darkred") +
    facet_grid(Class~MP) +
    geom_hline(yintercept = 0.1, linetype = "dashed",
               colour = "gray10") +
    coord_cartesian(ylim = c(0, 1)) +
    scale_fill_manual(values = my_pal[15]) +
    scale_colour_manual(values = my_pal[15]) +
    my_theme() +
    labs(x = "Performance Metric",
         y = expression(LRP[Safety])) +
    theme(legend.position = "none")
p02e

#####@> Exporting figures...
ggsave("05_Results/SCRS_Doc/Fig02_Violin_LRP_RefMPs_Reference_ver00.png",
       plot = p02a, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig02_Violin_LRP_IRMPs_Reference_ver00.png",
       plot = p02b, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig02_Violin_LRP_CEMPs_Reference_ver00.png",
       plot = p02c, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig02_Violin_LRP_SPMPs_Reference_ver00.png",
       plot = p02d, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig02_Violin_LRP_SPSSMPs_Reference_ver00.png",
       plot = p02e, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

#####@> AvC...
temp <- PM_output %>%
    filter(PM %in% c("AvC_long", "AvC_med", "AvC_short")) %>%
    mutate(PM = factor(PM, levels = c("AvC_short", "AvC_med",
                                      "AvC_long")),
           MP = factor(MP, levels = c("FMSYref", "FMSYref75",
                                      "FMSYref110", "IR_01", "IR_02",
                                      "IR_03", "CE_01", "CE_02",
                                      "CE_03", "SP_01", "SP_02",
                                      "SP_03", "SP_04"))) %>%
    mutate(Class = case_when(
               grepl("FMSY", MP) ~ "Reference",
               grepl("IR", MP) ~ "Index-based: Iratio",
               grepl("CE", MP) ~ "Index-based: Exploitation rate",
               grepl("SP_01|SP_03", MP) ~
                   "Model-based: Surplus Production",
               grepl("SP_02|SP_04", MP) ~
                   "Model-based: State-Space Surplus Production"))

p03 <- ggplot(data = temp,
              aes(x = PM, y = Values, fill = Class,
                  colour = Class)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5,
                scale = "width") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.2,
                 colour = "darkred") +
    facet_wrap(~MP, ncol = 3) +
    ## geom_hline(yintercept = 0, linetype = "dashed", colour = "gray10") +
    scale_fill_manual(values = my_pal[c(3, 5, 8, 13, 15)]) +
    scale_colour_manual(values = my_pal[c(3, 5, 8, 13, 15)]) +
    my_theme() +
    labs(x = "Performance Metric",
         y = expression(AvC[Yield])) +
    theme(legend.position = "none")
p03

p03a <- ggplot(data = filter(temp,
                             Class == "Reference"),
              aes(x = PM, y = Values, fill = Class,
                  colour = Class)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5,
                scale = "width") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.2,
                 colour = "darkred") +
    facet_grid(Class~MP) +
    ## geom_hline(yintercept = 0.1, linetype = "dashed",
    ##            colour = "gray10") +
    ## coord_cartesian(ylim = c(0, 1)) +
    scale_fill_manual(values = my_pal[3]) +
    scale_colour_manual(values = my_pal[3]) +
    my_theme() +
    labs(x = "Performance Metric",
         y = expression(AvC[Yield])) +
    theme(legend.position = "none")
p03a

p03b <- ggplot(data = filter(temp,
                             Class == "Index-based: Iratio"),
               aes(x = PM, y = Values, fill = Class,
                   colour = Class)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5,
                scale = "width") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.2,
                 colour = "darkred") +
    facet_grid(Class~MP) +
    ## geom_hline(yintercept = 0.1, linetype = "dashed",
    ##            colour = "gray10") +
    ## coord_cartesian(ylim = c(0, 1)) +
    scale_fill_manual(values = my_pal[5]) +
    scale_colour_manual(values = my_pal[5]) +
    my_theme() +
    labs(x = "Performance Metric",
         y = expression(AvC[Yield])) +
    theme(legend.position = "none")
p03b

p03c <- ggplot(data = filter(temp,
                             Class == "Index-based: Exploitation rate"),
               aes(x = PM, y = Values, fill = Class,
                   colour = Class)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5,
                scale = "width") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.2,
                 colour = "darkred") +
    facet_grid(Class~MP) +
    ## geom_hline(yintercept = 0.1, linetype = "dashed",
    ##            colour = "gray10") +
    ## coord_cartesian(ylim = c(0, 1)) +
    scale_fill_manual(values = my_pal[8]) +
    scale_colour_manual(values = my_pal[8]) +
    my_theme() +
    labs(x = "Performance Metric",
         y = expression(AvC[Yield])) +
    theme(legend.position = "none")
p03c

p03d <- ggplot(data = filter(temp,
                             Class == "Model-based: Surplus Production"),
               aes(x = PM, y = Values, fill = Class,
                   colour = Class)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5,
                scale = "width") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.2,
                 colour = "darkred") +
    facet_grid(Class~MP) +
    ## geom_hline(yintercept = 0.1, linetype = "dashed",
    ##            colour = "gray10") +
    ## coord_cartesian(ylim = c(0, 1)) +
    scale_fill_manual(values = my_pal[13]) +
    scale_colour_manual(values = my_pal[13]) +
    my_theme() +
    labs(x = "Performance Metric",
         y = expression(AvC[Yield])) +
    theme(legend.position = "none")
p03d

p03e <-
    ggplot(data =
               filter(temp,
                      Class == "Model-based: State-Space Surplus Production"),
               aes(x = PM, y = Values, fill = Class,
                   colour = Class)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5,
                scale = "width") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.2,
                 colour = "darkred") +
    facet_grid(Class~MP) +
    ## geom_hline(yintercept = 0.1, linetype = "dashed",
    ##            colour = "gray10") +
    ## coord_cartesian(ylim = c(0, 1)) +
    scale_fill_manual(values = my_pal[15]) +
    scale_colour_manual(values = my_pal[15]) +
    my_theme() +
    labs(x = "Performance Metric",
         y = expression(AvC[Yield])) +
    theme(legend.position = "none")
p03e

#####@> Exporting figures...
ggsave("05_Results/SCRS_Doc/Fig02_Violin_AvC_RefMPs_Reference_ver00.png",
       plot = p03a, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig02_Violin_AvC_IRMPs_Reference_ver00.png",
       plot = p03b, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig02_Violin_AvC_CEMPs_Reference_ver00.png",
       plot = p03c, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig02_Violin_AvC_SPMPs_Reference_ver00.png",
       plot = p03d, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig02_Violin_AvC_SPSSMPs_Reference_ver00.png",
       plot = p03e, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

#####@> Looking to the PMs - Estimating the average between all OMs...
tab01 <- PM_output %>%
    group_by(MP, PM, Name, Caption) %>%
    summarise(q50 = mean(Values, na.rm = TRUE)) %>%
    as.data.frame()

######@> Consolidating general probabilities for a heatmap table...
ProbTab <- tab01 %>%
    filter(PM %in% c("LRP", "LRP_long", "LRP_med", "LRP_short", "PGK",
                     "PGK_long", "PGK_med", "PGK_short", "PNOF", "POF",
                     "VarC","VarClong", "VarCmedium", "nLRP", "nLRP_long",
                     "nLRP_med", "nLRP_short", "rAvC_long", "rAvC_med",
                     "rAvC_short", "AvC_short", "AvC_med", "AvC_long")) %>%
    select(MP, PM, q50) %>%
    pivot_wider(names_from = "PM", values_from = "q50") %>%
    as.data.frame()

#####@> Guilty plot...
temp <- ProbTab %>%
    ## filter(MP %in% c("IR_01", "IR_02", "IR_03",
    ##                  "CE_01", "CE_02", "CE_03",
    ##                  "SP_01", "SP_02", "SP_03", "SP_04")) %>%
    select(MP, AvC_short, AvC_med, AvC_long, PGK_short, PGK_med,
           PGK_long, PGK, PNOF, nLRP_short, nLRP_med, nLRP_long, nLRP,
           VarCmedium, VarClong, VarC) %>%
    mutate(MP = factor(MP, levels = c("FMSYref", "FMSYref75",
                                      "FMSYref110", "IR_01", "IR_02",
                                      "IR_03", "CE_01", "CE_02",
                                      "CE_03", "SP_01", "SP_02",
                                      "SP_03", "SP_04"))) %>%
    mutate_if(is.numeric, ~replace_na(.x, NA)) %>%
    mutate(AvC_short = round(AvC_short),
           AvC_med = round(AvC_med),
           AvC_long = round(AvC_long)) %>%
    arrange(MP)

guiltplot <- temp %>%
    gt() %>%
    fmt_number(decimals = 2, columns = c(PGK_short, PGK_med,
               PGK_long, PGK, PNOF, nLRP_short, nLRP_med,
               nLRP_long, nLRP, VarCmedium, VarClong, VarC)) %>%
    data_color(
        columns = c(AvC_short, AvC_med, AvC_long, PGK_short, PGK_med,
                    PGK_long, PGK, PNOF, nLRP_short, nLRP_med,
                    nLRP_long, nLRP, VarCmedium, VarClong, VarC),
        colors = scales::col_bin(palette = c("white", "#8E9C97"),
                                 domain = NULL, bins = 4))

gtsave(data = guiltplot,
       filename = "Fig06_Guiltplot_Reference_ver01.html",
       path = "05_Results/SCRS_Doc/")

#####@> for SCRS SCI 115...
guiltplot <- temp %>%
    filter(!MP %in% c("FMSYref", "FMSYref75", "FMSYref110")) %>%
    gt() %>%
    fmt_number(decimals = 2,
               columns = c(PGK_short, PGK_med,
                           PGK_long, PGK, PNOF, nLRP_short, nLRP_med,
                           nLRP_long, nLRP, VarCmedium, VarClong, VarC)) %>%
    data_color(
        columns = c(AvC_short, AvC_med, AvC_long, PGK_short, PGK_med,
                    PGK_long, PGK, PNOF, nLRP_short, nLRP_med,
                    nLRP_long, nLRP, VarCmedium, VarClong, VarC),
        colors = scales::col_bin(palette = c("white", "#8E9C97"),
                                 domain = NULL, bins = 4))

gtsave(data = guiltplot,
       filename = "Fig06_Guiltplot_Reference_ver02.html",
       path = "05_Results/Response_Comm/")

######@> Heatmaps for ProbTab...

#####@> Figure for Yield...
tmp <- ProbTab %>%
    select(MP:AvC_short) %>%
    filter(MP %in% c("FMSYref", "FMSYref75", "FMSYref110",
                     "IR_01", "IR_02", "IR_03",
                     "CE_01", "CE_02", "CE_03",
                     "SP_01", "SP_02", "SP_03", "SP_04")) %>%
    pivot_longer(names_to = "PM", values_to = "value", 2:4) %>%
    mutate(class = "Yield",
           breaks = cut(value,
                        breaks = quantile(value, seq(0, 1, 0.2)),
                        include.lowest = TRUE, right = FALSE),
           PM = factor(PM, levels = c("AvC_short",
                                      "AvC_med",
                                      "AvC_long")),
           MP = factor(MP,
                       levels = rev(c("FMSYref", "FMSYref75",
                                      "FMSYref110",
                                      "IR_01", "IR_02", "IR_03",
                                      "CE_01", "CE_02", "CE_03",
                                      "SP_01", "SP_02", "SP_03",
                                      "SP_04"))),
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
    filter(MP %in% c("FMSYref", "FMSYref75", "FMSYref110",
                     "IR_01", "IR_02", "IR_03",
                     "CE_01", "CE_02", "CE_03",
                     "SP_01", "SP_02", "SP_03", "SP_04")) %>%
    pivot_longer(names_to = "PM", values_to = "value", 2:5) %>%
    mutate(class = "Status",
           ## breaks = cut(value,
           ##              breaks = quantile(value, seq(0, 1, 0.2)),
           ##              include.lowest = TRUE, right = FALSE),
           PM = factor(PM, levels = c("PGK_short",
                                      "PGK_med",
                                      "PGK_long",
                                      "PGK")),
           MP = factor(MP,
                       levels = rev(c("FMSYref", "FMSYref75",
                                      "FMSYref110", 
                                      "IR_01", "IR_02", "IR_03",
                                      "CE_01", "CE_02", "CE_03",
                                      "SP_01", "SP_02", "SP_03",
                                      "SP_04"))),
           class2 = ifelse(value >= 0.7, "Up", "Down"))

p01 <- ggplot(data = tmp, aes(x = PM, y = MP, fill = class2)) +
    geom_tile(colour = "black", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.2f", value)), size = 5,
              fontface = "bold") +
    scale_fill_brewer(palette = "Greys", direction = 1) +
    facet_wrap(~class) +
    labs(x = "Performance Metrics", y = "Management Procedures") +
    my_theme() +
    theme(legend.position = "none")
p01

tmp <- ProbTab %>%
    select(MP, PNOF:POF) %>%
    filter(MP %in% c("FMSYref", "FMSYref75", "FMSYref110",
                     "IR_01", "IR_02", "IR_03",
                     "CE_01", "CE_02", "CE_03",
                     "SP_01", "SP_02", "SP_03", "SP_04")) %>%
    pivot_longer(names_to = "PM", values_to = "value", 2:3) %>%
    mutate(class = "Status",
           ## breaks = cut(value,
           ##              breaks = quantile(value, seq(0, 1, 0.2)),
           ##              include.lowest = TRUE, right = FALSE),
           PM = factor(PM, levels = c("PNOF",
                                      "POF")),
           MP = factor(MP,
                       levels = rev(c("FMSYref", "FMSYref75",
                                      "FMSYref110",
                                      "IR_01", "IR_02", "IR_03",
                                      "CE_01", "CE_02", "CE_03",
                                      "SP_01", "SP_02", "SP_03",
                                      "SP_04"))), 
           class2 = ifelse(value >= 0.7, "Up", "Down"))
p01b <- ggplot(data = tmp, aes(x = PM, y = MP)) +
    geom_tile(colour = "black", fill = "white", alpha = 0.1) +
    geom_text(aes(label = sprintf("%.2f", value)), size = 5,
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
    filter(MP %in% c("FMSYref", "FMSYref75", "FMSYref110",
                     "IR_01", "IR_02", "IR_03",
                     "CE_01", "CE_02", "CE_03",
                     "SP_01", "SP_02", "SP_03", "SP_04")) %>%
    pivot_longer(names_to = "PM", values_to = "value", 2:5) %>%
    mutate(class = "Safety",
           ## breaks = cut(value,
           ##              breaks = quantile(value, seq(0, 1, 0.2)),
           ##              include.lowest = TRUE, right = FALSE),
           PM = factor(PM, levels = c("LRP_short",
                                      "LRP_med",
                                      "LRP_long",
                                      "LRP")),
           MP = factor(MP,
                       levels = rev(c("FMSYref", "FMSYref75", "FMSYref110",
                                      "IR_01", "IR_02", "IR_03",
                                      "CE_01", "CE_02", "CE_03",
                                      "SP_01", "SP_02", "SP_03", "SP_04"))),
           class2 = ifelse(value <= 0.1, "Up", "Down"))
p02 <- ggplot(data = tmp, aes(x = PM, y = MP, fill = class2)) +
    geom_tile(colour = "black", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.2f", value)), size = 5,
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
    filter(MP %in% c(## "FMSYref", "FMSYref75", "FMSYref110",
                     "IR_01", "IR_02", "IR_03",
                     "CE_01", "CE_02", "CE_03",
                     "SP_01", "SP_02", "SP_03", "SP_04")) %>%
    pivot_longer(names_to = "PM", values_to = "value", 2:4) %>%
    mutate(class = "Stability",
           ## breaks = cut(value,
           ##              breaks = quantile(value, seq(0, 1, 0.2)),
           ##              include.lowest = TRUE, right = FALSE),
           PM = factor(PM, levels = c("VarCmedium",
                                      "VarClong",
                                      "VarC")),
           MP = factor(MP,
                       levels = rev(c(## "FMSYref", "FMSYref75",
                           ## "FMSYref110",
                                      "IR_01", "IR_02", "IR_03",
                                      "CE_01", "CE_02", "CE_03",
                                      "SP_01", "SP_02", "SP_03",
                           "SP_04"))),
           class2 = ifelse(value <= 0.2, "Up", "Down"))
p03 <- ggplot(data = tmp, aes(x = PM, y = MP, fill = class2)) +
    geom_tile(colour = "black", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.2f", value)), size = 5,
              fontface = "bold") +
    scale_fill_brewer(palette = "Greys", direction = 1) +
    facet_wrap(~class) +
    labs(x = "Performance Metrics", y = "Management Procedures") +
    my_theme() +
    theme(legend.position = "none")
p03

#####@> Exporting figures...
ggsave("05_Results/SCRS_Doc/Fig03_Table_Yield_Reference_ver00.png",
       plot = p00, device = "png", units = "cm", w = 25, h = 40,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig03_Table_Status_Reference_ver00.png",
       plot = p01, device = "png", units = "cm", w = 25, h = 40,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig03_Table_Status2_Reference_ver00.png",
       plot = p01b, device = "png", units = "cm", w = 25, h = 40,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig03_Table_Safety_Reference_ver00.png",
       plot = p02, device = "png", units = "cm", w = 25, h = 40,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig03_Table_Stability_Reference_ver00.png",
       plot = p03, device = "png", units = "cm", w = 25, h = 40,
       dpi = 500, bg = "white")

######@> Extracting References Points from MSE objects...

#####@> Extracting SB_SBMSY values - SubMPs...

####@> Mode old school...
## tmp01 <- list()
## for(i in seq_along(MSEsSub)) {
##     om <- names(MSEsSub)[i]
##     output <- MSEsSub[[i]]@SB_SBMSY %>%
##         structure(
##             dimnames = list(
##                 Sim = 1:MSEsSub[[i]]@nsim,
##                 MPs = MSEsSub[[i]]@MPs,
##                 Year = 2020 + 1:MSEsSub[[i]]@proyears
##             )
##         ) %>%
##         reshape2::melt() %>%
##         mutate(OM = om)
##     tmp01[[i]] <- output
## }
## temp00 <- do.call(rbind.data.frame, tmp01)

####@> New mode - Sub MPs...
temp00 <- map2(MSEsSub, names(MSEsSub), function(mse, om) {
    output <- mse@SB_SBMSY %>%
        structure(
            dimnames = list(
                Sim = 1:mse@nsim,
                MPs = mse@MPs,
                Year = 2020 + 1:mse@proyears
            )
        ) %>%
        as.data.frame.table() %>%
        as_tibble() %>%
        ## { print(colnames(.)); . } %>%
        rename(Sim = 1, MPs = 2, Year = 3, Value = Freq) %>%
        mutate(OM = om) %>%
        as.data.frame
    return(output)
}) %>%
    bind_rows()

####@> New mode - SPSS...
temp01 <- map2(MSEsSPSS, names(MSEsSPSS), function(mse, om) {
    output <- mse@SB_SBMSY %>%
        structure(
            dimnames = list(
                Sim = 1:mse@nsim,
                MPs = mse@MPs,
                Year = 2020 + 1:mse@proyears
            )
        ) %>%
        as.data.frame.table() %>%
        as_tibble() %>%
        ## { print(colnames(.)); . } %>%
        rename(Sim = 1, MPs = 2, Year = 3, Value = Freq) %>%
        mutate(OM = om) %>%
        as.data.frame
    return(output)
}) %>%
    bind_rows()

####@> Binding everything...
SB_SBMSY_output <- bind_rows(temp00, temp01) %>%
    mutate(Year = as.numeric(as.character(Year)))

####@> Looking to the RP - Estimating the average between all MOMs...
tab02 <- SB_SBMSY_output %>%
    group_by(MPs, Year) %>%
    summarise(SB_SBMSY = median(Value, na.rm = TRUE),
              SB_SBMSY_Lwr = quantile(Value, probs = 0.025,
                                      na.rm = TRUE),
              SB_SBMSY_Upr = quantile(Value, probs = 0.975,
                                      na.rm = TRUE),
              SB_SBMSY_Lwr2 = quantile(Value, probs = 0.15,
                                      na.rm = TRUE),
              SB_SBMSY_Upr2 = quantile(Value, probs = 0.85,
                                       na.rm = TRUE)) %>%
    mutate(Class = case_when(
               grepl("FMSY", MPs) ~ "Reference",
               grepl("IR", MPs) ~ "Index-based: Iratio",
               grepl("CE", MPs) ~ "Index-based: Exploitation rate",
               grepl("SP_01|SP_03", MPs) ~
                   "Model-based: Surplus Production",
               grepl("SP_02|SP_04", MPs) ~
                   "Model-based: State-Space Surplus Production")) %>%
    as.data.frame()

###@> Figure Reference Points...
p01a <- ggplot(data =
                  filter(tab02,
                         Class == "Reference",
                         Year %in% 2021:2050)) +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr,
                    ymax = SB_SBMSY_Upr),
                alpha = 0.5, fill = my_pal[3]) +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr2,
                    ymax = SB_SBMSY_Upr2),
                alpha = 0.5, fill = my_pal[3]) +
    geom_line(aes(x = Year, y = SB_SBMSY), linewidth = 1,
              colour = "black") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_grid(Class~MPs) +
    scale_y_continuous(limits = c(0, 5), expand = c(0, 0)) +
    coord_cartesian(xlim = c(2021, 2050)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
    labs(x = "Year", y = expression(SSB/SSB[MSY])) +
    my_theme()
p01a

p01b <- ggplot(data =
                   filter(tab02,
                          Class == "Index-based: Iratio",
                          Year %in% 2021:2050)) +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr,
                    ymax = SB_SBMSY_Upr),
                alpha = 0.5, fill = my_pal[3]) +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr2,
                    ymax = SB_SBMSY_Upr2),
                alpha = 0.5, fill = my_pal[3]) +
    geom_line(aes(x = Year, y = SB_SBMSY), linewidth = 1,
              colour = "black") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_grid(Class~MPs) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 5)) +
    coord_cartesian(xlim = c(2021, 2050)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
    labs(x = "Year", y = expression(SSB/SSB[MSY])) +
    my_theme()
p01b

p01c <- ggplot(data =
                   filter(tab02,
                          Class == "Index-based: Exploitation rate",
                          Year %in% 2021:2050)) +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr,
                    ymax = SB_SBMSY_Upr),
                alpha = 0.5, fill = my_pal[3]) +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr2,
                    ymax = SB_SBMSY_Upr2),
                alpha = 0.5, fill = my_pal[3]) +
    geom_line(aes(x = Year, y = SB_SBMSY), linewidth = 1,
              colour = "black") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_grid(Class~MPs) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 5)) +
    coord_cartesian(xlim = c(2021, 2050)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
    labs(x = "Year", y = expression(SSB/SSB[MSY])) +
    my_theme()
p01c

p01d <- ggplot(data =
                   filter(tab02,
                          Class == "Model-based: Surplus Production",
                          Year %in% 2021:2050)) +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr,
                    ymax = SB_SBMSY_Upr),
                alpha = 0.5, fill = my_pal[3]) +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr2,
                    ymax = SB_SBMSY_Upr2),
                alpha = 0.5, fill = my_pal[3]) +
    geom_line(aes(x = Year, y = SB_SBMSY), linewidth = 1,
              colour = "black") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_grid(Class~MPs) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 5)) +
    coord_cartesian(xlim = c(2021, 2050)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
    labs(x = "Year", y = expression(SSB/SSB[MSY])) +
    my_theme()
p01d

p01e <- ggplot(data =
                   filter(tab02,
                          Class == "Model-based: State-Space Surplus Production",
                          Year %in% 2021:2050)) +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr,
                    ymax = SB_SBMSY_Upr),
                alpha = 0.5, fill = my_pal[3]) +
    geom_ribbon(aes(x = Year,
                    ymin = SB_SBMSY_Lwr2,
                    ymax = SB_SBMSY_Upr2),
                alpha = 0.5, fill = my_pal[3]) +
    geom_line(aes(x = Year, y = SB_SBMSY), linewidth = 1,
              colour = "black") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_grid(Class~MPs) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 8)) +
    coord_cartesian(xlim = c(2021, 2050)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
    labs(x = "Year", y = expression(SSB/SSB[MSY])) +
    my_theme()
p01e

#####@> Exporting figures...
ggsave("05_Results/SCRS_Doc/Fig04_TS_SSB_SSBMSY_RefsMPs_Reference_ver00.png",
       plot = p01a, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig04_TS_SSB_SSBMSY_IRMPs_Reference_ver00.png",
       plot = p01b, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig04_TS_SSB_SSBMSY_CEMPs_Reference_ver00.png",
       plot = p01c, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig04_TS_SSB_SSBMSY_SPMPs_Reference_ver00.png",
       plot = p01d, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig04_TS_SSB_SSBMSY_SPSSMPs_Reference_ver00.png",
       plot = p01e, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

#####@> Extracting F_FMSY values...
## tmp02 <- list()
## for(i in seq_along(MSEs)) {
##     om <- names(MSEs)[i]
##     output <- MSEs[[i]]@F_FMSY %>%
##         structure(
##             dimnames = list(
##                 Sim = 1:MSEs[[i]]@nsim,
##                 MPs = MSEs[[i]]@MPs,
##                 Year = 2020 + 1:MSEs[[i]]@proyears
##             )
##         ) %>%
##         reshape2::melt() %>%
##         mutate(OM = om) %>%
##         as.data.frame()
##     tmp02[[i]] <- output
## }
## F_FMSY_output <- do.call(rbind.data.frame, tmp02)

####@> New mode - Sub MPs...
temp00 <- map2(MSEsSub, names(MSEsSub), function(mse, om) {
    output <- mse@F_FMSY %>%
        structure(
            dimnames = list(
                Sim = 1:mse@nsim,
                MPs = mse@MPs,
                Year = 2020 + 1:mse@proyears
            )
        ) %>%
        as.data.frame.table() %>%
        as_tibble() %>%
        ## { print(colnames(.)); . } %>%
        rename(Sim = 1, MPs = 2, Year = 3, Value = Freq) %>%
        mutate(OM = om) %>%
        as.data.frame
    return(output)
}) %>%
    bind_rows()

####@> New mode - SPSS...
temp01 <- map2(MSEsSPSS, names(MSEsSPSS), function(mse, om) {
    output <- mse@F_FMSY %>%
        structure(
            dimnames = list(
                Sim = 1:mse@nsim,
                MPs = mse@MPs,
                Year = 2020 + 1:mse@proyears
            )
        ) %>%
        as.data.frame.table() %>%
        as_tibble() %>%
        ## { print(colnames(.)); . } %>%
        rename(Sim = 1, MPs = 2, Year = 3, Value = Freq) %>%
        mutate(OM = om) %>%
        as.data.frame
    return(output)
}) %>%
    bind_rows()

####@> Binding everything...
F_FMSY_output <- bind_rows(temp00, temp01) %>%
    mutate(Year = as.numeric(as.character(Year)))

####@> Looking to the RP - Estimating the average between all MOMs...
tab03 <- F_FMSY_output %>%
    group_by(MPs, Year) %>%
    summarise(F_FMSY = median(Value, na.rm = TRUE),
              F_FMSY_Lwr = quantile(Value, probs = 0.025,
                                    na.rm = TRUE),
              F_FMSY_Upr = quantile(Value, probs = 0.975,
                                    na.rm = TRUE),
              F_FMSY_Lwr2 = quantile(Value, probs = 0.15,
                                       na.rm = TRUE),
              F_FMSY_Upr2 = quantile(Value, probs = 0.85,
                                     na.rm = TRUE)) %>%
    mutate(Class = case_when(
               grepl("FMSY", MPs) ~ "Reference",
               grepl("IR", MPs) ~ "Index-based: Iratio",
               grepl("CE", MPs) ~ "Index-based: Exploitation rate",
               grepl("SP_01|SP_03", MPs) ~
                   "Model-based: Surplus Production",
               grepl("SP_02|SP_04", MPs) ~
                   "Model-based: State-Space Surplus Production")) %>%
    as.data.frame()

###@> Figure Reference Points...
p01a <- ggplot(data =
                  filter(tab03,
                         Class == "Reference",
                         Year %in% 2021:2050)) +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr,
                    ymax = F_FMSY_Upr),
                alpha = 0.5, fill = my_pal[3]) +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr2,
                    ymax = F_FMSY_Upr2),
                alpha = 0.5, fill = my_pal[3]) +
    geom_line(aes(x = Year, y = F_FMSY), linewidth = 1,
              colour = "black") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_grid(Class~MPs) +
    scale_y_continuous(limits = c(0, 3), expand = c(0, 0)) +
    coord_cartesian(xlim = c(2021, 2050)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
    labs(x = "Year", y = expression(F/F[MSY])) +
    my_theme()
p01a

p01b <- ggplot(data =
                   filter(tab03,
                          Class == "Index-based: Iratio",
                          Year %in% 2021:2050)) +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr,
                    ymax = F_FMSY_Upr),
                alpha = 0.5, fill = my_pal[3]) +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr2,
                    ymax = F_FMSY_Upr2),
                alpha = 0.5, fill = my_pal[3]) +
    geom_line(aes(x = Year, y = F_FMSY), linewidth = 1,
              colour = "black") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_grid(Class~MPs) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 8)) +
    coord_cartesian(xlim = c(2021, 2050), ylim = c(0, 6)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
    labs(x = "Year", y = expression(F/F[MSY])) +
    my_theme()
p01b

p01c <- ggplot(data =
                   filter(tab03,
                          Class == "Index-based: Exploitation rate",
                          Year %in% 2021:2050)) +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr,
                    ymax = F_FMSY_Upr),
                alpha = 0.5, fill = my_pal[3]) +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr2,
                    ymax = F_FMSY_Upr2),
                alpha = 0.5, fill = my_pal[3]) +
    geom_line(aes(x = Year, y = F_FMSY), linewidth = 1,
              colour = "black") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_grid(Class~MPs) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 8)) +
    coord_cartesian(xlim = c(2021, 2050), ylim = c(0, 6)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
    labs(x = "Year", y = expression(F/F[MSY])) +
    my_theme()
p01c

p01d <- ggplot(data =
                   filter(tab03,
                          Class == "Model-based: Surplus Production",
                          Year %in% 2021:2050)) +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr,
                    ymax = F_FMSY_Upr),
                alpha = 0.5, fill = my_pal[3]) +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr2,
                    ymax = F_FMSY_Upr2),
                alpha = 0.5, fill = my_pal[3]) +
    geom_line(aes(x = Year, y = F_FMSY), linewidth = 1,
              colour = "black") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_grid(Class~MPs) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 8)) +
    coord_cartesian(xlim = c(2021, 2050), ylim = c(0, 6)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
    labs(x = "Year", y = expression(F/F[MSY])) +
    my_theme()
p01d

p01e <- ggplot(data =
                   filter(tab03,
                          Class == "Model-based: State-Space Surplus Production",
                          Year %in% 2021:2050)) +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr,
                    ymax = F_FMSY_Upr),
                alpha = 0.5, fill = my_pal[3]) +
    geom_ribbon(aes(x = Year,
                    ymin = F_FMSY_Lwr2,
                    ymax = F_FMSY_Upr2),
                alpha = 0.5, fill = my_pal[3]) +
    geom_line(aes(x = Year, y = F_FMSY), linewidth = 1,
              colour = "black") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_grid(Class~MPs) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 8)) +
    coord_cartesian(xlim = c(2021, 2050), ylim = c(0, 6)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
    labs(x = "Year", y = expression(F/F[MSY])) +
    my_theme()
p01e

#####@> Exporting figures...
ggsave("05_Results/SCRS_Doc/Fig04_TS_F_FMSY_RefsMPs_Reference_ver00.png",
       plot = p01a, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig04_TS_F_FMSY_IRMPs_Reference_ver00.png",
       plot = p01b, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig04_TS_F_FMSY_CEMPs_Reference_ver00.png",
       plot = p01c, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig04_TS_F_FMSY_SPMPs_Reference_ver00.png",
       plot = p01d, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig04_TS_F_FMSY_SPSSMPs_Reference_ver00.png",
       plot = p01e, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

#####@> Extracting Catches values...
tmp03 <- list()
for(i in seq_along(MSEsSub)) {
    om <- names(MSEsSub)[i]
    output <- MSEsSub[[i]]@Catch %>%
        structure(
            dimnames = list(
                Sim = 1:MSEsSub[[i]]@nsim,
                MPs = MSEsSub[[i]]@MPs,
                Year = 2020 + 1:MSEsSub[[i]]@proyears
            )
        ) %>%
        reshape2::melt() %>%
        mutate(MOM = om) %>%
        as.data.frame()
    tmp03[[i]] <- output
}
temp00 <- do.call(rbind.data.frame, tmp03)

tmp03 <- list()
for(i in seq_along(MSEsSPSS)) {
    om <- names(MSEsSPSS)[i]
    output <- MSEsSPSS[[i]]@Catch %>%
        structure(
            dimnames = list(
                Sim = 1:MSEsSPSS[[i]]@nsim,
                MPs = MSEsSPSS[[i]]@MPs,
                Year = 2020 + 1:MSEsSPSS[[i]]@proyears
            )
        ) %>%
        reshape2::melt() %>%
        mutate(MOM = om) %>%
        as.data.frame()
    tmp03[[i]] <- output
}
temp01 <- do.call(rbind.data.frame, tmp03)
Catch_output <- bind_rows(temp00, temp01) %>%
    mutate(Year = as.numeric(as.character(Year)))


####@> New mode - Sub MPs...
temp00 <- map2(MSEsSub, names(MSEsSub), function(mse, om) {
    output <- mse@Catch %>%
        structure(
            dimnames = list(
                Sim = 1:mse@nsim,
                MPs = mse@MPs,
                Year = 2020 + 1:mse@proyears
            )
        ) %>%
        as.data.frame.table() %>%
        as_tibble() %>%
        ## { print(colnames(.)); . } %>%
        rename(Sim = 1, MPs = 2, Year = 3, Value = Freq) %>%
        mutate(OM = om) %>%
        as.data.frame
    return(output)
}) %>%
    bind_rows()

####@> New mode - SPSS...
temp01 <- map2(MSEsSPSS, names(MSEsSPSS), function(mse, om) {
    output <- mse@Catch %>%
        structure(
            dimnames = list(
                Sim = 1:mse@nsim,
                MPs = mse@MPs,
                Year = 2020 + 1:mse@proyears
            )
        ) %>%
        as.data.frame.table() %>%
        as_tibble() %>%
        ## { print(colnames(.)); . } %>%
        rename(Sim = 1, MPs = 2, Year = 3, Value = Freq) %>%
        mutate(OM = om) %>%
        as.data.frame
    return(output)
}) %>%
    bind_rows()

####@> Binding everything...
Catch_output <- bind_rows(temp00, temp01) %>%
    mutate(Year = as.numeric(as.character(Year)))

####@> Looking to the RP - Estimating the average between all MOMs...
tab04 <- Catch_output %>%
    group_by(MPs, Year) %>%
    summarise(Catch = median(Value, na.rm = TRUE),
              Catch_Lwr = quantile(Value, probs = 0.025,
                                   na.rm = TRUE),
              Catch_Upr = quantile(Value, probs = 0.975,
                                   na.rm = TRUE),
              Catch_Lwr2 = quantile(Value, probs = 0.15,
                                   na.rm = TRUE),
              Catch_Upr2 = quantile(Value, probs = 0.85,
                                    na.rm = TRUE)) %>%
    mutate(Class = case_when(
               grepl("FMSY", MPs) ~ "Reference",
               grepl("IR", MPs) ~ "Index-based: Iratio",
               grepl("CE", MPs) ~ "Index-based: Exploitation rate",
               grepl("SP_01|SP_03", MPs) ~
                   "Model-based: Surplus Production",
               grepl("SP_02|SP_04", MPs) ~
                   "Model-based: State-Space Surplus Production")) %>%
    as.data.frame()

###@> Figure Reference Points...
p01a <- ggplot(data =
                  filter(tab04,
                         Class == "Reference",
                         Year %in% 2021:2050)) +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr,
                    ymax = Catch_Upr),
                alpha = 0.5, fill = my_pal[3]) +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr2,
                    ymax = Catch_Upr2),
                alpha = 0.5, fill = my_pal[3]) +
    geom_line(aes(x = Year, y = Catch), linewidth = 1,
              colour = "black") +
    ## geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
    facet_grid(Class~MPs) +
    ## scale_y_continuous(limits = c(0, 3), expand = c(0, 0)) +
    coord_cartesian(xlim = c(2021, 2050)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
    labs(x = "Year", y = "Catches (t)") +
    my_theme()
p01a

p01b <- ggplot(data =
                   filter(tab04,
                          Class == "Index-based: Iratio",
                          Year %in% 2021:2050)) +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr,
                    ymax = Catch_Upr),
                alpha = 0.5, fill = my_pal[3]) +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr2,
                    ymax = Catch_Upr2),
                alpha = 0.5, fill = my_pal[3]) +
    geom_line(aes(x = Year, y = Catch), linewidth = 1,
              colour = "black") +
    geom_hline(yintercept = 21377, linetype = "dashed", colour = "gray30") +
    facet_grid(Class~MPs) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 60000)) +
    coord_cartesian(xlim = c(2021, 2050), ylim = c(0, 40000)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
    labs(x = "Year", y = "Catches (t)") +
    my_theme()
p01b

p01c <- ggplot(data =
                   filter(tab04,
                          Class == "Index-based: Exploitation rate",
                          Year %in% 2021:2050)) +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr,
                    ymax = Catch_Upr),
                alpha = 0.5, fill = my_pal[3]) +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr2,
                    ymax = Catch_Upr2),
                alpha = 0.5, fill = my_pal[3]) +
    geom_line(aes(x = Year, y = Catch), linewidth = 1,
              colour = "black") +
    geom_hline(yintercept = 21377, linetype = "dashed", colour = "gray30") +
    facet_grid(Class~MPs) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 60000)) +
    coord_cartesian(xlim = c(2021, 2050), ylim = c(0, 40000)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
    labs(x = "Year", y = "Catches (t)") +
    my_theme()
p01c

p01d <- ggplot(data =
                   filter(tab04,
                          Class == "Model-based: Surplus Production",
                          Year %in% 2021:2050)) +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr,
                    ymax = Catch_Upr),
                alpha = 0.5, fill = my_pal[3]) +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr2,
                    ymax = Catch_Upr2),
                alpha = 0.5, fill = my_pal[3]) +
    geom_line(aes(x = Year, y = Catch), linewidth = 1,
              colour = "black") +
    geom_hline(yintercept = 21377, linetype = "dashed", colour = "gray30") +
    facet_grid(Class~MPs) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 60000)) +
    coord_cartesian(xlim = c(2021, 2050), ylim = c(0, 40000)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
    labs(x = "Year", y = "Catches (t)") +
    my_theme()
p01d

p01e <- ggplot(data =
                   filter(tab04,
                          Class == "Model-based: State-Space Surplus Production",
                          Year %in% 2021:2050)) +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr,
                    ymax = Catch_Upr),
                alpha = 0.5, fill = my_pal[3]) +
    geom_ribbon(aes(x = Year,
                    ymin = Catch_Lwr2,
                    ymax = Catch_Upr2),
                alpha = 0.5, fill = my_pal[3]) +
    geom_line(aes(x = Year, y = Catch), linewidth = 1,
              colour = "black") +
    geom_hline(yintercept = 21377, linetype = "dashed", colour = "gray30") +
    facet_grid(Class~MPs) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 60000)) +
    coord_cartesian(xlim = c(2021, 2050), ylim = c(0, 40000)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
    labs(x = "Year", y = "Catches (t)") +
    my_theme()
p01e

#####@> Exporting figures...
ggsave("05_Results/SCRS_Doc/Fig04_TS_Catch_RefsMPs_Reference_ver00.png",
       plot = p01a, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig04_TS_Catch_IRMPs_Reference_ver00.png",
       plot = p01b, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig04_TS_Catch_CEMPs_Reference_ver00.png",
       plot = p01c, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig04_TS_Catch_SPMPs_Reference_ver00.png",
       plot = p01d, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig04_TS_Catch_SPSSMPs_Reference_ver00.png",
       plot = p01e, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

#####@> Extracting TAC values...

####@> New mode - Sub MPs...
temp00 <- map2(MSEsSub, names(MSEsSub), function(mse, om) {
    output <- mse@TAC %>%
        structure(
            dimnames = list(
                Sim = 1:mse@nsim,
                MPs = mse@MPs,
                Year = 2020 + 1:mse@proyears
            )
        ) %>%
        as.data.frame.table() %>%
        as_tibble() %>%
        ## { print(colnames(.)); . } %>%
        rename(Sim = 1, MPs = 2, Year = 3, Value = Freq) %>%
        mutate(OM = om) %>%
        as.data.frame
    return(output)
}) %>%
    bind_rows()

####@> New mode - SPSS...
temp01 <- map2(MSEsSPSS, names(MSEsSPSS), function(mse, om) {
    output <- mse@TAC %>%
        structure(
            dimnames = list(
                Sim = 1:mse@nsim,
                MPs = mse@MPs,
                Year = 2020 + 1:mse@proyears
            )
        ) %>%
        as.data.frame.table() %>%
        as_tibble() %>%
        ## { print(colnames(.)); . } %>%
        rename(Sim = 1, MPs = 2, Year = 3, Value = Freq) %>%
        mutate(OM = om) %>%
        as.data.frame
    return(output)
}) %>%
    bind_rows()

####@> Binding everything...
TAC_output <- bind_rows(temp00, temp01) %>%
    mutate(Year = as.numeric(as.character(Year)))

####@> Looking to the RP - Estimating the average between all MOMs...
tab05 <- TAC_output %>%
    group_by(MPs, Year) %>%
    summarise(TAC = median(Value, na.rm = TRUE),
              TAC_Lwr = quantile(Value, probs = 0.025,
                                 na.rm = TRUE),
              TAC_Upr = quantile(Value, probs = 0.975,
                                 na.rm = TRUE),
              TAC_Lwr2 = quantile(Value, probs = 0.1,
                                 na.rm = TRUE),
              TAC_Upr2 = quantile(Value, probs = 0.9,
                                  na.rm = TRUE)) %>%
    mutate(Class = case_when(
               grepl("FMSY", MPs) ~ "Reference",
               grepl("IR", MPs) ~ "Index-based: Iratio",
               grepl("CE", MPs) ~ "Index-based: Exploitation rate",
               grepl("SP_01|SP_03", MPs) ~
                   "Model-based: Surplus Production",
               grepl("SP_02|SP_04", MPs) ~
                   "Model-based: State-Space Surplus Production")) %>%
    as.data.frame()

out <- data.frame(MPs = NULL, Year = NULL, TAC = NULL, TAC_Lwr = NULL,
                  TAC_Upr = NULL, TAC_Lwr2 = NULL, TAC_Upr2 = NULL)
for(i in unique(tab05$MPs)) {
    tmp <- filter(tab05, MPs == i)
    o1 <- max_value_by_3years(tmp, "Year", "TAC")
    o2 <- max_value_by_3years(tmp, "Year", "TAC_Lwr")
    o3 <- max_value_by_3years(tmp, "Year", "TAC_Upr")
    o4 <- max_value_by_3years(tmp, "Year", "TAC_Lwr2")
    o5 <- max_value_by_3years(tmp, "Year", "TAC_Upr2")
    x <- o1 %>%
        left_join(o2, by = "year") %>%
        left_join(o3, by = "year") %>%
        left_join(o4, by = "year") %>%
        left_join(o5, by = "year") %>%
        rename("TAC" = max_value.x, "TAC_Lwr" = max_value.y,
               "TAC_Upr" = max_value.x.x, "TAC_Lwr2" = max_value.y.y,
               "TAC_Upr2" = max_value)
    out <- rbind(out, data.frame(MPs = i,
                                 Year = x$year,
                                 TAC = x$TAC,
                                 TAC_Lwr = x$TAC_Lwr,
                                 TAC_Upr = x$TAC_Upr,
                                 TAC_Lwr2 = x$TAC_Lwr2,
                                 TAC_Upr2 = x$TAC_Upr2))
}

## ###@> Figure Reference Points...
## p02a <- ggplot(data =
##                   filter(tab05,
##                          Class == "Reference",
##                          Year %in% 2021:2050)) +
##     geom_ribbon(aes(x = Year,
##                     ymin = TAC_Lwr,
##                     ymax = TAC_Upr),
##                 alpha = 0.5, fill = my_pal[3]) +
##     geom_ribbon(aes(x = Year,
##                     ymin = TAC_Lwr2,
##                     ymax = TAC_Upr2),
##                 alpha = 0.5, fill = my_pal[3]) +
##     geom_line(aes(x = Year, y = TAC), linewidth = 1,
##               colour = "black") +
##     ## geom_hline(yintercept = 1, linetype = "dashed", colour = "gray30") +
##     facet_grid(Class~MPs) +
##     ## scale_y_continuous(limits = c(0, 3), expand = c(0, 0)) +
##     coord_cartesian(xlim = c(2021, 2050)) +
##     scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
##     labs(x = "Year", y = "TAC (t)") +
##     my_theme()
## p02a

## p02b <- ggplot(data =
##                    filter(tab05,
##                           Class == "Index-based: Iratio",
##                           Year %in% 2021:2050)) +
##     geom_ribbon(aes(x = Year,
##                     ymin = TAC_Lwr,
##                     ymax = TAC_Upr),
##                 alpha = 0.5, fill = my_pal[3]) +
##     geom_ribbon(aes(x = Year,
##                     ymin = TAC_Lwr2,
##                     ymax = TAC_Upr2),
##                 alpha = 0.5, fill = my_pal[3]) +
##     geom_line(aes(x = Year, y = TAC), linewidth = 1,
##               colour = "black") +
##     geom_line(data = filter(out, MPs %in% c("IR_01", "IR_02", "IR_03")),
##               aes(x = Year, y = Value), linewidth = 1,
##               colour = "red") +
##     geom_hline(yintercept = 21377, linetype = "dashed", colour = "gray30") +
##     facet_grid(Class~MPs) +
##     scale_y_continuous(expand = c(0, 0), limits = c(0, 60000)) +
##     coord_cartesian(xlim = c(2021, 2050), ylim = c(0, 40000)) +
##     scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
##     labs(x = "Year", y = "TAC (t)") +
##     my_theme()
## p02b

## p02c <- ggplot(data =
##                    filter(tab05,
##                           Class == "Index-based: Exploitation rate",
##                           Year %in% 2021:2050)) +
##     geom_ribbon(aes(x = Year,
##                     ymin = TAC_Lwr,
##                     ymax = TAC_Upr),
##                 alpha = 0.5, fill = my_pal[3]) +
##     geom_ribbon(aes(x = Year,
##                     ymin = TAC_Lwr2,
##                     ymax = TAC_Upr2),
##                 alpha = 0.5, fill = my_pal[3]) +
##     ## geom_line(aes(x = Year, y = TAC), linewidth = 1,
##     ##           colour = "black") +
##     geom_line(data = filter(out, MPs %in% c("CE_01", "CE_02", "CE_03")),
##               aes(x = Year, y = Value), linewidth = 1,
##               colour = "red") +
##     geom_hline(yintercept = 21377, linetype = "dashed", colour = "gray30") +
##     facet_grid(Class~MPs) +
##     scale_y_continuous(expand = c(0, 0)) +
##     coord_cartesian(xlim = c(2021, 2050), ylim = c(0, 40000)) +
##     scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
##     labs(x = "Year", y = "TAC (t)") +
##     my_theme()
## p02c

## p02d <- ggplot(data =
##                    filter(tab05,
##                           Class == "Model-based: Surplus Production",
##                           Year %in% 2021:2050)) +
##     geom_ribbon(aes(x = Year,
##                     ymin = TAC_Lwr,
##                     ymax = TAC_Upr),
##                 alpha = 0.5, fill = my_pal[3]) +
##     geom_ribbon(aes(x = Year,
##                     ymin = TAC_Lwr2,
##                     ymax = TAC_Upr2),
##                 alpha = 0.5, fill = my_pal[3]) +
##     ## geom_line(aes(x = Year, y = TAC), linewidth = 1,
##     ##           colour = "black") +
##     geom_line(data = filter(out, MPs %in% c("SP_01", "SP_03")),
##               aes(x = Year, y = Value), linewidth = 1,
##               colour = "red") +
##     geom_hline(yintercept = 21377, linetype = "dashed", colour = "gray30") +
##     facet_grid(Class~MPs) +
##     scale_y_continuous(expand = c(0, 0)) +
##     coord_cartesian(xlim = c(2021, 2050), ylim = c(0, 40000)) +
##     scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
##     labs(x = "Year", y = "TAC (t)") +
##     my_theme()
## p02d

## p02e <- ggplot(data =
##                    filter(tab05,
##                           Class == "Model-based: State-Space Surplus Production",
##                           Year %in% 2021:2050)) +
##     geom_ribbon(aes(x = Year,
##                     ymin = TAC_Lwr,
##                     ymax = TAC_Upr),
##                 alpha = 0.5, fill = my_pal[3]) +
##     geom_ribbon(aes(x = Year,
##                     ymin = TAC_Lwr2,
##                     ymax = TAC_Upr2),
##                 alpha = 0.5, fill = my_pal[3]) +
##     geom_line(aes(x = Year, y = TAC), linewidth = 1,
##               colour = "black") +
##     geom_line(data = filter(out, MPs %in% c("SP_02", "SP_04")),
##               aes(x = Year, y = Value), linewidth = 1,
##               colour = "red") +
##     geom_hline(yintercept = 21377, linetype = "dashed", colour = "gray30") +
##     facet_grid(Class~MPs) +
##     scale_y_continuous(expand = c(0, 0)) +
##     coord_cartesian(xlim = c(2021, 2050), ylim = c(0, 40000)) +
##     scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
##     labs(x = "Year", y = "TACes (t)") +
##     my_theme()
## p02e

## #####@> Exporting figures...
## ggsave("05_Results/SCRS_Doc/Fig04_TS_TAC_RefsMPs_Reference_ver00.png",
##        plot = p01a, device = "png", units = "cm", w = 40, h = 25,
##        dpi = 500, bg = "white")

## ggsave("05_Results/SCRS_Doc/Fig04_TS_TAC_IRMPs_Reference_ver00.png",
##        plot = p01b, device = "png", units = "cm", w = 40, h = 25,
##        dpi = 500, bg = "white")

## ggsave("05_Results/SCRS_Doc/Fig04_TS_TAC_CEMPs_Reference_ver00.png",
##        plot = p01c, device = "png", units = "cm", w = 40, h = 25,
##        dpi = 500, bg = "white")

## ggsave("05_Results/SCRS_Doc/Fig04_TS_TAC_SPMPs_Reference_ver00.png",
##        plot = p01d, device = "png", units = "cm", w = 40, h = 25,
##        dpi = 500, bg = "white")

## ggsave("05_Results/SCRS_Doc/Fig04_TS_TAC_SPSSMPs_Reference_ver00.png",
##        plot = p01e, device = "png", units = "cm", w = 40, h = 25,
##        dpi = 500, bg = "white")

######@> Figures to the Response to the Comission...
tmp02 <- tab02 %>%
    select(MPs:SB_SBMSY_Upr2) %>%
    rename("m" = SB_SBMSY, "lwr" = SB_SBMSY_Lwr, "lwr2" = SB_SBMSY_Lwr2,
           "upr" = SB_SBMSY_Upr, "upr2" = SB_SBMSY_Upr2) %>%
    mutate(type = "SSB/SSBMSY",
           year2 = Year + 3)
tmp03 <- tab03 %>%
    select(MPs:F_FMSY_Upr2) %>%
    rename("m" = F_FMSY, "lwr" = F_FMSY_Lwr, "lwr2" = F_FMSY_Lwr2,
           "upr" = F_FMSY_Upr, "upr2" = F_FMSY_Upr2) %>%
    mutate(type = "F/FMSY",
           year2 = Year + 3)
tmp04 <- out %>%
    select(MPs:TAC_Upr2) %>%
    rename("m" = TAC, "lwr" = TAC_Lwr, "lwr2" = TAC_Lwr2,
           "upr" = TAC_Upr, "upr2" = TAC_Upr2) %>%
    mutate(type = "TAC",
           year2 = Year + 3)
TS_output <- bind_rows(tmp02, tmp03, tmp04)

periods <- bind_rows(
    data.frame(Year = unique(TS_output$Year),
               year2 = unique(TS_output$year2),
               y = 1,
               type = "F/FMSY",
               period = factor(c(rep(1, 3), rep(2, 7), rep(3, 20)))),
    data.frame(Year = unique(TS_output$Year),
               year2 = unique(TS_output$year2),
               y = 1,
               type = "SSB/SSBMSY",
               period = factor(c(rep(1, 4), rep(2, 7), rep(3, 19)))),
    data.frame(Year = unique(TS_output$Year),
               year2 = unique(TS_output$year2),
               y = 21377,
               type = "TAC",
               period = factor(c(rep(1, 4), rep(2, 7), rep(3, 19)))))

ggplot_colors <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

p00 <- ggplot(data = filter(TS_output, !MPs %in% c("FMSYref",
                                                   "FMSYref75",
                                                   "FMSYref110")),
              aes(x = Year)) +
    geom_ribbon(aes(ymin = lwr,
                    ymax = upr),
                alpha = 0.3, fill = my_pal[3]) +
    geom_ribbon(aes(ymin = lwr2,
                    ymax = upr2),
                alpha = 0.3, fill = my_pal[3]) +
    geom_line(aes(y = m), linewidth = 1,
              colour = "black") +
    geom_line(data = periods, aes(x = Year, y = y, colour = period),
              linewidth = 1) +
    ## geom_hline(yintercept = 21377, linetype = "dashed", colour = "gray30") +
    facet_grid(type~MPs, scales = "free_y") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_colour_manual(values = rev(ggplot_colors(3))) +
    ## coord_cartesian(xlim = c(2021, 2050), ylim = c(0, 40000)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 10)) +
    labs(x = "Year", y = "Median") +
    my_theme() +
    theme(legend.position = "none")
p00

ggsave("05_Results/Response_Comm/Fig04_TS_ALL_Together_ver01.png",
       plot = p00, device = "png", units = "cm", w = 50, h = 30,
       dpi = 500, bg = "white")

p01 <- ggplot(data = TS_output,
              aes(x = year2)) +
    geom_ribbon(aes(ymin = lwr,
                    ymax = upr),
                alpha = 0.3, fill = my_pal[3]) +
    geom_ribbon(aes(ymin = lwr2,
                    ymax = upr2),
                alpha = 0.3, fill = my_pal[3]) +
    geom_line(aes(y = m), linewidth = 1,
              colour = "black") +
    geom_line(data = periods, aes(x = year2, y = y, colour = period),
              linewidth = 1) +
    ## geom_hline(yintercept = 21377, linetype = "dashed", colour = "gray30") +
    facet_grid(type~MPs, scales = "free_y") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_colour_manual(values = rev(ggplot_colors(3))) +
    ## coord_cartesian(xlim = c(2021, 2050), ylim = c(0, 40000)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 10)) +
    labs(x = "Year", y = "Median") +
    my_theme() +
    theme(legend.position = "none")
p01

ggsave("05_Results/Response_Comm/Fig04_TS_ALL_Together_ver00.png",
       plot = p00, device = "png", units = "cm", w = 50, h = 30,
       dpi = 500, bg = "white")

######@> Creating the Kobe plots (by Year and general)...

#####@> Merging SB/SBMSY and F/FMSY...
df00 <- SB_SBMSY_output %>%
    left_join(F_FMSY_output, by = c("Sim", "MPs", "Year", "OM")) %>%
    rename("SB_SBMSY" = Value.x, "F_FMSY" = Value.y) %>%
    mutate(Class = case_when(
               grepl("FMSY", MPs) ~ "Reference",
               grepl("IR", MPs) ~ "Index-based: Iratio",
               grepl("CE", MPs) ~ "Index-based: Exploitation rate",
               grepl("SP_01|SP_03", MPs) ~
                   "Model-based: Surplus Production",
               grepl("SP_02|SP_04", MPs) ~
                   "Model-based: State-Space Surplus Production"),
           Class2 = case_when(
               grepl("FMSY", MPs) ~ "Reference",
               grepl("IR", MPs) ~ "Index-based: IR",
               grepl("CE", MPs) ~ "Index-based: CE",
               grepl("SP_01|SP_03", MPs) ~
                   "Model-based: SP",
               grepl("SP_02|SP_04", MPs) ~
                   "Model-based: SPSS"))

#####@> Estimating Kobe values for the last year...
tmp <- df00 %>%
    filter(Year == 2050) %>%
    group_by(Sim, MPs, Class) %>%
    summarise(SB_SBMSY = mean(SB_SBMSY, na.rm = TRUE),
              F_FMSY = mean(F_FMSY, na.rm = TRUE)) %>%
    mutate(OF = ifelse(F_FMSY < 1, FALSE, TRUE),
           OFD = ifelse(SB_SBMSY < 1, TRUE, FALSE))

#####@> Percentage of cases...
valdf <- tmp %>%
    group_by(MPs, Class) %>%
    summarise(BL = sum(OF == FALSE & OFD == TRUE)/300 * 100,
              BR = sum(OF == FALSE & OFD == FALSE)/300 * 100,
              TL = sum(OF == TRUE & OFD == TRUE)/300 * 100,
              TR = sum(OF == TRUE & OFD == FALSE)/300 * 100)
valdf <- valdf %>% tidyr::pivot_longer(., cols = 3:6)
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
                     data.frame(x = c(1, 1, 3, 3),
                                y = c(0, 1, 1, 0),
                                fill = "br"),
                     data.frame(x = c(0, 0, 1, 1),
                                y = c(1, 3, 3, 1),
                                fill = "tl"),
                     data.frame(x = c(1, 1, 3, 3),
                                y = c(1, 3, 3, 1),
                                fill = "tr"))
kobe_df$alpha <- 0.2

#####@> Figure dots...
p05 <- ggplot() +
    geom_polygon(data = kobe_df,
                 aes(x = x, y = y, fill = fill, alpha = alpha)) +
    scale_fill_manual(values = c("#F8DC7A", "#67C18B", "#D8775D",
                                 "#FDBD56")) +
    geom_point(data = filter(tmp, Class == "Reference"),
               aes(x = SB_SBMSY, y = F_FMSY), alpha = 0.5,
               size = 4) +
    geom_text(data = filter(valdf, Class == "Reference"),
              fontface = "bold", size = 4,
              aes(x = x, y = y, label = value,
                  hjust = hjustvar, vjust = vjustvar),
              colour = "black") +
    expand_limits(x = c(0, 3), y = c(0, 3)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 3)) +
    geom_hline(yintercept = 1, color = "darkgray", linetype = 2) +
    geom_vline(xintercept = 1, color = "darkgray", linetype = 2) +
    labs(x = expression(SSB/SSB[MSY]), y = expression(F/F[MSY])) +
    facet_grid(Class~MPs) +
    my_theme() +
    theme(legend.position = "none")
p05

p05b <- ggplot() +
    geom_polygon(data = kobe_df,
                 aes(x = x, y = y, fill = fill, alpha = alpha)) +
    scale_fill_manual(values = c("#F8DC7A", "#67C18B", "#D8775D",
                                 "#FDBD56")) +
    geom_point(data = filter(tmp, Class == "Index-based: Iratio"),
               aes(x = SB_SBMSY, y = F_FMSY), alpha = 0.5,
               size = 4) +
    geom_text(data = filter(valdf, Class == "Index-based: Iratio"),
              fontface = "bold", size = 4,
              aes(x = x, y = y, label = value,
                  hjust = hjustvar, vjust = vjustvar),
              colour = "black") +
    expand_limits(x = c(0, 3), y = c(0, 3)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 3)) +
    geom_hline(yintercept = 1, color = "darkgray", linetype = 2) +
    geom_vline(xintercept = 1, color = "darkgray", linetype = 2) +
    labs(x = expression(SSB/SSB[MSY]), y = expression(F/F[MSY])) +
    facet_grid(Class~MPs) +
    my_theme() +
    theme(legend.position = "none")
p05b

p05c <- ggplot() +
    geom_polygon(data = kobe_df,
                 aes(x = x, y = y, fill = fill, alpha = alpha)) +
    scale_fill_manual(values = c("#F8DC7A", "#67C18B", "#D8775D",
                                 "#FDBD56")) +
    geom_point(data = filter(tmp, Class == "Index-based: Exploitation rate"),
               aes(x = SB_SBMSY, y = F_FMSY), alpha = 0.5,
               size = 4) +
    geom_text(data = filter(valdf, Class == "Index-based: Exploitation rate"),
              fontface = "bold", size = 4,
              aes(x = x, y = y, label = value,
                  hjust = hjustvar, vjust = vjustvar),
              colour = "black") +
    expand_limits(x = c(0, 3), y = c(0, 3)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 3)) +
    geom_hline(yintercept = 1, color = "darkgray", linetype = 2) +
    geom_vline(xintercept = 1, color = "darkgray", linetype = 2) +
    labs(x = expression(SSB/SSB[MSY]), y = expression(F/F[MSY])) +
    facet_grid(Class~MPs) +
    my_theme() +
    theme(legend.position = "none")
p05c

p05d <- ggplot() +
    geom_polygon(data = kobe_df,
                 aes(x = x, y = y, fill = fill, alpha = alpha)) +
    scale_fill_manual(values = c("#F8DC7A", "#67C18B", "#D8775D",
                                 "#FDBD56")) +
    geom_point(data = filter(tmp,
                             Class == "Model-based: Surplus Production"),
               aes(x = SB_SBMSY, y = F_FMSY), alpha = 0.5,
               size = 4) +
    geom_text(data = filter(valdf,
                            Class == "Model-based: Surplus Production"),
              fontface = "bold", size = 4,
              aes(x = x, y = y, label = value,
                  hjust = hjustvar, vjust = vjustvar),
              colour = "black") +
    expand_limits(x = c(0, 3), y = c(0, 3)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 3)) +
    geom_hline(yintercept = 1, color = "darkgray", linetype = 2) +
    geom_vline(xintercept = 1, color = "darkgray", linetype = 2) +
    labs(x = expression(SSB/SSB[MSY]), y = expression(F/F[MSY])) +
    facet_grid(Class~MPs) +
    my_theme() +
    theme(legend.position = "none")
p05d

p05e <- ggplot() +
    geom_polygon(data = kobe_df,
                 aes(x = x, y = y, fill = fill, alpha = alpha)) +
    scale_fill_manual(values = c("#F8DC7A", "#67C18B", "#D8775D",
                                 "#FDBD56")) +
    geom_point(data =
                   filter(tmp,
                          Class == "Model-based: State-Space Surplus Production"),
               aes(x = SB_SBMSY, y = F_FMSY), alpha = 0.5,
               size = 4) +
    geom_text(data = filter(valdf,
                            Class == "Model-based: State-Space Surplus Production"),
              fontface = "bold", size = 4,
              aes(x = x, y = y, label = value,
                  hjust = hjustvar, vjust = vjustvar),
              colour = "black") +
    expand_limits(x = c(0, 3), y = c(0, 3)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 3)) +
    geom_hline(yintercept = 1, color = "darkgray", linetype = 2) +
    geom_vline(xintercept = 1, color = "darkgray", linetype = 2) +
    labs(x = expression(SSB/SSB[MSY]), y = expression(F/F[MSY])) +
    facet_grid(Class~MPs) +
    my_theme() +
    theme(legend.position = "none")
p05e

######@> Estimating Kobe values per year...
tmp <- df00 %>%
    ## group_by(Sim, Year, Class2, MPs) %>%
    ## summarise(SB_SBMSY = median(SB_SBMSY, na.rm = TRUE),
    ##           F_FMSY = median(F_FMSY, na.rm = TRUE)) %>%
    mutate(OF = ifelse(F_FMSY < 1, FALSE, TRUE),
           OFD = ifelse(SB_SBMSY < 1, TRUE, FALSE))

#####@> Proportions by year...
valdf <- tmp %>%
    group_by(Year, Class2, MPs) %>%
    summarise(Yellow = sum(OF == FALSE & OFD == TRUE)/2700 * 100,
              Green = sum(OF == FALSE & OFD == FALSE)/2700 * 100,
              Red = sum(OF == TRUE & OFD == TRUE)/2700 * 100,
              Orange = sum(OF == TRUE & OFD == FALSE)/2700 * 100) %>%
    pivot_longer(names_to = "Cond", values_to = "Perc", 4:7) %>%
    mutate(Cond = factor(Cond, levels = c("Green", "Yellow", "Orange",
                                          "Red")))

## p06 <- ggplot() +
##     geom_bar(data = filter(valdf, Year %in% 2021:2050),
##              aes(x = Year, y = Perc, fill = Cond),
##              stat = "identity", colour = "black",
##              width = 1) +
##     geom_hline(yintercept = 30, lty = 2, alpha = 0.5) +
##     facet_wrap(~MPs) +
##     scale_fill_manual(values = c("#67C18B", "#F8DC7A", "#FDBD56",
##                                  "#D8775D")) +
##     expand_limits(y = c(0, 100)) +
##     scale_y_continuous(expand = c(0, 0)) +
##     scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
##     labs(x = "Year", y = "%") +
##     my_theme() +
##     theme(legend.position = "none")
## p06

p07 <- ggplot() +
    geom_area(data = filter(valdf, Year %in% 2021:2050,
                            Class2 == "Reference"),
              aes(x = Year, y = Perc, fill = rev(Cond)),
              stat = "identity", colour = "black") +
    geom_hline(yintercept = 70, lty = 2, alpha = 0.5) +
    facet_grid(Class2~MPs) +
    scale_fill_manual(values = rev(c("#67C18B", "#F8DC7A", "#FDBD56",
                                     "#D8775D"))) +
    labs(x = "Year", y = "%") +
    expand_limits(y = c(0, 100)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 10)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
    annotate("rect", xmin = 2023, xmax = 2028, ymin = 0, ymax = 100,
             fill = "white", alpha = 0.5) +
    annotate("text", x = 2025.4, y = 95, label = "Kobe Matrix Period",
             alpha = 0.6, size = 2) +
    my_theme() +
    theme(legend.position = "none")
p07

p07b <- ggplot() +
    geom_area(data = filter(valdf, Year %in% 2021:2050,
                            Class2 == "Index-based: IR"),
              aes(x = Year, y = Perc, fill = rev(Cond)),
              stat = "identity", colour = "black") +
    geom_hline(yintercept = 70, lty = 2, alpha = 0.5) +
    facet_grid(Class2~MPs) +
    scale_fill_manual(values = rev(c("#67C18B", "#F8DC7A", "#FDBD56",
                                     "#D8775D"))) +
    labs(x = "Year", y = "%") +
    expand_limits(y = c(0, 100)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 10)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
    annotate("rect", xmin = 2023, xmax = 2028, ymin = 0, ymax = 100,
             fill = "white", alpha = 0.5) +
    annotate("text", x = 2025.4, y = 95, label = "Kobe Matrix Period",
             alpha = 0.6, size = 2) +
    my_theme() +
    theme(legend.position = "none")
p07b

p07c <- ggplot() +
    geom_area(data = filter(valdf, Year %in% 2021:2050,
                            Class2 == "Index-based: CE"),
              aes(x = Year, y = Perc, fill = rev(Cond)),
              stat = "identity", colour = "black") +
    geom_hline(yintercept = 70, lty = 2, alpha = 0.5) +
    facet_grid(Class2~MPs) +
    scale_fill_manual(values = rev(c("#67C18B", "#F8DC7A", "#FDBD56",
                                     "#D8775D"))) +
    labs(x = "Year", y = "%") +
    expand_limits(y = c(0, 100)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 10)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
    annotate("rect", xmin = 2023, xmax = 2028, ymin = 0, ymax = 100,
             fill = "white", alpha = 0.5) +
    annotate("text", x = 2025.4, y = 95, label = "Kobe Matrix Period",
             alpha = 0.6, size = 2) +
    my_theme() +
    theme(legend.position = "none")
p07c

p07d <- ggplot() +
    geom_area(data = filter(valdf, Year %in% 2021:2050,
                            Class2 == "Model-based: SP"),
              aes(x = Year, y = Perc, fill = rev(Cond)),
              stat = "identity", colour = "black") +
    geom_hline(yintercept = 70, lty = 2, alpha = 0.5) +
    facet_grid(Class2~MPs) +
    scale_fill_manual(values = rev(c("#67C18B", "#F8DC7A", "#FDBD56",
                                     "#D8775D"))) +
    labs(x = "Year", y = "%") +
    expand_limits(y = c(0, 100)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 10)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
    annotate("rect", xmin = 2023, xmax = 2028, ymin = 0, ymax = 100,
             fill = "white", alpha = 0.5) +
    annotate("text", x = 2025.4, y = 95, label = "Kobe Matrix Period",
             alpha = 0.6, size = 2) +
    my_theme() +
    theme(legend.position = "none")
p07d

p07e <- ggplot() +
    geom_area(data = filter(valdf, Year %in% 2021:2050,
                            Class2 == "Model-based: SPSS"),
              aes(x = Year, y = Perc, fill = rev(Cond)),
              stat = "identity", colour = "black") +
    geom_hline(yintercept = 70, lty = 2, alpha = 0.5) +
    facet_grid(Class2~MPs) +
    scale_fill_manual(values = rev(c("#67C18B", "#F8DC7A", "#FDBD56",
                                     "#D8775D"))) +
    labs(x = "Year", y = "%") +
    expand_limits(y = c(0, 100)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 10)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
    annotate("rect", xmin = 2023, xmax = 2028, ymin = 0, ymax = 100,
             fill = "white", alpha = 0.5) +
    annotate("text", x = 2025.4, y = 95, label = "Kobe Matrix Period",
             alpha = 0.6, size = 2) +
    my_theme() +
    theme(legend.position = "none")
p07e

p07f <- ggplot() +
    geom_area(data = filter(valdf, Year %in% 2021:2050),
              aes(x = Year, y = Perc, fill = rev(Cond)),
              stat = "identity", colour = "black") +
    geom_hline(yintercept = 70, lty = 2, alpha = 0.5) +
    facet_grid(~MPs) +
    scale_fill_manual(values = rev(c("#67C18B", "#F8DC7A", "#FDBD56",
                                     "#D8775D"))) +
    labs(x = "Year", y = "%") +
    expand_limits(y = c(0, 100)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 10)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(2020, 2050, 5)) +
    annotate("rect", xmin = 2023, xmax = 2028, ymin = 0, ymax = 100,
             fill = "white", alpha = 0.5) +
    annotate("text", x = 2025.4, y = 95, label = "Kobe Matrix Period",
             alpha = 0.6, size = 2) +
    my_theme() +
    theme(legend.position = "none")
p07f

plot07 <- p07/p07b/p07c/p07d/p07e

plot07b <- p07b/p07c/p07d/p07e

#####@> Exporting figures...
ggsave("05_Results/SCRS_Doc/Fig07_TS_Kobe_RefsMPs_Reference_ver00.png",
       plot = p07, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig07_TS_Kobe_IRMPs_Reference_ver00.png",
       plot = p07b, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig07_TS_Kobe_CEMPs_Reference_ver00.png",
       plot = p07c, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig07_TS_Kobe_SPMPs_Reference_ver00.png",
       plot = p07d, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig07_TS_Kobe_SPSSMPs_Reference_ver00.png",
       plot = p07e, device = "png", units = "cm", w = 40, h = 25,
       dpi = 500, bg = "white")

ggsave("05_Results/SCRS_Doc/Fig07_TS_Kobe_ALL_ver00.png",
       plot = plot07, device = "png", units = "cm", w = 50, h = 50,
       dpi = 300, bg = "white")

ggsave("05_Results/Response_Comm/Fig07_TS_Kobe_ALL_B_ver00.png",
       plot = plot07b, device = "png", units = "cm", w = 50, h = 40,
       dpi = 300, bg = "white")

######@> Boxplot...
temp <- PM_output %>%
    filter(PM %in% c("VarC")) %>%
    filter(MP %in% c("IR_01", "IR_02",
                     "IR_03", "CE_01", "CE_02",
                     "CE_03", "SP_01", "SP_02",
                     "SP_03", "SP_04")) %>%
    mutate(MP = factor(MP, levels = c("IR_01", "IR_02",
                                      "IR_03", "CE_01", "CE_02",
                                      "CE_03", "SP_01", "SP_02",
                                      "SP_03", "SP_04"))) %>%
    mutate(Class = case_when(
               grepl("IR", MP) ~ "Index-based: Iratio",
               grepl("CE", MP) ~ "Index-based: Exploitation rate",
               grepl("SP_01|SP_03", MP) ~
                   "Model-based: Surplus Production",
               grepl("SP_02|SP_04", MP) ~
                   "Model-based: State-Space Surplus Production"))

p01 <- ggplot(data = temp,
              aes(x = MP, y = Values, fill = MP,
                  colour = MP)) +
    geom_jitter(pch = 21, size = 4, alpha = 0.1) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.8,
                scale = "width", colour = "black") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.2,
                 colour = "darkred") +
    ## facet_wrap(~MP, ncol = 3) +
    geom_hline(yintercept = 0.2, linetype = "dashed", colour = "gray10") +
    ## scale_fill_manual(values = my_pal[c(3, 5, 8, 13, 15)]) +
    ## scale_colour_manual(values = my_pal[c(3, 5, 8, 13, 15)]) +
    coord_cartesian(ylim = c(0, 1)) +
    my_theme() +
    labs(x = "Candidate Management Procedures",
         y = "Absolute change in TAC (%)") +
    theme(legend.position = "none")
p01

ggsave("05_Results/Response_Comm/Absolute_Change_TAC_ver01.png",
       plot = p01,
       device = "png", dpi = 600, bg = "white", unit = "cm",
       w = 40, h = 35)

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
