#V3.30.18.00;_safe;_compile_date:_Sep 30 2021;_Stock_Synthesis_by_Richard_Methot_(NOAA)_using_ADMB_12.3
#_Stock_Synthesis_is_a_work_of_the_U.S._Government_and_is_not_subject_to_copyright_protection_in_the_United_States.
#_Foreign_copyrights_may_apply._See_copyright.txt_for_more_information.
#_User_support_available_at:NMFS.Stock.Synthesis@noaa.gov
#_User_info_available_at:https://vlab.noaa.gov/group/stock-synthesis
#_Source_code_at:_https://github.com/nmfs-stock-synthesis/stock-synthesis

#C file created using the SS_writectl function in the R package r4ss
#C file write time: 2020-03-30 14:33:05
#_data_and_control_files: data.ss_new // control_modified.ss
0  # 0 means do not read wtatage.ss; 1 means read and use wtatage.ss and also read and use growth parameters
1  #_N_Growth_Patterns (Growth Patterns, Morphs, Bio Patterns, GP are terms used interchangeably in SS)
1 #_N_platoons_Within_GrowthPattern 
#_Cond 1 #_Platoon_within/between_stdev_ratio (no read if N_platoons=1)
#_Cond  1 #vector_platoon_dist_(-1_in_first_val_gives_normal_approx)
#
4 # recr_dist_method for parameters:  2=main effects for GP, Area, Settle timing; 3=each Settle entity; 4=none (only when N_GP*Nsettle*pop==1)
1 # not yet implemented; Future usage: Spawner-Recruitment: 1=global; 2=by area
1 #  number of recruitment settlement assignments 
0 # unused option
#GPattern month  area  age (for each settlement assignment)
 1 1 1 0
#
#_Cond 0 # N_movement_definitions goes here if Nareas > 1
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) also cond on do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, morph=1, source=1 dest=2, age1=4, age2=10
#
1 #_Nblock_Patterns
 1 #_blocks_per_pattern 
# begin and end years of blocks
 2015 2020
#
# controls for all timevary parameters 
1 #_time-vary parm bound check (1=warn relative to base parm bounds; 3=no bound check); Also see env (3) and dev (5) options to constrain with base bounds
#
# AUTOGEN
 1 1 1 1 1 # autogen: 1st element for biology, 2nd for SR, 3rd for Q, 4th reserved, 5th for selex
# where: 0 = autogen time-varying parms of this category; 1 = read each time-varying parm line; 2 = read then autogen if parm min==-12345
#
#_Available timevary codes
#_Block types: 0: P_block=P_base*exp(TVP); 1: P_block=P_base+TVP; 2: P_block=TVP; 3: P_block=P_block(-1) + TVP
#_Block_trends: -1: trend bounded by base parm min-max and parms in transformed units (beware); -2: endtrend and infl_year direct values; -3: end and infl as fraction of base range
#_EnvLinks:  1: P(y)=P_base*exp(TVP*env(y));  2: P(y)=P_base+TVP*env(y);  3: P(y)=f(TVP,env_Zscore) w/ logit to stay in min-max;  4: P(y)=2.0/(1.0+exp(-TVP1*env(y) - TVP2))
#_DevLinks:  1: P(y)*=exp(dev(y)*dev_se;  2: P(y)+=dev(y)*dev_se;  3: random walk;  4: zero-reverting random walk with rho;  5: like 4 with logit transform to stay in base min-max
#_DevLinks(more):  21-25 keep last dev for rest of years
#
#_Prior_codes:  0=none; 6=normal; 1=symmetric beta; 2=CASAL's beta; 3=lognormal; 4=lognormal with biascorr; 5=gamma
#
# setup for M, growth, wt-len, maturity, fecundity, (hermaphro), recr_distr, cohort_grow, (movement), (age error), (catch_mult), sex ratio 
#_NATMORT
2 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate;_5=BETA:_Maunder_link_to_maturity
6 #_reference age for Lorenzen M; read 1P per morph
#
1 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr; 5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
1 #_Age(post-settlement)_for_L1;linear growth below this
999 #_Growth_Age_for_L2 (999 to use as Linf)
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0  #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
#
1 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
1 #_First_Mature_Age
2 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_parameter_offset_approach for M, G, CV_G:  1- direct, no offset**; 2- male=fem_parm*exp(male_parm); 3: male=female*exp(parm) then old=young*exp(parm)
#_** in option 1, any male parameter with value = 0.0 and phase <0 is set equal to female parameter
#
#_growth_parms
#_ LO HI INIT PRIOR PR_SD PR_type PHASE env_var&link dev_link dev_minyr dev_maxyr dev_PH Block Block_Fxn
# Sex: 1  BioPattern: 1  NatMort
 0.001 2 0.5 0 0 0 -3 0 0 0 0 0.5 0 0 # NatM_Lorenzen_Fem_GP_1
# Sex: 1  BioPattern: 1  Growth
 10 50 38 0 0 0 -3 0 0 0 0 0.5 0 0 # L_at_Amin_Fem_GP_1
 40 120 91.001 0 0 0 -3 0 0 0 0 0.5 0 0 # L_at_Amax_Fem_GP_1
 0.1 1.2 0.53 0 0 0 -2 0 0 0 0 0.5 0 0 # VonBert_K_Fem_GP_1
 0.1 1 0.2 0 0 0 -2 0 0 0 0 0.5 0 0 # CV_young_Fem_GP_1
 0.1 1 0.2 0 0 0 -3 0 0 0 0 0.5 0 0 # CV_old_Fem_GP_1
# Sex: 1  BioPattern: 1  WtLen
 0 3 7.48e-06 7.48e-06 99 0 -99 0 0 0 0 0 0 0 # Wtlen_1_Fem_GP_1
 2 4 3.253 3.253 99 0 -99 0 0 0 0 0 0 0 # Wtlen_2_Fem_GP_1
# Sex: 1  BioPattern: 1  Maturity&Fecundity
 0.0001 1000 42 42 99 0 -99 0 0 0 0 0 0 0 # Mat50%_Fem_GP_1
 -2 4 -0.226495 -0.226495 99 0 -99 0 0 0 0 0 0 0 # Mat_slope_Fem_GP_1
 0 3 7.48e-06 7.48e-06 0.8 0 -3 0 0 0 0 0 0 0 # Eggs_scalar_Fem_GP_1
 0 4 3.253 3.253 0.8 0 -3 0 0 0 0 0 0 0 # Eggs_exp_len_Fem_GP_1
# Hermaphroditism
#  Recruitment Distribution  
#  Cohort growth dev base
 0.1 10 1 1 1 0 -1 0 0 0 0 0 0 0 # CohortGrowDev
#  Movement
#  Age Error from parameters
#  catch multiplier
#  fraction female, by GP
 0.01 0.99 0.5 0 0 0 -99 0 0 0 0 0 0 0 # FracFemale_GP_1
#  M2 parameter for each predator fleet
#
#_no timevary MG parameters
#
#_seasonal_effects_on_biology_parms
 0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_ LO HI INIT PRIOR PR_SD PR_type PHASE
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
3 #_Spawner-Recruitment; Options: 1=NA; 2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=B-H_flattop; 7=survival_3Parm; 8=Shepherd_3Parm; 9=RickerPower_3parm
0  # 0/1 to use steepness in initial equ recruitment calculation
0  #  future feature:  0/1 to make realized sigmaR a function of SR curvature
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn #  parm_name
        0.0001            20       11.3715             7             0             0          1          0          0          0          0          0          0          0 # SR_LN(R0)
           0.2             1           0.7           0.7             0             0         -3          0          0          0          0          0          0          0 # SR_BH_steep
             0             2           0.4           0.4             0             0         -5          0          0          0          0          0          0          0 # SR_sigmaR
            -5             5             0             0             0             0         -6          0          0          0          0          0          0          0 # SR_regime
             0             2             0             1             0             0         -6          0          0          0          0          0          0          0 # SR_autocorr
#_no timevary SR parameters
1 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
1993 # first year of main recr_devs; early devs can preceed this era
2018 # last year of main recr_devs; forecast devs start in following year
3 #_recdev phase 
1 # (0/1) to read 13 advanced options
 0 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
 -5 #_recdev_early_phase
 -1 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
 1 #_lambda for Fcast_recr_like occurring before endyr+1
 1964 #_last_yr_nobias_adj_in_MPD; begin of ramp
 2002.8 #_first_yr_fullbias_adj_in_MPD; begin of plateau
 2018 #_last_yr_fullbias_adj_in_MPD
 2029.1 #_end_yr_for_ramp_in_MPD (can be in forecast to shape ramp, but SS sets bias_adj to 0.0 for fcast yrs)
 0.8109 #_max_bias_adj_in_MPD (typical ~0.8; -3 sets all years to 0.0; -2 sets all non-forecast yrs w/ estimated recdevs to 1.0; -1 sets biasadj=1.0 for all yrs w/ recdevs)
 0 #_period of cycles in recruitment (N parms read below)
 -5 #min rec_dev
 5 #max rec_dev
 0 #_read_recdevs
#_end of advanced SR options
#
#_placeholder for full parameter lines for recruitment cycles
# read specified recr devs
#_Yr Input_value
#
# all recruitment deviations
#  1993R 1994R 1995R 1996R 1997R 1998R 1999R 2000R 2001R 2002R 2003R 2004R 2005R 2006R 2007R 2008R 2009R 2010R 2011R 2012R 2013R 2014R 2015R 2016R 2017R 2018R 2019F 2020F 2021F
#  0.0915541 -0.387446 0.598208 0.0545103 0.41442 -0.293786 0.245295 0.0854021 0.0103837 -0.218039 0.147237 0.18182 0.349578 -0.0605376 0.011875 0.0824851 -0.266479 0.403661 0.136021 0.0660721 0.0208579 0.596466 -0.645563 -0.834614 -0.48459 -0.304792 0 0 0
#
#Fishing Mortality info 
0.5 # F ballpark value in units of annual_F
-1999 # F ballpark year (neg value to disable)
3 # F_Method:  1=Pope midseason rate; 2=F as parameter; 3=F as hybrid; 4=fleet-specific parm/hybrid (#4 is superset of #2 and #3 and is recommended)
4 # max F (methods 2-4) or harvest fraction (method 1)
4  # N iterations for tuning in hybrid mode; recommend 3 (faster) to 5 (more precise if many fleets)
#
#_initial_F_parms; for each fleet x season that has init_catch; nest season in fleet; count = 5
#_for unconstrained init_F, use an arbitrary initial catch and set lambda=0 for its logL
#_ LO HI INIT PRIOR PR_SD  PR_type  PHASE
 0 1000 1e-20 1 999 0 -1 # InitF_seas_1_flt_1PS_West
 0 1000 1e-20 1 999 0 -1 # InitF_seas_1_flt_2BB_West
 0 1000 1e-20 1 999 0 -1 # InitF_seas_1_flt_3LL_USMX
 0 1000 1e-20 1 999 0 -1 # InitF_seas_1_flt_4LL_OTH
 0 1000 1e-20 1 999 0 -1 # InitF_seas_1_flt_5HL_RR
#
# F rates by fleet x season
# Yr:  1952 1953 1954 1955 1956 1957 1958 1959 1960 1961 1962 1963 1964 1965 1966 1967 1968 1969 1970 1971 1972 1973 1974 1975 1976 1977 1978 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021
# seas:  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
# PS_West 0 0 0 0 0 0 0 0 0 0 0.0100923 0.0657856 0.0885312 0.0014113 0.000874304 0.00069957 0.00295268 0.00222341 0 0 0.0053139 0.000629091 0.000610051 0.0042869 0.0153806 0.00734325 0.0380837 0.0164412 0.0668199 0.116915 0.275659 0.305366 0.357366 0.342992 0.191874 0.189326 0.0873154 0.0939209 0.124082 0.282011 0.308407 0.559657 0.250098 0.126751 0.0946618 0.167053 0.109945 0.133394 0.104845 0.20249 0.0827516 0.111211 0.108732 0.0728885 0.0622333 0.0467415 0.0341641 0.074562 0.0947314 0.0579391 0.0575387 0.0371331 0.0460902 0.0434375 0.0508098 0.0903855 0.0610413 0.0436201 0.0389365 0.0389365
# BB_West 0.0207865 0.021774 0.0233511 0.0238395 0.0257091 0.0335406 0.0283473 0.0314664 0.0564769 0.0573993 0.0269435 0.0169486 0.0190157 0.0261 0.0287557 0.046234 0.041964 0.02876 0.0381735 0.0294466 0.0241981 0.0332104 0.0516145 0.049475 0.0505269 0.0454001 0.0434527 0.0752254 0.172769 0.362511 0.513744 0.506365 0.446859 0.854189 0.819381 0.582302 0.642256 0.719264 0.689526 0.788997 0.714798 0.726834 0.817257 0.866586 0.559358 0.808914 0.581137 0.838395 0.721783 0.769648 0.591041 0.845306 0.775357 0.780614 0.587356 0.725214 0.63259 0.694414 0.881537 0.760575 0.88918 1.05963 0.857375 0.351687 0.5724 0.857608 0.917158 0.872527 0.599617 0.599617
# LL_USMX 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7.13548e-05 6.80902e-05 0.0001765 0.00017574 0.000387701 5.70341e-05 3.32726e-05 0.0001129 9.18824e-06 1.44287e-05 0.000112151 6.65805e-05 0.00135663 0.000374446 0.000614653 0.000189003 0.000188628 0.000128633 0.000208618 0.000295833 0.000117004 0.000135869 0.000152969 0.000118713 0.000540017 0.000132589 6.72106e-05 0.00018397 0.000600032 0.000236195 0.000641046 0.00353239 0.000969153 0.00215888 0.00205893 0.00259116 0.000508251 0.000417474 0.000198334 0.00978371 0.00374084 0.000526871 0.00736635 0.00540038 0.00175765 0.00745477 0.00329271 0.00497308 0.00767161 0.00306512 0.00306512
# LL_OTH 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.000607949 0.000622184 0.000615552 0.000608432 0.00166023 0.00163127 0.0016613 0.00341503 0.00231692 0.00154795 0.00106914 0.000854052 0.00127589 0.00108524 0.000949079 0.00161878 0.00308286 0.00510357 0.00853285 0.00602169 0.0149464 0.00510524 0.00382435 0.00664722 0.00750941 0.0124228 0.00939442 0.00931527 0.0319379 0.0504307 0.00965205 0.00444305 0.00757241 0.00551172 0.00782454 0.0161063 0.00641991 0.00569954 0.006868 0.00642119 0.00609371 0.0034697 0.00509458 0.00132085 0.00371938 0.0044645 0.00305553 0.00473949 0.00609408 0.00370166 0.002931 0.00425106 0.00529046 0.00592185 0.00206502 0.00206502
# HL_RR 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.38547e-05 0 0.000228223 0.000760346 0.0035609 1.39017e-05 0.00248434 0.000348907 0.00194831 0.000688739 0.00134473 0.003354 0.00224987 0.00115777 0.000720477 0.00173428 0.00260674 0.00159437 0.00251037 0.00211014 0.00287217 0.00161327 0.00178352 0.0018692 0.00370578 0.000902565 0.00148017 0.00197808 0.00214115 0.00248375 0.00066795 0.00114977 0.0015641 0.00147022 0.00268687 0.00268877 0.00963481 0.00737708 0.0133313 0.014244 0.00884092 0.0310656 0.209799 0.204605 0.0957622 0.0750377 0.0750377
#
#_Q_setup for fleets with cpue or survey data
#_1:  fleet number
#_2:  link type: (1=simple q, 1 parm; 2=mirror simple q, 1 mirrored parm; 3=q and power, 2 parm; 4=mirror with offset, 2 parm)
#_3:  extra input for link, i.e. mirror fleet# or dev index number
#_4:  0/1 to select extra sd parameter
#_5:  0/1 for biasadj or not
#_6:  0/1 to float
#_   fleet      link link_info  extra_se   biasadj     float  #  fleetname
         1         1         0         0         0         1  #  PS_West
         2         1         0         0         0         1  #  BB_West
         3         1         0         0         0         1  #  LL_USMX
         5         1         0         0         0         1  #  HL_RR
         6         1         0         0         0         1  #  BRA_BB_hist
-9999 0 0 0 0 0
#
#_Q_parms(if_any);Qunits_are_ln(q)
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
           -15            15      -10.2145             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_PS_West(1)
           -15            15        -10.31             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_BB_West(2)
           -15            15      -9.50597             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_LL_USMX(3)
           -15            15      -11.5462             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_HL_RR(5)
           -15            15      -10.2781             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_BRA_BB_hist(6)
#_no timevary Q parameters
#
#_size_selex_patterns
#Pattern:_0;  parm=0; selex=1.0 for all sizes
#Pattern:_1;  parm=2; logistic; with 95% width specification
#Pattern:_2;  parm=6; modification of pattern 24 with improved sex-specific offset
#Pattern:_5;  parm=2; mirror another size selex; PARMS pick the min-max bin to mirror
#Pattern:_11; parm=2; selex=1.0  for specified min-max population length bin range
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_6;  parm=2+special; non-parm len selex
#Pattern:_43; parm=2+special+2;  like 6, with 2 additional param for scaling (average over bin range)
#Pattern:_8;  parm=8; double_logistic with smooth transitions and constant above Linf option
#Pattern:_9;  parm=6; simple 4-parm double logistic with starting length; parm 5 is first length; parm 6=1 does desc as offset
#Pattern:_21; parm=2+special; non-parm len selex, read as pairs of size, then selex
#Pattern:_22; parm=4; double_normal as in CASAL
#Pattern:_23; parm=6; double_normal where final value is directly equal to sp(6) so can be >1.0
#Pattern:_24; parm=6; double_normal with sel(minL) and sel(maxL), using joiners
#Pattern:_25; parm=3; exponential-logistic in length
#Pattern:_27; parm=special+3; cubic spline in length; parm1==1 resets knots; parm1==2 resets all 
#Pattern:_42; parm=special+3+2; cubic spline; like 27, with 2 additional param for scaling (average over bin range)
#_discard_options:_0=none;_1=define_retention;_2=retention&mortality;_3=all_discarded_dead;_4=define_dome-shaped_retention
#_Pattern Discard Male Special
 24 0 0 0 # 1 PS_West
 24 0 0 0 # 2 BB_West
 1 0 0 0 # 3 LL_USMX
 1 0 0 0 # 4 LL_OTH
 24 0 0 0 # 5 HL_RR
 15 0 0 2 # 6 BRA_BB_hist
#
#_age_selex_patterns
#Pattern:_0; parm=0; selex=1.0 for ages 0 to maxage
#Pattern:_10; parm=0; selex=1.0 for ages 1 to maxage
#Pattern:_11; parm=2; selex=1.0  for specified min-max age
#Pattern:_12; parm=2; age logistic
#Pattern:_13; parm=8; age double logistic
#Pattern:_14; parm=nages+1; age empirical
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_16; parm=2; Coleraine - Gaussian
#Pattern:_17; parm=nages+1; empirical as random walk  N parameters to read can be overridden by setting special to non-zero
#Pattern:_41; parm=2+nages+1; // like 17, with 2 additional param for scaling (average over bin range)
#Pattern:_18; parm=8; double logistic - smooth transition
#Pattern:_19; parm=6; simple 4-parm double logistic with starting age
#Pattern:_20; parm=6; double_normal,using joiners
#Pattern:_26; parm=3; exponential-logistic in age
#Pattern:_27; parm=3+special; cubic spline in age; parm1==1 resets knots; parm1==2 resets all 
#Pattern:_42; parm=2+special+3; // cubic spline; with 2 additional param for scaling (average over bin range)
#Age patterns entered with value >100 create Min_selage from first digit and pattern from remainder
#_Pattern Discard Male Special
 0 0 0 0 # 1 PS_West
 0 0 0 0 # 2 BB_West
 0 0 0 0 # 3 LL_USMX
 0 0 0 0 # 4 LL_OTH
 0 0 0 0 # 5 HL_RR
 0 0 0 0 # 6 BRA_BB_hist
#
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
# 1   PS_West LenSelex
            20            90       46.7746            47            99             0          2          0          0          0          0          0          1          2  #  Size_DblN_peak_PS_West(1)
           -15            15      -12.1294      -3.57422            99             0          2          0          0          0          0          0          1          2  #  Size_DblN_top_logit_PS_West(1)
            -4            12       4.16633        3.1391            99             0          3          0          0          0          0          0          1          2  #  Size_DblN_ascend_se_PS_West(1)
           -10             6       5.07148       4.49981            99             0          3          0          0          0          0          0          1          2  #  Size_DblN_descend_se_PS_West(1)
          -999            15           -15           -10            99             0         -4          0          0          0          0          0          1          2  #  Size_DblN_start_logit_PS_West(1)
           -20            20      -3.57935      -20.7233            99             0          3          0          0          0          0          0          1          2  #  Size_DblN_end_logit_PS_West(1)
# 2   BB_West LenSelex
            20            90       54.3119            52            99             0          2          0          0          0          0          0          0          0  #  Size_DblN_peak_BB_West(2)
           -15            15      -11.9419      -3.23868            99             0          2          0          0          0          0          0          0          0  #  Size_DblN_top_logit_BB_West(2)
            -4            12        4.9282       3.95003            99             0          3          0          0          0          0          0          0          0  #  Size_DblN_ascend_se_BB_West(2)
           -10             6       4.84398       4.49981            99             0          3          0          0          0          0          0          0          0  #  Size_DblN_descend_se_BB_West(2)
          -999            15           -15           -10            99             0         -4          0          0          0          0          0          0          0  #  Size_DblN_start_logit_BB_West(2)
           -20            20      -10.5019      -20.7233            99             0          3          0          0          0          0          0          0          0  #  Size_DblN_end_logit_BB_West(2)
# 3   LL_USMX LenSelex
            20           126       45.5077             0             0             0          2          0          0          0          0          0          0          0  #  Size_inflection_LL_USMX(3)
          0.01           100       8.13491             0             0             0          3          0          0          0          0          0          0          0  #  Size_95%width_LL_USMX(3)
# 4   LL_OTH LenSelex
            20           126       69.7592             0             0             0          2          0          0          0          0          0          0          0  #  Size_inflection_LL_OTH(4)
          0.01           100       9.40763             0             0             0          3          0          0          0          0          0          0          0  #  Size_95%width_LL_OTH(4)
# 5   HL_RR LenSelex
            20            90       49.0478            52            99             0          2          0          0          0          0          0          0          0  #  Size_DblN_peak_HL_RR(5)
           -15            15      -11.2587      -3.23868            99             0          2          0          0          0          0          0          0          0  #  Size_DblN_top_logit_HL_RR(5)
           -10            15       4.70265        4.5254            99             0          3          0          0          0          0          0          0          0  #  Size_DblN_ascend_se_HL_RR(5)
           -10            15       6.29462       4.49981            99             0          3          0          0          0          0          0          0          0  #  Size_DblN_descend_se_HL_RR(5)
          -999            15           -15           -10            99             0         -4          0          0          0          0          0          0          0  #  Size_DblN_start_logit_HL_RR(5)
           -20            20      -11.1847      -20.7233            99             0          3          0          0          0          0          0          0          0  #  Size_DblN_end_logit_HL_RR(5)
# 6   BRA_BB_hist LenSelex
# 1   PS_West AgeSelex
# 2   BB_West AgeSelex
# 3   LL_USMX AgeSelex
# 4   LL_OTH AgeSelex
# 5   HL_RR AgeSelex
# 6   BRA_BB_hist AgeSelex
#_No_Dirichlet parameters
# timevary selex parameters 
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type    PHASE  #  parm_name
            20            90       57.8113            47            99             0      2  # Size_DblN_peak_PS_West(1)_BLK1repl_2015
           -15            15      -3.88825      -3.57422            99             0      2  # Size_DblN_top_logit_PS_West(1)_BLK1repl_2015
            -4            12       4.50829        3.1391            99             0      3  # Size_DblN_ascend_se_PS_West(1)_BLK1repl_2015
           -10             6       3.86693       4.49981            99             0      3  # Size_DblN_descend_se_PS_West(1)_BLK1repl_2015
          -999            15           -15           -10            99             0      -4  # Size_DblN_start_logit_PS_West(1)_BLK1repl_2015
           -20            20      -1.90095      -20.7233            99             0      3  # Size_DblN_end_logit_PS_West(1)_BLK1repl_2015
# info on dev vectors created for selex parms are reported with other devs after tag parameter section 
#
0   #  use 2D_AR1 selectivity(0/1)
#_no 2D_AR1 selex offset used
#
# Tag loss and Tag reporting parameters go next
0  # TG_custom:  0=no read and autogen if tag data exist; 1=read
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# deviation vectors for timevary parameters
#  base   base first block   block  env  env   dev   dev   dev   dev   dev
#  type  index  parm trend pattern link  var  vectr link _mnyr  mxyr phase  dev_vector
#      5     1     1     1     2     0     0     0     0     0     0     0
#      5     2     2     1     2     0     0     0     0     0     0     0
#      5     3     3     1     2     0     0     0     0     0     0     0
#      5     4     4     1     2     0     0     0     0     0     0     0
#      5     5     5     1     2     0     0     0     0     0     0     0
#      5     6     6     1     2     0     0     0     0     0     0     0
     #
# Input variance adjustments factors: 
 #_1=add_to_survey_CV
 #_2=add_to_discard_stddev
 #_3=add_to_bodywt_CV
 #_4=mult_by_lencomp_N
 #_5=mult_by_agecomp_N
 #_6=mult_by_size-at-age_N
 #_7=mult_by_generalized_sizecomp
#_Factor  Fleet  Value
      4      1    2.5943
      4      2     1.594
      4      3    2.5338
      4      4    0.3018
      4      5    1.6986
 -9999   1    0  # terminator
#
15 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 12 changes to default Lambdas (default value is 1.0)
# Like_comp codes:  1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch; 9=init_equ_catch; 
# 10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin; 17=F_ballpark; 18=initEQregime
#like_comp fleet  phase  value  sizefreq_method
 8 1 1 1 1
 8 2 1 1 1
 8 3 1 1 1
 8 4 1 1 1
 8 5 1 1 1
 8 6 1 1 1
 9 1 1 0 1
 9 2 1 0 1
 9 3 1 0 1
 9 4 1 0 1
 9 5 1 0 1
 9 6 1 0 1
-9999  1  1  1  1  #  terminator
#
# lambdas (for info only; columns are phases)
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_CPUE/survey:_1
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_CPUE/survey:_2
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_CPUE/survey:_3
#  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 #_CPUE/survey:_4
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_CPUE/survey:_5
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_CPUE/survey:_6
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_lencomp:_1
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_lencomp:_2
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_lencomp:_3
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_lencomp:_4
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_lencomp:_5
#  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 #_lencomp:_6
#  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 #_init_equ_catch1
#  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 #_init_equ_catch2
#  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 #_init_equ_catch3
#  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 #_init_equ_catch4
#  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 #_init_equ_catch5
#  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 #_init_equ_catch6
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_recruitments
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_parameter-priors
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_parameter-dev-vectors
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_crashPenLambda
#  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 # F_ballpark_lambda
0 # (0/1/2) read specs for more stddev reporting: 0 = skip, 1 = read specs for reporting stdev for selectivity, size, and numbers, 2 = add options for M,Dyn. Bzero, SmryBio
 # 0 2 0 0 # Selectivity: (1) fleet, (2) 1=len/2=age/3=both, (3) year, (4) N selex bins
 # 0 0 # Growth: (1) growth pattern, (2) growth ages
 # 0 0 0 # Numbers-at-age: (1) area(-1 for all), (2) year, (3) N ages
 # -1 # list of bin #'s for selex std (-1 in first bin to self-generate)
 # -1 # list of ages for growth std (-1 in first bin to self-generate)
 # -1 # list of ages for NatAge std (-1 in first bin to self-generate)
999

