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
 0.001 2 0.6001 0 0 0 -3 0 0 0 0 0.5 0 0 # NatM_Lorenzen_Fem_GP_1
# Sex: 1  BioPattern: 1  Growth
 10 50 38 0 0 0 -3 0 0 0 0 0.5 0 0 # L_at_Amin_Fem_GP_1
 40 120 76 0 0 0 -3 0 0 0 0 0.5 0 0 # L_at_Amax_Fem_GP_1
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
        0.0001            20       11.7486             7             0             0          1          0          0          0          0          0          0          0 # SR_LN(R0)
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
#  -0.00770444 -0.341757 0.560855 0.0419172 0.407497 -0.264213 0.32223 0.0529883 0.0233359 -0.246585 0.221689 0.228951 0.405882 -0.131302 0.0705392 0.100113 -0.288059 0.584886 -0.0201245 0.0659815 0.00874107 0.734433 -0.914016 -0.898487 -0.414352 -0.30344 0 0 0
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
# PS_West 0 0 0 0 0 0 0 0 0 0 0.00690894 0.0449202 0.0603487 0.000966395 0.000598716 0.000478334 0.00201856 0.00152126 0 0 0.00363703 0.000430514 0.000417095 0.00293072 0.0105085 0.00502032 0.0259937 0.0112107 0.0451797 0.0778359 0.179592 0.19768 0.229932 0.211776 0.117949 0.116133 0.0527816 0.0557748 0.0730766 0.160908 0.174641 0.303679 0.147284 0.064516 0.0594393 0.0899332 0.0653161 0.0676974 0.0573308 0.109869 0.0456292 0.0594348 0.0589371 0.0387423 0.0338671 0.0251024 0.018495 0.0411005 0.0499563 0.0298844 0.0317863 0.0198849 0.0250452 0.0280636 0.0288495 0.0512835 0.0359438 0.0263533 0.0249315 0.0249315
# BB_West 0.0155471 0.0162818 0.0174638 0.0178312 0.0192291 0.0250757 0.0211974 0.0235283 0.0421679 0.0428436 0.0201385 0.01264 0.0141536 0.0195169 0.0215248 0.0345504 0.0313431 0.021497 0.0285261 0.0220125 0.0180913 0.0248241 0.038543 0.0369378 0.0377049 0.0339023 0.0324041 0.056039 0.127672 0.263737 0.365941 0.359069 0.316336 0.582809 0.557249 0.396482 0.431901 0.475065 0.451099 0.500253 0.448937 0.438121 0.530744 0.489157 0.39775 0.487946 0.385735 0.470263 0.44024 0.460256 0.359764 0.495979 0.470211 0.462426 0.357307 0.426288 0.37928 0.423717 0.508326 0.444816 0.533679 0.628808 0.514381 0.220416 0.30821 0.453428 0.510186 0.509592 0.375435 0.375435
# LL_USMX 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00010365 9.89571e-05 0.0002567 0.00025577 0.00056417 8.29501e-05 4.83369e-05 0.000163928 1.33308e-05 2.09302e-05 0.000162209 9.50513e-05 0.00187945 0.000498757 0.000783696 0.000227572 0.000216824 0.000144353 0.000230333 0.000321407 0.000124365 0.000139943 0.000151531 0.000117797 0.00051672 0.000130739 6.58852e-05 0.000182718 0.000591333 0.000227631 0.000616506 0.00342226 0.000957505 0.00208313 0.00192323 0.00239592 0.000484747 0.000406643 0.000196131 0.00973603 0.00344186 0.000489962 0.00696303 0.005092 0.00160918 0.00708165 0.00328085 0.00486113 0.00721585 0.00284444 0.00284444
# LL_OTH 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00179478 0.00183727 0.00181735 0.00179549 0.00490083 0.00481589 0.00490542 0.010089 0.00685151 0.00457901 0.00316321 0.00252532 0.00377342 0.00320959 0.0028153 0.00482493 0.00921895 0.0151387 0.024858 0.017151 0.0407048 0.013303 0.00966034 0.0163853 0.0180611 0.0292506 0.0213085 0.020432 0.0673678 0.105557 0.0202292 0.0100874 0.0163135 0.0119795 0.0160525 0.0332943 0.0130549 0.0117215 0.0140698 0.0132185 0.01239 0.00706812 0.0102675 0.00272989 0.00770249 0.00911742 0.00613442 0.00935502 0.0119831 0.00755062 0.00610547 0.00812489 0.00994107 0.0113232 0.00407157 0.00407157
# HL_RR 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.23324e-05 0 0.000213627 0.000711347 0.00332618 1.29015e-05 0.00226917 0.000311366 0.00171792 0.00060622 0.0011479 0.00281927 0.00189828 0.000971208 0.000594363 0.00141163 0.00205876 0.00123554 0.00187285 0.00166524 0.00204544 0.00135601 0.00139178 0.0015088 0.00270825 0.000663722 0.00110133 0.00147891 0.00157049 0.00181736 0.000487791 0.000857524 0.00116852 0.00108053 0.00201428 0.00194545 0.00674732 0.00556256 0.00971398 0.0104275 0.00665235 0.0227905 0.138742 0.133569 0.0665654 0.056314 0.056314
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
           -15            15      -10.8023             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_PS_West(1)
           -15            15      -10.8329             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_BB_West(2)
           -15            15      -9.77876             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_LL_USMX(3)
           -15            15      -11.8377             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_HL_RR(5)
           -15            15      -10.7745             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_BRA_BB_hist(6)
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
            20            90       48.8713            47            99             0          2          0          0          0          0          0          1          2  #  Size_DblN_peak_PS_West(1)
           -15            15      -12.1695      -3.57422            99             0          2          0          0          0          0          0          1          2  #  Size_DblN_top_logit_PS_West(1)
            -4            12       4.39367        3.1391            99             0          3          0          0          0          0          0          1          2  #  Size_DblN_ascend_se_PS_West(1)
           -10             6       4.79317       4.49981            99             0          3          0          0          0          0          0          1          2  #  Size_DblN_descend_se_PS_West(1)
          -999            15           -15           -10            99             0         -4          0          0          0          0          0          1          2  #  Size_DblN_start_logit_PS_West(1)
           -20            20      -2.45953      -20.7233            99             0          3          0          0          0          0          0          1          2  #  Size_DblN_end_logit_PS_West(1)
# 2   BB_West LenSelex
            20            90       55.8122            52            99             0          2          0          0          0          0          0          0          0  #  Size_DblN_peak_BB_West(2)
           -15            15      -11.9043      -3.23868            99             0          2          0          0          0          0          0          0          0  #  Size_DblN_top_logit_BB_West(2)
            -4            12       4.91162       3.95003            99             0          3          0          0          0          0          0          0          0  #  Size_DblN_ascend_se_BB_West(2)
           -10             6       4.68551       4.49981            99             0          3          0          0          0          0          0          0          0  #  Size_DblN_descend_se_BB_West(2)
          -999            15           -15           -10            99             0         -4          0          0          0          0          0          0          0  #  Size_DblN_start_logit_BB_West(2)
           -20            20      -3.83445      -20.7233            99             0          3          0          0          0          0          0          0          0  #  Size_DblN_end_logit_BB_West(2)
# 3   LL_USMX LenSelex
            20           126       48.1446             0             0             0          2          0          0          0          0          0          0          0  #  Size_inflection_LL_USMX(3)
          0.01           100       8.85973             0             0             0          3          0          0          0          0          0          0          0  #  Size_95%width_LL_USMX(3)
# 4   LL_OTH LenSelex
            20           126       77.1716             0             0             0          2          0          0          0          0          0          0          0  #  Size_inflection_LL_OTH(4)
          0.01           100        13.561             0             0             0          3          0          0          0          0          0          0          0  #  Size_95%width_LL_OTH(4)
# 5   HL_RR LenSelex
            20            90        53.011            52            99             0          2          0          0          0          0          0          0          0  #  Size_DblN_peak_HL_RR(5)
           -15            15      -10.8117      -3.23868            99             0          2          0          0          0          0          0          0          0  #  Size_DblN_top_logit_HL_RR(5)
           -10            15       4.93568        4.5254            99             0          3          0          0          0          0          0          0          0  #  Size_DblN_ascend_se_HL_RR(5)
           -10            15        3.0283       4.49981            99             0          3          0          0          0          0          0          0          0  #  Size_DblN_descend_se_HL_RR(5)
          -999            15           -15           -10            99             0         -4          0          0          0          0          0          0          0  #  Size_DblN_start_logit_HL_RR(5)
           -20            20     -0.693945      -20.7233            99             0          3          0          0          0          0          0          0          0  #  Size_DblN_end_logit_HL_RR(5)
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
            20            90       57.6449            47            99             0      2  # Size_DblN_peak_PS_West(1)_BLK1repl_2015
           -15            15      -3.06442      -3.57422            99             0      2  # Size_DblN_top_logit_PS_West(1)_BLK1repl_2015
            -4            12       4.38309        3.1391            99             0      3  # Size_DblN_ascend_se_PS_West(1)_BLK1repl_2015
           -10             6       3.55683       4.49981            99             0      3  # Size_DblN_descend_se_PS_West(1)_BLK1repl_2015
          -999            15           -15           -10            99             0      -4  # Size_DblN_start_logit_PS_West(1)_BLK1repl_2015
           -20            20      -1.05796      -20.7233            99             0      3  # Size_DblN_end_logit_PS_West(1)_BLK1repl_2015
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

