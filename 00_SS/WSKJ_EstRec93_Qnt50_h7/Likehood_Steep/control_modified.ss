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
        0.0001            20       11.2537             7             0             0          1          0          0          0          0          0          0          0 # SR_LN(R0)
 0.2 1 0.9101 0.7 0 0 -3 0 0 0 0 0 0 0 # SR_BH_steep
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
#  0.00716653 -0.339071 0.523067 0.0379306 0.402648 -0.249943 0.346826 0.0638619 0.0286024 -0.248217 0.228641 0.230573 0.394795 -0.136897 0.0760913 0.0947499 -0.270558 0.598783 -0.0233276 0.103285 0.0365544 0.705095 -0.989681 -0.886079 -0.380703 -0.354194 0 0 0
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
# PS_West 0 0 0 0 0 0 0 0 0 0 0.008131 0.0528826 0.0712012 0.00113985 0.000704317 0.000562354 0.00237386 0.00178848 0 0 0.00427223 0.00050571 0.000490304 0.0034484 0.0123737 0.00591322 0.0306304 0.0132311 0.0536025 0.0937151 0.221887 0.249913 0.292876 0.274357 0.153954 0.148213 0.0658553 0.069047 0.0903584 0.20045 0.21913 0.383674 0.184473 0.0799079 0.0742813 0.111009 0.0816916 0.0836419 0.0705806 0.135344 0.0562821 0.0731638 0.0727216 0.0475963 0.0419456 0.0309539 0.0228502 0.0509854 0.0614534 0.0368144 0.0392898 0.0245304 0.030724 0.0348046 0.035116 0.0624525 0.0437155 0.0325281 0.0296417 0.0296417
# BB_West 0.0180357 0.0189078 0.0202923 0.0207249 0.0223529 0.0291585 0.0246554 0.0273665 0.0490842 0.0499352 0.0234722 0.0147357 0.0165369 0.0228027 0.0250781 0.0402171 0.0364944 0.0250241 0.0331933 0.0256105 0.0210411 0.0288716 0.0448607 0.0430371 0.0439662 0.0395463 0.0378139 0.0654996 0.150023 0.31467 0.4487 0.451813 0.40187 0.753357 0.727787 0.506501 0.537824 0.585504 0.555024 0.620096 0.56134 0.551993 0.664704 0.605104 0.494372 0.600165 0.479507 0.578962 0.539567 0.564031 0.441403 0.607151 0.577506 0.565596 0.439616 0.5224 0.465674 0.522225 0.623137 0.54548 0.655929 0.773272 0.629085 0.27351 0.375513 0.553514 0.619448 0.630455 0.441971 0.441971
# LL_USMX 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.000107656 0.000102725 0.000266388 0.000265529 0.000586294 8.63021e-05 5.03385e-05 0.000170856 1.39186e-05 2.19642e-05 0.000172783 0.000104423 0.00215095 0.000589472 0.00095776 0.000287266 0.000272223 0.000176505 0.000276204 0.000381984 0.000148174 0.000168498 0.000184445 0.000144104 0.000625366 0.00015798 7.86038e-05 0.000218683 0.000699582 0.000268371 0.000724033 0.00401075 0.00111407 0.00244273 0.00225727 0.00282013 0.000565359 0.00047299 0.000228464 0.0113189 0.0040469 0.000575765 0.0082168 0.00600276 0.00191585 0.00820369 0.00375347 0.00558556 0.0084444 0.00326712 0.00326712
# LL_OTH 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00166431 0.0017024 0.00168249 0.00166084 0.00452981 0.00444858 0.00452852 0.00930959 0.00632293 0.00422893 0.00292477 0.00233762 0.00349699 0.00297942 0.0026254 0.00455823 0.00896196 0.0153823 0.0264752 0.0191972 0.0481017 0.016028 0.0113959 0.0187785 0.0203011 0.0327072 0.0240496 0.0233873 0.0784947 0.122739 0.0230415 0.0113674 0.0182578 0.0133404 0.0177445 0.0365708 0.0142761 0.0127305 0.0152955 0.0144116 0.0134798 0.00767012 0.0110765 0.00293234 0.00832746 0.00989073 0.00670335 0.0103409 0.0132995 0.00826716 0.00661195 0.0087466 0.0107045 0.0122381 0.0044143 0.0044143
# HL_RR 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.56338e-05 0 0.000245579 0.00081825 0.00383232 1.49415e-05 0.00266686 0.000376247 0.00213933 0.000767967 0.00148838 0.00371309 0.00245914 0.00122933 0.000743672 0.001759 0.00257981 0.00156041 0.00238492 0.00211199 0.00257266 0.00170223 0.00173452 0.00189476 0.00337712 0.000821347 0.00136399 0.00183174 0.00194087 0.00224748 0.000601719 0.00106388 0.00144654 0.00133643 0.00249778 0.00240301 0.00832406 0.00690778 0.0120639 0.0128949 0.00831389 0.0281485 0.170865 0.163152 0.0827824 0.0669034 0.0669034
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
           -15            15      -10.5905             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_PS_West(1)
           -15            15      -10.6306             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_BB_West(2)
           -15            15      -9.60172             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_LL_USMX(3)
           -15            15      -11.6173             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_HL_RR(5)
           -15            15      -10.5589             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_BRA_BB_hist(6)
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
            20            90       48.7977            47            99             0          2          0          0          0          0          0          1          2  #  Size_DblN_peak_PS_West(1)
           -15            15      -12.1872      -3.57422            99             0          2          0          0          0          0          0          1          2  #  Size_DblN_top_logit_PS_West(1)
            -4            12       4.39271        3.1391            99             0          3          0          0          0          0          0          1          2  #  Size_DblN_ascend_se_PS_West(1)
           -10             6       4.80059       4.49981            99             0          3          0          0          0          0          0          1          2  #  Size_DblN_descend_se_PS_West(1)
          -999            15           -15           -10            99             0         -4          0          0          0          0          0          1          2  #  Size_DblN_start_logit_PS_West(1)
           -20            20      -2.57044      -20.7233            99             0          3          0          0          0          0          0          1          2  #  Size_DblN_end_logit_PS_West(1)
# 2   BB_West LenSelex
            20            90        55.561            52            99             0          2          0          0          0          0          0          0          0  #  Size_DblN_peak_BB_West(2)
           -15            15      -11.9476      -3.23868            99             0          2          0          0          0          0          0          0          0  #  Size_DblN_top_logit_BB_West(2)
            -4            12       4.89013       3.95003            99             0          3          0          0          0          0          0          0          0  #  Size_DblN_ascend_se_BB_West(2)
           -10             6       4.69562       4.49981            99             0          3          0          0          0          0          0          0          0  #  Size_DblN_descend_se_BB_West(2)
          -999            15           -15           -10            99             0         -4          0          0          0          0          0          0          0  #  Size_DblN_start_logit_BB_West(2)
           -20            20      -4.18907      -20.7233            99             0          3          0          0          0          0          0          0          0  #  Size_DblN_end_logit_BB_West(2)
# 3   LL_USMX LenSelex
            20           126       47.9544             0             0             0          2          0          0          0          0          0          0          0  #  Size_inflection_LL_USMX(3)
          0.01           100       8.76354             0             0             0          3          0          0          0          0          0          0          0  #  Size_95%width_LL_USMX(3)
# 4   LL_OTH LenSelex
            20           126       76.6849             0             0             0          2          0          0          0          0          0          0          0  #  Size_inflection_LL_OTH(4)
          0.01           100        13.486             0             0             0          3          0          0          0          0          0          0          0  #  Size_95%width_LL_OTH(4)
# 5   HL_RR LenSelex
            20            90       52.8994            52            99             0          2          0          0          0          0          0          0          0  #  Size_DblN_peak_HL_RR(5)
           -15            15      -10.8655      -3.23868            99             0          2          0          0          0          0          0          0          0  #  Size_DblN_top_logit_HL_RR(5)
           -10            15       4.94186        4.5254            99             0          3          0          0          0          0          0          0          0  #  Size_DblN_ascend_se_HL_RR(5)
           -10            15       3.16224       4.49981            99             0          3          0          0          0          0          0          0          0  #  Size_DblN_descend_se_HL_RR(5)
          -999            15           -15           -10            99             0         -4          0          0          0          0          0          0          0  #  Size_DblN_start_logit_HL_RR(5)
           -20            20     -0.839894      -20.7233            99             0          3          0          0          0          0          0          0          0  #  Size_DblN_end_logit_HL_RR(5)
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
            20            90       57.6419            47            99             0      2  # Size_DblN_peak_PS_West(1)_BLK1repl_2015
           -15            15      -3.13691      -3.57422            99             0      2  # Size_DblN_top_logit_PS_West(1)_BLK1repl_2015
            -4            12       4.39019        3.1391            99             0      3  # Size_DblN_ascend_se_PS_West(1)_BLK1repl_2015
           -10             6       3.59178       4.49981            99             0      3  # Size_DblN_descend_se_PS_West(1)_BLK1repl_2015
          -999            15           -15           -10            99             0      -4  # Size_DblN_start_logit_PS_West(1)_BLK1repl_2015
           -20            20      -1.21623      -20.7233            99             0      3  # Size_DblN_end_logit_PS_West(1)_BLK1repl_2015
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

