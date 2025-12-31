#C file created using an r4ss function
#C file write time: 2024-11-03  12:35:36
#
0 # 0 means do not read wtatage.ss; 1 means read and usewtatage.ss and also read and use growth parameters
1 #_N_Growth_Patterns
1 #_N_platoons_Within_GrowthPattern
4 # recr_dist_method for parameters
1 # not yet implemented; Future usage:Spawner-Recruitment; 1=global; 2=by area
1 # number of recruitment settlement assignments 
0 # unused option
# for each settlement assignment:
#_GPattern	month	area	age
1	1	1	0	#_recr_dist_pattern1
#
#_Cond 0 # N_movement_definitions goes here if N_areas > 1
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) also cond on do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, morph=1, source=1 dest=2, age1=4, age2=10
#
6 #_Nblock_Patterns
2 4 1 1 1 1 #_blocks_per_pattern
#_begin and end years of blocks
1996 2005 2006 2024
1990 2004 2005 2006 2007 2016 2017 2024
2017 2024
2014 2016
1976 1976
1976 2006
#
# controls for all timevary parameters 
1 #_env/block/dev_adjust_method for all time-vary parms (1=warn relative to base parm bounds; 3=no bound check)
#
# AUTOGEN
1 1 1 1 1 # autogen: 1st element for biology, 2nd for SR, 3rd for Q, 4th reserved, 5th for selex
# where: 0 = autogen all time-varying parms; 1 = read each time-varying parm line; 2 = read then autogen if parm min==-12345
#
# setup for M, growth, maturity, fecundity, recruitment distibution, movement
#
0 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate;_5=Maunder_M;_6=Age-range_Lorenzen
#_no additional input for selected M option; read 1P per morph
1 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr;5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
1.5 #_Age(post-settlement)_for_L1;linear growth below this
999 #_Growth_Age_for_L2 (999 to use as Linf)
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0 #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
2 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
1 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
1 #_First_Mature_Age
1 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
#
#_growth_parms
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env_var&link	dev_link	dev_minyr	dev_maxyr	dev_PH	Block	Block_Fxn
  0.1	     1.5	  0.499189	 -0.81	 0.41	3	  5	0	0	0	0	0	4	2	#_NatM_p_1_Fem_GP_1  
    0	      50	    6.3923	6.1252	   99	0	  1	0	0	0	0	0	0	0	#_L_at_Amin_Fem_GP_1 
   70	     130	   99.4617	 99.46	0.015	6	  1	0	0	0	0	0	0	0	#_L_at_Amax_Fem_GP_1 
    0	       1	  0.188842	0.1966	 0.03	6	  1	0	0	0	0	0	0	0	#_VonBert_K_Fem_GP_1 
    0	      10	   3.82037	     0	    0	0	 10	0	0	0	0	0	0	0	#_CV_young_Fem_GP_1  
    0	      20	   7.42895	     0	    0	0	 10	0	0	0	0	0	0	0	#_CV_old_Fem_GP_1    
  -99	      99	3.4297e-06	     0	    0	0	 -3	0	0	0	0	0	0	0	#_Wtlen_1_Fem_GP_1   
  -99	      99	   3.27469	     0	    0	0	 -3	0	0	0	0	0	0	0	#_Wtlen_2_Fem_GP_1   
  -99	      99	      53.7	     0	    0	0	 -1	0	0	0	0	0	0	0	#_Mat50%_Fem_GP_1    
  -99	      99	 -0.273657	     0	    0	0	 -1	0	0	0	0	0	0	0	#_Mat_slope_Fem_GP_1 
  -99	      99	         1	     0	    0	0	 -1	0	0	0	0	0	0	0	#_Eggs_alpha_Fem_GP_1
  -99	      99	         0	     0	    0	0	 -1	0	0	0	0	0	0	0	#_Eggs_beta_Fem_GP_1 
  0.1	      10	         1	     1	    1	0	 -1	0	0	0	0	0	0	0	#_CohortGrowDev      
1e-06	0.999999	       0.5	   0.5	  0.5	0	-99	0	0	0	0	0	0	0	#_FracFemale_GP_1    
#_timevary MG parameters
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE
0.1	2	0.874837	-0.81	0.41	3	6	#_NatM_p_1_Fem_GP_1_BLK4repl_2014
# info on dev vectors created for MGparms are reported with other devs after tag parameter section
#
#_seasonal_effects_on_biology_parms
0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_ LO HI INIT PRIOR PR_SD PR_type PHASE
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
3 #_Spawner-Recruitment; 2=Ricker; 3=std_B-H; 4=SCAA;5=Hockey; 6=B-H_flattop; 7=survival_3Parm;8=Shepard_3Parm
0 # 0/1 to use steepness in initial equ recruitment calculation
0 # future feature: 0/1 to make realized sigmaR a function of SR curvature
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn # parm_name
 10	20	13.1552	   0	0	0	 1	0	0	0	0	0	0	0	#_SR_LN(R0)  
  0	 1	      1	   1	0	0	-1	0	0	0	0	0	0	0	#_SR_BH_steep
  0	10	   0.44	0.44	0	0	-4	0	0	0	0	0	0	0	#_SR_sigmaR  
 -5	 5	      0	   0	0	0	-3	0	0	0	0	0	5	1	#_SR_regime  
-99	99	      0	   0	0	0	-1	0	0	0	0	0	0	0	#_SR_autocorr
# timevary SR parameters
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE
-10	10	-0.523818	0	0	0	1	#_SR_regime_BLK5add_1976
2 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
1978 # first year of main recr_devs; early devs can preceed this era
2024 # last year of main recr_devs; forecast devs start in following year
1 #_recdev phase
1 # (0/1) to read 13 advanced options
1967 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
2 #_recdev_early_phase
-1 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
1 #_lambda for Fcast_recr_like occurring before endyr+1
1960.73622545148 #_last_yr_nobias_adj_in_MPD; begin of ramp
1988.86106210699 #_first_yr_fullbias_adj_in_MPD; begin of plateau
2020.97137754315 #_last_yr_fullbias_adj_in_MPD
2023.57715575014 #_end_yr_for_ramp_in_MPD (can be in forecast to shape ramp, but SS sets bias_adj to 0.0 for fcast yrs)
0.913810328379987 #_max_bias_adj_in_MPD (-1 to override ramp and set biasadj=1.0 for all estimated recdevs)
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
#Fishing Mortality info
0 # F ballpark
-1999 # F ballpark year (neg value to disable)
3 # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
5 # max F or harvest rate, depends on F_Method
5 # N iterations for tuning F in hybrid method (recommend 3 to 7)
#
#_initial_F_parms; count = 0
#
#_Q_setup for fleets with cpue or survey data
#_fleet	link	link_info	extra_se	biasadj	float  #  fleetname
    4	1	0	0	0	0	#_Srv       
    5	1	0	0	0	0	#_LLSrv     
    6	0	0	0	0	0	#_IPHCLL    
    7	0	0	0	0	0	#_ADFG      
    8	3	0	0	0	0	#_SPAWN     
    9	3	0	0	0	0	#_Seine     
-9999	0	0	0	0	0	#_terminator
#_Q_parms(if_any);Qunits_are_ln(q)
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn  #  parm_name
-10	10	0.101348	0	99	0	 1	  0	0	0	0	0	0	0	#_LnQ_base_Srv(4)   
-10	10	0.115124	0	99	0	 1	101	0	0	0	0	0	0	#_LnQ_base_LLSrv(5) 
-25	25	       0	0	 1	0	-1	  0	0	0	0	0	0	0	#_LnQ_base_IPHCLL(6)
-25	25	       0	0	 1	0	-1	  0	0	0	0	0	0	0	#_LnQ_base_ADFG(7)  
-25	25	       0	0	99	0	-1	  0	0	0	0	0	0	0	#_LnQ_base_SPAWN(8) 
 -1	 5	       0	0	99	0	-1	  0	0	0	0	0	0	0	#_Q_power_SPAWN(8)  
-50	25	       0	0	99	0	-1	  0	0	0	0	0	0	0	#_LnQ_base_Seine(9) 
 -1	 5	       0	0	99	0	-1	  0	0	0	0	0	0	0	#_Q_power_Seine(9)  
# timevary Q parameters
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE
-10	10	1.12292	0	99	0	5	#_LnQ_base_LLSrv(5)_ENV_add
# info on dev vectors created for Q parms are reported with other devs after tag parameter section
#
#_size_selex_patterns
#_Pattern	Discard	Male	Special
24	0	0	0	#_1 FshTrawl
24	0	0	0	#_2 FshLL   
24	0	0	0	#_3 FshPot  
24	0	0	0	#_4 Srv     
24	0	0	0	#_5 LLSrv   
 0	0	0	0	#_6 IPHCLL  
 0	0	0	0	#_7 ADFG    
 0	0	0	0	#_8 SPAWN   
 0	0	0	0	#_9 Seine   
#
#_age_selex_patterns
#_Pattern	Discard	Male	Special
10	0	0	0	#_1 FshTrawl
10	0	0	0	#_2 FshLL   
10	0	0	0	#_3 FshPot  
10	0	0	0	#_4 Srv     
10	0	0	0	#_5 LLSrv   
 0	0	0	0	#_6 IPHCLL  
 0	0	0	0	#_7 ADFG    
 0	0	0	0	#_8 SPAWN   
 0	0	0	0	#_9 Seine   
#
#_SizeSelex
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn  #  parm_name
   10	    100	  57.6551	 50	0	0	 1	0	1	1977	1989	3	2	2	#_SizeSel_P_1_FshTrawl(1)
  -20	     10	 -4.57915	  0	0	0	 2	0	0	1977	1989	0	2	2	#_SizeSel_P_2_FshTrawl(1)
  -10	     10	  5.09148	  0	0	0	 2	0	1	1977	1989	3	2	2	#_SizeSel_P_3_FshTrawl(1)
  -10	     10	-0.301675	 10	0	0	 2	0	1	1977	1989	3	2	2	#_SizeSel_P_4_FshTrawl(1)
-1000	2.71828	     -999	-10	0	0	-2	0	0	   0	   0	0	0	0	#_SizeSel_P_5_FshTrawl(1)
  -10	     10	       10	 10	0	0	-2	0	0	   0	   0	0	0	0	#_SizeSel_P_6_FshTrawl(1)
   10	     90	  66.4253	 50	0	0	 1	0	1	1978	1989	3	2	2	#_SizeSel_P_1_FshLL(2)   
  -20	     10	 -5.06337	  0	0	0	 2	0	0	1978	1989	0	2	2	#_SizeSel_P_2_FshLL(2)   
  -10	     10	  5.14687	  0	0	0	 2	0	1	1978	1989	3	2	2	#_SizeSel_P_3_FshLL(2)   
    0	     10	       10	 10	0	0	-2	0	0	   0	   0	0	2	2	#_SizeSel_P_4_FshLL(2)   
-1000	2.71828	     -999	-10	0	0	-2	0	0	   0	   0	0	0	0	#_SizeSel_P_5_FshLL(2)   
  -10	     10	       10	 10	0	0	-2	0	0	   0	   0	0	0	0	#_SizeSel_P_6_FshLL(2)   
   10	     90	  70.7282	 50	0	0	 1	0	0	   0	   0	0	3	2	#_SizeSel_P_1_FshPot(3)  
  -20	     10	 -12.0643	  0	0	0	 2	0	0	   0	   0	0	3	2	#_SizeSel_P_2_FshPot(3)  
  -10	     10	  5.01625	  0	0	0	 2	0	0	   0	   0	0	3	2	#_SizeSel_P_3_FshPot(3)  
    0	     10	   4.0558	 10	0	0	 2	0	0	   0	   0	0	0	0	#_SizeSel_P_4_FshPot(3)  
-1000	2.71828	     -999	-10	0	0	-2	0	0	   0	   0	0	0	0	#_SizeSel_P_5_FshPot(3)  
  -10	     10	 0.284913	 10	0	0	 2	0	0	   0	   0	0	0	0	#_SizeSel_P_6_FshPot(3)  
   10	     90	  60.2246	 50	0	0	 1	0	0	   0	   0	0	1	2	#_SizeSel_P_1_Srv(4)     
  -20	     10	 -11.9482	  0	0	0	 2	0	0	   0	   0	0	1	2	#_SizeSel_P_2_Srv(4)     
  -10	     10	  5.54831	  0	0	0	 2	0	0	   0	   0	0	1	2	#_SizeSel_P_3_Srv(4)     
    0	     10	   3.9945	 10	0	0	 5	0	0	   0	   0	0	1	2	#_SizeSel_P_4_Srv(4)     
  -10	2.71828	 -3.10266	-10	0	0	 2	0	0	   0	   0	0	1	2	#_SizeSel_P_5_Srv(4)     
  -10	     10	 -0.66555	 10	0	0	 5	0	0	   0	   0	0	1	2	#_SizeSel_P_6_Srv(4)     
   10	     90	  65.5449	 50	0	0	 1	0	0	   0	   0	0	0	0	#_SizeSel_P_1_LLSrv(5)   
  -20	     10	  -12.478	  0	0	0	 2	0	0	   0	   0	0	0	0	#_SizeSel_P_2_LLSrv(5)   
  -10	     10	  4.68122	  0	0	0	 2	0	0	   0	   0	0	0	0	#_SizeSel_P_3_LLSrv(5)   
    0	     10	  4.64978	 10	0	0	 2	0	0	   0	   0	0	0	0	#_SizeSel_P_4_LLSrv(5)   
-1000	2.71828	     -999	-10	0	0	-2	0	0	   0	   0	0	0	0	#_SizeSel_P_5_LLSrv(5)   
  -10	     10	-0.466947	 10	0	0	 2	0	0	   0	   0	0	0	0	#_SizeSel_P_6_LLSrv(5)   
#_AgeSelex
#_No age_selex_parm
# timevary selex parameters 
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE
   10	     90	  78.6311	 50	  0	0	 1	#_SizeSel_P_1_FshTrawl(1)_BLK2repl_1990
   10	    100	  94.5791	 50	  0	0	 1	#_SizeSel_P_1_FshTrawl(1)_BLK2repl_2005
   10	     90	  79.9048	 50	  0	0	 1	#_SizeSel_P_1_FshTrawl(1)_BLK2repl_2007
   10	     90	  77.2556	 50	  0	0	-1	#_SizeSel_P_1_FshTrawl(1)_BLK2repl_2017
1e-04	      2	      0.2	0.2	0.5	6	-5	#_SizeSel_P_1_FshTrawl(1)_dev_se       
-0.99	   0.99	        0	  0	0.5	6	-6	#_SizeSel_P_1_FshTrawl(1)_dev_autocorr 
  -20	     10	 0.434816	  0	  0	0	 2	#_SizeSel_P_2_FshTrawl(1)_BLK2repl_1990
  -20	     10	 -4.93158	  0	  0	0	 2	#_SizeSel_P_2_FshTrawl(1)_BLK2repl_2005
  -20	     10	 -4.75885	  0	  0	0	 2	#_SizeSel_P_2_FshTrawl(1)_BLK2repl_2007
  -20	     10	 -4.88944	  0	  0	0	 2	#_SizeSel_P_2_FshTrawl(1)_BLK2repl_2017
  -10	     10	  6.02471	  0	  0	0	 2	#_SizeSel_P_3_FshTrawl(1)_BLK2repl_1990
  -10	     10	  6.26352	  0	  0	0	 2	#_SizeSel_P_3_FshTrawl(1)_BLK2repl_2005
  -10	     10	  6.24936	  0	  0	0	 2	#_SizeSel_P_3_FshTrawl(1)_BLK2repl_2007
  -10	     10	  6.29422	  0	  0	0	 2	#_SizeSel_P_3_FshTrawl(1)_BLK2repl_2017
1e-04	      2	      0.2	0.2	0.5	6	-5	#_SizeSel_P_3_FshTrawl(1)_dev_se       
-0.99	   0.99	        0	  0	0.5	6	-6	#_SizeSel_P_3_FshTrawl(1)_dev_autocorr 
    0	     10	  5.10604	 10	  0	0	 2	#_SizeSel_P_4_FshTrawl(1)_BLK2repl_1990
    0	     10	  5.15973	 10	  0	0	 2	#_SizeSel_P_4_FshTrawl(1)_BLK2repl_2005
    0	     10	  5.79906	 10	  0	0	 2	#_SizeSel_P_4_FshTrawl(1)_BLK2repl_2007
    0	     10	  5.26328	 10	  0	0	 2	#_SizeSel_P_4_FshTrawl(1)_BLK2repl_2017
1e-04	      2	      0.2	0.2	0.5	6	-5	#_SizeSel_P_4_FshTrawl(1)_dev_se       
-0.99	   0.99	        0	  0	0.5	6	-6	#_SizeSel_P_4_FshTrawl(1)_dev_autocorr 
   10	     90	  71.8699	 50	  0	0	 1	#_SizeSel_P_1_FshLL(2)_BLK2repl_1990   
   10	     90	  68.4629	 50	  0	0	 1	#_SizeSel_P_1_FshLL(2)_BLK2repl_2005   
   10	     90	  72.5472	 50	  0	0	 1	#_SizeSel_P_1_FshLL(2)_BLK2repl_2007   
   10	     90	  71.9494	 50	  0	0	 1	#_SizeSel_P_1_FshLL(2)_BLK2repl_2017   
1e-04	      2	      0.2	0.2	0.5	6	-5	#_SizeSel_P_1_FshLL(2)_dev_se          
-0.99	   0.99	        0	  0	0.5	6	-6	#_SizeSel_P_1_FshLL(2)_dev_autocorr    
  -20	     10	-0.224231	  0	  0	0	 2	#_SizeSel_P_2_FshLL(2)_BLK2repl_1990   
  -20	     10	 -5.10633	  0	  0	0	 2	#_SizeSel_P_2_FshLL(2)_BLK2repl_2005   
  -20	     10	  -5.1146	  0	  0	0	 2	#_SizeSel_P_2_FshLL(2)_BLK2repl_2007   
  -20	     10	 -4.69045	  0	  0	0	 2	#_SizeSel_P_2_FshLL(2)_BLK2repl_2017   
  -10	     10	  5.36118	  0	  0	0	 2	#_SizeSel_P_3_FshLL(2)_BLK2repl_1990   
  -10	     10	  5.25191	  0	  0	0	 2	#_SizeSel_P_3_FshLL(2)_BLK2repl_2005   
  -10	     10	  5.44255	  0	  0	0	 2	#_SizeSel_P_3_FshLL(2)_BLK2repl_2007   
  -10	     10	  5.35388	  0	  0	0	 2	#_SizeSel_P_3_FshLL(2)_BLK2repl_2017   
1e-04	      2	      0.2	0.2	0.5	6	-5	#_SizeSel_P_3_FshLL(2)_dev_se          
-0.99	   0.99	        0	  0	0.5	6	-6	#_SizeSel_P_3_FshLL(2)_dev_autocorr    
    0	     10	       10	 10	  0	0	-2	#_SizeSel_P_4_FshLL(2)_BLK2repl_1990   
    0	     10	       10	 10	  0	0	-2	#_SizeSel_P_4_FshLL(2)_BLK2repl_2005   
    0	     10	       10	 10	  0	0	-2	#_SizeSel_P_4_FshLL(2)_BLK2repl_2007   
    0	     10	       10	 10	  0	0	-2	#_SizeSel_P_4_FshLL(2)_BLK2repl_2017   
   10	     90	  75.9646	 50	  0	0	 1	#_SizeSel_P_1_FshPot(3)_BLK3repl_2017  
  -20	     10	  2.80531	  0	  0	0	 2	#_SizeSel_P_2_FshPot(3)_BLK3repl_2017  
  -10	     10	  5.47946	  0	  0	0	 2	#_SizeSel_P_3_FshPot(3)_BLK3repl_2017  
   10	     90	   59.965	 50	  0	0	 1	#_SizeSel_P_1_Srv(4)_BLK1repl_1996     
   10	     90	  56.5473	 50	  0	0	 1	#_SizeSel_P_1_Srv(4)_BLK1repl_2006     
  -20	     10	 -4.15407	  0	  0	0	 2	#_SizeSel_P_2_Srv(4)_BLK1repl_1996     
  -20	     10	 -3.73717	  0	  0	0	 2	#_SizeSel_P_2_Srv(4)_BLK1repl_2006     
  -10	     10	   5.4719	  0	  0	0	 2	#_SizeSel_P_3_Srv(4)_BLK1repl_1996     
  -10	     10	  5.09389	  0	  0	0	 2	#_SizeSel_P_3_Srv(4)_BLK1repl_2006     
    0	     10	  4.21246	 10	  0	0	 5	#_SizeSel_P_4_Srv(4)_BLK1repl_1996     
    0	     10	    4.195	 10	  0	0	 5	#_SizeSel_P_4_Srv(4)_BLK1repl_2006     
  -10	2.71828	 -2.09021	-10	  0	0	 2	#_SizeSel_P_5_Srv(4)_BLK1repl_1996     
  -10	2.71828	 -2.01521	-10	  0	0	 2	#_SizeSel_P_5_Srv(4)_BLK1repl_2006     
    0	     10	       10	 10	  0	0	-5	#_SizeSel_P_6_Srv(4)_BLK1repl_1996     
    0	     10	       10	 10	  0	0	-5	#_SizeSel_P_6_Srv(4)_BLK1repl_2006     
# info on dev vectors created for selex parms are reported with other devs after tag parameter section
#
0 #  use 2D_AR1 selectivity(0/1):  experimental feature
#_no 2D_AR1 selex offset used
# Tag loss and Tag reporting parameters go next
0 # TG_custom:  0=no read; 1=read if tags exist
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# Input variance adjustments factors: 
#_Data_type	Fleet	Value
    1	8	0	#_Variance_adjustment_list1 
    1	9	0	#_Variance_adjustment_list2 
    1	6	0	#_Variance_adjustment_list3 
    1	7	0	#_Variance_adjustment_list4 
    1	4	0	#_Variance_adjustment_list5 
    1	5	0	#_Variance_adjustment_list6 
    4	1	1	#_Variance_adjustment_list7 
    4	2	1	#_Variance_adjustment_list8 
    4	3	1	#_Variance_adjustment_list9 
    4	4	1	#_Variance_adjustment_list10
    4	5	1	#_Variance_adjustment_list11
    6	3	1	#_Variance_adjustment_list12
    6	4	1	#_Variance_adjustment_list13
    6	5	0	#_Variance_adjustment_list14
    4	6	0	#_Variance_adjustment_list15
    5	6	0	#_Variance_adjustment_list16
    6	6	0	#_Variance_adjustment_list17
    4	7	0	#_Variance_adjustment_list18
    5	7	0	#_Variance_adjustment_list19
    6	7	0	#_Variance_adjustment_list20
-9999	0	0	#_terminator                
#
1 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 4 changes to default Lambdas (default value is 1.0)
#_like_comp	fleet	phase	value	sizefreq_method
    1	8	1	0	1	#_Surv_SPAWN_Phz1 
    1	6	1	0	1	#_Surv_IPHCLL_Phz1
    1	7	1	0	1	#_Surv_ADFG_Phz1  
    1	9	1	0	1	#_Surv_Seine_Phz1 
-9999	0	0	0	0	#_terminator      
#
0 # 0/1 read specs for more stddev reporting
#
999
