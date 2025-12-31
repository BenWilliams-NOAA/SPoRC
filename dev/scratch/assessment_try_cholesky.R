# Purpose: To create an RTMB model for the sablefish assessment using SPoRC to compare likelihoods (iid-Logistic-Normal)
# Creator: Matthew LH. Cheng (UAF CFOS)
# Date 5/28/25

unloadNamespace('SPoRC')
library(here)
library(SPoRC)
library(ggplot2)
library(RTMB)
library(tidyr)
library(tidyverse)

mod_name <- "25_12_iidLN"

# Load in assessment data
SPoRC_one_reg <- readRDS(here("/Users/matthewcheng/Desktop/Side projects/SexCompLikelihoods/inputs/SPoRC_one_reg_age_joint_no_len_kp_JPN.rds"))

#------------------ Initialize model dimensions and data list
input_list <- Setup_Mod_Dim(years = 1:SPoRC_one_reg$nyrs,                      # vector of years
                            ages = 1:SPoRC_one_reg$nages,                      # vector of ages
                            lens = SPoRC_one_reg$length_bins,                  # length bins
                            n_regions = SPoRC_one_reg$nreg,                    # number of regions, if nreg==1 then AK-wide, if n==5 then 1==BS, 2==AI, 3==WG, 4==CG, 5==EG
                            n_sexes = SPoRC_one_reg$nsex,                      # number of sexes, 1==female, 2==male
                            n_fish_fleets = SPoRC_one_reg$nflt,                # number of fishery fleet, 1==fixed gear, 2==trawl gear
                            n_srv_fleets = SPoRC_one_reg$nsrv_flt,             # number of survey fleets, 1==LLS, 2==GOA Trawl Survey, 3==JPN LLS
                            verbose = TRUE
)

#------------------- Setup recruitment stuff (using defaults for other stuff)
input_list <- Setup_Mod_Rec(input_list = input_list,                           # input data list from above
                            do_rec_bias_ramp = 1,                              # do bias ramp (0 == don't do bias ramp, 1 == do bias ramp)
                            bias_year = c(length(SPoRC_one_reg$styr:1980),     # breakpoints for bias ramp--NOTE these are flexibly coded and don't need to be changed as add more years to model (first entry == years no bias ramp - 1960 to 1980--years 1 to 21;
                                          length(SPoRC_one_reg$styr:1990),     # second entry == years ascending limb of bias ramp - 1981 to 1990--years 21 to 30;
                                          (length(SPoRC_one_reg$styr:          # third entry == years full bias correction - 1990 to 5 years from terminal year (currently 2023)--years 31 to 60 currently--start descending limb in this year
                                                    SPoRC_one_reg$endyr)-5),
                                          (length(SPoRC_one_reg$styr:          # fourth entry == years no bias correction at end of model - terminal year of recruitment estimate (2023)-- year 64 currently--(and again in terminal year of model), which is terminal year of model-1 since terminal year is fixed at average recruitment)
                                                    SPoRC_one_reg$endyr) - 1)) ,
                            sigmaR_switch = as.integer(length(                 # when to switch from early to late sigmaR
                              SPoRC_one_reg$styr:1975)),
                            dont_est_recdev_last = 1,                          # don't estimate last recruitment deviate (set equal to mean recruitment)
                            ln_sigmaR = log(c(0.4, 0.9)),                      # early and late sigma_R starting values
                            rec_model = "mean_rec",                            # recruitment model, use mean recruitment (not BH SRR)
                            sigmaR_spec = "fix",                               # fix early sigmaR, estimate late sigmaR
                            InitDevs_spec = NULL,                              # estimate all initial deviations
                            RecDevs_spec = NULL,                               # estimate all recruitment deviations
                            init_age_strc = 1,                                 # Initial age structure is derived by assuming a geometric series (init_age_strc = 1; the alternative is iterating age structure to some equilibrium init_age_strc = 0),
                            init_F_prop = 0.1                                  # A 10% fraction of the mean fishing mortality rate is assumed for initializing the age structure.
)

#------------------- Setup Biological Inputs
M_prior <- data.frame(
  regionblk = 1,
  yearblk = 1,
  ageblk = c(1),
  sexblk = c(1,2),
  mu = 0.085,
  sd = 0.1
)

input_list <- Setup_Mod_Biologicals(input_list = input_list,
                                    # Data inputs
                                    WAA = SPoRC_one_reg$WAA, # weight-at-age
                                    MatAA = SPoRC_one_reg$MatAA, # maturity at age
                                    AgeingError = as.matrix(SPoRC_one_reg$AgeError), # ageing error
                                    SizeAgeTrans = SPoRC_one_reg$SizeAgeTrans, # size age transition matrix
                                    # Model options
                                    Use_M_prior = 1, # use natural mortality prior
                                    M_prior = M_prior, # mean and sd for M prior
                                    fit_lengths = 1, # fitting length compositions
                                    M_spec = "est_ln_M",
                                    M_sexblk_spec = list(1,2)
)

#------------------- Setup Movement (No Movement in Operational Assessment)
input_list <- Setup_Mod_Movement(input_list = input_list,
                                 use_fixed_movement = 1,                       # don't est move
                                 Fixed_Movement = NA,                          # use identity move
                                 do_recruits_move = 0                          # recruits don't move
)

#------------------- Setup Tagging (No Tagging Fit in Operational Assessment)
input_list <- Setup_Mod_Tagging(input_list = input_list,
                                UseTagging = 0                                 # don't fit/use tagging data
)

#------------------- Setup Catch Specifications
input_list <- Setup_Mod_Catch_and_F(input_list = input_list,
                                    ObsCatch = SPoRC_one_reg$ObsCatch,         # Observed Catch in mt
                                    Catch_Type =                               # just a switch to determine whether to fit region-specific catch, ==1 aggregate across regions (does nothing for single region model)
                                      array(1, dim = c(length(input_list$data$years),
                                                       input_list$data$n_fish_fleets)),
                                    UseCatch = SPoRC_one_reg$UseCatch,         # Years in which catch data is available and should be fit in the model

                                    # Model options
                                    Use_F_pen = 1,                             # whether to use f penalty, == 0 don't use, == 1 use
                                    ln_sigmaF = array(log(1),
                                                      dim = c(SPoRC_one_reg$nreg,
                                                              SPoRC_one_reg$nflt)),          # the sd for the F prior (high value gives low 'weight' to the F prior...1.0 is the default value)
                                    sigmaC_spec = 'fix',                       # Whether the variance term for fitting catch data is estimated or fixed
                                    ln_sigmaC = array(log(.05),
                                                      dim = c(SPoRC_one_reg$nreg,
                                                              SPoRC_one_reg$nyrs,
                                                              SPoRC_one_reg$nflt)),        # 0.05 is common value to fit catch relatively closely, but allow for some error
                                    Catch_Constant = c(0.0, 0.0)               # The small constant to add to the catch term to avoid log(0)--Fournier robust likelihood approach; fixed gear then trawl fishery
)

#------------------- Setup fishery indices and compositions
input_list <- Setup_Mod_FishIdx_and_Comps(input_list = input_list,
                                          ObsFishIdx =                         # Observed fishery indices (CPUE)
                                            SPoRC_one_reg$ObsFishIdx,
                                          ObsFishIdx_SE =                      # Fishery index standard errors
                                            (SPoRC_one_reg$ObsFishIdx_SE/
                                               SPoRC_one_reg$ObsFishIdx) * 2,  # converting to true variance for use in the lognormal likelihood; *2 is scaling fishery CPUE to mean CV of 0.2 (calc historically as mean of 0.1 which is too small)
                                          UseFishIdx =                         # Switch to identify years in which fishery indices exist and whether to fit them in the model
                                            SPoRC_one_reg$UseFishIdx,
                                          ObsFishAgeComps =                    # Observed fishery age comps
                                            SPoRC_one_reg$ObsFishAgeComps,
                                          UseFishAgeComps =                    # Switch identifying whether fishery age comps exist and whether to fit them
                                            SPoRC_one_reg$UseFishAgeComps,
                                          ISS_FishAgeComps =                   # The input sample sizes for fishery age comps
                                            SPoRC_one_reg$ISS_FishAgeComps,
                                          ObsFishLenComps =                    # Observed fishery length comps
                                            SPoRC_one_reg$ObsFishLenComps,
                                          UseFishLenComps =                    # Whether to fit fishery length comps in a given year
                                            SPoRC_one_reg$UseFishLenComps,
                                          ISS_FishLenComps =                   # Input sample size for fishery length comps
                                            SPoRC_one_reg$ISS_FishLenComps,

                                          # Model options
                                          fish_idx_type =                      # Define units (e.g. biomass vs. abundance/numbers) for fishery indices (none is used if no index for a given fleet)
                                            c("biom", "none"),
                                          FishAgeComps_LikeType =              # age comp likelihoods by fleet (none if not fit or not available for given fleet)
                                            c("Multinomial", "none"),
                                          FishLenComps_LikeType =              # length comp likelihoods for fishery fleet 1 and 2
                                            c("none", "Multinomial"),
                                          FishAgeComps_Type =                  # age comp structure for fishery fleet 1 and 2 (agg_Year_1-terminal_Fleet_1 indicates that age comps are aggregated across sex for all years for fleet 1--fixed gear)
                                            c("spltRjntS_Year_1-terminal_Fleet_1",   # NOTE: these entries are needed for all years and all fleets, so if want to agg for some years then disagg for remaining years, would need two entries for that fleet
                                              "none_Year_1-terminal_Fleet_2"), # NOTE: if want joint by sex would use spltRjntS_Year_1-terminal_Fleet_1
                                          FishLenComps_Type =                  # length comp structure for fishery fleet 1 and 2 (spltRspltS_Year_1-terminal_Fleet_1 indicates that length comps are disaggregated across sex and fit using 'split' approach for all years for fleet 1--fixed gear)
                                            c("none_Year_1-terminal_Fleet_1",
                                              "spltRjntS_Year_1-terminal_Fleet_2")
)

#------------------- Setup survey indices and compositions
SPoRC_one_reg$UseSrvIdx[,,2] <- 0                                              # set use_trawl_index to 0, so it is not fit in the model
SPoRC_one_reg$UseSrvLenComps[,,2] <- 0                                         # same for length comps

input_list <- Setup_Mod_SrvIdx_and_Comps(input_list = input_list,
                                         ObsSrvIdx =                           # Observed fishery independent survey indices
                                           SPoRC_one_reg$ObsSrvIdx,
                                         ObsSrvIdx_SE =                        # Observed fishery independent survey standard errors
                                           (SPoRC_one_reg$ObsSrvIdx_SE/
                                              SPoRC_one_reg$ObsSrvIdx)*2,        # since moving to variance weighted likelihoods, need to convert to appropriate variance (not SE); generally LLS (JPN and DOM) CV ~ 0.05, TS ~ 0.1
                                         UseSrvIdx =                           # Whether to fit fishery independent surveys in each year
                                           SPoRC_one_reg$UseSrvIdx,
                                         ObsSrvAgeComps =                      # Observed age comps for fishery independent surveys
                                           SPoRC_one_reg$ObsSrvAgeComps,
                                         ISS_SrvAgeComps =                     # Input sample size for fishery independent survey age comps
                                           SPoRC_one_reg$ISS_SrvAgeComps,
                                         UseSrvAgeComps =                      # Whether to fit fishery independent survey age comps in a given year
                                           SPoRC_one_reg$UseSrvAgeComps,
                                         ObsSrvLenComps =                      # Observed fishery independent survey length comps
                                           SPoRC_one_reg$ObsSrvLenComps,
                                         UseSrvLenComps =                      # Whether to fit fishery independent survey length comps in a given year
                                           SPoRC_one_reg$UseSrvLenComps,
                                         ISS_SrvLenComps =                     # Input sample size for fishery independent survey length comps
                                           SPoRC_one_reg$ISS_SrvLenComps,

                                         # Model options
                                         srv_idx_type =                        # Define units (e.g. biomass vs. abundance/numbers) for fishery independent survey indices (none is used if no index for a given fleet)
                                           c("abd", "biom", "abd"),            # NOAA LLS uses RPNs, Trawl Survey uses biomass, and JPN LLS uses RPNs
                                         SrvAgeComps_LikeType =                # age comp likelihoods by survey fleet (none if not fit or not available for given fleet)
                                           c("Multinomial", "none", "Multinomial"),
                                         SrvLenComps_LikeType =                # length comp likelihoods for fishery survey fleets
                                           c("none", "none", "none"),
                                         SrvAgeComps_Type =                    # age comp structure for survey fleets (agg_Year_1-terminal_Fleet_1 indicates that age comps are aggregated across sex for all years for fleet 1--NOAA LLS)
                                           c("spltRjntS_Year_1-terminal_Fleet_1",
                                             "none_Year_1-terminal_Fleet_2",
                                             "spltRjntS_Year_1-terminal_Fleet_3"),
                                         SrvLenComps_Type =                    # length comp structure for survey fleets (spltRspltS_Year_1-terminal_Fleet_1 indicates that length comps are disaggregated across sex and fit using 'split' approach for all years for fleet 1--NOAA LLS)
                                           c("none_Year_1-terminal_Fleet_1",
                                             "none_Year_1-terminal_Fleet_2",
                                             "none_Year_1-terminal_Fleet_3")
)


#------------------- Setup up fishery selectivity
# Setup priors for selectivity
# Define valid fleet-block combinations
fleet_blocks <- data.frame(
  fleet = c(1, 1, 2),                                                       # fleets corresponding to time blocks (3 fixed gear time blocks, no trawl gear blocks)
  block = c(1, 2, 1)                                                        # corresponding time blocks
)

# Define sex and parameter combinations
sex_par <- expand.grid(sex = 1:2, par = 1:2)

# Merge to get all valid combinations
fish_selex_structure <- merge(fleet_blocks, sex_par) %>%
  dplyr::filter(!(fleet == 2 & block == 1 & sex == 2 & par == 1))              # remove priors for any unestimated pars -- par1=a50, par2=delta; NEEDS TO MATCH PARAMETER MAPPING

# Add the lognormal prior values - creates a dataframe, each row is a unique parameter combination to apply the prior to
fish_selex_prior <- cbind(
  region = 1,
  fish_selex_structure,
  mu = 2,                                                                      # All selex means = 1 (means should be defined in normal space)
  sd = 2                                                                       # All selex sd = 5
)

fish_selex_prior_tf <- fish_selex_prior %>%                                    # set tighter selex prior for TF
  dplyr::filter((fleet == 2 & par == 1)) %>%
  dplyr::mutate(mu = 5, sd = 1) %>%
  dplyr::full_join(fish_selex_prior %>%  dplyr::filter(!(fleet == 2 & par == 1 )))

fish_selex_prior_tf <- fish_selex_prior_tf %>%                                    # set tighter selex prior for TF
  dplyr::filter((fleet == 2 & par == 2)) %>%
  dplyr::mutate(mu = 8, sd = 2) %>%
  dplyr::full_join(fish_selex_prior_tf %>%  dplyr::filter(!(fleet == 2 & par == 2)))

input_list <- Setup_Mod_Fishsel_and_Q(input_list = input_list,
                                      cont_tv_fish_sel =                       # whether continuous time-varying selectivity is used for either fishery fleet
                                        c("none_Fleet_1", "none_Fleet_2"),
                                      fish_sel_blocks =                        # fishery selectivity time blocks if not TV specified above for a given fleet
                                        c("Block_1_Year_1-56_Fleet_1",         # 1960 - 2015 time block
                                          "Block_2_Year_57-terminal_Fleet_1",  # Recent time block for fixed gear fishery--2016 to terminal year
                                          "none_Fleet_2"),                     # no blocks for trawl fishery
                                      fish_sel_model =                         # fishery selectivity form
                                        c("logist1_Fleet_1",                   # logistic selectivity for fixed gear fleet
                                          "gamma_Fleet_2"),                    # gamma function selectivity for trawl gear fleet
                                      fish_q_blocks =                          # fishery indices catchability time blocks
                                        c("Block_1_Year_1-35_Fleet_1",         # pre-IFQ time block for fixed gear fishery 1994 and before
                                          "Block_2_Year_36-56_Fleet_1",        # IFQ time block for fixed gear fishery-- 1995 to 2015
                                          "Block_3_Year_57-terminal_Fleet_1",  # Recent time block for fixed gear fishery--2016 to terminal year
                                          "none_Fleet_2"),                     # no blocks for trawl fishery since no CPUE index used
                                      fish_fixed_sel_pars =                    # Whether to estimate all fixed effects for fishery selectivity, but can later modify to fix and share parameters via the parameter mapping
                                        c("est_all", "est_all"),               # Estimate all selectivity parameters for both fishery fleets
                                      fish_q_spec =                            # Whether to estimate all fixed effects for fishery catchability
                                        c("est_all", "fix"),                   # Estimate fishery q for fixed gear, not for trawl (no CPUE index)
                                      Use_fish_selex_prior = 1,                # Using selex priors
                                      fish_selex_prior = fish_selex_prior_tf
)


#------------------- Setup survey selectivity
# Setup priors for selectivity
# Define valid fleet-block combinations
fleet_blocks <- data.frame(
  fleet = c(1, 1, 2, 3),                                                       # fleets for each time block (LLS has 2 timeb blocks TS and JPN have 1 each)
  block = c(1, 2, 1, 1)                                                        # corresponding blocks
)

# Define sex and parameter combinations
sex_par <- expand.grid(sex = 1:2, par = 1:2)

# Merge to get all valid combinations
srv_selex_structure <- merge(fleet_blocks, sex_par) %>%
  dplyr::filter(!(fleet == 3 & block == 1 & sex == 2 & par == 2)) %>%              # remove priors for any unestimated pars -- par1=a50, par2=delta; NEEDS TO MATCH PARAMETER MAPPING
  dplyr::filter(!(fleet == 2))                                   # remove priors for any unestimated pars -- par1=a50, par2=delta; NEEDS TO MATCH PARAMETER MAPPING

# Add the lognormal prior values - creates a dataframe, each row is a unique parameter combination to apply the prior to
srv_selex_prior <- cbind(
  region = 1,
  srv_selex_structure,
  mu = 2,                                                                      # All selex means = 1 (means should be defined in normal space)
  sd = 2                                                                       # All selex sd = 5
)

input_list <- Setup_Mod_Srvsel_and_Q(input_list = input_list,
                                     cont_tv_srv_sel =                         # whether continuous time-varying selectivity is used for either survey fleet
                                       c("none_Fleet_1",
                                         "none_Fleet_2",
                                         "none_Fleet_3"),
                                     srv_sel_blocks =                          # survey selectivity time blocks if not TV specified above for a given fleet
                                       c("Block_1_Year_1-56_Fleet_1",          # Early time block LLS-- 1960 to 2016
                                         "Block_2_Year_57-terminal_Fleet_1",   # Recent time block for LLS--2017 to terminal year
                                         "none_Fleet_2",                       # No blocks for trawl survey
                                         "none_Fleet_3"                        # No blocks for JPN LLS
                                       ),
                                     srv_sel_model =                           # Survey selectivity form
                                       c("logist1_Fleet_1",                    # logistic selectivity for LLS
                                         "exponential_Fleet_2",                # exponential declining selectivity for trawl survey
                                         "logist1_Fleet_3"),                   # logistic selectivity for JPN LLS
                                     srv_q_blocks =                            # survey indices catchability time blocks (no q time blocks for any fleets)
                                       c("none_Fleet_1",
                                         "none_Fleet_2",
                                         "none_Fleet_3"),
                                     srv_fixed_sel_pars_spec =                 # Whether to estimate all fixed effects for survey selectivity, but can later modify to fix and share parameters via the parameter mapping
                                       c("est_all",                            # Estimate all selectivity parameters for all survey fleets
                                         "fix",
                                         "est_all"),
                                     srv_q_spec =                              # Whether to estimate all fixed effects for survey catchability
                                       c("est_all",                            # Estimate q for all survey fleets
                                         "fix",
                                         "est_all"),
                                     Use_srv_selex_prior = 1,                 # Using selex priors
                                     srv_selex_prior = srv_selex_prior
)

Wt_FishAgeComps <- array(1, dim = c(SPoRC_one_reg$nreg,
                                    SPoRC_one_reg$nyrs,
                                    SPoRC_one_reg$nsex,
                                    SPoRC_one_reg$nflt))

Wt_FishLenComps <- array(1, dim = c(SPoRC_one_reg$nreg,
                                    SPoRC_one_reg$nyrs,
                                    SPoRC_one_reg$nsex,
                                    SPoRC_one_reg$nflt))

Wt_SrvAgeComps <- array(1, dim = c(SPoRC_one_reg$nreg,
                                   SPoRC_one_reg$nyrs,
                                   SPoRC_one_reg$nsex,
                                   SPoRC_one_reg$nsrv_flt))

Wt_SrvLenComps <- array(1, dim = c(SPoRC_one_reg$nreg,
                                   SPoRC_one_reg$nyrs,
                                   SPoRC_one_reg$nsex,
                                   SPoRC_one_reg$nsrv_flt))

#------------------- Setup data weighting
input_list <- Setup_Mod_Weighting(input_list = input_list,
                                  Wt_Catch = 1,                                # Weight (lambda) for catch data
                                  Wt_FishIdx = 1,                              # Weight (lambda) for fishery indices data
                                  Wt_SrvIdx = 1,                               # Weight (lambda) for survey indices data
                                  Wt_Rec = 1,                                  # Weight (lambda) for stock-recruit penalty term
                                  Wt_F = 1,                                    # Weight (lambda) for fishing mortality penalty term
                                  Wt_Tagging = 0,                              # Weight (lambda) for tagging data (not fit in single region model)
                                  Wt_FishAgeComps = Wt_FishAgeComps,           # Weights (lambda) for fishery age comp data (note these are updated by Francis reweighting and input at start of this script)
                                  Wt_FishLenComps = Wt_FishLenComps,           # Weights (lambda) for fishery length comp data (note these are updated by Francis reweighting and input at start of this script)
                                  Wt_SrvAgeComps = Wt_SrvAgeComps,             # Weights (lambda) for survey age comp data (note these are updated by Francis reweighting and input at start of this script)
                                  Wt_SrvLenComps = Wt_SrvLenComps              # Weights (lambda) for survey length comp data (note these are updated by Francis reweighting and input at start of this script)
)

input_list$map$ln_fish_fixed_sel_pars <- factor(c(1:8,                         # sharing logistic delta across sexes from early LLF fishery (first time block), all other pars estimated for LLF across all time blocks
                                                  rep(9:10,2),                # for trawl fishery estimate both female gamma pars, but no time blocks (repeat 3 times to fill the 3 time blocks associated with the LLF)
                                                  rep(c(9,11),2)))            # for trawl fishery males, share the bmax par of the gamma function with the female estimate par (repeat 3 times to fill the 3 time blocks associated with the LLF)

# Share JPN sex-delta
input_list$map$ln_srv_fixed_sel_pars <-
  factor(c(1:8,                                                                # domestic ll survey (2 time blocks, 2 sexes),
           rep(NA,4),                                                          # domestic trawl survey (single time block) and only one parameter (NOT ESTIMATED) (exponential selex fxn), so only one parameter estimated per sex across blocks (so repeat value 4 times for each sex)
           rep(NA, 4),                                                         # domestic trawl survey (single time block) and only one parameter (NOT ESTIMATED) (exponential selex fxn), so only one parameter estimated per sex across blocks (so repeat value 4 times for each sex)
           rep(c(9:10),2),                                                     # coop JPN LL survey (single time block, females)
           rep(c(11:10),2)))                                                   # coop JPN LL survey (single time block, males), share delta

# Selex starting values
# input_list$par$ln_srv_fixed_sel_pars[1,,,,] <- log(2)
# input_list$par$ln_fish_fixed_sel_pars[1,,,,1] <- log(2)
input_list$par$ln_fish_fixed_sel_pars[1,,2,,2] <- log(2)

# Fit Model ---------------------------------------------------------------

data <- input_list$data
parameters <- input_list$par
mapping <- input_list$map
parameters$ln_global_R0 <- log(20)
data$ISS_FishLenComps[] <- 5

# run model
# twod_obj <- fit_model(data,
#                       parameters,
#                       mapping,
#                       random = NULL,
#                       newton_loops = 3,
#                       silent = FALSE
#
# )

# twod_obj$sd_rep <- RTMB::sdreport(twod_obj)

devtools::load_all(here("R"))
# data$FishAgeComps_LikeType <- c(5,999)
data$SrvAgeComps_LikeType <- c(5,999,2)

data$B <- 60
n_fact <- 1:10
chol <- list()

for(i in 1:length(n_fact)) {

  data$n_factors <- n_fact[i]
  n_L_params <- data$n_factors * (data$B - 1) - data$n_factors * (data$n_factors - 1) / 2
  parameters$Lvec <- rep(0.01, n_L_params)
  parameters$dvec <- rep(log(0.05), 1)
  # mapping$ln_FishAge_theta <- factor(rep(NA, 4))
  # mapping$FishAge_corr_pars <- factor(rep(NA, 8))
  # mapping$Lvec <- factor(rep(seq(1,n_L_params/3, by = 1), each = 3))
  # mapping$dvec <- factor(NA)
  # mapping$ln_SrvAge_theta <- factor(c(rep(NA,4), 1, NA))
  # mapping$SrvAge_corr_pars <- factor(rep(NA, 12))

  devtools::load_all(here("R"))

  # run model
  chol_obj <- fit_model(data,
                        parameters,
                        mapping,
                        random = NULL,
                        newton_loops = 3,
                        silent = FALSE

  )

  chol_obj$sd_rep <- RTMB::sdreport(chol_obj)
  chol[[i]] <- chol_obj
}

reshape2::melt(chol_obj$rep$Sigma) %>%
  mutate(type = 'chol') %>%
  # bind_rows(
  #   reshape2::melt(twod_obj$rep$Sigma) %>%
  #     mutate(type = '2d')
  # ) %>%
  filter(
        # Var1 != Var2,
         Var1 <= 60, Var2 <= 60
         ) %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  facet_wrap(~type) +
  scale_fill_viridis_c()

plot(as.vector(chol_obj$rep$Sigma), type = 'l')
lines(as.vector(twod_obj$rep$Sigma), type = 'l')

plot(as.vector(chol_obj$rep$Sigma) , as.vector(twod_obj$rep$Sigma))

marg_AIC(chol_obj$optim)
marg_AIC(twod_obj$optim)



# c(1,2,3,4,6,8,9)
marg_AIC(chol[[1]]$optim) # 7904.747
marg_AIC(chol[[2]]$optim) # 7557.935
marg_AIC(chol[[3]]$optim) # 7548.762
marg_AIC(chol[[4]]$optim) # 7570.157
marg_AIC(chol[[5]]$optim) # 7570.157

marg_AIC(chol[[6]]$optim) # 7532.734
marg_AIC(chol[[7]]$optim) # 7532.734

marg_AIC(chol[[8]]$optim) # 7460.594
marg_AIC(chol[[9]]$optim) # 7448
marg_AIC(chol[[10]]$optim) # 7448
marg_AIC(chol[[11]]$optim) # 7448
marg_AIC(chol[[12]]$optim) # 7448
marg_AIC(chol[[13]]$optim) # 7448
marg_AIC(chol[[14]]$optim) # 7448
marg_AIC(chol[[15]]$optim) # 7448

# OSA ---------------------------------------------------------------------

# extract out parameters for logistic normal
sigma <- exp(twod_obj$sd_rep$par.fixed[names(twod_obj$sd_rep$par.fixed) %in% "ln_FishAge_theta"][1])
corr_b <- rho_trans(twod_obj$sd_rep$par.fixed[names(twod_obj$sd_rep$par.fixed) %in% "FishAge_corr_pars"][1])
corr_s <- rho_trans(twod_obj$sd_rep$par.fixed[names(twod_obj$sd_rep$par.fixed) %in% "FishAge_corr_pars"][2])
ar1_mat <- SPoRC:::get_AR1_CorrMat(30, corr_b)
sex_mat <- SPoRC:::get_Constant_CorrMat(2, corr_s)
Sigma <- kronecker(sex_mat, ar1_mat) * sigma^2


nf <- 8
L = matrix(0, nrow = 59, ncol = nf)
dmat = matrix(0,59,59) # diagonal matrix

idx = 1
Lvec <- chol[[nf]]$optim$par[names(chol[[nf]]$optim$par) == 'Lvec']
for (j in 1:1) {
  for (i in j:(60-1)) { # fill only lower diag
    L[i, j] = Lvec[idx]
    idx = idx + 1
  } # end i
} # end j

dvec <- chol[[nf]]$optim$par[names(chol[[nf]]$optim$par) == 'dvec']
diag(dmat) = exp(dvec)^2
chol_sig = (L %*% t(L)) + dmat # lower triangular factorization

models <- list(chol_obj, twod_obj)
sigma <- list(chol_sig, Sigma)
for(i in 1:length(models)) {

  tmp_rep <- models[[i]]$rep
  # get composition proportions
  tmp <- get_comp_prop(data, tmp_rep, 2:31, seq(41,99,2), year_labels = 1960:2024)

  # get osas
  tmp_osas <- SPoRC::get_osa(obs_mat = tmp$Obs_FishAge_mat,
                             exp_mat = tmp$Pred_FishAge_mat,
                             N = NULL,
                             years = which(data$UseFishAgeComps[,,1] == 1),
                             fleet = 1,
                             bins = 2:31,
                             comp_type = 2,
                             bin_label = "Ages",
                             comp_like = 4,
                             LN_Sigma = array(sigma[[i]], c(1, dim(sigma[[i]]))))

  tmp_plots <- plot_resids(tmp_osas)

  # QQ plot
  tmp_qq <- tmp_plots[[1]] +
    facet_null()

  # Residual plot
  tmp_resids <- tmp_plots[[2]] +
    facet_wrap(~sex, labeller = labeller(sex = c("1" = "Female", "2" = "Male")))

  # aggregated fits
  tmp_agg <- tmp$Fishery_Ages %>%
    group_by(Fleet, Age, Sex) %>%
    summarize(pred = mean(pred), obs = mean(obs)) %>%
    filter(Fleet == 1) %>%
    ggplot() +
    geom_col(aes(x = Age, y = obs), fill = 'darkgreen', alpha = 0.5, color = 'black') +
    geom_line(aes(x = Age, y = pred), col = 'black', lwd = 1.3, lty = 2) +
    facet_wrap(~Sex, labeller = labeller(Sex = c("1" = "Female", "2" = "Male"))) +
    theme_sablefish() +
    labs(x = 'Age', y = 'Proportion')

  tmp_bottom <- cowplot::plot_grid(tmp_qq, tmp_agg, labels = c("B", "C"), label_size = 25)

  print(
    cowplot::plot_grid(tmp_resids, tmp_bottom, ncol = 1, rel_heights = c(0.6, 0.4),
                       labels = c("A", ""), label_size = 25)
  )

}

# if(Likelihood_Type == 5) {
#
#   # Build L matrix for B-1 (unconstrained part only)
#   L = matrix(0, nrow = B, ncol = n_factors)
#   dmat = matrix(0, B, B)
#   idx = 1
#   for (j in 1:n_factors) {
#     for (i in j:(B-1)) {
#       L[i, j] = exp( Lvec[idx])
#       idx = idx + 1
#     }
#   }
#   diag(dmat) = exp(dvec)
#   Sigma = (L %*% t(L)) + dmat # (B-1) x (B-1)
#
#   RTMB::REPORT(Sigma)
#
#   # Dealing with zeros - normalize full vector
#   tmp_Obs = Obs[r,,] / sum(Obs[r,,])
#   zeros = which(tmp_Obs == 0)
#   if(length(zeros) > 0) {
#     tmp_Obs = tmp_Obs[-zeros]
#     tmp_Exp = tmp_Exp[-zeros]
#     tmp_Obs = tmp_Obs / sum(tmp_Obs)
#     tmp_Exp = tmp_Exp / sum(tmp_Exp)
#     # Adjust Sigma
#     Sigma = Sigma[-zeros, -zeros]
#   }
#
#   Sigma = Sigma[-nrow(Sigma), -ncol(Sigma)] # remove last bins
#
#   comp_nLL[r,1] = -1 * dlogistnormal(obs = tmp_Obs, pred = tmp_Exp, Sigma_or_Q = Sigma, type = 'dmvnorm', TRUE) # Logistic Normal likelihood (1dar1 by age, constant corr by sex)
# }

# Lvec = Lvec,
# dvec = dvec,
# n_factors = n_factors,
# B = B
#
# sel_nLL = sel_nLL + sum(Lvec^2)
# sel_nLL = sel_nLL -RTMB::dnorm(dvec, log(0.2), 0.1, TRUE)
# sel_nLL = sel_nLL -sum(RTMB::dnorm(Lvec, log(0.1), 0.5, TRUE))

