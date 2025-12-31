# Purpose: To bridge to the 2024 GOA Pcod model (24.0)
# Creator: Matthew LH. Cheng
# 7/3/25

devtools::load_all(here('R'))

# Function to construct ALK given mean length and sd
get_al_trans_matrix = function(age_bins, len_bins, mean_length, sd) {
  # container array
  age_length = matrix(0.0, nrow = length(age_bins), ncol = length(len_bins))

  for(a in 1:length(age_bins)) {
    # Use actual bin lower limits as per SS3 specification
    # Assume len_bins contains the lower limits of each bin
    bin_lower_limits = len_bins

    # Calculate cumulative probabilities at bin lower limits
    AL = pnorm(bin_lower_limits, mean_length[a], sd[a])

    # Calculate bin probabilities according to the SS3 specification
    for(l in 1:length(len_bins)) {
      if(l == 1) {
        # First bin: from -infinity to lower limit of bin 2
        age_length[a, l] = pnorm(bin_lower_limits[2], mean_length[a], sd[a])
      } else if(l == length(len_bins)) {
        # Last bin: from lower limit of last bin to +infinity
        age_length[a, l] = 1 - pnorm(bin_lower_limits[l], mean_length[a], sd[a])
      } else {
        # Middle bins: from lower limit of bin l to lower limit of bin l+1
        age_length[a, l] = pnorm(bin_lower_limits[l+1], mean_length[a], sd[a]) -
          pnorm(bin_lower_limits[l], mean_length[a], sd[a])
      }
    }
  }
  return(age_length)
}

# Installation
# install.packages("devtools") # install dev tools
# install.packages("TMB") # install TMB
# install.packages("RTMB") # install RTMB
# TMB:::install.contrib("https://github.com/vtrijoulet/OSA_multivariate_dists/archive/main.zip") # get multivariate OSA distributions
#
# # optional packages to install
# devtools::install_github("fishfollower/compResidual/compResidual")
# devtools::install_github("noaa-afsc/afscOSA", dependencies = TRUE)
#
# # install SPoRC
# devtools::install_github("chengmatt/SPoRC", dependencies = TRUE)

# Setup -------------------------------------------------------------------

library(here)
library(SPoRC)
library(r4ss)
library(tidyverse)

rep <- r4ss::SS_read(here("dev", "dev_data", "goapcod"))
ss3_pars <- read.csv(here("dev", "dev_data", "goapcod", "model_parameters.csv"))
out <- r4ss::SS_output(here("dev", "dev_data", "goapcod"), verbose = TRUE)

# Extract Data ------------------------------------------------------------

# Setup dimensions
n_regions <- 1 # number of regions
n_sexes <- 1 # number of sexes
n_fish_fleets <- 3 # number of fishery fleets
n_srv_fleets <- 2 # number of survey fleets
years <- rep$dat$styr:rep$dat$endyr # years
n_yrs <- length(years) # number of years
ages <- 0:10 # vector of ages
n_ages <- length(ages) # number of ages
lens <- rep$dat$lbin_vector # vector of lengths
n_lens <- length(lens) # number of length bins
spwn_month <- rep$dat$spawn_month - 1 # spawning month (start of year)

# Fishery Data ------------------------------------------------------------

# fishery catches (fleets 1 - 3)
ss3_catch <- rep$dat$catch %>% dplyr::filter(year != -999) # get ss3 catches
ObsFishCatch <- array(0, dim = c(n_regions, n_yrs, n_fish_fleets))

# input catches
ObsFishCatch[1,,1] <- as.numeric(gsub(",", "", ss3_catch$catch[ss3_catch$fleet == 1])) # converting characters to numbers (Fleet 1)
ObsFishCatch[1,,2] <- as.numeric(gsub(",", "", ss3_catch$catch[ss3_catch$fleet == 2])) # converting characters to numbers (Fleet 2)
ObsFishCatch[1,,3] <- as.numeric(gsub(",", "", ss3_catch$catch[ss3_catch$fleet == 3])) # converting characters to numbers (Fleet 3)
Catch_Type <- array(1, dim = c(n_yrs, n_fish_fleets)) # Catch type (ignoreed in a single area case)

# Use catch indicator
UseCatch <- array(1, dim = c(n_regions, n_yrs, n_fish_fleets)) # Use catch indicator
UseCatch[1,which(ObsFishCatch[1,,1] == 0),1] <- 0 # don't fit if 0 catch
UseCatch[1,which(ObsFishCatch[1,,2] == 0),2] <- 0 # don't fit if 0 catch
UseCatch[1,which(ObsFishCatch[1,,3] == 0),3] <- 0 # don't fit if 0 catch

# fishery index (none used)
ObsFishIdx <- array(NA, dim = c(n_regions, n_yrs, n_fish_fleets))
ObsFishIdx_SE <- array(NA, dim = c(n_regions, n_yrs, n_fish_fleets))
UseFishIdx <- array(0, dim = c(n_regions, n_yrs, n_fish_fleets))

# fishery ages
ObsFishAgeComps <- array(0, dim = c(n_regions, n_yrs, n_ages-1, n_sexes, n_fish_fleets))
# note n_ages - 1 (the ageing error collapses this so the first age (age-0) gets ignored, similar to rectangular rockfish case)
UseFishAgeComps <- array(0, dim = c(n_regions, n_yrs, n_fish_fleets))
ISS_FishAgeComps <- array(0, dim = c(n_regions, n_yrs, n_sexes, n_fish_fleets))

# fleet 1 ages
ages_fleet1 <- rep$dat$agecomp %>% dplyr::filter(fleet == -1)
ages_fleet1_yrs <- ages_fleet1$year

# loop through to populate fleet 1
for(i in ages_fleet1_yrs) {
  ObsFishAgeComps[1,which(years == i),,1,1] <- unlist((ages_fleet1 %>% dplyr::filter(year == i))[,-c(1:9)]) # observed ages
  UseFishAgeComps[1,which(years == i),1] <- 1 # whether to use fishery ages
  ISS_FishAgeComps[1,which(years == i),1,1] <- (ages_fleet1 %>% dplyr::filter(year == i))$Nsamp
} # end i loop

# fleet 2 ages
ages_fleet2 <- rep$dat$agecomp %>% dplyr::filter(fleet == -2)
ages_fleet2_yrs <- ages_fleet2$year

# loop through to populate fleet 2
for(i in ages_fleet2_yrs) {
  ObsFishAgeComps[1,which(years == i),,1,2] <- unlist((ages_fleet2 %>% dplyr::filter(year == i))[,-c(1:9)]) # observed ages
  UseFishAgeComps[1,which(years == i),2] <- 1 # whether to use fishery ages
  ISS_FishAgeComps[1,which(years == i),1,2] <- (ages_fleet2 %>% dplyr::filter(year == i))$Nsamp
} # end i loop

# fleet 3 ages
ages_fleet3 <- rep$dat$agecomp %>% dplyr::filter(fleet == -3)
ages_fleet3_yrs <- ages_fleet3$year

# loop through to populate fleet 3
for(i in ages_fleet3_yrs) {
  ObsFishAgeComps[1,which(years == i),,1,3] <- unlist((ages_fleet3 %>% dplyr::filter(year == i))[,-c(1:9)]) # observed ages
  UseFishAgeComps[1,which(years == i),3] <- 1 # whether to use fishery ages
  ISS_FishAgeComps[1,which(years == i),1,3] <- (ages_fleet3 %>% dplyr::filter(year == i))$Nsamp
}

# fishery lengths
ObsFishLenComps <- array(0, dim = c(n_regions, n_yrs, n_lens, n_sexes, n_fish_fleets))
UseFishLenComps <- array(0, dim = c(n_regions, n_yrs, n_fish_fleets))
ISS_FishLenComps <- array(0, dim = c(n_regions, n_yrs, n_sexes, n_fish_fleets))

# fleet 1 lengths
lens_fleet1 <- rep$dat$lencomp %>% dplyr::filter(fleet == 1)
lens_fleet1_yrs <- lens_fleet1$year

# loop through to populate fleet 1
for(i in lens_fleet1_yrs) {
  ObsFishLenComps[1,which(years == i),,1,1] <- unlist((lens_fleet1 %>% dplyr::filter(year == i))[,-c(1:6)]) # observed lens
  UseFishLenComps[1,which(years == i),1] <- 1 # whether to use fishery lens
  ISS_FishLenComps[1,which(years == i),1,1] <- (lens_fleet1 %>% dplyr::filter(year == i))$Nsamp
}

# fleet 2 lengths
lens_fleet2 <- rep$dat$lencomp %>% dplyr::filter(fleet == 2)
lens_fleet2_yrs <- lens_fleet2$year

# loop through to populate fleet 2
for(i in lens_fleet2_yrs) {
  ObsFishLenComps[1,which(years == i),,1,2] <- unlist((lens_fleet2 %>% dplyr::filter(year == i))[,-c(1:6)]) # observed lens
  UseFishLenComps[1,which(years == i),2] <- 1 # whether to use fishery lens
  ISS_FishLenComps[1,which(years == i),1,2] <- (lens_fleet2 %>% dplyr::filter(year == i))$Nsamp
}

# fleet 3 lengths
lens_fleet3 <- rep$dat$lencomp %>% dplyr::filter(fleet == 3)
lens_fleet3_yrs <- lens_fleet3$year

# loop through to populate fleet 3
for(i in lens_fleet3_yrs) {
  ObsFishLenComps[1,which(years == i),,1,3] <- unlist((lens_fleet3 %>% dplyr::filter(year == i))[,-c(1:6)]) # observed lens
  UseFishLenComps[1,which(years == i),3] <- 1 # whether to use fishery lens
  ISS_FishLenComps[1,which(years == i),1,3] <- (lens_fleet3 %>% dplyr::filter(year == i))$Nsamp
}

# Survey Data -------------------------------------------------------------

# survey index
ss3_idx <- rep$dat$CPUE # get ss3 cpue
ObsSrvIdx <- array(0, dim = c(n_regions, n_yrs, n_srv_fleets))
ObsSrvIdx_SE <- array(0, dim = c(n_regions, n_yrs, n_srv_fleets))
UseSrvIdx <- array(0, dim = c(n_regions, n_yrs, n_srv_fleets))

# fleet 4
ObsSrvIdx[1,which(years %in% c(ss3_idx$year[ss3_idx$index == 4])),1] <- ss3_idx$obs[ss3_idx$index == 4] # survey index
ObsSrvIdx_SE[1,which(years %in% c(ss3_idx$year[ss3_idx$index == 4])),1] <- ss3_idx$se_log[ss3_idx$index == 4] # survey index se
UseSrvIdx[1,which(years %in% c(ss3_idx$year[ss3_idx$index == 4])),1] <- 1 # when to use survey index

# fleet 5
ObsSrvIdx[1,which(years %in% c(ss3_idx$year[ss3_idx$index == 5])),2] <- ss3_idx$obs[ss3_idx$index == 5] # survey index
ObsSrvIdx_SE[1,which(years %in% c(ss3_idx$year[ss3_idx$index == 5])),2] <- ss3_idx$se_log[ss3_idx$index == 5] # survey index se
UseSrvIdx[1,which(years %in% c(ss3_idx$year[ss3_idx$index == 5])),2] <- 1 # when to use survey index

# survey ages
ObsSrvAgeComps <- array(0, dim = c(n_regions, n_yrs, n_ages-1, n_sexes, n_srv_fleets))
# note n_ages - 1 (the ageing error collapses this so the first age (age-0) gets ignored, similar to rectangular rockfish case)
UseSrvAgeComps <- array(0, dim = c(n_regions, n_yrs, n_srv_fleets))
ISS_SrvAgeComps <- array(0, dim = c(n_regions, n_yrs, n_sexes, n_srv_fleets))

# fleet 4 ages
ages_srv4 <- rep$dat$agecomp %>% dplyr::filter(fleet == -4) # get ss3 ages
ages_srv4_yrs <- ages_srv4$year

for(i in ages_srv4_yrs) {
  ObsSrvAgeComps[1,which(years == i),,1,1] <- unlist((ages_srv4 %>% dplyr::filter(year == i))[,-c(1:9)]) # observed ages
  UseSrvAgeComps[1,which(years == i),1] <- 1 # whether to use fishery ages
  ISS_SrvAgeComps[1,which(years == i),1,1] <- (ages_srv4 %>% dplyr::filter(year == i))$Nsamp
} # end i loop

# fleet 5 (no ages for survey fleet 5)

# survey lengths
ObsSrvLenComps <- array(0, dim = c(n_regions, n_yrs, n_lens, n_sexes, n_srv_fleets))
UseSrvLenComps <- array(0, dim = c(n_regions, n_yrs, n_srv_fleets))
ISS_SrvLenComps <- array(0, dim = c(n_regions, n_yrs, n_sexes, n_srv_fleets))

# fleet 4 lengths
lens_fleet4 <- rep$dat$lencomp %>% dplyr::filter(fleet == 4)
lens_fleet4_yrs <- lens_fleet4$year

# loop through to populate fleet 4
for(i in lens_fleet4_yrs) {
  ObsSrvLenComps[1,which(years == i),,1,2] <- unlist((lens_fleet4 %>% dplyr::filter(year == i))[,-c(1:6)]) # observed lens
  UseSrvLenComps[1,which(years == i),2] <- 1 # whether to use fishery lens
  ISS_SrvLenComps[1,which(years == i),1,2] <- (lens_fleet4 %>% dplyr::filter(year == i))$Nsamp
}

# fleet 5
lens_fleet5 <- rep$dat$lencomp %>% dplyr::filter(fleet == 5)
lens_fleet5_yrs <- lens_fleet5$year

# loop through to populate fleet 5
for(i in lens_fleet5_yrs) {
  ObsSrvLenComps[1,which(years == i),,1,1] <- unlist((lens_fleet5 %>% dplyr::filter(year == i))[,-c(1:6)]) # observed lens
  UseSrvLenComps[1,which(years == i),1] <- 1 # whether to use fishery lens
  ISS_SrvLenComps[1,which(years == i),1,1] <- (lens_fleet5 %>% dplyr::filter(year == i))$Nsamp
}

# Demographics ------------------------------------------------------------

### Weight at age -----------------------------------------------------------
waa_arr <- array(0, dim = c(n_regions, n_yrs, n_ages, n_sexes))
for(i in 1:length(years))  waa_arr[1,i,,1] <- out$endgrowth$Wt_Beg

### Maturity at age ---------------------------------------------------------
mataa_arr <- array(rep(out$endgrowth$Len_Mat, each = n_yrs),
                   dim = c(n_regions, n_yrs, n_ages, n_sexes))

### Size Age Transition -----------------------------------------------------
laa <- out$endgrowth$Len_Mid # length at age
sd <- out$endgrowth$SD_Mid # sd of length at age

# construct alk
alk <- get_al_trans_matrix(age_bins = ages,
                           len_bins  = lens,
                           mean_length =  as.numeric(unlist(laa)),
                           sd = as.numeric(unlist(sd))
)

# normalize ALK in case it doesnt sum to 1
alk <- t(apply(alk, 1, function(x) x / sum(x)))

# populate size age transition matrix
sizeage <- array(0, dim = c(n_regions, n_yrs, n_lens, n_ages, n_sexes))
for(i in 1:n_yrs) sizeage[1,i,,,1] <- t(alk)


### Ageing Error ------------------------------------------------------------
ageerror <- out$AAK[2,,] # Is this ageing error + bias? (right now, no time-varying ageing stuff allowed)
ageerror <- t(ageerror[nrow(ageerror):1,]) # need to flip the matrix and also transpose

ageerror2 <- out$AAK[1,,] # Is this ageing error + bias? (right now, no time-varying ageing stuff allowed)
ageerror2 <- t(ageerror2[nrow(ageerror2):1,]) # need to flip the matrix and also transpose

# time-varying ageing error
ageerror_t <- array(0, dim = c(n_yrs, dim(ageerror2)))
for(i in 1:n_yrs) {
  if(i < length(1977:2007)) ageerror_t[i,,] <- ageerror2
  else ageerror_t[i,,] <- ageerror
}

### Natural Mortality -------------------------------------------------------
# Setup fixed natural mortality array first
fixed_natmort <- array(0, dim = c(n_regions, n_yrs, n_ages, n_sexes))
fixed_natmort[] <- 0.493858 # natural mortality in regular period
fixed_natmort[,which(years %in% c(2014:2016)),,] <- 0.819832 # time blocked natural mortality

# Setup Model -------------------------------------------------------------

# Setup model dimensions
input_list <- Setup_Mod_Dim(
  years = years,
  # vector of years
  ages = ages,
  # vector of ages
  lens = lens,
  # number of lengths
  n_regions = n_regions,
  # number of regions
  n_sexes = n_sexes,
  # number of sexes
  n_fish_fleets = n_fish_fleets,
  # number of fishery fleets
  n_srv_fleets = n_srv_fleets, # number of survey fleets
  verbose = TRUE # whether to output messages
)

# Setup recruitment stuff (using defaults for other stuff)
input_list <- Setup_Mod_Rec(
  input_list = input_list,

  # Model options
  do_rec_bias_ramp = 1, # do bias ramp (0 == don't do bias ramp, 1 == do bias ramp)
  bias_year = c(-16, length(years[1]:1989), length(years[1]:2021), length(years[1]:2024)), # recruitment bias ramp (according to SS3 specifications, full bias correction until 2018) # out$recruit$biasadjuster
  max_bias_ramp_fct = max(out$recruit$biasadjuster), # maximum bias factor to apply to ramp
  sigmaR_switch = 1, # when to switch from early to late sigmaR (set at 1 to switch in first yearl i.e., not used)
  rec_model = "mean_rec",   # recruitment model
  sigmaR_spec = "fix", # fix early sigmaR and late sigmaR
  init_age_strc = 1, # geometric series to derive initial age structure
  equil_init_age_strc = 2, # estimate all age deviations inlcluding for plus group

  # Fixed values
  ln_sigmaR = c(log(0.44), log(0.44)), # 2 values for early and late sigma (early sigma is not used); starting values / fixed valuesfor early and late sigmaR
  t_spawn = spwn_month # spawning month

)

# Setup biologicals
input_list <- Setup_Mod_Biologicals(
  input_list = input_list,

  # Model options
  fit_lengths = 1, # fit lengths
  Selex_Type = "length", # length-based selectivity
  M_spec = "fix", # fixing natural mortality
  addtocomp = 1e-4,

  # Data inputs
  WAA = waa_arr,
  MatAA = mataa_arr,
  SizeAgeTrans = sizeage, # size age transition
  AgeingError = ageerror_t, # ageing error matrix
  Fixed_natmort = fixed_natmort
)

# Setup movement stuff (using defaults for other stuff) (not used)
input_list <- Setup_Mod_Movement(
  input_list = input_list,
  use_fixed_movement = 1,
  Fixed_Movement = NA,
  do_recruits_move = 0
)

# set up tagging stuff (not used)
input_list <- Setup_Mod_Tagging(input_list = input_list, UseTagging = 0)

# setup catch and fishing mortality stuff
input_list <- Setup_Mod_Catch_and_F(
  input_list = input_list,

  # Data inputs
  ObsCatch = ObsFishCatch,
  Catch_Type = Catch_Type,
  UseCatch = UseCatch,

  # Model options
  Use_F_pen = 1, # whether to use f penalty, == 0 don't use, == 1 use
  sigmaC_spec = "fix", # fixing sigma catch

  # starting values / fixed values
  ln_sigmaC = array(log(0.01), dim = c(n_regions, n_yrs, n_fish_fleets)), # sigma catch (basically known)
  ln_sigmaF = array(log(1), dim = c(n_regions, n_fish_fleets)) # sigma of f penalty(set large because catch is known)
)

# setup fishery index and comps stuff
input_list <- Setup_Mod_FishIdx_and_Comps(
  input_list = input_list,

  # data inputs
  ObsFishIdx = ObsFishIdx,
  ObsFishIdx_SE = ObsFishIdx_SE,
  UseFishIdx = UseFishIdx,
  ObsFishAgeComps = ObsFishAgeComps,
  UseFishAgeComps = UseFishAgeComps,
  ISS_FishAgeComps = ISS_FishAgeComps,
  ObsFishLenComps = ObsFishLenComps,
  UseFishLenComps = UseFishLenComps,
  ISS_FishLenComps = ISS_FishLenComps,

  # Model options
  # indices for fishery
  fish_idx_type = c("none",
                    "none",
                    "none"),

  # age comp likelihoods for fishery fleet
  FishAgeComps_LikeType = c("Multinomial",
                            "Multinomial",
                            "Multinomial"),

  # len comp likelihoods for fishery fleet
  FishLenComps_LikeType = c("Multinomial",
                            "Multinomial",
                            "Multinomial"),

  # aggregated age comps (b/c only single sex)
  FishAgeComps_Type = c("agg_Year_1-terminal_Fleet_1",
                        "agg_Year_1-terminal_Fleet_2",
                        "agg_Year_1-terminal_Fleet_3"),

  # aggregated length comps (b/c only single sex)
  FishLenComps_Type = c("agg_Year_1-terminal_Fleet_1",
                        "agg_Year_1-terminal_Fleet_2",
                        "agg_Year_1-terminal_Fleet_3"),
  FishAge_comp_agg_type = c(1, 1, 1),
  # ADMB aggregation quirks, ideally get rid of this (just order of operations from apply ageing error, normalize, etc)
  FishLen_comp_agg_type = c(0, 0, 0)
  # ADMB aggregation quirks, ideally get rid of this
)

# setup survey index and comps stuff
input_list <- Setup_Mod_SrvIdx_and_Comps(
  input_list = input_list,

  # data inputs
  ObsSrvIdx = ObsSrvIdx,
  ObsSrvIdx_SE = ObsSrvIdx_SE,
  UseSrvIdx = UseSrvIdx,
  ObsSrvAgeComps = ObsSrvAgeComps,
  ISS_SrvAgeComps = ISS_SrvAgeComps,
  UseSrvAgeComps = UseSrvAgeComps,
  ObsSrvLenComps = ObsSrvLenComps,
  UseSrvLenComps = UseSrvLenComps,
  ISS_SrvLenComps = ISS_SrvLenComps,

  # Model options
  # abundance and biomass for survey fleet 1 and 2, respectively
  srv_idx_type = c("biom",
                   "abd"),

  # survey age composition likelihood for survey fleet 1, none used for survey fleet 2
  SrvAgeComps_LikeType = c("Multinomial",
                           "none"),

  # survey len composition likelihood for survey fleet 1, none used for survey fleet 2
  SrvLenComps_LikeType = c("Multinomial",
                           "Multinomial"),

  # survey age aggregation type
  SrvAgeComps_Type = c("agg_Year_1-terminal_Fleet_1",
                       "none_Year_1-terminal_Fleet_2"),

  # survey length comp type
  SrvLenComps_Type = c("agg_Year_1-terminal_Fleet_1",
                       "agg_Year_1-terminal_Fleet_2"),

  SrvAge_comp_agg_type = c(1, 1),
  # ADMB aggregation quirks, ideally get rid of this
  SrvLen_comp_agg_type = c(0, 0)
  # ADMB aggregation quirks, ideally get rid of this
)


# setup fishery selectivity and catchability
input_list <- Setup_Mod_Fishsel_and_Q(

  input_list = input_list,

  # Model options
  # fishery selectivity, whether continuous time-varying
  cont_tv_fish_sel = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3"),

  # fishery selectivity blocks (setup to mimic assessment)
  fish_sel_blocks = c("Block_1_Year_1-13_Fleet_1", # longline (1977 - 1989)
                      "Block_2_Year_14-terminal_Fleet_1", # longline (1990 - temrinal)
                      "Block_1_Year_1-13_Fleet_2", # trawl (1977 - 1989)
                      "Block_2_Year_14-28_Fleet_2", # trawl (1990 - 2004)
                      "Block_3_Year_29-30_Fleet_2", # trawl (2005 - 2006); this should probably be removed and subsumed into another block
                      "Block_4_Year_31-40_Fleet_2", # trawl (2007 - 2016)
                      "Block_5_Year_41-terminal_Fleet_2", # trawl (2017 - terminal)
                      "Block_1_Year_1-36_Fleet_3", # pot (1977 - 2012)
                      "Block_2_Year_37-terminal_Fleet_3" # pot (2012 - terminal)
  ),

  # fishery selectivity form
  fish_sel_model = c("logist1_Fleet_1", # logistic a50 and slope
                     "logist1_Fleet_2",
                     "gamma_Fleet_3"),

  # fishery catchability blocks (none b/c no fishery idx)
  fish_q_blocks = c("none_Fleet_1",
                    "none_Fleet_2",
                    "none_Fleet_3"),

  # whether to estiamte all fixed effects for fishery selectivity
  fish_fixed_sel_pars_spec = c("est_all",
                                "est_all",
                                "est_all"),

  # whether to estiamte all fixed effects for fishery catchability (fix b/c no fishery idx)
  fish_q_spec = c("fix",
                  "fix",
                  "fix")
)

# Setup survey selectivity and catchability
input_list <- Setup_Mod_Srvsel_and_Q(

  input_list = input_list,

  # Model options
  # survey selectivity, whether continuous time-varying
  cont_tv_srv_sel = c("none_Fleet_1",
                      "none_Fleet_2"),

  # survey selectivity blocks
  srv_sel_blocks = c("Block_1_Year_1-19_Fleet_1",
                     "Block_1_Year_20-30_Fleet_1",
                     "Block_1_Year_31-terminal_Fleet_1",
                     "none_Fleet_2"),

  # survey selectivity form
  srv_sel_model = c("logist1_Fleet_1",
                    "logist1_Fleet_2"),

  # survey catchability blocks
  srv_q_blocks = c("none_Fleet_1",
                   "none_Fleet_2"),

  # whether to estiamte all fixed effects for survey selectivity
  srv_fixed_sel_pars_spec = c("est_all",
                              "est_all"),

  # whether to estiamte all fixed effects for survey catchability
  srv_q_spec = c("est_all",
                 "est_all"),

  # Set up formula list
  srv_q_formula = list(
    "Region_1_Fleet_1" = ~ NULL,
    "Region_1_Fleet_2" = ~ 0 + temp
  ),

  srv_q_cov_dat = list(
    temp = c(0, 0, rep$dat$envdat$value)
  )

)

# Setup model weighting (empahsis factors / lambdas for SS3)
input_list <- Setup_Mod_Weighting(
  input_list = input_list,
  Wt_Catch = 1,
  Wt_FishIdx = 0,
  Wt_SrvIdx = 1,
  Wt_Rec = 1,
  Wt_F = 1,
  Wt_Tagging = 0,
  Wt_FishAgeComps = array(1, dim = c(input_list$data$n_regions,
                                     length(input_list$data$years),
                                     input_list$data$n_sexes,
                                     input_list$data$n_fish_fleets)),
  Wt_FishLenComps = array(1, dim = c(input_list$data$n_regions,
                                     length(input_list$data$years),
                                     input_list$data$n_sexes,
                                     input_list$data$n_fish_fleets)),
  Wt_SrvAgeComps = array(1, dim = c(input_list$data$n_regions,
                                    length(input_list$data$years),
                                    input_list$data$n_sexes,
                                    input_list$data$n_srv_fleets)),
  Wt_SrvLenComps = array(1, dim = c(input_list$data$n_regions,
                                    length(input_list$data$years),
                                    input_list$data$n_sexes,
                                    input_list$data$n_srv_fleets))
)

# extract data, parameters, and mapping from input_list (via the setup functions defined above)
data <- input_list$data
parameters <- input_list$par
mapping <- input_list$map


# Setup more starting values ----------------------------------------------

# Starting values to help out the model
parameters$ln_global_R0 <- ss3_pars$Value[ss3_pars$Label == 'SR_LN(R0)'] # mean recruitment
parameters$ln_srv_q[] <- c(ss3_pars$Value[ss3_pars$Label == 'LnQ_base_Srv(4)'], ss3_pars$Value[ss3_pars$Label == 'LnQ_base_LLSrv(5)']) # survey qs
parameters$ln_InitDevs[] <- rev(ss3_pars$Value[str_detect(ss3_pars$Label, "Early_Init")]) # early initial age deviations
parameters$ln_RecDevs[] <- ss3_pars$Value[str_detect(ss3_pars$Label, "RecrDev")] # recruitment age deviations

# setup starting values for mean fishing mortality
ss3_f <- out$timeseries %>% filter(Era == 'TIME')
parameters$ln_F_mean[] <- c(mean(log(ss3_f$`F:_1`[ss3_f$`F:_1` != 0])),
                            mean(log(ss3_f$`F:_2`[ss3_f$`F:_2` != 0])),
                            mean(log(ss3_f$`F:_3`[ss3_f$`F:_3` != 0]))) # log mean F by fleet
parameters$ln_F_devs[,,1] <- log(ss3_f$`F:_1`) - parameters$ln_F_mean[,1] # f devs for fleet 1
parameters$ln_F_devs[,,2] <- log(ss3_f$`F:_2`) - parameters$ln_F_mean[,2] # f devs for fleet 2
parameters$ln_F_devs[,,3] <- log(ss3_f$`F:_3`) - parameters$ln_F_mean[,3] # f devs for fleet 3
parameters$ln_F_devs[!is.finite(parameters$ln_F_devs)] <- 0 # don't estimate f devs when catch = 0

# Setup starting values for selectivity

# Notes: For these fixed effects selex arrays,
# its dimensioned by: n_regions, n_max_sel_pars (maximum selectivity parameters for a given functional form defined across fleets,
# n_max_sel_blocks, n_sexes, n_fleets)
# n_max_sel_pars = If you have 3 fleets, fleet 1 = logistic, fleet 2 = logistic, fleet 3 = double normal, this dimension will always be 6
# n_max_sel_blocks = If you have 5 blocks for fleet 1, 3 blocks for fleet 2, and 1 block for fleet 3, this dimension will be 5.

# parameters$ln_fish_fixed_sel_pars[,1,,,1:2] <- log(60) # l50 parameter for logistic, fleet 1 and 2
# parameters$ln_fish_fixed_sel_pars[,2,,,1:2] <- log(0.5) # slope parameter for logistic, fleet 1 and 2
# parameters$ln_fish_fixed_sel_pars[,1,,,3] <- log(50) # lmax parameter for gamma, fleet 3
# parameters$ln_fish_fixed_sel_pars[,2,,,3] <- log(10) # slope for gamma fleet 3
#
# parameters$ln_srv_fixed_sel_pars[,1,,,] <- log(50) # lmax parameter for gamma, fleet 3
# parameters$ln_srv_fixed_sel_pars[,2,,,] <- log(0.5) # slope for gamma fleet 3


# Fit Model ---------------------------------------------------------------

# Fit model here
goapcod_obj <- SPoRC::fit_model(data = data,
                                parameters = parameters,
                                mapping = mapping,
                                random = NULL,
                                newton_loops = 3,
                                silent = FALSE)

goapcod_obj$sdrep <- RTMB::sdreport(goapcod_obj) # get standard error report


# Comparisons -------------------------------------------------------------

# Compare bias ramp differences
plot(goapcod_obj$rep$bias_ramp, main = 'Points = SPoRC, SS3 = Red', ylab = 'Bias Correction')
lines((out$recruit %>% filter(era %in% c("Early", "Main")))$biasadjuster, col = 'red')

# Think SS3 does some more bias ramp stuff for the initial age deviations, which
# we can't do right now.

# compare initial age structure (think differences is due to the bias ramp stuff to the early age deviations)
plot(goapcod_obj$rep$NAA[1,1,-1,1], main = 'Points = SPoRC, SS3 = Red', ylab = 'Initial Numbers at Age')
lines(unlist((out$natage %>% filter(Era == 'TIME', `Beg/Mid` == 'B', Time == 1977))[,-c(1:13)]), type = 'l', col = 'red')

# Fishing Mortality Rates
fmort_ts <- data.frame(
  ss3 = c(ss3_f$`F:_1`, ss3_f$`F:_2`, ss3_f$`F:_3`),
  SPoRC = c(goapcod_obj$rep$Fmort[,,1], goapcod_obj$rep$Fmort[,,2], goapcod_obj$rep$Fmort[,,3]),
  Type = rep(c("Fleet 1", "Fleet 2", "Fleet 3"), each = length(data$years)),
  Year = rep(data$years,3)
)

# Pretty different ... likely due to selex differences
ggplot(fmort_ts) +
  geom_line(aes(x = Year, y = ss3, lty = 'SS3', col = 'SS3'), lwd = 1.3, alpha = 0.8) +
  geom_line(aes(x = Year, y = SPoRC, lty = 'SPoRC', col = 'SPoRC'), lwd = 1.3, alpha = 0.8) +
  theme_sablefish() +
  facet_wrap(~Type, scales = 'free') +
  labs(x = 'Year', y = 'Instantaneous Fishing Mortality', color = 'Model', lty = 'Model')

# Numbers at age
naa_ts <- data.frame(
  ss3 = rowSums((out$natage %>% filter(Era == 'TIME', `Beg/Mid` == 'B'))[,-c(1:12)]),
  SPoRC = rowSums(goapcod_obj$rep$NAA[1,-49,,1]),
  Type = 'Numbers at Age',
  Year = data$years
)

ggplot(naa_ts) +
  geom_line(aes(x = Year, y = ss3, lty = 'SS3', col = 'SS3'), lwd = 1.3, alpha = 0.8) +
  geom_line(aes(x = Year, y = SPoRC, lty = 'SPoRC', col = 'SPoRC'), lwd = 1.3, alpha = 0.8) +
  theme_sablefish() +
  labs(x = 'Year', y = 'Numbers at age', color = 'Model', lty = 'Model')

# Total Biomass
totbiom_ts <- data.frame(
  ss3 = rowSums((out$natage %>% filter(Era == 'TIME', `Beg/Mid` == 'B'))[,-c(1:12)] * waa_arr),
  SPoRC = as.vector(goapcod_obj$rep$Total_Biom),
  Type = 'Total Biomass',
  Year = data$years
)

ggplot(totbiom_ts) +
  geom_line(aes(x = Year, y = ss3, lty = 'SS3', col = 'SS3'), lwd = 1.3, alpha = 0.8) +
  geom_line(aes(x = Year, y = SPoRC, lty = 'SPoRC', col = 'SPoRC'), lwd = 1.3, alpha = 0.8) +
  theme_sablefish() +
  labs(x = 'Year', y = 'Total Biomass', color = 'Model', lty = 'Model')

# Spawning Biomass
ssb_ts <- data.frame(
  ss3 = rowSums((out$natage %>% filter(Era == 'TIME', `Beg/Mid` == 'B'))[,-c(1:12)] * waa_arr * mataa_arr) * 0.5,
  SPoRC = as.vector(goapcod_obj$rep$SSB),
  Type = 'Spawning Stock Biomass',
  Year = data$years
)

ggplot(ssb_ts) +
  geom_line(aes(x = Year, y = ss3, lty = 'SS3', col = 'SS3'), lwd = 1.3, alpha = 0.8) +
  geom_line(aes(x = Year, y = SPoRC, lty = 'SPoRC', col = 'SPoRC'), lwd = 1.3, alpha = 0.8) +
  theme_sablefish() +
  labs(x = 'Year', y = 'Spawning Stock Biomass', color = 'Model', lty = 'Model')

# Compare Terminal Age Based Fishery and Survey Selectivity
sel_a_ts <- data.frame(
  ss3 = c(
    unlist((out$ageselex %>% filter(Label == '2024_1_Asel2'))[,-c(1:7)]),
    unlist((out$ageselex %>% filter(Label == '2024_2_Asel2'))[,-c(1:7)]),
    unlist((out$ageselex %>% filter(Label == '2024_3_Asel2'))[,-c(1:7)]),
    unlist((out$ageselex %>% filter(Label == '2024_4_Asel2'))[,-c(1:7)]),
    unlist((out$ageselex %>% filter(Label == '2024_5_Asel2'))[,-c(1:7)])
  ),
  SPoRC = c(goapcod_obj$rep$fish_sel[1,48,,1,1],
            goapcod_obj$rep$fish_sel[1,48,,1,2],
            goapcod_obj$rep$fish_sel[1,48,,1,3],
            goapcod_obj$rep$srv_sel[1,48,,1,1],
            goapcod_obj$rep$srv_sel[1,48,,1,2]),
  Type = rep(c("Fishery Fleet 1",
               "Fishery Fleet 2",
               "Fishery Fleet 3",
               "Survey Fleet 1",
               "Survey Fleet 2"), each = length(data$ages)),
  Ages = rep(data$ages, 5)
)


# think SS3 does finer length-bins when doing length-based selectivtiy so we could be missing something there
ggplot(sel_a_ts) +
  geom_line(aes(x = Ages, y = ss3, lty = 'SS3', col = 'SS3'), lwd = 1.3, alpha = 0.8) +
  geom_line(aes(x = Ages, y = SPoRC, lty = 'SPoRC', col = 'SPoRC'), lwd = 1.3, alpha = 0.8) +
  facet_wrap(~Type) +
  theme_sablefish() +
  labs(x = 'Age', y = 'Selectivity', color = 'Model', lty = 'Model')

# Recruitment
rec_ts <- data.frame(
  ss3 = (out$timeseries %>% filter(Era == 'TIME'))$Recruit_0,
  SPoRC = as.vector(goapcod_obj$rep$Rec),
  Type = 'Recruitment',
  Year = data$years
)

ggplot(rec_ts) +
  geom_line(aes(x = Year, y = ss3, lty = 'SS3', col = 'SS3'), lwd = 1.3, alpha = 0.8) +
  geom_line(aes(x = Year, y = SPoRC, lty = 'SPoRC', col = 'SPoRC'), lwd = 1.3, alpha = 0.8) +
  theme_sablefish() +
  labs(x = 'Year', y = 'Recruitment', color = 'Model', lty = 'Model')

# Some notes on estimation
# 1) could be different weights with CAAL vs. marginal ages and lengths being used,
# 2) Seems to not be able to predict that early recruitment spike
# 3) Could be due to time-varying ageing error + ageing bias

# # Deterministic -----------------------------------------------------------
#
# # we can also do some deterministic comparison. Using triple colon because these functions are hidden function
# obj <- RTMB::MakeADFun(SPoRC:::cmb(SPoRC:::SPoRC_rtmb, data), parameters = parameters, map = mapping, random = NULL, silent = F)
# obj$rep <- obj$report(obj$env$last.par.best) # get report of unoptimized values (deterministic comparison)
#
# # remember we setup the starting values above, based on ss3 values, so just using those for comparison
#
# # compare initial age structure (not quite sure how they are so different ... pretty sure some weird SS3 age0 stuff;
# # it estimates it back into the right scale which is good.)
# plot(obj$rep$NAA[1,1,-1,1], main = 'Points = SPoRC, SS3 = Red', ylab = 'Initial Numbers at Age')
# lines(unlist((out$natage %>% filter(Era == 'TIME', `Beg/Mid` == 'B', Time == 1977))[,-c(1:13)]), type = 'l', col = 'red')
#
# # Fishing Mortality Rates
# fmort_ts <- data.frame(
#   ss3 = c(ss3_f$`F:_1`, ss3_f$`F:_2`, ss3_f$`F:_3`),
#   SPoRC = c(obj$rep$Fmort[,,1], obj$rep$Fmort[,,2], obj$rep$Fmort[,,3]),
#   Type = rep(c("Fleet 1", "Fleet 2", "Fleet 3"), each = length(data$years)),
#   Year = rep(data$years,3)
# )
#
# # Same F if feeding in the same values
# ggplot(fmort_ts) +
#   geom_line(aes(x = Year, y = ss3, lty = 'SS3', col = 'SS3'), lwd = 1.3, alpha = 0.8) +
#   geom_line(aes(x = Year, y = SPoRC, lty = 'SPoRC', col = 'SPoRC'), lwd = 1.3, alpha = 0.8) +
#   theme_sablefish() +
#   facet_wrap(~Type, scales = 'free') +
#   labs(x = 'Year', y = 'Instantaneous Fishing Mortality', color = 'Model', lty = 'Model')
#
# # Numbers at age
# naa_ts <- data.frame(
#   ss3 = rowSums((out$natage %>% filter(Era == 'TIME', `Beg/Mid` == 'B'))[,-c(1:12)]),
#   SPoRC = rowSums(obj$rep$NAA[1,-49,,1]),
#   Type = 'Numbers at Age',
#   Year = data$years
# )
#
# # Numbers at age are pretty similar except for early period, b/c of initial age deviations
# # other differences are probably because selectivity is parameterized different
# ggplot(naa_ts) +
#   geom_line(aes(x = Year, y = ss3, lty = 'SS3', col = 'SS3'), lwd = 1.3, alpha = 0.8) +
#   geom_line(aes(x = Year, y = SPoRC, lty = 'SPoRC', col = 'SPoRC'), lwd = 1.3, alpha = 0.8) +
#   theme_sablefish() +
#   labs(x = 'Year', y = 'Numbers at age', color = 'Model', lty = 'Model')
#
# # Total Biomass
# totbiom_ts <- data.frame(
#   ss3 = rowSums((out$natage %>% filter(Era == 'TIME', `Beg/Mid` == 'B'))[,-c(1:12)] * waa_arr),
#   SPoRC = as.vector(obj$rep$Total_Biom),
#   Type = 'Total Biomass',
#   Year = data$years
# )
#
# # again, pretty close for total biomass (all these can be chalked up to selex not being the same;
# # if we feed in the same selectivity values from ss3 into SPoRC, these will be identical)
# ggplot(totbiom_ts) +
#   geom_line(aes(x = Year, y = ss3, lty = 'SS3', col = 'SS3'), lwd = 1.3, alpha = 0.8) +
#   geom_line(aes(x = Year, y = SPoRC, lty = 'SPoRC', col = 'SPoRC'), lwd = 1.3, alpha = 0.8) +
#   theme_sablefish() +
#   labs(x = 'Year', y = 'Total Biomass', color = 'Model', lty = 'Model')
#
# # Spawning Biomass
# ssb_ts <- data.frame(
#   ss3 = rowSums((out$natage %>% filter(Era == 'TIME', `Beg/Mid` == 'B'))[,-c(1:12)] * waa_arr * mataa_arr) * 0.5,
#   SPoRC = as.vector(obj$rep$SSB),
#   Type = 'Spawning Stock Biomass',
#   Year = data$years
# )
#
# ggplot(ssb_ts) +
#   geom_line(aes(x = Year, y = ss3, lty = 'SS3', col = 'SS3'), lwd = 1.3, alpha = 0.8) +
#   geom_line(aes(x = Year, y = SPoRC, lty = 'SPoRC', col = 'SPoRC'), lwd = 1.3, alpha = 0.8) +
#   theme_sablefish() +
#   labs(x = 'Year', y = 'Spawning Stock Biomass', color = 'Model', lty = 'Model')
#
# # Compare Terminal Age Based Fishery and Survey Selectivity
# sel_a_ts <- data.frame(
#   ss3 = c(
#     unlist((out$ageselex %>% filter(Label == '2024_1_Asel2'))[,-c(1:7)]),
#     unlist((out$ageselex %>% filter(Label == '2024_2_Asel2'))[,-c(1:7)]),
#     unlist((out$ageselex %>% filter(Label == '2024_3_Asel2'))[,-c(1:7)])
#   ),
#   SPoRC = c(obj$rep$fish_sel[1,48,,1,1],
#             obj$rep$fish_sel[1,48,,1,2],
#             obj$rep$fish_sel[1,48,,1,3]),
#   Type = rep(c("Fishery Fleet 1",
#                "Fishery Fleet 2",
#                "Fishery Fleet 3"
#   ), each = length(data$ages)),
#   Ages = rep(data$ages, data$n_fish_fleets)
# )
#
#
# # think SS3 does finer length-bins when doing length-based selectivtiy so we could be missing something there
# # also, our selectivity from SPoRC here is different from SS3 specifications, ours is simplified and these are
# # unoptimized values with different functional forms
# ggplot(sel_a_ts) +
#   geom_line(aes(x = Ages, y = ss3, lty = 'SS3', col = 'SS3'), lwd = 1.3, alpha = 0.8) +
#   geom_line(aes(x = Ages, y = SPoRC, lty = 'SPoRC', col = 'SPoRC'), lwd = 1.3, alpha = 0.8) +
#   facet_wrap(~Type) +
#   theme_sablefish() +
#   labs(x = 'Age', y = 'Selectivity', color = 'Model', lty = 'Model')
#
# # Recruitment
# rec_ts <- data.frame(
#   ss3 = (out$timeseries %>% filter(Era == 'TIME'))$Recruit_0,
#   SPoRC = as.vector(obj$rep$Rec),
#   Type = 'Recruitment',
#   Year = data$years
# )
#
# # Recruitment values are basicially the same if feeding back the same values estimated from SS3
# # remaining differences due to bias ramp
# ggplot(rec_ts) +
#   geom_line(aes(x = Year, y = ss3, lty = 'SS3', col = 'SS3'), lwd = 1.3, alpha = 0.8) +
#   geom_line(aes(x = Year, y = SPoRC, lty = 'SPoRC', col = 'SPoRC'), lwd = 1.3, alpha = 0.8) +
#   theme_sablefish() +
#   labs(x = 'Year', y = 'Recruitment', color = 'Model', lty = 'Model')
#
#
#
# # Dig SS3 WAA -------------------------------------------------------------
#
# # Trying to figure out some weight-at-age stuff in SS3?
# # Is this time-varying and fleet-specific??
#
# # weight at age
# waa_arr <- array(0, dim = c(n_yrs, n_ages, 5))
# fleet <- 1:5
#
# # get ss3 weights
# ss3_wts <- out$ageselex %>% filter(str_detect(Label, "bodywt"))
#
# for(f in 1:length(fleet)) {
#   tmp <- ss3_wts %>% filter(Fleet == f)
#   for(i in tmp$Yr) waa_arr[which(years == i),,f] <- unlist((tmp %>% filter(Yr == i))[,-c(1:7)])
# }
#
# reshape2::melt(waa_arr) %>%
#   ggplot(aes(x = Var1, y = value,  group = Var3, color = factor(Var3))) +
#   geom_line() +
#   facet_wrap(~Var2, scales = 'free') +
#   labs(x = 'Year', y = 'Weight', color = 'Fleet') +
#   theme_sablefish()
#
# # Scratch -----------------------------------------------------------------
#
# # starting values if double normal is specified
# # double normal selectivity for fleet 1
# # parameters$ln_fish_fixed_sel_pars[1,1:6,1,1,1] <- c(
# #   log((out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_peak_FshTrawl(1)'] - min(lens)) / (max(lens) - out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_peak_FshTrawl(1)'])),
# #   out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_top_logit_FshTrawl(1)'],
# #   out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_ascend_se_FshTrawl(1)'],
# #   out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_descend_se_FshTrawl(1)'],
# #   -999,
# #   10
# # )
#
# # double normal selectivity for fleet 2
# # parameters$ln_fish_fixed_sel_pars[1,1:6,1,1,2] <- c(
# #   log((out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_peak_FshLL(2)'] - min(lens)) / (max(lens) - out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_peak_FshLL(2)'])),
# #   out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_top_logit_FshLL(2)'],
# #   out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_ascend_se_FshLL(2)'],
# #   10,
# #   -999,
# #   10
# # )
#
# # double normal selectivity for fleet 3
# # parameters$ln_fish_fixed_sel_pars[1,1:6,1,1,3] <- c(
# #   log((out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_peak_FshPot(3)'] - min(lens)) / (max(lens) - out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_peak_FshPot(3)'])),
# #   -999,
# #   out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_ascend_se_FshPot(3)'],
# #   out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_descend_se_FshPot(3)'],
# #   -999,
# #   out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_end_logit_FshPot(3)']
# # )
#
# # double normal selectivity for LLS fleet 4
# # parameters$ln_srv_fixed_sel_pars[1,1:6,1,1,1] <- c(
# #   log((out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_peak_Srv(4)'] - min(lens)) / (max(lens) - out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_peak_Srv(4)'])),
# #   -999,
# #   out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_ascend_se_Srv(4)'],
# #   out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_descend_se_Srv(4)'],
# #   -999,
# #   out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_end_logit_Srv(4)']
# # )
#
# # double normal selectivity for LLS fleet 5
# # parameters$ln_srv_fixed_sel_pars[1,1:6,1,1,2] <- c(
# #   log((out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_peak_LLSrv(5)'] - min(lens)) / (max(lens) - out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_peak_LLSrv(5)'])),
# #   -999,
# #   6,
# #   1,
# #   -999,
# #  0.9
# # )
