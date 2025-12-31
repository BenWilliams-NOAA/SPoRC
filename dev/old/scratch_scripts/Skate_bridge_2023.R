# Purpose: To bridge the 2023 BSAI Skate Model from SS3 to SPoRC
# Creator: Matthew LH. Cheng
# Date: 7/3/25

get_al_trans_matrix = function(age_bins, len_bins, mean_length, sd) {
  # container array
  age_length = matrix(0.0, nrow = length(age_bins), ncol = length(len_bins))

  for(a in 1:length(age_bins)) {
    # Use actual bin lower limits as per specification
    # Assume len_bins contains the lower limits of each bin
    bin_lower_limits = len_bins

    # Calculate cumulative probabilities at bin lower limits
    AL = pnorm(bin_lower_limits, mean_length[a], sd[a])

    # Calculate bin probabilities according to the specification
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

# Setup -------------------------------------------------------------------

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

library(here)
library(SPoRC)
library(r4ss)

# Read in parameters, report, and data from SS3 files
rep <- r4ss::SS_read(here("dev", "dev_data", "bsaiskate"))
waa <- r4ss::SS_readwtatage(here("dev", "dev_data", "bsaiskate", "wtatage.ss_new"))
out <- r4ss::SS_output(here("dev", "dev_data", "bsaiskate"), verbose = TRUE)

# Extract Data ------------------------------------------------------------

# define dimensions
n_regions <- 1 # number of regions
n_sexes <- 1 # number of sexes
n_fish_fleets <- 2 # number of fishery fleets
n_srv_fleets <- 1 # number of survey fleets
years <- rep$dat$styr:rep$dat$endyr # years
n_yrs <- length(years) # number of years
ages <- 0:25 # vector of ages
n_ages <- length(ages) # number of ages
lens <- rep$dat$lbin_vector # vector of lengths
n_lens <- length(lens) # number of length bins
spwn_month <- rep$dat$spawn_month - 1 # spawning month (start of year)

### Demographics ------------------------------------------------------------

# weight at age
waa <- out$ageselex %>% filter(Label == '1952_1_bodywt')
waa_arr <- array(0, dim = c(n_regions, n_yrs, n_ages, n_sexes))
waa_arr[] <- rep(out$endgrowth$Wt_Mid, each = n_yrs)

# get length at age and sd to construct ALK
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

# maturity array
mataa_arr <- array(rep(out$endgrowth$Len_Mat, each = n_yrs), dim = c(n_regions, n_yrs, n_ages, n_sexes))

# Natural mortality array
fixed_natmort <- array(0.13, dim = c(n_regions, n_yrs, n_ages, n_sexes))


### Fishery Stuff -----------------------------------------------------------

# fishery catches
ss3_catch <- rep$dat$catch %>% dplyr::filter(year != -999) # get ss3 catches

ObsFishCatch <- array(0, dim = c(n_regions, n_yrs, n_fish_fleets)) # observed array
ObsFishCatch[1,,1] <- as.numeric(gsub(",", "", ss3_catch$catch[ss3_catch$fleet == 1])) # converting characters to numbers (Fleet 1)
ObsFishCatch[1,,2] <- as.numeric(gsub(",", "", ss3_catch$catch[ss3_catch$fleet == 2])) # converting characters to numbers (Fleet 2)
Catch_Type <- array(1, dim = c(n_yrs, n_fish_fleets)) # Catch type (set at 1 since its a single region model)
UseCatch <- array(1, dim = c(n_regions, n_yrs, n_fish_fleets)) # Use catch indicator
UseCatch[1,which(ObsFishCatch[1,,1] == 0),1] <- 0 # don't fit if 0 catch
UseCatch[1,which(ObsFishCatch[1,,2] == 0),2] <- 0 # don't fit if 0 catch

# ObsFishCatch[1,74,1] <- 15 # NOTE: This needs to be corrected!

# fishery index (none used)
ObsFishIdx <- array(NA, dim = c(n_regions, n_yrs, n_fish_fleets))
ObsFishIdx_SE <- array(NA, dim = c(n_regions, n_yrs, n_fish_fleets))
UseFishIdx <- array(0, dim = c(n_regions, n_yrs, n_fish_fleets))

# fishery ages (none used)
ObsFishAgeComps <- array(NA, dim = c(n_regions, n_yrs, n_ages, n_sexes, n_fish_fleets))
UseFishAgeComps <- array(0, dim = c(n_regions, n_yrs, n_fish_fleets))
ISS_FishAgeComps <- array(0, dim = c(n_regions, n_yrs, n_sexes, n_fish_fleets))

# fishery lengths
ObsFishLenComps <- array(0, dim = c(n_regions, n_yrs, n_lens, n_sexes, n_fish_fleets))
UseFishLenComps <- array(0, dim = c(n_regions, n_yrs, n_fish_fleets))
ISS_FishLenComps <- array(0, dim = c(n_regions, n_yrs, n_sexes, n_fish_fleets))

# format ss3 data (fleet 1)
ss3_fish1_lens <- rep$dat$lencomp %>% dplyr::filter(fleet == 1) # get ss3 lenghts
ss3_fish1_yrs <- ss3_fish1_lens$year

# loop through to populate
for(i in ss3_fish1_yrs) {
  ObsFishLenComps[1,which(years == i),,1,1] <- unlist((ss3_fish1_lens %>% dplyr::filter(year == i))[,-c(1:6)]) # input fish1 lengths
  UseFishLenComps[1,which(years == i),1] <- 1 # use indicator
  ISS_FishLenComps[1,which(years == i),1,1] <- unlist((ss3_fish1_lens %>% dplyr::filter(year == i))[,6]) # input sample size for fish1 lengths
} # end i loop

# format ss3 data (fleet 2)
ss3_fish2_lens <- rep$dat$lencomp %>% dplyr::filter(fleet == 2) # get ss3 lenghts
ss3_fish2_yrs <- ss3_fish2_lens$year

# loop through to populate
for(i in ss3_fish2_yrs) {
  ObsFishLenComps[1,which(years == i),,1,2] <- unlist((ss3_fish2_lens %>% dplyr::filter(year == i))[,-c(1:6)]) # input fish1 lengths
  UseFishLenComps[1,which(years == i),2] <- 1 # use indicator
  ISS_FishLenComps[1,which(years == i),1,2] <- unlist((ss3_fish2_lens %>% dplyr::filter(year == i))[,6]) # input sample size for fish1 lengths
} # end i loop


## Survey Stuff ------------------------------------------------------------
# survey index
ss3_idx <- rep$dat$CPUE # get ss3 cpue
ss3_idx$year <- rep$dat$CPUE$year
ObsSrvIdx <- array(0, dim = c(n_regions, n_yrs, n_srv_fleets))
ObsSrvIdx_SE <- array(0, dim = c(n_regions, n_yrs, n_srv_fleets))
UseSrvIdx <- array(0, dim = c(n_regions, n_yrs, n_srv_fleets))
ObsSrvIdx[1,which(years %in% c(ss3_idx$year)),1] <- rep$dat$CPUE$obs
ObsSrvIdx_SE[1,which(years %in% c(ss3_idx$year)),1] <- rep$dat$CPUE$se_log
UseSrvIdx[1,which(years %in% c(ss3_idx$year)),1] <- 1

# survey ages (none used)
ObsSrvAgeComps <- array(NA, dim = c(n_regions, n_yrs, n_ages, n_sexes, n_srv_fleets))
UseSrvAgeComps <- array(0, dim = c(n_regions, n_yrs, n_srv_fleets))
ISS_SrvAgeComps <- array(0, dim = c(n_regions, n_yrs, n_sexes, n_srv_fleets))

# survey lengths
ObsSrvLenComps <- array(0, dim = c(n_regions, n_yrs, n_lens, n_sexes, n_srv_fleets))
UseSrvLenComps <- array(0, dim = c(n_regions, n_yrs, n_srv_fleets))
ISS_SrvLenComps <- array(0, dim = c(n_regions, n_yrs, n_sexes, n_srv_fleets))

# format ss3 data
ss3_srv_lens <- rep$dat$lencomp %>% dplyr::filter(fleet == 3) # get ss3 lenghts
ss3_srv_yrs <- ss3_srv_lens$year

# loop through to populate
for(i in ss3_srv_yrs) {
  ObsSrvLenComps[1,which(years == i),,1,1] <- unlist((ss3_srv_lens %>% dplyr::filter(year == i))[,-c(1:6)]) # input srv lengths
  UseSrvLenComps[1,which(years == i),1] <- 1 # use indicator
  ISS_SrvLenComps[1,which(years == i),1,1] <- unlist((ss3_srv_lens %>% dplyr::filter(year == i))[,6]) # input sample size for srv lengths
} # end i loop


# Setup Model -------------------------------------------------------------

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
  do_rec_bias_ramp = 1,
  # recruitment bias ramp (according to SS3 specifications, full bias correction until 2018)
  bias_year = c(0, 0, length(1950:2018),length(1950:2024)),
  # do bias ramp (0 == don't do bias ramp, 1 == do bias ramp)
  sigmaR_switch = 1,
  # when to switch from early to late sigmaR (switch in first year)
  ln_sigmaR = c(NA, log(0.4)), # 2 values for early and late sigma (early sigma is not used); starting values for early and late sigmaR

  # recruitment model
  rec_model = "bh_rec",
  h_spec = 'fix',
  steepness_h = 100, # mean reruitment (because steepness is set at 1)
  sigmaR_spec = "fix",
  # fix early sigmaR and late sigmaR
  sexratio = as.vector(c(1.0)),
  # recruitment sex ratio
  init_age_strc = 1, # geometric series to derive initial age structure
  t_spawn = spwn_month,
  equil_init_age_strc = 0
)

input_list <- Setup_Mod_Biologicals(
  input_list = input_list,

  # Data inputs
  WAA = waa_arr,
  MatAA = mataa_arr,

  # Model options
  fit_lengths = 1, # fit lengths
  Selex_Type = "length", # length-based selectivity
  SizeAgeTrans = sizeage,
  AgeingError = NULL,
  M_spec = "fix",
  # fixing natural mortality
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
  sigmaC_spec = "fix",
  Catch_Constant = rep(0.00001, 2),
  ln_sigmaC = array(log(0.01), dim = c(n_regions, n_fish_fleets)), # sigma catch
  ln_sigmaF = array(log(10), dim = c(n_regions, n_fish_fleets)) # sigma of f penalty(set large because catch is known)
)

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
  fish_idx_type = c("none", "none"),
  # indices for fishery
  FishAgeComps_LikeType = c("none", "none"),
  # age comp likelihoods for fishery fleet
  FishLenComps_LikeType = c("Multinomial", "Multinomial"),
  # length comp likelihoods for fishery
  FishAgeComps_Type = c("none_Year_1-terminal_Fleet_1", "none_Year_1-terminal_Fleet_2"),
  # age comp structure for fishery
  FishLenComps_Type = c("agg_Year_1-terminal_Fleet_1", "agg_Year_1-terminal_Fleet_2"),
  # length comp structure for fishery
  FishAge_comp_agg_type = c(NA, NA),
  # ADMB aggregation quirks, ideally get rid of this
  FishLen_comp_agg_type = c(0, 0)
  # ADMB aggregation quirks, ideally get rid of this
)

# Setup survey indices and compositions
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
  srv_idx_type = c("biom"),
  # abundance and biomass for survey fleet 1
  SrvAgeComps_LikeType = c("none"),
  # survey age composition likelihood for survey fleet 1
  SrvLenComps_LikeType = c("Multinomial"),
  #  survey length composition likelihood for survey fleet 1
  SrvAgeComps_Type = c("none_Year_1-terminal_Fleet_1"),
  # survey age comp type
  SrvLenComps_Type = c( "agg_Year_1-terminal_Fleet_1"),
  # survey length comp type
  SrvAge_comp_agg_type = NA,
  # ADMB aggregation quirks, ideally get rid of this
  SrvLen_comp_agg_type = 0
  # ADMB aggregation quirks, ideally get rid of this
)


# setup fishery selectivity and catchability
input_list <- Setup_Mod_Fishsel_and_Q(

  input_list = input_list,

  # Model options
  # fishery selectivity, whether continuous time-varying
  cont_tv_fish_sel = c("none_Fleet_1", "none_Fleet_2"),
  # fishery selectivity blocks
  fish_sel_blocks = c("none_Fleet_1", "none_Fleet_2"),
  # fishery selectivity form
  fish_sel_model = c("dbnrml_Fleet_1", "dbnrml_Fleet_2"),
  # fishery catchability blocks
  fish_q_blocks = c("none_Fleet_1", "none_Fleet_2"),
  # whether to estiamte all fixed effects for fishery selectivity
  fish_fixed_sel_pars_spec = c("est_all", "est_all"),
  # whether to estiamte all fixed effects for fishery catchability
  fish_q_spec = c("fix", "fix")
)

# Setup survey selectivity and catchability
input_list <- Setup_Mod_Srvsel_and_Q(
  input_list = input_list,

  # Model options
  # survey selectivity, whether continuous time-varying
  cont_tv_srv_sel = c("none_Fleet_1"),
  # survey selectivity blocks
  srv_sel_blocks = c("none_Fleet_1"),
  # survey selectivity form
  srv_sel_model = c("dbnrml_Fleet_1"),
  # survey catchability blocks
  srv_q_blocks = c("none_Fleet_1"),
  # whether to estiamte all fixed effects for survey selectivity
  srv_fixed_sel_pars_spec = c("est_all"),
  # whether to estiamte all fixed effects for survey catchability
  srv_q_spec = c("fix")
)

input_list$par$ln_srv_q[] <- log(1) # set starting value / fix parameter at catchabilty of 1

# Setup model weighting (empahsis factors for SS3)
input_list <- Setup_Mod_Weighting(
  input_list = input_list,
  sablefish_ADMB = 0,
  likelihoods = 1, # using TMB likelihoods
  Wt_Catch = 1,
  Wt_FishIdx = 0,
  Wt_SrvIdx = 1,
  Wt_Rec = 1,
  Wt_F = 1,
  Wt_Tagging = 0,
  Wt_FishAgeComps = array(0, dim = c(input_list$data$n_regions,
                                     length(input_list$data$years),
                                     input_list$data$n_sexes,
                                     input_list$data$n_fish_fleets)),
  Wt_FishLenComps = array(1, dim = c(input_list$data$n_regions,
                                     length(input_list$data$years),
                                     input_list$data$n_sexes,
                                     input_list$data$n_fish_fleets)),
  Wt_SrvAgeComps = array(0, dim = c(input_list$data$n_regions,
                                    length(input_list$data$years),
                                    input_list$data$n_sexes,
                                    input_list$data$n_srv_fleets)),
  Wt_SrvLenComps = array(1, dim = c(input_list$data$n_regions,
                                    length(input_list$data$years),
                                    input_list$data$n_sexes,
                                    input_list$data$n_srv_fleets))
)

# Deterministic Comparisons -----------------------------------------------

# extract data, parameters, and mapping from input_list (via the setup functions defined above)
data <- input_list$data
parameters <- input_list$par
mapping <- input_list$map

## Setup starting values ---------------------------------------------------

# We are setting up starting values here to make sure when the model is unoptimized, we can get the same predictions
# of SSB, biomass, etc back, just to ensure that the general equations and dynamics are identicial (or close to identical)

# Setting additional starting vaues
parameters$ln_global_R0 <- out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'SR_LN(R0)'] # mean recruitment
parameters$ln_RecDevs[] <- out$recruitpars$Value[1:74] # recruitment deviations

# setup starting values for mean fishing mortality
ss3_f <- out$timeseries %>% filter(Era == 'TIME')
parameters$ln_F_mean[] <- c(mean(log(ss3_f$`F:_1`[ss3_f$`F:_1` != 0])), mean(log(ss3_f$`F:_2`[ss3_f$`F:_2` != 0]))) # log mean F by fleet
parameters$ln_F_devs[,,1] <- log(ss3_f$`F:_1`) - parameters$ln_F_mean[,1] # f devs for fleet 1
parameters$ln_F_devs[,,2] <- log(ss3_f$`F:_2`) - parameters$ln_F_mean[,2] # f devs for fleet 2
parameters$ln_F_devs[!is.finite(parameters$ln_F_devs)] <- 0 # don't estimate f devs when catch = 0

# double normal selectivity for fleet 1
parameters$ln_fish_fixed_sel_pars[1,1:6,1,1,1] <- c(
  log((out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_peak_LGL(1)'] - min(lens)) / (max(lens) - out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_peak_LGL(1)'])),
  out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_top_logit_LGL(1)'],
  out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_ascend_se_LGL(1)'],
  out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_descend_se_LGL(1)'],
  out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_start_logit_LGL(1)'],
  out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_end_logit_LGL(1)']
)

# double normal selectivity for fleet 2
parameters$ln_fish_fixed_sel_pars[1,1:6,1,1,2] <- c(
  log((out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_peak_TWL(2)'] - min(lens)) / (max(lens) - out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_peak_TWL(2)'])),
  out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_top_logit_TWL(2)'],
  out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_ascend_se_TWL(2)'],
  out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_descend_se_TWL(2)'],
  out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_start_logit_TWL(2)'],
  out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_end_logit_TWL(2)']
)

# double normal selectivity parameters for survey fleet
parameters$ln_srv_fixed_sel_pars[] <- c(
  log((out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_peak_SURV(3)'] - min(lens)) / (max(lens) - out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_peak_SURV(3)'])),
  out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_top_logit_SURV(3)'],
  out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_ascend_se_SURV(3)'],
  out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_descend_se_SURV(3)'],
  out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_start_logit_SURV(3)'],
  out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'Size_DblN_end_logit_SURV(3)']
)


# Compare Deterministic Computations (Pre Optimization) -------------------

# make AD model function (using triple colon because the RTMB model is hidden in the package)
obj <- RTMB::MakeADFun(SPoRC:::cmb(SPoRC:::SPoRC_rtmb, data), parameters = parameters, map = mapping, random = NULL, silent = F)
obj$rep <- obj$report(obj$env$last.par.best) # get report of unoptimized values (deterministic comparison)

# Fishing Mortality Rates
fmort_ts <- data.frame(
  ss3 = c(ss3_f$`F:_1`, ss3_f$`F:_2`),
  SPoRC = c(obj$rep$Fmort[,,1], obj$rep$Fmort[,,2]),
  Type = rep(c("Fleet 1", "Fleet 2"), each = length(data$years)),
  Year = rep(data$years,2)
)

ggplot(fmort_ts) +
  geom_line(aes(x = Year, y = ss3, lty = 'SS3', col = 'SS3'), lwd = 1.3, alpha = 0.8) +
  geom_line(aes(x = Year, y = SPoRC, lty = 'SPoRC', col = 'SPoRC'), lwd = 1.3, alpha = 0.8) +
  theme_sablefish() +
  facet_wrap(~Type) +
  labs(x = 'Year', y = 'Instantaneous Fishing Mortality', color = 'Model', lty = 'Model')

# Numbers at age
naa_ts <- data.frame(
  ss3 = rowSums((out$natage %>% filter(Era == 'TIME', `Beg/Mid` == 'B'))[,-c(1:12)]),
  SPoRC = rowSums(obj$rep$NAA[1,-75,,1]),
  Type = 'Numbers at Age',
  Year = data$years
)

# Slight differences likely due to how length-based selectivity is computed
ggplot(naa_ts) +
  geom_line(aes(x = Year, y = ss3, lty = 'SS3', col = 'SS3'), lwd = 1.3, alpha = 0.8) +
  geom_line(aes(x = Year, y = SPoRC, lty = 'SPoRC', col = 'SPoRC'), lwd = 1.3, alpha = 0.8) +
  theme_sablefish() +
  labs(x = 'Year', y = 'Numbers at age', color = 'Model', lty = 'Model')

# Total Biomass
totbiom_ts <- data.frame(
  ss3 = rowSums((out$natage %>% filter(Era == 'TIME', `Beg/Mid` == 'B'))[,-c(1:12)] * waa_arr),
  SPoRC = as.vector(obj$rep$Total_Biom),
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
  SPoRC = as.vector(obj$rep$SSB),
  Type = 'Spawning Stock Biomass',
  Year = data$years
)

ggplot(ssb_ts) +
  geom_line(aes(x = Year, y = ss3, lty = 'SS3', col = 'SS3'), lwd = 1.3, alpha = 0.8) +
  geom_line(aes(x = Year, y = SPoRC, lty = 'SPoRC', col = 'SPoRC'), lwd = 1.3, alpha = 0.8) +
  theme_sablefish() +
  labs(x = 'Year', y = 'Spawning Stock Biomass', color = 'Model', lty = 'Model')

# Compare Age Based Fishery and Survey Selectivity
sel_a_ts <- data.frame(
  ss3 = c(
    unlist((out$ageselex %>% filter(Label == '1980_1_Asel2'))[,-c(1:7)]),
    unlist((out$ageselex %>% filter(Label == '1980_2_Asel2'))[,-c(1:7)]),
    unlist((out$ageselex %>% filter(Label == '1980_3_Asel2'))[,-c(1:7)])
  ),
  SPoRC = c(obj$rep$fish_sel[1,1,,1,1], obj$rep$fish_sel[1,1,,1,2], obj$rep$srv_sel[1,1,,1,1]),
  Type = rep(c("Fishery Fleet 1", "Fishery Fleet 2", "Survey Fleet 1"), each = length(data$ages)),
  Ages = rep(data$ages, 3)
)

# Basically, the differences are likely because SS3 (i think ...) uses finer scale length bins
# to compute the ALK and selectivity compared to the modelled number of lengths, and also seemingly
# appears to set the age-0 selectivity at 0.

ggplot(sel_a_ts) +
  geom_line(aes(x = Ages, y = ss3, lty = 'SS3', col = 'SS3'), lwd = 1.3, alpha = 0.8) +
  geom_line(aes(x = Ages, y = SPoRC, lty = 'SPoRC', col = 'SPoRC'), lwd = 1.3, alpha = 0.8) +
  facet_wrap(~Type) +
  theme_sablefish() +
  labs(x = 'Age', y = 'Selectivity', color = 'Model', lty = 'Model')


# Optimized Comparisons v1 ---------------------------------------------------

# In this optimized comparison, we will attempt to do a direct bridge comparison
# by optimizing values

# Fit model here
skate_obj <- SPoRC::fit_model(data = data,
                              parameters = parameters,
                              mapping = mapping,
                              random = NULL,
                              newton_loops = 3,
                              silent = FALSE)

# save sdreport
skate_obj$sdrep <- RTMB::sdreport(skate_obj) # get standard error report

post_optim_sanity_checks(skate_obj$sdrep,
                         skate_obj$rep,
                         gradient_tol = 1e-3,
                         se_tol = 10,
                         corr_tol = 1)

# NOTE: Model does not converge with double normal selectivity ... (it is
# hitting paramter bounds )

# we can inspect what this is estimated at:
get_selex_plot(list(skate_obj$rep), model_names = 1, Selex_Type = 'length')

# we can also look at SS3 values to see what's going on
out$estimated_non_dev_parameters
# inspecting some of the double normal parameters, it looks like they are estimated
# all very close to the initial starting values, potentially suggesting that
# the model is very unstable and barely moved from those values. Likely suggests there
# is not enough information to estimate these parameters and that we need to simplifiy


# Optimized Comparison v2 -------------------------------------------------

# In light of the results from v1, we will simplify the model to have
# logistic selectivity for the fishery fleets and the survey fleet.

# We can iteratively modify the input list already defined above to do so.

# setup fishery selectivity and catchability
input_list_v2 <- Setup_Mod_Fishsel_and_Q(

  input_list = input_list,

  # Model options
  # fishery selectivity, whether continuous time-varying
  cont_tv_fish_sel = c("none_Fleet_1", "none_Fleet_2"),
  # fishery selectivity blocks
  fish_sel_blocks = c("none_Fleet_1", "none_Fleet_2"),
  # fishery selectivity form
  fish_sel_model = c("logist1_Fleet_1", "logist1_Fleet_2"), # changing this to logistic!
  # fishery catchability blocks
  fish_q_blocks = c("none_Fleet_1", "none_Fleet_2"),
  # whether to estiamte all fixed effects for fishery selectivity
  fish_fixed_sel_pars_spec = c("est_all", "est_all"),
  # whether to estiamte all fixed effects for fishery catchability
  fish_q_spec = c("fix", "fix")
)

# Setup survey selectivity and catchability
input_list_v2 <- Setup_Mod_Srvsel_and_Q(
  input_list = input_list_v2,

  # Model options
  # survey selectivity, whether continuous time-varying
  cont_tv_srv_sel = c("none_Fleet_1"),
  # survey selectivity blocks
  srv_sel_blocks = c("none_Fleet_1"),
  # survey selectivity form
  srv_sel_model = c("logist1_Fleet_1"), # changing this to logistic!
  # survey catchability blocks
  srv_q_blocks = c("none_Fleet_1"),
  # whether to estiamte all fixed effects for survey selectivity
  srv_fixed_sel_pars_spec = c("est_all"),
  # whether to estiamte all fixed effects for survey catchability
  srv_q_spec = c("fix")
)

input_list_v2$par$ln_srv_q[] <- log(1) # set starting value / fix parameter at catchabilty of 1

# extract data, parameters, and mapping from input_list (via the setup functions defined above)
data <- input_list_v2$data
parameters <- input_list_v2$par
mapping <- input_list_v2$map

# Setup some starting values to help the model
parameters$ln_fish_fixed_sel_pars[,1,,,] <- log(60) # l50 for fishery fleet 1 and fleet 2, respectively
parameters$ln_fish_fixed_sel_pars[,2,,,] <- log(1) # slope for fishery fleet 1 and fleet 2, respectively
parameters$ln_srv_fixed_sel_pars[,1,,,] <- log(60) # l50 for survey fleet 1
parameters$ln_srv_fixed_sel_pars[,2,,,] <- log(1) # l95 for survey fleet 1

# Setting additional starting values
parameters$ln_global_R0 <- out$estimated_non_dev_parameters$Value[rownames(out$estimated_non_dev_parameters) == 'SR_LN(R0)'] # mean recruitment
parameters$ln_RecDevs[] <- out$recruitpars$Value[1:74] # recruitment deviations

# setup starting values for mean fishing mortality
ss3_f <- out$timeseries %>% filter(Era == 'TIME')
parameters$ln_F_mean[] <- c(mean(log(ss3_f$`F:_1`[ss3_f$`F:_1` != 0])), mean(log(ss3_f$`F:_2`[ss3_f$`F:_2` != 0]))) # log mean F by fleet
parameters$ln_F_devs[,,1] <- log(ss3_f$`F:_1`) - parameters$ln_F_mean[,1] # f devs for fleet 1
parameters$ln_F_devs[,,2] <- log(ss3_f$`F:_2`) - parameters$ln_F_mean[,2] # f devs for fleet 2
parameters$ln_F_devs[!is.finite(parameters$ln_F_devs)] <- 0 # don't estimate f devs when catch = 0

# Fit model here
skate_obj_v2 <- SPoRC::fit_model(data = data,
                              parameters = parameters,
                              mapping = mapping,
                              random = NULL,
                              newton_loops = 3,
                              silent = FALSE)

# save sdreport
skate_obj_v2$sdrep <- RTMB::sdreport(skate_obj_v2) # get standard error report

post_optim_sanity_checks(skate_obj_v2$sdrep,
                         skate_obj_v2$rep,
                         gradient_tol = 1e-4,
                         se_tol = 5,
                         corr_tol = 1)


### Compare Model Results -----------------------------------------------------------

# Fishing Mortality Rates
fmort_ts <- data.frame(
  ss3 = c(ss3_f$`F:_1`, ss3_f$`F:_2`),
  SPoRC = c(skate_obj_v2$rep$Fmort[,,1], skate_obj_v2$rep$Fmort[,,2]),
  Type = rep(c("Fleet 1", "Fleet 2"), each = length(data$years)),
  Year = rep(data$years,2)
)

ggplot(fmort_ts) +
  geom_line(aes(x = Year, y = ss3, lty = 'SS3', col = 'SS3'), lwd = 1.3, alpha = 0.8) +
  geom_line(aes(x = Year, y = SPoRC, lty = 'SPoRC', col = 'SPoRC'), lwd = 1.3, alpha = 0.8) +
  theme_sablefish() +
  facet_wrap(~Type) +
  labs(x = 'Year', y = 'Instantaneous Fishing Mortality', color = 'Model', lty = 'Model')

# Numbers at age
naa_ts <- data.frame(
  ss3 = rowSums((out$natage %>% filter(Era == 'TIME', `Beg/Mid` == 'B'))[,-c(1:12)]),
  SPoRC = rowSums(skate_obj_v2$rep$NAA[1,-75,,1]),
  Type = 'Numbers at Age',
  Year = data$years
)

# Slight differences likely due to how length-based selectivity is computed
ggplot(naa_ts) +
  geom_line(aes(x = Year, y = ss3, lty = 'SS3', col = 'SS3'), lwd = 1.3, alpha = 0.8) +
  geom_line(aes(x = Year, y = SPoRC, lty = 'SPoRC', col = 'SPoRC'), lwd = 1.3, alpha = 0.8) +
  theme_sablefish() +
  labs(x = 'Year', y = 'Numbers at age', color = 'Model', lty = 'Model')

# Total Biomass
totbiom_ts <- data.frame(
  ss3 = rowSums((out$natage %>% filter(Era == 'TIME', `Beg/Mid` == 'B'))[,-c(1:12)] * waa_arr),
  SPoRC = as.vector(skate_obj_v2$rep$Total_Biom),
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
  SPoRC = as.vector(skate_obj_v2$rep$SSB),
  Type = 'Spawning Stock Biomass',
  Year = data$years
)

ggplot(ssb_ts) +
  geom_line(aes(x = Year, y = ss3, lty = 'SS3', col = 'SS3'), lwd = 1.3, alpha = 0.8) +
  geom_line(aes(x = Year, y = SPoRC, lty = 'SPoRC', col = 'SPoRC'), lwd = 1.3, alpha = 0.8) +
  theme_sablefish() +
  labs(x = 'Year', y = 'Spawning Stock Biomass', color = 'Model', lty = 'Model')

# Compare Age Based Fishery and Survey Selectivity
sel_a_ts <- data.frame(
  ss3 = c(
    unlist((out$ageselex %>% filter(Label == '1980_1_Asel2'))[,-c(1:7)]),
    unlist((out$ageselex %>% filter(Label == '1980_2_Asel2'))[,-c(1:7)]),
    unlist((out$ageselex %>% filter(Label == '1980_3_Asel2'))[,-c(1:7)])
  ),
  SPoRC = c(skate_obj_v2$rep$fish_sel[1,1,,1,1], skate_obj_v2$rep$fish_sel[1,1,,1,2], skate_obj_v2$rep$srv_sel[1,1,,1,1]),
  Type = rep(c("Fishery Fleet 1", "Fishery Fleet 2", "Survey Fleet 1"), each = length(data$ages)),
  Ages = rep(data$ages, 3)
)

ggplot(sel_a_ts) +
  geom_line(aes(x = Ages, y = ss3, lty = 'SS3', col = 'SS3'), lwd = 1.3, alpha = 0.8) +
  geom_line(aes(x = Ages, y = SPoRC, lty = 'SPoRC', col = 'SPoRC'), lwd = 1.3, alpha = 0.8) +
  facet_wrap(~Type) +
  theme_sablefish() +
  labs(x = 'Age', y = 'Selectivity', color = 'Model', lty = 'Model')



# Compare Reference Points ------------------------------------------------

# get quick look at ABC and reference points
# setup projection model options
proj_model_opt <- list(
  n_proj_yrs = 2,
  n_avg_yrs = 1,
  HCR_function = function(x, frp, brp, alpha = 0.05) {
    stock_status <- x / brp
    if (stock_status >= 1) f <- frp
    if (stock_status > alpha && stock_status < 1) f <- frp * (stock_status - alpha) / (1 - alpha)
    if (stock_status < alpha) f <- 0
    return(f)
  },
  recruitment_opt = "mean_rec",
  fmort_opt = "HCR"
)

# setup reference point options
reference_points_opt <- list(
  SPR_x = 0.4,
  t_spwn = 0.5,
  sex_ratio_f = 0.5,
  calc_rec_st_yr = which(data$years == 1979),
  rec_age = 0,
  type = "single_region",
  what = 'SPR'
)

ref_table <- get_key_quants(data = list(data),
                            rep = list(skate_obj_v2$rep),
                            reference_points_opt =  reference_points_opt,
                            proj_model_opt = proj_model_opt,
                            model_names = "SPoRC")

ref_table[[2]]

# Compare to 2023 assessment (not too far off)
# B40 = 69152 vs. 67536.77
# F40 = 0.08 vs. 0.08075
# ABC = 27950 vs. 28743.54
# Terminal SSB = 106549 vs. 108450.3

# Retrospective -----------------------------------------------------------

ret <- do_retrospective(10,
                        data,
                        parameters,
                        mapping,
                        random = NULL,
                        do_par = TRUE,
                        n_cores = 8,
                        do_francis = FALSE)

ret1 <- get_retrospective_relative_difference(ret)
get_retrospective_plot(retro_output = ret, Rec_Age = 1)
