# Purpose: To bridge to the spatial assessment in Cheng and Marsh et al. 2025  using SPoRC
# Creator: Matthew LH. Cheng
# Date Created: 2/5/25

# Set up ------------------------------------------------------------------

library(here)
library(tidyverse)
library(RTMB)
library(SPoRC)

# Initialize model dimensions and data list
input_list <- Setup_Mod_Dim(years = 1:length(mlt_rg_sable_data$years), # vector of years (1 - 62)
                            ages = mlt_rg_sable_data$ages, # vector of ages (1 - 30)
                            lens = mlt_rg_sable_data$lens, # number of lengths (41 - 99)
                            n_regions = mlt_rg_sable_data$n_regions, # number of regions (5)
                            n_sexes = mlt_rg_sable_data$n_sexes, # number of sexes (2)
                            n_fish_fleets = mlt_rg_sable_data$n_fish_fleets, # number of fishery fleet (2)
                            n_srv_fleets = mlt_rg_sable_data$n_srv_fleets, # number of survey fleets (2)
                            verbose = TRUE
)

# Setup recruitment stuff (using defaults for other stuff)
input_list <- Setup_Mod_Rec(input_list = input_list, # input data list from above
                            do_rec_bias_ramp = 0, # not using bias ramp
                            sigmaR_switch = 16, # switch to using late sigma in year 16
                            dont_est_recdev_last = 1, # don't estimate last rec dev
                            sexratio = c(0.5, 0.5), # fix sex ratio at 0.5

                            # Model options
                            rec_model = "mean_rec", # recruitment model
                            sigmaR_spec = "fix", # fixing
                            InitDevs_spec = "est_shared_r", # initial deviations are shared across regions,
                            # but recruitment deviations are region specific
                            ln_sigmaR = log(c(0.4, 1.2)), # values to fix sigmaR at, or starting values
                            ln_global_R0 = log(20),
                            # starting value for global R0
                            R0_prop = array(c(0.2, 0.2, 0.2, 0.2),
                                            dim = c(input_list$data$n_regions - 1))
                            # starting value for R0 proportions in multinomial logit space
)

# Setup biological stuff (using defaults for other stuff)
input_list <- Setup_Mod_Biologicals(input_list = input_list,
                                    WAA = mlt_rg_sable_data$WAA, # weight at age
                                    MatAA = mlt_rg_sable_data$MatAA, # maturity at age
                                    AgeingError = mlt_rg_sable_data$AgeingError, # ageing error matrix
                                    fit_lengths = 1, # fitting lengths
                                    SizeAgeTrans = mlt_rg_sable_data$SizeAgeTrans, # size age transition matrix
                                    M_spec = "fix", # fix natural mortality
                                    ln_M = log(0.104884) # value to fix natural mortality at
)

# setting up movement parameterization
input_list <- Setup_Mod_Movement(input_list = input_list,
                                 # Model options
                                 Movement_ageblk_spec = list(c(1:6), c(7:15), c(16:30)), # estimating movement in 3 age blocks
                                 # (ages 1-6, ages 7-15, ages 16-30)
                                 Movement_yearblk_spec = "constant", # time-invariant for movement
                                 Movement_sexblk_spec = "constant", # sex-invariant movement
                                 do_recruits_move = 0, # recruits do not move
                                 use_fixed_movement = 0, # estimating movement
                                 Use_Movement_Prior = 1, # priors used for movement
                                 Movement_prior = 1.75 # vague prior to penalize movement away from the extremes
)

# setting up tagging parameterization
input_list <- Setup_Mod_Tagging(input_list = input_list,
                                UseTagging = 1, # using tagging data
                                max_tag_liberty = 15, # maximum number of years to track a cohort

                                # Data Inputs
                                tag_release_indicator = mlt_rg_sable_data$tag_release_indicator,
                                # tag release indicator (first col = tag region, second col = tag year),
                                # total number of rows = number of tagged cohorts
                                Tagged_Fish = mlt_rg_sable_data$Tagged_Fish, # Released fish
                                # dimensioned by total number of tagged cohorts, (implicitly
                                # tracks the release year and region), age, and sex
                                Obs_Tag_Recap = mlt_rg_sable_data$Obs_Tag_Recap,
                                # dimensioned by max tag liberty, tagged cohorts, regions,
                                # ages, and sexes

                                # Model options
                                Tag_LikeType = "NegBin", # Negative Binomial
                                mixing_period = 2, # Don't fit tagging until release year + 1
                                t_tagging = 0.5, # tagging happens midway through the year,
                                # movement does not occur within that year
                                tag_selex = "SexSp_AllFleet", # tagging recapture selectivity
                                # is a weighted average of fishery selectivity of two fleets
                                tag_natmort = "AgeSp_SexSp", # tagging natural mortality is
                                # age and sex-specific
                                Use_TagRep_Prior = 1, # tag reporting rate priors are used
                                TagRep_PriorType = 0, # symmetric beta is used for tag reporting prior
                                TagRep_mu = NA, # Tag reporting mean
                                TagRep_sd = 1, # tag reporting sd controlling symmetric beta
                                move_age_tag_pool = list(c(1:6), c(7:15), c(16:30)), # whether or
                                # not to pool tagging data when fitting (for computational cost)
                                move_sex_tag_pool = list(c(1:2)), # whether or not to pool
                                # sex-specific data when fitting
                                Init_Tag_Mort_spec = "fix", # fixing initial tag mortality
                                Tag_Shed_spec = "fix", # fixing chronic shedding
                                TagRep_spec = "est_shared_r", # tag reporting rates are not region specific
                                # Time blocks for tag reporting rates
                                Tag_Reporting_blocks = c(
                                  paste("Block_1_Year_1-35_Region_", c(1:input_list$data$n_regions), sep = ''),
                                  paste("Block_2_Year_36-62_Region_", c(1:input_list$data$n_regions), sep = '')
                                ),

                                # Specify starting values or fixing values
                                ln_Init_Tag_Mort = log(0.1), # fixing initial tag mortality
                                ln_Tag_Shed = log(0.02),  # fixing tag shedding
                                ln_tag_theta = log(0.5), # starting value for tagging overdispersion
                                Tag_Reporting_Pars = array(log(0.2 / (1-0.2)), # starting values for tag reporting pars
                                                           dim = c(input_list$data$n_regions, 3))
)

# Catch munging
mlt_rg_sable_data$ObsCatch[1,1:17,1] <- colSums(mlt_rg_sable_data$ObsCatch[,1:17,1])
mlt_rg_sable_data$ObsCatch[-1,1:17,1] <- 0
mlt_rg_sable_data$ObsCatch[1,1:17,2] <- colSums(mlt_rg_sable_data$ObsCatch[,1:17,2])
mlt_rg_sable_data$ObsCatch[-1,1:17,2] <- 0

# Set catch type to be aggregated
mlt_rg_sable_data$Catch_Type[1:17,] <- 0

# From Sasaki 1985, Table 16
ub_catch_censor <- array(NA, dim = c(5, 17, 2))
lb_catch_censor <- array(NA, dim = c(5, 17, 2))

# lower bound is from Sasaki
lb_catch_censor[1,5:17,1] <- c(2231, 1043, 2270, 1839, 310, 589, 1253, 1156, 268, 294, 799, 487, 673) / 1000
lb_catch_censor[2,5:17,1] <- c(888, 217, 1123, 783, 438, 320, 627, 602, 485, 798, 969, 1041, 1136) / 1000

ub_catch_censor[1,5:17,1] <- mlt_rg_sable_data$ObsCatch[1,5:17,1] # upper bound is aggregated catch
ub_catch_censor[2,5:17,1] <- mlt_rg_sable_data$ObsCatch[1,5:17,1] # upper bound is aggregated catch

lb_catch_censor[1,1:3,1] <- c(1861, 15627, 25989) / 1000
lb_catch_censor[2,1:3,1] <- c(1193, 451, 390) / 1000
ub_catch_censor[1,1:3,1] <- mlt_rg_sable_data$ObsCatch[1,1:3,1]
ub_catch_censor[2,1:3,1] <- mlt_rg_sable_data$ObsCatch[2,1:3,1]

# lb_catch_censor[3:5,8:17,1] <- 0 # lower bound here is 0
# ub_catch_censor[3,8:17,1] <- c(570, 10351, 14414, 20478, 20903, 26914, 21266, 20061, 18103, 17905) / 1000
# ub_catch_censor[4,8:17,1] <- c(570, 10351, 14414, 20478, 20903, 26914, 21266, 20061, 18103, 17905) / 1000
# ub_catch_censor[5,8:17,1] <- c(570, 10351, 14414, 20478, 20903, 26914, 21266, 20061, 18103, 17905) / 1000

# setting up catch data
input_list <- Setup_Mod_Catch_and_F(input_list = input_list,
                                    # Data inputs
                                    ObsCatch = mlt_rg_sable_data$ObsCatch,
                                    Catch_Type = mlt_rg_sable_data$Catch_Type,
                                    UseCatch = mlt_rg_sable_data$UseCatch,
                                    # Model options
                                    Use_F_pen = 1,
                                    est_all_regional_F = 1,
                                    # whether to use f penalty, == 0 don't use, == 1 use
                                    sigmaC_spec = 'fix',
                                    sigmaF_spec = "fix",
                                    ln_sigmaC =
                                      array(log(0.02), dim = c(input_list$data$n_regions,
                                                               input_list$data$n_fish_fleets)),
                                    # fixing catch sd at small value
                                    ln_F_mean = array(-2, dim = c(input_list$data$n_regions,
                                                                  input_list$data$n_fish_fleets)),
                                    # some starting values for fishing mortality
                                    censor_regional_catch = 0,
                                    ub_catch_censor = ub_catch_censor,
                                    lb_catch_censor = lb_catch_censor,
                                    catch_censor_sd = 0.025
                                    )

# Fishery Indices and Compositions
input_list <- Setup_Mod_FishIdx_and_Comps(input_list = input_list,
                                          # data inputs
                                          ObsFishIdx = mlt_rg_sable_data$ObsFishIdx,
                                          ObsFishIdx_SE = mlt_rg_sable_data$ObsFishIdx_SE,
                                          UseFishIdx =  mlt_rg_sable_data$UseFishIdx,
                                          ObsFishAgeComps = mlt_rg_sable_data$ObsFishAgeComps,
                                          UseFishAgeComps = mlt_rg_sable_data$UseFishAgeComps,
                                          ISS_FishAgeComps = mlt_rg_sable_data$ISS_FishAgeComps,
                                          ObsFishLenComps = mlt_rg_sable_data$ObsFishLenComps,
                                          UseFishLenComps = mlt_rg_sable_data$UseFishLenComps,
                                          ISS_FishLenComps = mlt_rg_sable_data$ISS_FishLenComps,

                                          # Model options
                                          fish_idx_type = c("none", "none"),
                                          # fishery indices not used
                                          FishAgeComps_LikeType = c("1d-Logistic-Normal", "none"),
                                          # age comp likelihoods for fishery fleet 1 and 2
                                          FishLenComps_LikeType = c("Multinomial", "Multinomial"),
                                          # length comp likelihoods for fishery fleet 1 and 2
                                          FishAgeComps_Type =
                                            c("spltRjntS_Year_1-62_Fleet_1",
                                              "none_Year_1-62_Fleet_2"),
                                          # age comp structure for fishery fleet 1 and 2
                                          FishLenComps_Type =
                                            c("spltRjntS_Year_1-62_Fleet_1",
                                              "spltRjntS_Year_1-62_Fleet_2"),
                                          # length comp structure for fishery fleet 1 and 2
                                          FishAge_comp_agg_type = c(NA,NA),
                                          # ADMB aggregation quirks, ideally get rid of this
                                          FishLen_comp_agg_type = c(NA,NA)
                                          # ADMB aggregation quirks, ideally get rid of this
)

# Survey Indices and Compositions
input_list <- Setup_Mod_SrvIdx_and_Comps(input_list = input_list,
                                         # data inputs
                                         ObsSrvIdx = mlt_rg_sable_data$ObsSrvIdx,
                                         ObsSrvIdx_SE = mlt_rg_sable_data$ObsSrvIdx_SE,
                                         UseSrvIdx =  mlt_rg_sable_data$UseSrvIdx,
                                         ObsSrvAgeComps = mlt_rg_sable_data$ObsSrvAgeComps,
                                         ISS_SrvAgeComps = mlt_rg_sable_data$ISS_SrvAgeComps,
                                         UseSrvAgeComps = mlt_rg_sable_data$UseSrvAgeComps,
                                         ObsSrvLenComps = mlt_rg_sable_data$ObsSrvLenComps,
                                         UseSrvLenComps = mlt_rg_sable_data$UseSrvLenComps,
                                         ISS_SrvLenComps = mlt_rg_sable_data$ISS_SrvLenComps,

                                         # Model options
                                         srv_idx_type = c("abd", "abd"),
                                         # abundance and biomass for survey fleet 1 and 2
                                         SrvAgeComps_LikeType =
                                           c("Multinomial", "Multinomial"),
                                         # survey age composition likelihood for survey fleet
                                         # 1, and 2
                                         SrvLenComps_LikeType =
                                           c("none", "none"),
                                         #  no length compositions used for survey
                                         SrvAgeComps_Type = c("spltRjntS_Year_1-62_Fleet_1",
                                                              "spltRjntS_Year_1-62_Fleet_2"),
                                         # survey age comp type
                                         SrvLenComps_Type = c("none_Year_1-62_Fleet_1",
                                                              "none_Year_1-62_Fleet_2"),
                                         # survey length comp type
                                         SrvAge_comp_agg_type = c(NA, NA),
                                         # ADMB aggregation quirks, ideally get rid of this
                                         SrvLen_comp_agg_type = c(NA, NA)
                                         # ADMB aggregation quirks, ideally get rid of this
)

# Fishery Selectivity and Catchability
input_list <- Setup_Mod_Fishsel_and_Q(input_list = input_list,

                                      # Model options
                                      cont_tv_fish_sel = c("none_Fleet_1", "none_Fleet_2"),
                                      # fishery selectivity, whether continuous time-varying

                                      # fishery selectivity blocks
                                      fish_sel_blocks =
                                        c("Block_1_Year_1-56_Fleet_1",
                                          # block 1, fishery ll selex
                                          "Block_2_Year_57-62_Fleet_1",
                                          # block 3 fishery ll selex
                                          "none_Fleet_2"),
                                      # no blocks for trawl fishery

                                      # fishery selectivity form
                                      fish_sel_model =
                                        c("logist1_Fleet_1",
                                          "gamma_Fleet_2"),

                                      # fishery catchability blocks
                                      fish_q_blocks =
                                        c("none_Fleet_1",
                                          "none_Fleet_2"),
                                      # no blocks since q is not estimated

                                      # whether to estimate all fixed effects
                                      # for fishery selectivity and later modify
                                      # to fix and share parameters
                                      fish_fixed_sel_pars =
                                        c("est_all", "est_all"),

                                      # whether to estimate all fixed effects
                                      # for fishery catchability
                                      fish_q_spec =
                                        c("fix", "fix")
                                      # fix fishery q since not used
)

# Custom parameter sharing for fishery selectivity
map_ln_fish_fixed_sel_pars <- input_list$par$ln_fish_fixed_sel_pars # mapping fishery selectivity

# Fixed gear fleet, unique parameters for each sex (time block 1)
map_ln_fish_fixed_sel_pars[,1,1,1,1] <- 1 # a50, female, time block 1, fixed gear
map_ln_fish_fixed_sel_pars[,2,1,1,1] <- 2 # delta, female, time block 1, fixed gear (shared with time block 2 and sex)
map_ln_fish_fixed_sel_pars[,1,1,2,1] <- 3 # a50, male, time block 1, fixed gear
map_ln_fish_fixed_sel_pars[,2,1,2,1] <- 2 # delta, male, time block 1, fixed gear (shared with time block 2 and sex)

# time block 2, fixed gear fishery
map_ln_fish_fixed_sel_pars[,1,2,1,1] <- 4 # a50, female, time block 2, fixed gear
map_ln_fish_fixed_sel_pars[,2,2,1,1] <- 2 # delta, female, time block 2, fixed gear (shared with time block 1 and sex)
map_ln_fish_fixed_sel_pars[,1,2,2,1] <- 5 # a50, male, time block 2, fixed gear
map_ln_fish_fixed_sel_pars[,2,2,2,1] <- 2 # delta, male, time block 2, fixed gear (shared with time block 1 and sex)

# time block 1 and 2, trawl gear fishery
map_ln_fish_fixed_sel_pars[,1,1,1,2] <- 6 # amax, female, time block 1, trawl gear
map_ln_fish_fixed_sel_pars[,2,1,1,2] <- 7 # delta, female, time block 1, trawl gear (shared by sex)
map_ln_fish_fixed_sel_pars[,1,1,2,2] <- 8 # amax, male, time block 1, trawl gear
map_ln_fish_fixed_sel_pars[,2,1,2,2] <- 7 # delta, male, time block 1, trawl gear (shared by sex)
map_ln_fish_fixed_sel_pars[,,2,,2] <- NA # no parameters estimated for time block 2 trawl gear

input_list$map$ln_fish_fixed_sel_pars <- factor(map_ln_fish_fixed_sel_pars) # input into map list
input_list$par$ln_fish_fixed_sel_pars[] <- log(5) # some more inforamtive starting values

# Survey Selectivity and Catchability
input_list <- Setup_Mod_Srvsel_and_Q(input_list = input_list,

                                     # Model options
                                     # survey selectivity, whether continuous time-varying
                                     cont_tv_srv_sel =
                                       c("none_Fleet_1",
                                         "none_Fleet_2"),

                                     # survey selectivity blocks
                                     srv_sel_blocks =
                                       c("none_Fleet_1",
                                         "none_Fleet_2"
                                       ), # no blocks for jp and domestic survey

                                     # survey selectivity form
                                     srv_sel_model =
                                       c("logist1_Fleet_1",
                                         "logist1_Fleet_2"),

                                     # survey catchability blocks
                                     srv_q_blocks =
                                       c("none_Fleet_1",
                                         "none_Fleet_2"),

                                     # whether to estiamte all fixed effects
                                     # for survey selectivity and later
                                     # modify to fix/share parameters
                                     srv_fixed_sel_pars_spec =
                                       c("est_all",
                                         "est_all"),

                                     # whether to estiamte all
                                     # fixed effects for survey catchability
                                     # spatially-invariant q
                                     srv_q_spec =
                                       c("est_shared_r",
                                         "est_shared_r"),

                                     # Starting values for survey catchability
                                     ln_srv_q = array(9,
                                                      dim = c(input_list$data$n_regions, 1,
                                                              input_list$data$n_srv_fleets))
)

# Custom mapping survey selectivity stuff
map_ln_srv_fixed_sel_pars <- input_list$par$ln_srv_fixed_sel_pars # set up mapping factor stuff

# Coop survey (japanese)
map_ln_srv_fixed_sel_pars[,1,1,1,1] <- 1 # a50, coop survey, time block 1, female
map_ln_srv_fixed_sel_pars[,2,1,1,1] <- 2 # delta, coop survey, time block 1, female (sharing with domestic survey)
map_ln_srv_fixed_sel_pars[,1,1,2,1] <- 3 # a50, coop survey, time block 1, male
map_ln_srv_fixed_sel_pars[,2,1,2,1] <- 4 # delta, coop survey, time block 1, male (sharing with domestic survey)

# domestic survey
map_ln_srv_fixed_sel_pars[,1,1,1,2] <- 5 # a50, domestic survey, time block 1, female
map_ln_srv_fixed_sel_pars[,2,1,1,2] <- 2 # delta, domestic survey, time block 1, female (sharing with coop survey)
map_ln_srv_fixed_sel_pars[,1,1,2,2] <- 6 # a50, domestic survey, time block 1, male
map_ln_srv_fixed_sel_pars[,2,1,2,2] <- 4 # delta, domestic survey, time block 1, male (sharing with coop survey)

input_list$map$ln_srv_fixed_sel_pars <- factor(map_ln_srv_fixed_sel_pars)  # input into map list
input_list$par$ln_srv_fixed_sel_pars[] <- log(5) # some more informative starting values


# set up model weighting stuff
input_list <- Setup_Mod_Weighting(input_list = input_list,
                                  sablefish_ADMB = 0, # don't use sablefish singl region ADMB quirks
                                  likelihoods = 1, # using TMB-style likelihoods,
                                  # and weight using sigmas, instead of lambdas and sigmas together
                                  Wt_Catch = 1,
                                  Wt_FishIdx = 1,
                                  Wt_SrvIdx = 1,
                                  Wt_Rec = 1,
                                  Wt_F = 2,
                                  # Composition model weighting
                                  Wt_FishAgeComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_fish_fleets)),
                                  Wt_FishLenComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_fish_fleets)),
                                  Wt_SrvAgeComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_srv_fleets)),
                                  Wt_SrvLenComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_srv_fleets))
)

# extract out lists updated with helper functions
data <- input_list$data
parameters <- input_list$par
mapping <- input_list$map

# parameters$ln_sigmaF[] <- log(0.2)

# Fit model
st <- Sys.time()
sabie_rtmb_model <- fit_model(data,
                              parameters,
                              mapping,
                              # random = NULL,
                              newton_loops = 3,
                              silent = FALSE
)
en <- Sys.time()
print(en - st)

sabie_rtmb_model$sd_rep <- sdreport(sabie_rtmb_model)

bind_rows(reshape2::melt(sabie_rtmb_model$rep$SSB) %>%
              mutate(type = 'interp') ) %>%
  ggplot(aes(x = Var2 + 1959, y = value, color = type)) +
  facet_wrap(~Var1, scales = 'free_y') +
  geom_line(lwd = 1.3)

bind_rows(reshape2::melt(sabie_rtmb_model$rep$Rec) %>%
            mutate(type = 'interp') ) %>%
  ggplot(aes(x = Var2 + 1959, y = value, color = type)) +
  facet_wrap(~Var1, scales = 'free_y') +
  geom_line(lwd = 1.3)

bind_rows(reshape2::melt(apply(sabie_rtmb_model$rep$FAA, c(1,2,5),max)) %>%
              mutate(type = 'interp') ) %>%
  ggplot(aes(x = Var2 + 1959, y = value)) +
  facet_grid(Var3~Var1, scales = 'free_y') +
  geom_line(lwd = 1.3) +
  geom_vline(xintercept = 1977, lty = 2)

bind_rows(reshape2::melt(sabie_rtmb_model$rep$PredCatch) %>%
              mutate(type = 'interp') ) %>%
  ggplot(aes(x = Var2 + 1959, y = value * 1000)) +
  facet_grid(Var3~Var1, scales = 'free_y') +
  geom_line(lwd = 1.3) +
  geom_vline(xintercept = 1977, lty = 2) +
  scale_y_continuous(labels = scales::comma)

reshape2::melt(sabie_rtmb_model$rep$fish_sel) %>%
  filter(Var2 == 60) %>%
  ggplot(aes(x = Var3, y = value)) +
  facet_grid(Var5~Var4, scales = 'free_y') +
  geom_line(lwd = 1.3) +
  coord_cartesian(ylim = c(0,1))

