library(SPoRC)
library(testthat)
data("three_rg_sable_data")

test_that("Three-region Sablefish RTMB model produces expected results", {

  # Initialize model dimensions and data list
  input_list <- Setup_Mod_Dim(years = 1:length(three_rg_sable_data$years), # vector of years (1 - 62)
                              ages = 1:30, # vector of ages (1 - 30)
                              lens = three_rg_sable_data$lens, # number of lengths (41 - 99)
                              n_regions = three_rg_sable_data$n_regions, # number of regions (5)
                              n_sexes = three_rg_sable_data$n_sexes, # number of sexes (2)
                              n_fish_fleets = three_rg_sable_data$n_fish_fleets, # number of fishery fleet (2)
                              n_srv_fleets = three_rg_sable_data$n_srv_fleets, # number of survey fleets (2)
                              verbose = FALSE
  )

  # Setup recruitment stuff (using defaults for other stuff)
  input_list <- Setup_Mod_Rec(input_list = input_list, # input data list from above
                              do_rec_bias_ramp = 0, # not using bias ramp
                              sigmaR_switch = 16, # switch to using late sigma in year 16
                              dont_est_recdev_last = 1, # don't estimate last rec dev

                              # Model options
                              rec_model = "mean_rec", # recruitment model
                              h_spec = 'fix',
                              sigmaR_spec = "fix", # fixing
                              InitDevs_spec = "est_shared_r", # initial deviations are shared across regions,
                              # but recruitment deviations are region specific
                              ln_sigmaR = log(c(0.4, 1.2)), # values to fix sigmaR at, or starting values
                              ln_global_R0 = log(20),
                              # starting value for global R0
                              R0_prop = array(c(0.2, 0.2),
                                              dim = c(input_list$data$n_regions - 1))
                              # starting value for R0 proportions in multinomial logit space
  )

  # Setup biological stuff (using defaults for other stuff)
  input_list <- Setup_Mod_Biologicals(input_list = input_list,
                                      WAA = three_rg_sable_data$WAA, # weight at age
                                      MatAA = three_rg_sable_data$MatAA, # maturity at age
                                      AgeingError = three_rg_sable_data$AgeingError, # ageing error matrix
                                      fit_lengths = 1, # fitting lengths
                                      SizeAgeTrans = three_rg_sable_data$SizeAgeTrans, # size age transition matrix
                                      M_spec = "fix", # fix natural mortality
                                      Fixed_natmort = array(0.104884, dim = c(three_rg_sable_data$n_regions,
                                                                              length(three_rg_sable_data$years),
                                                                              length(three_rg_sable_data$ages),
                                                                              three_rg_sable_data$n_sexes))
                                      # value to fix natural mortality at
  )

  # setting up movement parameterization
  Movement_prior <- expand.grid(
    region_from = 1:3, # regions
    age = c(6,7,16), # age blocks
    sex = 1, # sex
    alpha = I(list(rep(2.5, 3))) # prior alpha to each row
  )

  input_list <- Setup_Mod_Movement(input_list = input_list,
                                   # Model options
                                   Movement_ageblk_spec = list(c(1:6), c(7:15), c(16:30)), # estimating movement in 3 age blocks
                                   # (ages 1-6, ages 7-15, ages 16-30)
                                   Movement_yearblk_spec = "constant", # time-invariant for movement
                                   Movement_sexblk_spec = "constant", # sex-invariant movement
                                   do_recruits_move = 0, # recruits do not move
                                   use_fixed_movement = 0, # estimating movement
                                   Use_Movement_Prior = 1, # priors used for movement
                                   Movement_prior = Movement_prior, # vague prior to penalize movement away from the extremes
                                   cont_vary_movement = 'none'
  )

  # setting up tagging parameterization
  # setup tagging priors
  tag_prior <- data.frame(
    region = 1,
    block = c(1,2),
    mu = NA, # no mean, since symmetric beta
    sd = 5, # sd = 5
    type = 0 # symmetric beta
  )

  input_list <- Setup_Mod_Tagging(input_list = input_list,
                                  UseTagging = 1, # using tagging data
                                  max_tag_liberty = 15, # maximum number of years to track a cohort

                                  # Data Inputs
                                  tag_release_indicator = three_rg_sable_data$tag_release_indicator,
                                  # tag release indicator (first col = tag region, second col = tag year),
                                  # total number of rows = number of tagged cohorts
                                  Tagged_Fish = three_rg_sable_data$Tagged_Fish, # Released fish
                                  # dimensioned by total number of tagged cohorts, (implicitly
                                  # tracks the release year and region), age, and sex
                                  Obs_Tag_Recap = three_rg_sable_data$Obs_Tag_Recap,
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
                                  TagRep_Prior = tag_prior,
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
                                    paste("Block_2_Year_36-terminal_Region_", c(1:input_list$data$n_regions), sep = '')
                                  ),

                                  # Specify starting values or fixing values
                                  ln_Init_Tag_Mort = log(0.1), # fixing initial tag mortality
                                  ln_Tag_Shed = log(0.02),  # fixing tag shedding
                                  ln_tag_theta = log(0.5), # starting value for tagging overdispersion
                                  Tag_Reporting_Pars = array(log(0.2 / (1-0.2)), # starting values for tag reporting pars
                                                             dim = c(input_list$data$n_regions, 2))
  )


  # setting up catch data
  input_list <- Setup_Mod_Catch_and_F(input_list = input_list,
                                      # Data inputs
                                      ObsCatch = three_rg_sable_data$ObsCatch,
                                      Catch_Type = three_rg_sable_data$Catch_Type,
                                      UseCatch = three_rg_sable_data$UseCatch,
                                      # Model options
                                      Use_F_pen = 1,
                                      # whether to use f penalty, == 0 don't use, == 1 use
                                      sigmaC_spec = 'fix',
                                      ln_sigmaC =
                                        array(log(0.05), dim = c(input_list$data$n_regions,
                                                                 length(input_list$data$years),
                                                                 input_list$data$n_fish_fleets)),
                                      # fixing catch sd at small value
                                      ln_F_mean = array(-2, dim = c(input_list$data$n_regions,
                                                                    input_list$data$n_fish_fleets))
                                      # some starting values for fishing mortality
  )

  # Fishery Indices and Compositions
  input_list <- Setup_Mod_FishIdx_and_Comps(input_list = input_list,
                                            # data inputs
                                            ObsFishIdx = three_rg_sable_data$ObsFishIdx,
                                            ObsFishIdx_SE = three_rg_sable_data$ObsFishIdx_SE,
                                            UseFishIdx =  three_rg_sable_data$UseFishIdx,
                                            ObsFishAgeComps = three_rg_sable_data$ObsFishAgeComps,
                                            UseFishAgeComps = three_rg_sable_data$UseFishAgeComps,
                                            ISS_FishAgeComps = three_rg_sable_data$ISS_FishAgeComps,
                                            ObsFishLenComps = three_rg_sable_data$ObsFishLenComps,
                                            UseFishLenComps = three_rg_sable_data$UseFishLenComps,
                                            ISS_FishLenComps = three_rg_sable_data$ISS_FishLenComps,

                                            # Model options
                                            fish_idx_type = c("none", "none"),
                                            # fishery indices not used
                                            FishAgeComps_LikeType = c("Multinomial", "none"),
                                            # age comp likelihoods for fishery fleet 1 and 2
                                            FishLenComps_LikeType = c("Multinomial", "Multinomial"),
                                            # length comp likelihoods for fishery fleet 1 and 2
                                            FishAgeComps_Type =
                                              c("spltRjntS_Year_1-terminal_Fleet_1",
                                                "none_Year_1-terminal_Fleet_2"),
                                            # age comp structure for fishery fleet 1 and 2
                                            FishLenComps_Type =
                                              c("spltRjntS_Year_1-terminal_Fleet_1",
                                                "spltRjntS_Year_1-terminal_Fleet_2")
                                            # length comp structure for fishery fleet 1 and 2
  )

  # Survey Indices and Compositions
  input_list <- Setup_Mod_SrvIdx_and_Comps(input_list = input_list,
                                           # data inputs
                                           ObsSrvIdx = three_rg_sable_data$ObsSrvIdx,
                                           ObsSrvIdx_SE = three_rg_sable_data$ObsSrvIdx_SE,
                                           UseSrvIdx =  three_rg_sable_data$UseSrvIdx,
                                           ObsSrvAgeComps = three_rg_sable_data$ObsSrvAgeComps,
                                           ISS_SrvAgeComps = three_rg_sable_data$ISS_SrvAgeComps,
                                           UseSrvAgeComps = three_rg_sable_data$UseSrvAgeComps,
                                           ObsSrvLenComps = three_rg_sable_data$ObsSrvLenComps,
                                           UseSrvLenComps = three_rg_sable_data$UseSrvLenComps,
                                           ISS_SrvLenComps = three_rg_sable_data$ISS_SrvLenComps,

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
                                           SrvAgeComps_Type = c("spltRjntS_Year_1-terminal_Fleet_1",
                                                                "spltRjntS_Year_1-terminal_Fleet_2"),
                                           # survey age comp type
                                           SrvLenComps_Type = c("none_Year_1-terminal_Fleet_1",
                                                                "none_Year_1-terminal_Fleet_2")
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
                                            "Block_2_Year_57-terminal_Fleet_1",
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
                                        fish_fixed_sel_pars_spec =
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
  map_ln_fish_fixed_sel_pars[,2,1,2,1] <- 4 # delta, male, time block 1, fixed gear (shared with time block 2 and sex)

  # time block 2, fixed gear fishery
  map_ln_fish_fixed_sel_pars[,1,2,1,1] <- 5 # a50, female, time block 2, fixed gear
  map_ln_fish_fixed_sel_pars[,2,2,1,1] <- 2 # delta, female, time block 2, fixed gear (shared with time block 1 and sex)
  map_ln_fish_fixed_sel_pars[,1,2,2,1] <- 6 # a50, male, time block 2, fixed gear
  map_ln_fish_fixed_sel_pars[,2,2,2,1] <- 4 # delta, male, time block 2, fixed gear (shared with time block 1 and sex)

  # time block 1 and 2, trawl gear fishery
  map_ln_fish_fixed_sel_pars[,1,1,1,2] <- 7 # amax, female, time block 1, trawl gear
  map_ln_fish_fixed_sel_pars[,2,1,1,2] <- 8 # delta, female, time block 1, trawl gear (shared by sex)
  map_ln_fish_fixed_sel_pars[,1,1,2,2] <- 9 # amax, male, time block 1, trawl gear
  map_ln_fish_fixed_sel_pars[,2,1,2,2] <- 8 # delta, male, time block 1, trawl gear (shared by sex)
  map_ln_fish_fixed_sel_pars[,,2,,2] <- NA # no parameters estimated for time block 2 trawl gear

  input_list$map$ln_fish_fixed_sel_pars <- factor(map_ln_fish_fixed_sel_pars) # input into map list
  input_list$par$ln_fish_fixed_sel_pars[] <- log(0.1) # some more inforamtive starting values

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
                                         c("est_shared_r", "est_shared_r"),

                                       # Starting values for survey catchability
                                       ln_srv_q = array(8.75,
                                                        dim = c(input_list$data$n_regions, 1,
                                                                input_list$data$n_srv_fleets))
  )

  # Custom mapping survey selectivity stuff
  map_ln_srv_fixed_sel_pars <- input_list$par$ln_srv_fixed_sel_pars # set up mapping factor stuff

  # Coop survey (japanese)
  map_ln_srv_fixed_sel_pars[,1,1,1,1] <- 1 # a50, coop survey, time block 1, female
  map_ln_srv_fixed_sel_pars[,2,1,1,1] <- 2 # delta, coop survey, time block 1, female (sharing with domestic survey)
  map_ln_srv_fixed_sel_pars[,1,1,2,1] <- 3 # a50, coop survey, time block 1, male
  map_ln_srv_fixed_sel_pars[,2,1,2,1] <- 2 # delta, coop survey, time block 1, male (sharing with domestic survey)

  # domestic survey
  map_ln_srv_fixed_sel_pars[,1,1,1,2] <- 5 # a50, domestic survey, time block 1, female
  map_ln_srv_fixed_sel_pars[,2,1,1,2] <- 2 # delta, domestic survey, time block 1, female (sharing with coop survey)
  map_ln_srv_fixed_sel_pars[,1,1,2,2] <- 6 # a50, domestic survey, time block 1, male
  map_ln_srv_fixed_sel_pars[,2,1,2,2] <- 2 # delta, domestic survey, time block 1, male (sharing with coop survey)

  input_list$map$ln_srv_fixed_sel_pars <- factor(map_ln_srv_fixed_sel_pars)  # input into map list
  input_list$par$ln_srv_fixed_sel_pars[] <- log(0.1) # some more informative starting values


  # set up model weighting stuff
  input_list <- Setup_Mod_Weighting(input_list = input_list,
                                    Wt_Catch = 1,
                                    Wt_FishIdx = 1,
                                    Wt_SrvIdx = 1,
                                    Wt_Rec = 1,
                                    Wt_F = 1,
                                    Wt_Tagging = 1,
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

  # Pre-specify ISS
  data$ISS_FishAgeComps[] <- 50
  data$ISS_FishLenComps[] <- 40
  data$ISS_SrvAgeComps[] <- 60

  sabie_rtmb_model <- fit_model(data,
                                parameters,
                                mapping,
                                random = NULL,
                                newton_loops = 3,
                                silent = TRUE)

  # Save model results
  sabie_rtmb_model$sd_rep <- RTMB::sdreport(sabie_rtmb_model)

  ssb_expected_vec <- c(
    41.31334, 112.03276, 92.32823, 40.95897,
    111.17527, 91.43794, 39.41037, 107.68528,
    87.73502, 37.01079, 102.02302, 81.96514,
    35.59069, 98.41208, 79.47796, 35.30646,
    96.82819, 78.69146, 35.08223, 95.40563,
    78.10767, 34.45520, 93.11114, 76.56564,
    33.76643, 90.59166, 74.79041, 32.20515,
    86.27259, 71.18591, 30.23384, 81.10761,
    66.69770, 28.24654, 75.62611, 61.23012,
    25.56714, 69.23012, 55.67898, 21.85469,
    60.98703, 48.72356, 19.64360, 55.22301,
    43.56300, 17.59444, 49.60068, 38.29287,
    15.89943, 44.65376, 33.81197, 13.95575,
    39.21270, 29.02038, 12.41632, 36.41166,
    25.31809, 12.35180, 34.33274, 24.69786,
    13.02970, 32.33558, 23.34568, 14.34337,
    31.31804, 22.98218, 17.31755, 31.31101,
    22.64579, 22.98803, 33.96925, 23.50783,
    30.45777, 39.31391, 25.18314, 36.85791,
    46.66771, 29.37425, 39.92687, 55.55482,
    36.58013, 42.86800, 62.90875, 41.45861,
    44.81373, 69.04212, 44.73558, 41.98176,
    73.91971, 48.76300, 37.19320, 76.67455,
    52.17302, 33.05445, 76.00453, 52.62515,
    28.44692, 73.87178, 52.61725, 24.62396,
    70.64786, 51.15683, 22.01163, 66.01084,
    48.05534, 20.33407, 62.26051, 44.56067,
    18.84145, 58.99149, 41.78472, 17.72453,
    56.10474, 39.84279, 16.76462, 53.31921,
    38.73848, 16.19008, 50.73804, 37.21608,
    15.74297, 48.52365, 35.86079, 15.68860,
    46.53411, 33.79828, 16.15540, 46.18020,
    32.47626, 16.69367, 46.91151, 32.26707,
    17.68851, 48.11924, 32.16655, 18.50490,
    49.56069, 32.68358, 19.00383, 50.75666,
    33.94823, 19.53500, 51.69528, 34.80850,
    18.62709, 51.69105, 35.77683, 17.52899,
    50.92796, 36.11574, 16.54400, 49.59006,
    35.90766, 15.84797, 48.24508, 35.61607,
    15.49644, 46.56341, 34.61231, 15.23844,
    44.69829, 33.34705, 15.21171, 42.45868,
    31.39539, 14.99489, 40.38498, 30.09393,
    14.78778, 38.60958, 28.85505, 14.61202,
    37.17293, 27.90553, 15.71122, 36.84435,
    26.85475, 18.27793, 38.18578, 26.18840,
    23.12911, 42.46053, 26.94233, 30.07849,
    50.58409, 29.60340
  )

  rec_expected_vec <- c(
    9.6148867, 4.7628042, 2.6240565, 9.8152715,
    4.8122609, 2.6392956, 10.3084691, 4.9062933,
    2.6692695, 11.2873663, 5.0545166, 2.7171714,
    12.5870065, 5.2117803, 2.7658632, 14.0837548,
    5.3356052, 2.8075168, 15.3125159, 5.3874395,
    2.8289695, 15.2767541, 5.3153180, 2.8166022,
    14.2783119, 5.1464349, 2.7769095, 12.9401836,
    4.9052797, 2.7141051, 11.6751789, 4.6184366,
    2.6348140, 10.5300742, 4.2853360, 2.5350707,
    9.4087368, 3.8892240, 2.4059369, 8.4794735,
    3.4628713, 2.2636538, 8.2162219, 3.0666051,
    2.1299754, 3.2581110, 0.6988721, 0.6573834,
    4.9556619, 0.5842927, 0.6181444, 5.1883982,
    0.4982642, 0.6035069, 7.2912508, 0.5755407,
    0.7320428, 32.8433942, 0.8953928, 16.0971020,
    64.5291811, 1.0576272, 2.1431233, 3.9182895,
    0.8848496, 0.7160814, 19.2758958, 1.0128301,
    0.6689465, 85.2569618, 1.0189279, 0.7389304,
    33.3994656, 1.1411407, 0.9323679, 3.2782052,
    1.6011670, 0.9420896, 22.8645395, 10.8880136,
    0.9495787, 3.5951812, 6.1150578, 0.8634550,
    1.2750286, 2.8875786, 0.9550809, 1.0018263,
    12.4500175, 2.4370427, 0.8746270, 2.2597593,
    8.6965781, 1.2789112, 3.9514055, 14.4258385,
    4.0177584, 7.5509891, 1.5711284, 0.7996369,
    5.4826663, 2.9384787, 0.6768236, 10.4631048,
    3.0536955, 0.9335852, 2.3183874, 1.7352583,
    3.9893154, 2.1945989, 4.6834187, 7.7965838,
    6.2460951, 2.1524309, 3.1770034, 6.1593914,
    1.0547154, 9.7888633, 18.0359826, 0.8012966,
    1.8300217, 23.3820131, 1.4059043, 2.4862392,
    3.1247399, 3.3237795, 17.2191738, 17.4358114,
    1.2562539, 5.6582553, 6.0157993, 0.7718985,
    2.5049310, 3.1933669, 2.9156549, 3.2812893,
    6.0596172, 2.7902627, 3.0436321, 2.8755509,
    1.3909855, 2.6698167, 5.5923299, 3.8386324,
    1.4219191, 2.0031776, 2.6153778, 6.4917936,
    2.3274801, 6.8831119, 8.3406690, 4.0942803,
    1.2690962, 3.8437368, 1.2687131, 0.9711480,
    4.9356866, 1.1537111, 1.4220850, 0.8784587,
    0.9127894, 1.9236057, 3.6254646, 1.5864582,
    2.1099395, 6.2957189, 3.8180067, 0.6894480,
    32.4643659, 22.1930391, 1.6658426, 6.3974593,
    5.1580625, 2.6236741, 78.3922492, 54.7443987,
    4.7630196, 15.6643769, 22.3552314, 0.786718

  )

  expect_equal(as.vector(sabie_rtmb_model$rep$SSB), ssb_expected_vec, tolerance = 1e-2)
  expect_equal(as.vector(sabie_rtmb_model$rep$Rec), rec_expected_vec, tolerance = 1e-2)
  expect_true(sabie_rtmb_model$sd_rep$pdHess)

})

