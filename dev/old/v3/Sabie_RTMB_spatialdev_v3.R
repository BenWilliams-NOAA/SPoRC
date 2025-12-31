# Purpose: To do model development for a spatial model
# Creator: Matthew LH. Cheng
# Date Created: 12/24/24

# Set up ------------------------------------------------------------------

library(here)
library(R2admb)
library(tidyverse)
library(RTMB)
library(SPoRC)

sim_out = readRDS(here("sim_out.RDS"))

# Setting up with simulated data ---------------------------------------------
n_sims <- 500
ssb_mat <- array(NA, dim = c(length(1:10), 2, n_sims))
rec_mat <- array(NA, dim = c(length(1:10), 2, n_sims))

r0_mat <- matrix(NA, nrow = n_sims, ncol = 2)
tagrep <- vector()
status <- vector()
dir <- vector()
pd = vector()
grad = vector()
move_cv_mat = matrix(NA, nrow = 20, ncol = n_sims)

# Testing with simulated data ---------------------------------------------
for(sim in 1:n_sims) {
  # Prepare Data ------------------------------------------------------------

  # Initialize model dimensions and data list
  input_list <- Setup_Mod_Dim(years = 1:10, # vector of years
                              ages = 1:8, # vector of ages
                              lens = 1, # number of lengths
                              n_regions = 2, # number of regions
                              n_sexes = 1, # number of sexes
                              n_fish_fleets = 1, # number of fishery fleet
                              n_srv_fleets = 1, verbose = F # number of survey fleets
  )

  inv_steepness <- function(s) {
    qlogis((s - 0.2) / 0.8)
  }

  # Setup recruitment stuff (using defaults for other stuff)
  input_list <- Setup_Mod_Rec(input_list = input_list, # input data list from above

                              # Model options
                              ln_sigmaR = log(c(0.5, 0.5)),
                              rec_model = "bh_rec", # recruitment model
                              sigmaR_spec = "fix",
                              InitDevs_spec = "est_shared_r",
                              RecDevs_spec = "est_shared_r",
                              rec_dd = "global",
                              rec_lag = 1,
                              h_spec = "fix",
                              steepness_h = inv_steepness(c(0.8, 0.8))
  )

  # Setup biological stuff (using defaults for other stuff)
  input_list <- Setup_Mod_Biologicals(input_list = input_list,
                                      WAA = array(sim_out$WAA[,,,,sim, drop = FALSE], dim = c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes)),
                                      MatAA = array(sim_out$Maturity_AA[,,,,sim, drop = FALSE], dim = c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes)),
                                      ln_M = log(0.5)
  )

  # Setup movement stuff (using defaults for other stuff)
  input_list <- Setup_Mod_Movement(input_list = input_list,

                                   # Model options
                                   Movement_ageblk_spec = 'constant',
                                   Movement_yearblk_spec = 'constant',
                                   Movement_sexblk_spec = "constant",
                                   do_recruits_move = 0,
                                   cont_vary_movement = 'none'
  )

  # Setup catch and fishing mortality stuff
  input_list <- Setup_Mod_Catch_and_F(input_list = input_list,
                                      ObsCatch = array(sim_out$Obs_Catch[,,,sim, drop = FALSE],dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)),
                                      Catch_Type = array(1, dim = c(length(input_list$data$years), input_list$data$n_fish_fleets)),
                                      UseCatch = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)),

                                      # Model options
                                      sigmaC_spec = 'fix',
                                      sigmaF_spec = "fix"
                                      )

  # Setup fishery indices and compositions
  input_list <- Setup_Mod_FishIdx_and_Comps(input_list = input_list,
                                            ObsFishIdx = array(NA, c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)),
                                            ObsFishIdx_SE = array(NA, c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)),
                                            UseFishIdx =  array(0, c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)),
                                            ObsFishAgeComps = array(sim_out$Obs_FishAgeComps[,,,,,sim,drop = FALSE], dim = c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                            UseFishAgeComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)),
                                            ObsFishLenComps = array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$lens), input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                            UseFishLenComps = array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)),

                                            # Model options
                                            fish_idx_type = "none",
                                            FishAgeComps_LikeType = "Multinomial",
                                            FishLenComps_LikeType = "none",
                                            FishAgeComps_Type =  "spltRjntS_Year_1-terminal_Fleet_1",
                                            FishLenComps_Type =  "none_Year_1-terminal_Fleet_1"
  )


  # Setup survey indices and compositions
  input_list <- Setup_Mod_SrvIdx_and_Comps(input_list = input_list,
                                           ObsSrvIdx = array(sim_out$Obs_SrvIdx[,,,sim,drop = FALSE], dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets)),
                                           ObsSrvIdx_SE = array(0.2, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets)),
                                           UseSrvIdx =  array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)),
                                           ObsSrvAgeComps = array(sim_out$Obs_SrvAgeComps[,,,,,sim,drop = FALSE],  dim = c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes, input_list$data$n_srv_fleets)),
                                           UseSrvAgeComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)),
                                           ObsSrvLenComps = array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$lens), input_list$data$n_sexes, input_list$data$n_srv_fleets)),
                                           UseSrvLenComps = array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets)),

                                           # Model options
                                           srv_idx_type = "abd",
                                           SrvAgeComps_LikeType = "Multinomial",
                                           SrvLenComps_LikeType = "none",
                                           SrvAgeComps_Type = "spltRjntS_Year_1-terminal_Fleet_1",
                                           SrvLenComps_Type = "none_Year_1-terminal_Fleet_1"
  )

  # Setup fishery selectivity and catchability
  input_list <- Setup_Mod_Fishsel_and_Q(input_list = input_list,
                                        # Model options
                                        cont_tv_fish_sel = "none_Fleet_1",
                                        fish_sel_blocks = "none_Fleet_1",
                                        fish_sel_model = "logist1_Fleet_1",
                                        fish_q_blocks = "none_Fleet_1",
                                        fish_fixed_sel_pars = "est_shared_r"
  )

  # Setup survey selectivity and catchability
  input_list <- Setup_Mod_Srvsel_and_Q(input_list = input_list,
                                       cont_tv_srv_sel = "none_Fleet_1",
                                       srv_sel_blocks = "none_Fleet_1",
                                       srv_sel_model = "logist1_Fleet_1",
                                       srv_q_blocks = "none_Fleet_1",
                                       srv_fixed_sel_pars_spec = "est_shared_r",
                                       srv_q_spec = "est_shared_r"
  )
  # Setup tagging stuff
  input_list <- Setup_Mod_Tagging(input_list = input_list,
                                  tag_release_indicator = sim_out$Tag_Release_Ind,
                                  Tagged_Fish = array(sim_out$Tag_Fish[,,,sim], dim = c(dim(sim_out$Tag_Fish[,,,sim]), input_list$data$n_sexes)),
                                  Obs_Tag_Recap = array(sim_out$Obs_Tag_Recap[,,,,,sim], dim = c(dim(sim_out$Obs_Tag_Recap[,,,,,sim]), input_list$data$n_sexes)),

                                  # Model options
                                  UseTagging = 1,
                                  t_tagging = 0.5,
                                  max_tag_liberty = 30,
                                  Tag_LikeType = 'Poisson',
                                  tag_selex = "SexSp_DomFleet",
                                  tag_natmort = "AgeSp_SexSp",
                                  move_age_tag_pool = "all",
                                  move_sex_tag_pool = "all", # pooling all sex data
                                  Init_Tag_Mort_spec = "fix",
                                  Tag_Shed_spec = "fix",
                                  TagRep_spec = "est_shared_r",
                                  Tag_Reporting_blocks = c("none_Region_1", "none_Region_2")
                                  )

  input_list <- Setup_Mod_Weighting(input_list = input_list,
                                    sablefish_ADMB = 0,
                                    likelihoods = 1,
                                    Wt_Catch = 1,
                                    Wt_FishIdx = 1,
                                    Wt_SrvIdx = 1,
                                    Wt_Rec = 1,
                                    Wt_F = 1,
                                    Wt_FishAgeComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                    Wt_FishLenComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                    Wt_SrvAgeComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes, input_list$data$n_srv_fleets)),
                                    Wt_SrvLenComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes, input_list$data$n_srv_fleets))
  )

  data <- input_list$data
  parameters <- input_list$par
  mapping <- input_list$map

  # Do 3D Correlations
  # data$cont_vary_movement <- 0
  # map_pe <- parameters$move_pe_pars
  # map_pe[] <- NA
  # # map_pe[1,1:4] <- 1:4
  # mapping$move_pe_pars <- factor(map_pe)

  # Mapping -----------------------------------------------------------------
  # rho_trans = function(x) 2 / (1 + exp(-2 * x)) - 1
  # inv_rho_trans = function(y) -0.5 * log((1 - y) / (1 + y))
  # parameters$move_pe_pars[1,1:4] <- c(inv_rho_trans(0), inv_rho_trans(0), inv_rho_trans(0), log(0.05))


  sabie_rtmb_model <- fit_model(data,
                                parameters,
                                mapping,
                                random = NULL,
                                newton_loops = 3,
                                silent = T
                                # random = 'logit_move_devs'
                                )

  # for(j in 1:5) {
  #
  #   if(j == 1) { # reset weights at 1
  #     data$Wt_FishAgeComps[] <- 1
  #     data$Wt_FishLenComps[] <- 1
  #     data$Wt_SrvAgeComps[] <- 1
  #     data$Wt_SrvLenComps[] <- 1
  #   } else {
  #     data$Wt_FishAgeComps[] <- wts$new_fish_age_wts
  #     data$Wt_FishLenComps[] <- wts$new_fish_len_wts
  #     data$Wt_SrvAgeComps[] <- wts$new_srv_age_wts
  #     data$Wt_SrvLenComps[] <- wts$new_srv_len_wts
  #   }
  #
  #   sabie_rtmb_model <- fit_model(data,
  #                                 parameters,
  #                                 mapping,
  #                                 random = NULL,
  #                                 newton_loops = 3,
  #                                 silent = TRUE
  #   )
  #
  #   rep <- sabie_rtmb_model$report(sabie_rtmb_model$env$last.par.best) # Get report
  #   wts <- do_francis_reweighting(data = data, rep = rep, age_labels = 1:8,
  #                                 len_labels = NA, year_labels = 1:10)
  # }

  sabie_rtmb_model$sd_rep <- RTMB::sdreport(sabie_rtmb_model)

  prop <- sabie_rtmb_model$sd_rep$par.fixed[names(sabie_rtmb_model$sd_rep$par.fixed) == 'Rec_prop']
  Rec_prop <- exp(c(0, prop)) / sum(exp(c(0, prop)) )
  r0_mat[sim,] <- sabie_rtmb_model$rep$R0 * Rec_prop
  tagrep[sim] <- mean(sabie_rtmb_model$rep$Tag_Reporting)
  PD = sabie_rtmb_model$sd_rep$pdHess
  gradient = max(sabie_rtmb_model$sd_rep$gradient.fixed)
  pd[sim] = PD
  grad[sim] = gradient
  if(PD == TRUE && gradient < 0.005) status[sim] = TRUE else status[sim] = FALSE
  ssb_mat[,1,sim] <- (sabie_rtmb_model$rep$SSB[1,] - sim_out$SSB[1,,sim]) / sim_out$SSB[1,,sim]
  ssb_mat[,2,sim] <- (sabie_rtmb_model$rep$SSB[2,] - sim_out$SSB[2,,sim]) / sim_out$SSB[2,,sim]

  print(PD)
  print(gradient)

  par(mfrow = c(2,4))
  # plot(sim_out$CAA[20,2,,1,1,sim], type = 'l', col = 'red')
  # lines(sabie_rtmb_model$rep$CAA[2,20,,1,1], type = 'l')

  plot(sim_out$Total_Biom[1,,sim], col = 'red', type = 'l', ylab = 'Biomass region 1 (red = simulation, black = est)')
  lines(sabie_rtmb_model$rep$Total_Biom[1,], type = 'l')

  plot(sim_out$Total_Biom[2,,sim], col = 'red', type = 'l', ylab = 'Biomass region 2 (red = simulation, black = est)')
  lines(sabie_rtmb_model$rep$Total_Biom[2,], type = 'l')

  plot(sabie_rtmb_model$rep$Fmort[1,,1], type = 'l', ylab = 'Fmort Region 1 (red = simulation, black = est)', xlab = 'Year')
  lines(sim_out$Fmort[1,-51,1,sim], col = 'red')

  plot(sim_out$fish_sel[1,1,,1,1,sim], col = 'red')
  lines(sabie_rtmb_model$rep$fish_sel[1,1,,1,1])

  # par(mfrow = c(1,3))
  hist((r0_mat[,1] - 50) / 50, main = round(median((r0_mat[,1] - 50) / 50, na.rm = T),2), xlab = 'R0 bias region 1')
  hist((r0_mat[,2] - 50) / 50, main = round(median((r0_mat[,2] - 50) / 50, na.rm = T),2), xlab = 'R0 bias region 2')
  hist((tagrep - 0.2) / 0.2, main = round(median((tagrep - 0.2) / 0.2, na.rm = T),2), xlab = 'R0 bias region 2')

}

reshape2::melt(ssb_mat) %>%
  drop_na() %>%
  rename(Year = Var1, Region = Var2, Sim = Var3) %>%
  ggplot(aes(x = Year, y = value, group = Sim)) +
  geom_line() +
  # geom_ribbon(alpha = 0.2) +
  geom_hline(yintercept = 0, lty = 2) +
  coord_cartesian(ylim = c(-1, 1)) +
  facet_wrap(~Region) +
  theme_bw(base_size = 18) +
  labs(title = 'NegBin Tag DM Age')
