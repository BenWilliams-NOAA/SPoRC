# Purpose: To demonstrate how to run a closed loop simulation
# Creator: Matthew LH. Cheng
# Date: 5/21/25

# unloadNamespace("SPoRC")
# library(SPoRC)
library(ggplot2)

# Define Operating Model (Simulation Options) -----------------------------

# First, we need to define the operating model we want to condition on
# Set up model dimensions
sim_list <- Setup_Sim_Dim(n_sims = 5,
                          n_yrs = 50,
                          n_regions = 2,
                          n_ages = 8,
                          n_sexes = 1,
                          n_fish_fleets = 1,
                          n_srv_fleets = 1,
                          run_feedback = T,
                          feedback_start_yr = 25
                          )

# set up containers
sim_list <- Setup_Sim_Containers(sim_list)

# Setup fishing mortality
sim_list <- Setup_Sim_FishMort(sim_list = sim_list,
                               sigmaC = 1e-3,
                               init_F_val = 0,
                               Fmort_pattern = matrix(c('one-way', "one-way"), nrow = sim_list$n_regions, ncol = sim_list$n_fish_fleets),
                               Fmort_start = matrix(c(0.01, 0.01), nrow = sim_list$n_regions, ncol = sim_list$n_fish_fleets),
                               Fmort_fct = matrix(c(50, 50), nrow = sim_list$n_regions, ncol = sim_list$n_fish_fleets),
                               proc_error = TRUE,
                               proc_error_sd = 0
)

# Setup fishery selectivity
sim_list <- Setup_Sim_FishSel(sel_model = matrix(c('logistic', "logistic"), nrow = sim_list$n_regions, ncol = sim_list$n_fish_fleets),
                              # a50, k for logistic shared across regions
                              fixed_fish_sel_pars = array(c(3,3,1,1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
                              sim_list = sim_list
)

# Setup survey catchability and selectivity
sim_list <- Setup_Sim_Survey(sim_list = sim_list,
                             ObsSrvIdx_SE = array(0.05, dim = c(sim_list$n_regions, sim_list$n_srv_fleets)), # survey observation error
                             base_srv_q = array(1, dim = c(sim_list$n_regions, sim_list$n_srv_fleets)), # base survey catchability value
                             srv_q_pattern = matrix(c('constant', "constant"), nrow = sim_list$n_regions, ncol = sim_list$n_srv_fleets), # catchability pattern
                             sel_model = matrix(c('logistic', "logistic"), nrow = sim_list$n_regions, ncol = sim_list$n_srv_fleets), # selectivity model
                             # a50, k, for logistic shared across regions
                             fixed_srv_sel_pars = array(c(3,3,1,1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets, 2))
)

# Setup recruitment stuff
sim_list <- Setup_Sim_Rec(
  sim_list = sim_list,
  do_recruits_move = "dont_move", # == 0, recruits don't move , == 1 recruits move
  base_sexratio = 1, # single sex
  sexratio_vary = "constant",
  base_r0 = c(100, 50),
  r0_vary = "constant",
  base_h = c(0.8, 0.8),
  init_sigmaR = 0.5,
  sigmaR = 0.5,
  recruitment_opt = "bh_rec",
  rec_dd = "global",
  init_dd = "global",
  rec_lag = 1
)

# Setup biologicals
sim_list <- Setup_Sim_Biologicals(
  sim_list = sim_list,
  base_natmort_value = array(0.5, dim = c(sim_list$n_regions, sim_list$n_ages, sim_list$n_sims)),
  natmort_pattern = "constant",
  base_WAA_values = array(rep(5 / (1 + exp(-1 * (1:sim_list$n_ages - 5))), each = sim_list$n_regions * sim_list$n_sexes),
                          dim = c(sim_list$n_regions, sim_list$n_ages, sim_list$n_sexes)),
  base_WAA_fish_values = array(rep(5 / (1 + exp(-1 * (1:sim_list$n_ages - 5))), each = sim_list$n_regions * sim_list$n_sexes),
                               dim = c(sim_list$n_regions, sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets)),
  WAA_pattern = "constant",
  base_MatAA_values = array(rep(1 / (1 + exp(-2 * (1:sim_list$n_ages - 5))), each = sim_list$n_regions * sim_list$n_sexes),
                                  dim = c(sim_list$n_regions, sim_list$n_ages, sim_list$n_sexes)),
  MatAA_pattern = "constant"
)

ref <- 1
Movement <- array(0, dim = c(sim_list$n_regions, sim_list$n_regions, sim_list$n_yrs, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims)) # From, To
base <- matrix(c(0, 1), sim_list$n_regions, sim_list$n_regions, byrow = TRUE)

# Plug in movement process error
for(sim in 1:sim_list$n_sims) {
  for(a in 1:sim_list$n_ages) {
    for(s in 1:sim_list$n_sexes) {
      for(y in 1:sim_list$n_yrs) {
        for(r in 1:sim_list$n_regions) {
          tmp_move <- base[r,]
          tmp_move[-ref] <- tmp_move[-ref]
          Movement[r,,y,a,s,sim] <- exp(tmp_move) / sum(exp(tmp_move))
        } # end r loop
      } # end y loop
    } # end s loop
  } # end a loop
} # end sim loop

sim_list$Movement <- Movement

# Setup tagging stuff
sim_list <- Setup_Sim_Tagging(
  sim_list = sim_list,
  UseTagging = 0,
  n_tags = 5000,
  max_liberty = 30,
  tag_years = seq(1, sim_list$n_yrs, 1),
  t_tagging = 0.5,
  base_Tag_Reporting = c(0.2, 0.2),
  Tag_Reporting_pattern = "constant",
  Tag_Ind_Mort = 0,
  Tag_Shed = 0
)

# Setup observation processes
sim_list <- Setup_Sim_Observation_Proc(
  sim_list = sim_list,
  Comp_Structure = "spltR_jntS",
  Comp_Srv_Like = "Multinomial",
  Comp_Fish_Like = "Multinomial",
  ISS_FishAge_Pattern = 'constant',
  FishAgeTheta = 3,
  SrvAgeTheta = 2,
  Srv_Like_Pars = NA,
  base_ISS_FishAge = 500,
  base_ISS_SrvAge = 500,
  Tag_Like = "Poisson",
  Tag_Like_Pars = NA
)

# Define Assessment Model (Assessment Options) ----------------------------
# Next, we need to define assessment options (assessment model specifications)
# We will be defining skeleton data, parameter, and mapping lists using the maximum
# number of years in the simulation, which will get adjusted within the simulation function.

# Initialize model dimensions and data list
input_list <- Setup_Mod_Dim(years = 1:sim_list$n_yrs, # vector of years
                            ages = 1:sim_list$n_ages, # vector of ages
                            lens = 1, # number of lengths
                            n_regions = sim_list$n_regions, # number of regions
                            n_sexes = sim_list$n_sexes, # number of sexes
                            n_fish_fleets = sim_list$n_fish_fleets, # number of fishery fleet
                            n_srv_fleets = sim_list$n_srv_fleets, # number of survey fleets
                            verbose = F
)

# helper to define steepness starting values
inv_steepness <- function(s) qlogis((s - 0.2) / 0.8)

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
                                    WAA = array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes)),
                                    WAA_fish = array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                    MatAA = array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes)),
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
                                    ObsCatch = array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)),
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
                                          UseFishIdx = array(0, c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)),
                                          ObsFishAgeComps = array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                          UseFishAgeComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)),
                                          ObsFishLenComps = array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$lens), input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                          UseFishLenComps = array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)),
                                          ISS_FishAgeComps = array(500, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes, input_list$data$n_fish_fleets)),

                                          # Model options
                                          fish_idx_type = "none",
                                          FishAgeComps_LikeType = "Multinomial",
                                          FishLenComps_LikeType = "none",
                                          FishAgeComps_Type =  "spltRjntS_Year_1-terminal_Fleet_1",
                                          FishLenComps_Type =  "none_Year_1-terminal_Fleet_1"
                                          )


# Setup survey indices and compositions
input_list <- Setup_Mod_SrvIdx_and_Comps(input_list = input_list,
                                         ObsSrvIdx = array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets)),
                                         ObsSrvIdx_SE = array(0.05, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets)),
                                         UseSrvIdx =  array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)),
                                         ObsSrvAgeComps = array(NA,  dim = c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes, input_list$data$n_srv_fleets)),
                                         UseSrvAgeComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)),
                                         ObsSrvLenComps = array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$lens), input_list$data$n_sexes, input_list$data$n_srv_fleets)),
                                         UseSrvLenComps = array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets)),
                                         ISS_SrvAgeComps = array(500, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes, input_list$data$n_srv_fleets)),

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
                                      fish_fixed_sel_pars_spec = "est_shared_r",
                                      fish_q_spec = "fix"
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
                                tag_release_indicator = as.matrix(sim_list$tag_rel_indicator),
                                Tagged_Fish = array(NA, dim = c(dim(sim_list$Tag_Fish[1:nrow(as.matrix(sim_list$tag_rel_indicator)),,,1]), input_list$data$n_sexes)),
                                Obs_Tag_Recap = array(NA, dim = c(dim(sim_list$Obs_Tag_Recap[,1:nrow(as.matrix(sim_list$tag_rel_indicator)),,,,1]), input_list$data$n_sexes)),

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

# Now, define skeleton data, parameters, and mapping, which will get adjusted later
skeleton_data <- input_list$data
skeleton_parameters <- input_list$par
skeleton_mapping <- input_list$map


# Simulate Feedback Loop --------------------------------------------------

# Define HCR to use in simulation
HCR_function <- function(x, frp, brp, alpha = 0.05) {
  stock_status <- x / brp # define stock status
  # If stock status is > 1
  if(stock_status >= 1) f <- frp
  # If stock status is between brp and alpha
  if(stock_status > alpha && stock_status < 1) f <- frp * (stock_status - alpha) / (1 - alpha)
  # If stock status is less than alpha
  if(stock_status < alpha) f <- 0
  return(f)
}

# Setup simulation environment
sim_env <- Setup_sim_env(sim_list)

# Start Simulation
for (sim in 1:sim_env$n_sims) {
  for (y in 1:sim_env$n_yrs) {

    ### Run Annual Cycle --------------------------------------------------------
    run_annual_cycle(y, sim, sim_env)

    # Start Feedback Loop
    if(y >= sim_env$feedback_start_yr) {

      # Get feedback data (aligning year index)
      feedback <- Get_Feedback_Data(sim_env = sim_env, sim_list = sim_list,
                                    y = y, sim = sim,
                                    skeleton_data = skeleton_data,
                                    skeleton_parameters = skeleton_parameters,
                                    skeleton_mapping = skeleton_mapping
                                    )

      # Extract out data, parameters, and mapping
      data <- feedback$retro_data
      parameters <- feedback$retro_parameters
      mapping <- feedback$retro_mapping

      ### Run Assessment ----------------------------------------------------------
      obj <- fit_model(data,
                       parameters,
                       mapping,
                       random = NULL,
                       newton_loops = 3,
                       silent = F
                       )

      # Get necessary data inputs for later functions
      n_proj_yrs <- 2
      tmp_terminal_NAA <- array(obj$rep$NAA[,y,,], dim = c(sim_env$n_regions, sim_env$n_ages, sim_env$n_sexes)) # terminal numbers at age
      tmp_WAA <- array(rep(data$WAA[,y,,], each = n_proj_yrs), dim = c(sim_env$n_regions, n_proj_yrs, sim_env$n_ages, sim_env$n_sexes)) # weight at age
      tmp_WAA_fish <- array(rep(data$WAA[,y,,], each = n_proj_yrs), dim = c(sim_env$n_regions, n_proj_yrs, sim_env$n_ages, sim_env$n_sexes, sim_env$n_fish_fleets)) # weight at age fishery
      tmp_MatAA <- array(rep(data$MatAA[,y,,], each = n_proj_yrs), dim = c(sim_env$n_regions, n_proj_yrs, sim_env$n_ages, sim_env$n_sexes)) # maturity at age
      tmp_fish_sel <- array(rep(obj$rep$fish_sel[,y,,,], each = n_proj_yrs), dim = c(sim_env$n_regions, n_proj_yrs, sim_env$n_ages, sim_env$n_sexes, sim_env$n_fish_fleets)) # selectivity
      tmp_Movement <- array(rep(obj$rep$Movement[,,y,,,drop=FALSE], n_proj_yrs), dim = c(sim_env$n_regions, sim_env$n_regions, n_proj_yrs, sim_env$n_ages, sim_env$n_sexes))
      tmp_terminal_F <- array(obj$rep$Fmort[,y,], dim = c(sim_env$n_regions, sim_env$n_fish_fleets)) # terminal fishing mortality
      tmp_natmort <- array(obj$rep$natmort[,y,,], dim = c(sim_env$n_regions, n_proj_yrs, sim_env$n_ages, sim_env$n_sexes)) # natural mortality
      tmp_recruitment <- array(obj$rep$Rec[,1:(y - 1)], dim = c(sim_env$n_regions, length(1:(y - 1)))) # recruitment to use for projections
      tmp_sexratio <- array(1, dim = c(sim_env$n_regions, n_proj_yrs, sim_env$n_sexes)) # recruitment sex ratio

      ### Run reference point module ----------------------------------------------
      reference_points <- Get_Reference_Points(data = data,
                                               rep = obj$rep,
                                               SPR_x = 0.4,
                                               t_spwn = 0,
                                               sex_ratio_f = 0.5,
                                               calc_rec_st_yr = 1,
                                               rec_age = 1,
                                               type = 'multi_region',
                                               what = 'independent_SPR'
                                               )

      # extract fishery and biological reference points
      tmp_f_ref_pt <- array(reference_points$f_ref_pt, dim = c(sim_env$n_regions, n_proj_yrs)) # fishery reference points
      tmp_b_ref_pt <- array(reference_points$b_ref_pt, dim = c(sim_env$n_regions, n_proj_yrs)) # biological reference points

     ### Run Projection Module ---------------------------------------------------

     # Do projection to get TAC
     proj <- Do_Population_Projection(n_proj_yrs = 2,
                                      n_regions = sim_env$n_regions,
                                      n_ages = sim_env$n_ages,
                                      n_sexes = sim_env$n_sexes,
                                      sexratio = tmp_sexratio,
                                      n_fish_fleets = sim_env$n_fish_fleets,
                                      do_recruits_move = 0,
                                      recruitment = tmp_recruitment,
                                      terminal_NAA = tmp_terminal_NAA,
                                      terminal_F = tmp_terminal_F,
                                      natmort = tmp_natmort,
                                      WAA = tmp_WAA,
                                      WAA_fish = tmp_WAA_fish,
                                      MatAA = tmp_MatAA,
                                      fish_sel = tmp_fish_sel,
                                      Movement = tmp_Movement,
                                      f_ref_pt = tmp_f_ref_pt,
                                      b_ref_pt = tmp_b_ref_pt,
                                      HCR_function = HCR_function,
                                      recruitment_opt = "mean_rec",
                                      fmort_opt = "HCR",
                                      t_spawn = 0
                                      )

      # Get TAC
      tmp_TAC <- proj$proj_Catch[,2,,drop = FALSE] # get TAC from projected year

      ### TAC to Fishing Mortality ------------------------------------------------
      # Only run till last year
      if(y < sim_env$n_yrs) {
        for(r in 1:sim_env$n_regions) {
          for(f in 1:sim_env$n_fish_fleets) {

            # Go from TAC to Fishing mortality (using true values from simulation)
            tmp_F <- bisection_F(f_guess = 0.05,
                                 catch = tmp_TAC[r,1,f],
                                 NAA = sim_env$NAA[r,y+1,,,sim],
                                 WAA = sim_env$WAA[r,y+1,,,sim],
                                 natmort = sim_env$M[r,y+1,,,sim],
                                 fish_sel = sim_env$fish_sel[r,y+1,,,f,sim]
                                 )

            sim_env$Fmort[r,y+1,f,sim] <- tmp_F

          } # end r loop
        } # end f loop
      } # end if

      plot(sim_env$SSB[1,1:y,sim], type = 'l')

    } # end if for feedback loop

  } # end y loop
} # end sim loop

png(here("vignettes", "figures", "h_ssb_closed_loop.png"), width = 800, height = 400)
reshape2::melt(sim_env$SSB) %>%
  dplyr::rename(Region = Var1, Year = Var2, Sim = Var3) %>%
  ggplot(aes(x = Year, y = value, group = Sim)) +
  geom_line() +
  geom_vline(xintercept = sim_env$feedback_start_yr, lty = 2) +
  facet_wrap(~Region) +
  theme_bw(base_size = 15) +
  labs(y = 'SSB') +
  ylim(0, NA)
dev.off()

png(here("vignettes", "figures", "h_catch_closed_loop.png"), width = 800, height = 400)
reshape2::melt(sim_env$Obs_Catch) %>%
  dplyr::rename(Region = Var1, Year = Var2, Fleet = Var3, Sim = Var4) %>%
  dplyr::filter(Year != 36) %>%
  ggplot(aes(x = Year, y = value, group = Sim)) +
  geom_line() +
  geom_vline(xintercept = sim_env$feedback_start_yr, lty = 2) +
  facet_wrap(~Region) +
  theme_bw(base_size = 15) +
  labs(y = 'Catch') +
  ylim(0, NA)
dev.off()


