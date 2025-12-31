  # Purpose: To set up simulations for a spatially-structured population
  # Author: Matthew LH. Cheng (UAF-CFOS)

  # load in libraries
  library(here)
  library(SPoRC)
  library(tidyverse)

  # Setup Simulation Parameters --------------------------------------------

  # Set up model dimensions
  sim_list <- Setup_Sim_Dim(n_sims = 500,
                            n_yrs = 10,
                            n_regions = 2,
                            n_ages = 8,
                            n_sexes = 1,
                            n_fish_fleets = 1,
                            n_srv_fleets = 1,
                            run_feedback = F,
                            feedback_start_yr = NULL
                            )

  # set up containers
  sim_list <- Setup_Sim_Containers(sim_list)

  # Setup fishing mortality
  sim_list <- Setup_Sim_FishMort(sim_list = sim_list,
                                 sigmaC = 1e-3,
                                 init_F_vals = matrix(0, nrow = sim_list$n_regions, ncol = sim_list$n_fish_fleets),
                                 Fmort_pattern = matrix(c('one-way', "one-way"), nrow = sim_list$n_regions, ncol = sim_list$n_fish_fleets),
                                 Fmort_start = matrix(c(0.01, 0.01), nrow = sim_list$n_regions, ncol = sim_list$n_fish_fleets),
                                 Fmort_fct = matrix(c(15, 15), nrow = sim_list$n_regions, ncol = sim_list$n_fish_fleets),
                                 proc_error = TRUE,
                                 proc_error_sd = 0.05
                                 )

  # Setup fishery selectivity
  sim_list <- Setup_Sim_FishSel(sel_model = matrix(c('logistic', "logistic"), nrow = sim_list$n_regions, ncol = sim_list$n_fish_fleets),
                                # a50, k for logistic shared across regions
                                fixed_fish_sel_pars = array(c(3,3,1,1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
                                sim_list = sim_list
                                )

  # Setup survey catchability and selectivity
  sim_list <- Setup_Sim_Survey(sim_list = sim_list,
                               sigmaSrvIdx = array(0.2, dim = c(sim_list$n_regions, sim_list$n_srv_fleets)), # survey observation error
                               base_srv_q = array(1, dim = c(sim_list$n_regions, sim_list$n_srv_fleets)), # base survey catchability value
                               srv_q_pattern = matrix(c('constant', "constant"), nrow = sim_list$n_regions, ncol = sim_list$n_srv_fleets), # catchability pattern
                               sel_model = matrix(c('logistic', "logistic"), nrow = sim_list$n_regions, ncol = sim_list$n_srv_fleets), # selectivity model
                               # a50, k, for logistic shared across regions
                               fixed_srv_sel_pars = array(c(3,3,1,1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets, 2))
                             )

  # Setup recruitment stuff
  sim_list <- Setup_Sim_Rec(sim_list = sim_list,
                            do_recruits_move = "dont_move", # == 0, recruits don't move , == 1 recruits move
                            base_rec_sexratio = 1, # single sex
                            rec_sexratio_vary = "constant",
                            base_r0 = c(50, 50),
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
  sim_list <- Setup_Sim_Biologicals(sim_list = sim_list,
                                    base_M_value = array(0.5, dim = c(sim_list$n_regions, sim_list$n_ages, sim_list$n_sims)),
                                    M_pattern = "constant",
                                    base_WAA_values = array(rep(5 / (1 + exp(-1 * (1:sim_list$n_ages - 5))), each = sim_list$n_regions * sim_list$n_sexes),
                                                            dim = c(sim_list$n_regions, sim_list$n_ages, sim_list$n_sexes)),
                                    WAA_pattern = "constant",
                                    base_Maturity_AA_values = array(rep(1 / (1 + exp(-2 * (1:sim_list$n_ages - 5))), each = sim_list$n_regions * sim_list$n_sexes),
                                         dim = c(sim_list$n_regions, sim_list$n_ages, sim_list$n_sexes)),
                                    Maturity_AA_pattern = "constant"
                                  )

  # Setup tagging stuff
  sim_list <- Setup_Sim_Tagging(sim_list = sim_list,
                                n_tags = 5e3,
                                max_liberty = 30,
                                tag_years = seq(1, sim_list$n_yrs, 1),
                                t_tagging = 0.5,
                                base_Tag_Reporting = c(0.2, 0.2),
                                Tag_Reporting_pattern = "constant",
                                Tag_Ind_Mort = 0,
                                Tag_Shed = 0
                              )

  # Setup observation processes
  sim_list <- Setup_Sim_Observation_Proc(sim_list = sim_list,
                                         Comp_Structure = "spltR_jntS",
                                         Comp_Srv_Like = "Multinomial",
                                         Comp_Fish_Like = "Multinomial",
                                         ISS_FishAge_Pattern = 'constant',
                                         FishAgeTheta = 3,
                                         SrvAgeTheta = 2,
                                         Srv_Like_Pars = NA,
                                         base_ISS_FishAge = 100,
                                         base_ISS_SrvAge = 100,
                                         Tag_Like = "Poisson",
                                         Tag_Like_Pars = NA
                                         )

  # IID Movement Matrix across years and ages
  ref <- 1
  movement_matrix <- array(0, dim = c(sim_list$n_regions, sim_list$n_regions, sim_list$n_yrs, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims)) # From, To
  base <- matrix(0, sim_list$n_regions, sim_list$n_regions)

  # Plug in movement process error
  for(sim in 1:sim_list$n_sims) {
    for(a in 1:sim_list$n_ages) {
      for(s in 1:sim_list$n_sexes) {
        for(y in 1:sim_list$n_yrs) {
          for(r in 1:sim_list$n_regions) {
            if(a > 1) pe_err <- rnorm(length(tmp_move[-ref]), 0, 0) # logit proces error
            else pe_err <- 0
            tmp_move <- base[r,]
            tmp_move[-ref] <- tmp_move[-ref] + pe_err
            movement_matrix[r,,y,a,s,sim] <- exp(tmp_move) / sum(exp(tmp_move))
          }
        } # end y loop
      } # end s loop
    } # end a loop
  } # end sim loop

  sim_list$movement_matrix <- movement_matrix

  # Run Simulation ----------------------------------------------------------
  Simulate_Pop_Static(sim_list = sim_list, output_path = here("sim_out.RDS"))

































  # Scracth Movement Simulation Code
  # # Within a given region and across years, all ages have the same deviation (if pcorr_age => 1)
  # age_cor <- -0.95 # Correlation between adjacent ages
  # year_cor <-  0.95 # Correlation between adjacent years
  # region_cor <- 0
  # Sdev_at = array(NA, dim=c(n_regions, n_yrs, n_ages))
  # Cov_aa = age_cor ^ as.matrix(dist(1:n_ages,diag=TRUE,upper=TRUE)) / (1 - age_cor^2)
  # Cov_rr = region_cor ^ as.matrix(dist(1:n_regions, diag = TRUE, upper = TRUE))  / (1 - region_cor^2)
  # Cov_yy = year_cor ^ as.matrix(dist(1:n_yrs, diag = TRUE, upper = TRUE)) / (1 - year_cor^2)
  # Cov_full = 0 * kronecker(Cov_rr, kronecker(Cov_yy, Cov_aa))
  #
  # movement_matrix <- array(0, dim = c(n_regions, n_regions, n_yrs, n_ages, n_sexes, n_sims))  # From, To
  # base <- matrix(rev(c(0.8, 2, 0, 0)), n_regions, n_regions)
  #
  # for(sim in 1:n_sims) {
  #
  #   set.seed(123)
  #   pe_array <- array(0, dim = c(n_regions, n_regions - 1, n_yrs, n_ages))
  #
  #   for(i in 1:(n_regions - 1)) {
  #     devs = mvtnorm::rmvnorm(n = 1, mean = rep(0, n_yrs * n_regions * n_ages), sigma = Cov_full)[1, ]
  #     Sdev_at <- array(devs, dim = c(n_ages, n_yrs, n_regions))
  #     for(r in 1:n_regions) for(y in 1:n_yrs) for(a in 1:n_ages) pe_array[r,i,y,a] <- Sdev_at[a,y,r]
  #   }
  #
  #   for(r in 1:n_regions) {
  #     for(s in 1:n_sexes) {
  #       for(y in 1:n_yrs) {
  #         for(a in 1:n_ages) {
  #           if(a > 1) pe_err <- pe_array[r,,y, a] # logit proces error
  #           else pe_err <- 0
  #           tmp_move <- base[r, ]
  #           tmp_move[-ref] <- tmp_move[-ref] + pe_err
  #           movement_matrix[r, , y, a, s, sim] <- exp(tmp_move) / sum(exp(tmp_move))
  #         }
  #       }
  #     }
  #   }
  # }

  # reshape2::melt(Sdev_at) %>%
  #   rename(Region = Var3, Year = Var2, Age = Var1) %>%
  #   ggplot(aes(x = Year, y = value, color = Age, group = Age)) +
  #   geom_line(lwd = 1.3) +
  #   facet_wrap(~Region, scales = 'free') +
  #   scale_color_viridis_c()
  #
  # reshape2::melt(Sdev_at) %>%
  #   rename(Region = Var3, Year = Var2, Age = Var1) %>%
  #   ggplot(aes(x = Year, y = value, color = factor(Region))) +
  #   geom_line() +
  #   facet_grid(Region~Age, scales = 'free')
