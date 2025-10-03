#' Constructs simulation objects in a new simulation environment for use in simulation functions
#'
#' @param sim_list Simulation list objects
#'
#' @returns A new simulation environment with objects from sim_list
#' @export Setup_sim_env
#'
#' @examples
#' \dontrun{
#' sim_env <- Setup_sim_env(sim_list)
#' }
#' @family Simulation Setup
Setup_sim_env <- function(sim_list) {

  sim_env <- new.env(parent = parent.frame()) # define new environment for simulation

  # Get SPoRC functions in simulation environment
  sim_env$Get_Det_Recruitment <- Get_Det_Recruitment
  sim_env$Get_Tagging_Mortality <- Get_Tagging_Mortality
  sim_env$rdirM <- rdirM
  sim_env$predict_sim_fish_iss_fmort <- predict_sim_fish_iss_fmort
  sim_env$rho_trans <- rho_trans

  # output into simulation environment
  list2env(sim_list, envir = sim_env)

  return(sim_env)
}

#' Simulate Age or Length Compositions
#'
#' Generates Observed fish compositions by age or length for a given region, year, fleet, and simulation iteration.
#' Supports multinomial, Dirichlet-multinomial, and logistic-normal likelihoods, with optional ageing error applied.
#'
#' @param r Integer. Region index.
#' @param y Integer. Year index.
#' @param f Integer. Fleet index.
#' @param sim Integer. Simulation iteration index.
#' @param Exp Array. Expected compositions (age or length) with dimensions [region, year, category, sex, fleet, sim].
#' @param ISS Array. Sample size (integer) for the Observed compositions with dimensions [region, year, sex, fleet, sim].
#' @param AgeingError Array. Ageing error matrix for each year, dimensions [year, category, category, sim].
#' @param comp_like Integer vector. Composition likelihood type per fleet: 0 = multinomial, 1 = Dirichlet-multinomial, 2-4 = logistic-normal.
#' @param ln_theta Array. Log-variance parameter for compositions per region, sex, and fleet, dimensions [region, sex, fleet].
#' @param corr_pars Array. Correlation parameters for logistic-normal likelihood, dimensions [region, sex, fleet, ?].
#' @param ln_theta_agg Numeric vector. Log-variance parameter for aggregated compositions per fleet.
#' @param corr_pars_agg Numeric vector. Correlation parameters for aggregated logistic-normal likelihood per fleet.
#' @param comp_type Integer array. Composition type: 0 = aggregated across regions, 1 = split by sex, 2 = joint across sexes, dimensions [year, fleet].
#' @param n_sexes Integer. Number of sexes.
#' @param n_regions Integer. Number of regions.
#' @param n_cat Integer. Number of categories (ages or lengths).
#' @param Obs Array. Observed compositions array to fill, same dimensions as `Exp`.
#' @param age_or_len Integer. Flag to indicate if ageing error should be applied: 0 = apply ageing error (for ages), 1 = do not apply (for lengths).
#'
#' @return Array of Observed compositions with the same dimensions as `Obs`, updated with simulated Observations.
#'
#' @details
#' The function handles three cases based on `comp_type`:
#' 1. Split by sex (comp_type = 1): compositions are simulated separately for each sex.
#' 2. Joint compositions across sexes (comp_type = 2): compositions simulated jointly and multiplied by a kronecker matrix for logistic-normal or Dirichlet-multinomial likelihoods.
#' 3. Aggregated across regions (comp_type = 0): only applied in the last region and averages across regions and sexes.
#'
#' The function normalizes expected compositions, applies the selected likelihood (`comp_like`), and multiplies by `AgeingError` when applicable.
#'
#' @keywords internal
simulate_comps <- function(r, y, f, sim,
                           Exp, ISS,
                           AgeingError,
                           comp_like,
                           ln_theta,
                           corr_pars,
                           ln_theta_agg,
                           corr_pars_agg,
                           comp_type,
                           n_sexes,
                           n_regions,
                           n_cat,
                           Obs,
                           age_or_len = 0) {

  if(comp_type[y,f] == 999 || comp_like[f] == 999) return(Obs)

  # helper functions
  get_expected <- function(prob_vec) prob_vec / sum(prob_vec)
  apply_error <- function(mat, age_or_len, AgeingError) {
    if(age_or_len == 0) return(mat %*% AgeingError)
    if(age_or_len == 1) return(mat)
  }

  if(age_or_len == 0) {
    if(comp_type[y,f] %in% c(0,1)) age_error_mat <- AgeingError[y,,,sim] # aggregated or split
    if(comp_type[y,f] == 2) age_error_mat <- kronecker(diag(n_sexes), AgeingError[y,,,sim]) # joint
  } else if(age_or_len == 1) age_error_mat <- NULL # length compositions

  # Split by sex
  if(comp_type[y,f] == 1) {
    for(s in 1:n_sexes) {
      tmp_prob <- Exp[r,y,,s,f,sim] # extract compositions

      # multinomial
      if(comp_like[f] == 0) {
        Obs[r,y,,s,f,sim] <- array(
          apply_error(as.vector(
            stats::rmultinom(n = 1, ISS[r,y,s,f,sim], get_expected(tmp_prob))), age_or_len, age_error_mat),
          dim = dim(Obs[r,y,,s,f,sim, drop = FALSE])
        )

        # dirichlet-multinomial
      } else if(comp_like[f] == 1) {
        Obs[r,y,,s,f,sim] <- array(
          apply_error(as.vector(
            rdirM(
              n = 1,
              N = ISS[r,y,s,f,sim],
              alpha = (exp(ln_theta[r,s,f]) * ISS[r,y,s,f,sim]) * get_expected(tmp_prob)
            )
          ), age_or_len, age_error_mat),
          dim = dim(Obs[r,y,,s,f,sim, drop = FALSE])
        )

        # logistic normal
      } else if(comp_like[f] %in% 2:4) {
        Obs[r,y,,s,f,sim] <- array(
          apply_error(as.vector(
            rlogistnormal(
              exp = get_expected(tmp_prob),
              pars = c(exp(ln_theta[r,s,f]), rho_trans(corr_pars[r,s,f,])),
              comp_like = comp_like[f]
            )
          ), age_or_len, age_error_mat),
          dim = dim(Obs[r,y,,s,f,sim, drop = FALSE])
        )
      }

    } # end s loop
  } # end split by sex

  # Joint compositions
  if(comp_type[y,f] == 2) {
    tmp_prob <- Exp[r,y,,,f,sim] # extract compositions

    # multinomial
    if(comp_like[f] == 0) {
      Obs[r,y,,,f,sim] <- array(
        apply_error(as.vector(stats::rmultinom(1, ISS[r,y,1,f,sim], get_expected(tmp_prob))),
                    age_or_len, age_error_mat),
        dim = dim(Obs[r,y,,,f,sim, drop = FALSE])
      )

      # dirichlet-multinomial
    } else if(comp_like[f] == 1) {
      Obs[r,y,,,f,sim] <- array(
        apply_error(as.vector(
          rdirM(
            n = 1,
            N = ISS[r,y,1,f,sim],
            alpha = (exp(ln_theta[r,1,f]) * ISS[r,y,1,f,sim]) * get_expected(tmp_prob)
          )
        ), age_or_len, age_error_mat),
        dim = dim(Obs[r,y,,,f,sim, drop = FALSE])
      )

      # logistic normal
    } else if(comp_like[f] %in% 2:4) {
      Obs[r,y,,,f,sim] <- array(
        apply_error(as.vector(
          rlogistnormal(
            exp = get_expected(tmp_prob),
            pars = c(exp(ln_theta[r,1,f]), rho_trans(corr_pars[r,1,f,])),
            comp_like = comp_like[f]
          )
        ), age_or_len, age_error_mat),
        dim = dim(Obs[r,y,,,f,sim, drop = FALSE])
      )
    }

  } # end joint compositions

  # Aggregated comps across regions
  if(r == n_regions && comp_type[y,f] == 0) {

    # extract compositions
    tmp_prob <- Exp[,y,,,f,sim]
    tmp_prob <- matrix(rowSums(matrix(tmp_prob, nrow = n_cat)) / (n_sexes * n_regions), nrow = 1)
    tmp_prob <- tmp_prob / sum(tmp_prob)

    # multinomial
    if(comp_like[f] == 0) {
      Obs[1,y,,1,f,sim] <- array(
        apply_error(as.vector(stats::rmultinom(1, ISS[1,y,1,f,sim], get_expected(tmp_prob))), age_or_len, age_error_mat),
        dim = dim(Obs[1,y,,1,f,sim, drop = FALSE])
      )

      # dirichlet-multinomial
    } else if(comp_like[f] == 1) {
      Obs[1,y,,1,f,sim] <- array(
        apply_error(as.vector(
          rdirM(
            n = 1,
            N = ISS[1,y,1,f,sim],
            alpha = (exp(ln_theta_agg[f]) * ISS[1,y,1,f,sim]) * get_expected(tmp_prob)
          )
        ), age_or_len, age_error_mat),
        dim = dim(Obs[1,y,,1,f,sim, drop = FALSE])
      )

      # logistic normal
    } else if(comp_like[f] %in% 2:4) {
      Obs[1,y,,1,f,sim] <- array(
        apply_error(as.vector(
          rlogistnormal(
            exp = get_expected(tmp_prob),
            pars = c(exp(ln_theta_agg[f]), rho_trans(corr_pars_agg[f])),
            comp_like = comp_like[f]
          )
        )),
        dim = dim(Obs[1,y,,1,f,sim, drop = FALSE])
      )
    }
  }

  return(Obs)
}

#' Run Annual Cycle in Simulation Environment
#'
#' @param y Year index
#' @param sim Simulation index
#' @param sim_env Simulation environment will all the necessary elements to run the annual cycle
#' @export run_annual_cycle
#' @importFrom stats rnorm rmultinom
#' @family Simulation Setup
run_annual_cycle <- function(y,
                             sim,
                             sim_env) {

  # Assign y and sim into simulation environment
  sim_env$y <- y
  sim_env$sim <- sim

  with(sim_env, {
    # Initialize Age Structure ------------------------------------------------
    if(y == 1) {

     if(init_age_strc == 0) {

       # Set up initial equilibrium age structure
       for(r in 1:n_regions) {
         for(s in 1:n_sexes) {
           tmp_cumsum_Z <- cumsum(natmort[r,1,1:(n_ages-1),s,sim] + init_F * fish_sel[r,1,1:(n_ages-1),s,1,sim]) # cumulative sum of total mortality
           Init_NAA[r,,s,sim] <- c(R0[r,1,sim] * sexratio[r,1,s,sim], R0[r,1,sim] * sexratio[r,1,s,sim] * exp(-tmp_cumsum_Z)) # exponential mortality model
         } # end s loop
       } # end r loop

       # Apply annual cycle and iterate to equilibrium
       for(i in 1:init_iter) {
         for(s in 1:n_sexes) {
           Init_NAA_next_year[,1,s,sim] <- R0[,1,sim] * sexratio[,1,s,sim] # recruitment
           # movement
           if(do_recruits_move == 0) for(a in 2:n_ages) Init_NAA[,a,s,sim] <- t(Init_NAA[,a,s,sim]) %*% Movement[,,1,a,s,sim] # recruits dont move
           if(do_recruits_move == 1) for(a in 1:n_ages) Init_NAA[,a,s,sim] <- t(Init_NAA[,a,s,sim]) %*% Movement[,,1,a,s,sim] # recruits move
           # ageing and mortality
           Init_NAA_next_year[,2:n_ages,s,sim] <- Init_NAA[,1:(n_ages-1),s,sim] * exp(-(natmort[,1,1:(n_ages-1),s,sim] + (init_F * fish_sel[,1,1:(n_ages-1),s,1,sim])))
           # accumulate plus group
           Init_NAA_next_year[,n_ages,s,sim] <- (Init_NAA_next_year[,n_ages,s,sim] * exp(-(natmort[,1,n_ages,s,sim] + (init_F * fish_sel[,1,n_ages,s,1,sim])))) +
                                                (Init_NAA[,n_ages,s,sim] * exp(-(natmort[,1,n_ages,s,sim] + (init_F * fish_sel[,1,n_ages,s,1,sim]))))
           Init_NAA <- Init_NAA_next_year # iterate to next cycle
         } # end s loop
       } # end i loop

       # Set up initial age deviations
       tmp_ln_init_devs <- NULL
       for(r in 1:n_regions) {
         if(exists("ln_InitDevs_input")) { # if exists in environment, then use input
           tmp_ln_init_devs <- ln_InitDevs_input[r,,sim]
         } else { # simulate new initial age devs otherwise
           if(init_dd == 0 && is.null(tmp_ln_init_devs)) tmp_ln_init_devs <- stats::rnorm(n_ages-1, 0, exp(ln_sigmaR[1])) # global initial age deviations
           if(init_dd == 1) tmp_ln_init_devs <- stats::rnorm(n_ages-1, 0, exp(ln_sigmaR[1])) # local initial age deviations
         }
         ln_InitDevs[r,,sim] <- tmp_ln_init_devs # save ln_init_devs
         for(s in 1:n_sexes) NAA[r,1,2:n_ages,s,sim] <- Init_NAA[r,2:n_ages,s,sim] * exp(ln_InitDevs[r,,sim]) # apply deviations
       }

     } else if(init_age_strc == 1) {

       # Set up initial age deviations
       tmp_ln_init_devs <- NULL
       init_age_idx <- 1:(n_ages - 2) # Get initial age indexing
       for(r in 1:n_regions) {

         if(exists("ln_InitDevs_input")) { # if exists in environment, then use input
           tmp_ln_init_devs <- ln_InitDevs_input[r,,sim]
         } else { # simulate new initial age devs otherwise
           if(init_dd == 0 && is.null(tmp_ln_init_devs)) tmp_ln_init_devs <- stats::rnorm(n_ages-1, 0, exp(ln_sigmaR[1])) # global initial age deviations
           if(init_dd == 1) tmp_ln_init_devs <- stats::rnorm(n_ages-1, 0, exp(ln_sigmaR[1])) # local initial age deviations
         }

         ln_InitDevs[r,,sim] <- tmp_ln_init_devs # save ln_init_devs

         for(s in 1:n_sexes) {
           NAA[r,1,init_age_idx + 1,s,sim] <- R0[r,1,sim] * exp(ln_InitDevs[r,init_age_idx,sim] - # not plus group
                                              (init_age_idx * (natmort[r,1, init_age_idx + 1, s,sim] +
                                              (init_F * fish_sel[r,1, init_age_idx + 1, s, 1,sim])))) * sexratio[r,1,s,sim]

           NAA[r,1,n_ages,s,sim] <- R0[r,1,sim] * exp(ln_InitDevs[r,n_ages - 1,sim] - # plus group calculations (geometric series)
                                    ((n_ages - 1) * (natmort[r,1, n_ages, s,sim] + (init_F * fish_sel[r,1, n_ages, s, 1,sim])))) /
                                    (1 - exp(-(natmort[r,1, n_ages, s,sim] + (init_F *
                                    fish_sel[r,1, n_ages, s, 1,sim])))) * sexratio[r,1,s,sim]

         } # end s loop
       } # end r loop

     }
      NAA0[,1,,,sim] <- NAA[,1,,,sim] # Initialize Unfished NAA
    } # end initializing age structure

    # Run Annual Cycle --------------------------------------------------------
    ### Recruitment (Year 1 only) ----------------------------------------------
    if(y == 1) {
      tmp_ln_RecDevs <- NULL # Initialize container vector for global recruitment deviations (remains NULL within a given year)
      for(r in 1:n_regions) {

        # if recruitment input exists, use input, if not, simulate new deviates
        if(exists("Rec_input")) for(s in 1:n_sexes) NAA[r,1,1,s,sim] <- Rec_input[r,y,sim] * sexratio[r,y,s,sim] else {

          if(rec_dd == 0 && is.null(tmp_ln_RecDevs)) tmp_ln_RecDevs <- stats::rnorm(1, 0, exp(ln_sigmaR[2])) # Global Recruitment Deviations
          if(rec_dd == 1) tmp_ln_RecDevs <- ln_RecDevs[r,y,sim] <- stats::rnorm(1, 0, exp(ln_sigmaR[2])) # Local Recruitment Deviations
          ln_RecDevs[r,y,sim] <- tmp_ln_RecDevs # Input recruitment deviations into vector

          # Get deterministic recruitment
          tmp_det_rec <- Get_Det_Recruitment(recruitment_model = recruitment_opt,
                                             recruitment_dd = rec_dd,
                                             y = y,
                                             rec_lag = rec_lag,
                                             R0 = sum(R0[,y,sim]), # sum to get global R0
                                             Rec_Prop = R0[,y,sim] / sum(R0[,y,sim]), # get R0 proportion
                                             h = h[,y,sim],
                                             n_regions = n_regions,
                                             n_ages = n_ages,
                                             WAA = array(WAA[,y,,1,sim], dim = c(n_regions, n_ages)),
                                             MatAA = array(MatAA[,y,,1,sim], dim = c(n_regions, n_ages)),
                                             natmort = array(natmort[,y,,1,sim], dim = c(n_regions, n_ages)),
                                             SSB_vals = array(SSB[,,sim], dim = c(n_regions, n_yrs))
          )

          # input deviates
          for(s in 1:n_sexes) NAA[r,1,1,s,sim] <- tmp_det_rec[r] * exp(ln_RecDevs[r,y,sim] - exp(ln_sigmaR[2])^2/2) * sexratio[r,y,s,sim]
        }

        Rec[r,1,sim] <- sum(NAA[r,1,1,,sim]) # Save recruitment estimates
        NAA0[r,y,1,,sim] = NAA[r,y,1,,sim] # populate unfished NAA

      } # end r loop
    } # end if first year recruitment

    #### Movement ----------------------------------------------------------------
    if(n_regions > 1) {
      if(do_recruits_move == 0) { # Recruits don't move
        for(a in 2:n_ages) { # apply movement after ageing processes - start movement at age 2
          for(s in 1:n_sexes) {
            NAA[,y,a,s,sim] <- t(NAA[,y,a,s,sim]) %*% Movement[,,y,a,s,sim] # Fished
            NAA0[,y,a,s,sim] <- t(NAA0[,y,a,s,sim]) %*% Movement[,,y,a,s,sim] # Unfished
          } # end s loop
        } # end a loop
      } # end if recruits don't move

      if(do_recruits_move == 1) { # Recruits move here
        for(a in 1:n_ages) {
          for(s in 1:n_sexes) {
            NAA[,y,a,s,sim] <- t(NAA[,y,a,s,sim]) %*% Movement[,,y,a,s,sim] # Fished
            NAA0[,y,a,s,sim] <- t(NAA0[,y,a,s,sim]) %*% Movement[,,y,a,s,sim] # Unfished
          } # end s loop
        } # end a loop
      } # end if
    } # only compute if spatial

    ### Mortality and Ageing ----------------------------------------------------
    for(r in 1:n_regions) for(a in 1:n_ages) for(s in 1:n_sexes) ZAA[r,y,a,s,sim] <- natmort[r,y,a,s,sim] + sum(Fmort[r,y,,sim] * fish_sel[r,y,a,s,,sim]) # compute total mortality
    NAA[,y+1,2:n_ages,,sim] = NAA[,y,1:(n_ages-1),,sim] * exp(-ZAA[,y,1:(n_ages-1),,sim]) # Exponential mortality for individuals not in plus group (fished)
    NAA[,y+1,n_ages,,sim] = NAA[,y+1,n_ages,,sim] + NAA[,y,n_ages,,sim] * exp(-ZAA[,y,n_ages,,sim]) # Acuumulate plus group (fished)
    NAA0[,y+1,2:n_ages,,sim] <- NAA0[,y,1:(n_ages-1),,sim] * exp(-natmort[,y,1:(n_ages-1),,sim]) # Exponential mortality for individuals not in plus group (unfished)
    NAA0[,y+1,n_ages,,sim] <- NAA0[,y+1,n_ages,,sim] + NAA0[,y,n_ages,,sim] * exp(-natmort[,y,n_ages,,sim]) # Acuumulate plus group (unfished)

    ### Compute Biomass Quantities ----------------------------------------------
    Total_Biom[, y, sim] <- apply(NAA[, y, , , sim,drop = FALSE] * WAA[, y, , , sim,drop = FALSE], 1, sum) # Total Biomass
    SSB[, y, sim] <- apply(NAA[, y, , 1, sim,drop = FALSE] * WAA[, y, , 1, sim,drop = FALSE] * MatAA[, y, , 1, sim,drop = FALSE] * exp(-ZAA[, y, , 1, sim,drop = FALSE] * t_spawn), 1, sum) # Spawning Stock Biomass
    Dynamic_SSB0[,y,sim] <- apply(NAA0[, y, , 1, sim,drop = FALSE] * WAA[, y, , 1, sim,drop = FALSE] * MatAA[, y, , 1, sim,drop = FALSE] * exp(-natmort[, y, , 1, sim,drop = FALSE] * t_spawn), 1, sum) # Dynamic B0
    if(n_sexes == 1) { # If single sex model, multiply SSB calculations by 0.5
      SSB[,y,sim] <- SSB[,y,sim] * 0.5
      Dynamic_SSB0[,y,sim] <- Dynamic_SSB0[,y,sim] * 0.5
    }

    for(r in 1:n_regions) {

      ### Generate Fishery Observations ----------------------------------------------------------
      for(f in 1:n_fish_fleets) {

        # Baranov's catch equation
        CAA[r,y,,,f,sim] <- (Fmort[r,y,f,sim] * fish_sel[r,y,,,f,sim]) / ZAA[r,y,,,sim] *  NAA[r,y,,,sim] * (1 - exp(-ZAA[r,y,,,sim]))
        if(exists("SizeAgeTrans") && !is.null(SizeAgeTrans)) for(s in 1:n_sexes) CAL[r,y,,s,f,sim] <- SizeAgeTrans[r,y,,,s,sim] %*% CAA[r,y,,s,f,sim] # Catch at length

        #### Catch -------------------------------------------------------------------
        TrueCatch[r,y,f,sim] <- sum(CAA[r,y,,,f,sim] * WAA_fish[r,y,,,f,sim]) # True Catch
        ObsCatch[r,y,f,sim] <- TrueCatch[r,y,f,sim] * exp(stats::rnorm(1, 0, exp(ln_sigmaC[r,y,f]))) # Observed Catch w/ lognormal deviations

        #### Fishery Index -------------------------------------------------------------------
        if(fish_idx_type[f] == 0) TrueFishIdx[r,y,f,sim] <- fish_q[r,y,f,sim] * sum(NAA[r,y,,,sim] * fish_sel[r,y,,,f,sim]) # True Fishery Index (abundance)
        if(fish_idx_type[f] == 1) TrueFishIdx[r,y,f,sim] <- fish_q[r,y,f,sim] * sum(NAA[r,y,,,sim] * fish_sel[r,y,,,f,sim] * WAA_fish[r,y,,,f,sim]) # True Fishery Index (biomass)
        ObsFishIdx[r,y,f,sim] <- fish_q[r,y,f,sim] * TrueFishIdx[r,y,f,sim] * exp(stats::rnorm(1, 0, ObsFishIdx_SE[r,y,f])) # Observed Fishery index w/ lognormal deviations

        #### Fishery Compositions ----------------------------------------------------
        if(Fmort[r,y,f,sim] > 0) { # only simulate if Fishing Mortality > 0

          # Age Compositions (Dynamic ISS based on feedback fishing mortality)
          if(exists("ISS_FishAgeComps_fill") && isTRUE(ISS_FishAgeComps_fill == "F_pattern") && isTRUE(run_feedback) && y >= feedback_start_yr + 1 && r == 1 && f == 1) {
            ISS_FishAgeComps[,1:y,,,sim] <- predict_sim_fish_iss_fmort(ISS_FishComps = ISS_FishAgeComps, Fmort = Fmort, y = y, sim = sim)
          }

          # Length Compositions (Dynamic ISS based on feedback fishing mortality)
          if(exists("ISS_FishLenComps_fill") && isTRUE(ISS_FishLenComps_fill == "F_pattern") && isTRUE(run_feedback) && y >= feedback_start_yr + 1 && r == 1 && f == 1) {
            ISS_FishLenComps[,1:y,,,sim] <- predict_sim_fish_iss_fmort(ISS_FishComps = ISS_FishLenComps, Fmort = Fmort, y = y, sim = sim)
          }

          # Sample fishery ages
          ObsFishAgeComps <- simulate_comps(r = r,
                                            y = y,
                                            f = f,
                                            sim = sim,
                                            Exp = CAA,
                                            ISS = ISS_FishAgeComps,
                                            AgeingError = AgeingError,
                                            comp_like = comp_fishage_like,
                                            ln_theta = ln_FishAge_theta,
                                            ln_theta_agg = ln_FishAge_theta_agg,
                                            corr_pars = FishAge_corr_pars,
                                            corr_pars_agg = FishAge_corr_pars_agg,
                                            comp_type = FishAgeComps_Type,
                                            n_sexes = n_sexes,
                                            n_regions = n_regions,
                                            n_cat = n_ages,
                                            Obs = ObsFishAgeComps,
                                            age_or_len = 0)

          # Sample fishery lengths
          if(exists("SizeAgeTrans") && !is.null(SizeAgeTrans)) {
            ObsFishLenComps <- simulate_comps(r = r,
                                              y = y,
                                              f = f,
                                              sim = sim,
                                              Exp = CAL,
                                              ISS = ISS_FishLenComps,
                                              AgeingError = NULL,
                                              comp_like = comp_fishlen_like,
                                              ln_theta = ln_FishLen_theta,
                                              ln_theta_agg = ln_FishLen_theta_agg,
                                              corr_pars = FishLen_corr_pars,
                                              corr_pars_agg = FishLen_corr_pars_agg,
                                              comp_type = FishLenComps_Type,
                                              n_sexes = n_sexes,
                                              n_regions = n_regions,
                                              n_cat = n_lens,
                                              Obs = ObsFishLenComps,
                                              age_or_len = 1)
          } # end if size age transition if availiable
        } # end if Fmort > 0

      } # end f loop
    } # end r loop

      ### Generate Survey Observations ----------------------------------------------------------
      for(r in 1:n_regions) {
        for(sf in 1:n_srv_fleets) {

          # Survey Ages Indexed (midpoint year)
          SrvIAA[r,y,,,sf,sim] <- NAA[r,y,,,sim] * srv_sel[r,y,,,sf,sim] * exp(-t_srv[r,sf] * ZAA[r,y,,,sim])
          if(exists("SizeAgeTrans") && !is.null(SizeAgeTrans)) for(s in 1:n_sexes) SrvIAL[r,y,,s,sf,sim] <- SizeAgeTrans[r,y,,,s,sim] %*% SrvIAA[r,y,,s,sf,sim] # Survey index at length

          #### Survey Index ------------------------------------------------------------
          if(srv_idx_type[sf] == 0) TrueSrvIdx[r,y,sf,sim] <- srv_q[r,y,sf,sim] * sum(SrvIAA[r,y,,,sf,sim]) # True Survey Index (abundance)
          if(srv_idx_type[sf] == 1) TrueSrvIdx[r,y,sf,sim] <- srv_q[r,y,sf,sim] * sum(SrvIAA[r,y,,,sf,sim] * WAA_srv[r,y,,,sf,sim]) # True Survey Index (biomass)
          ObsSrvIdx[r,y,sf,sim] <- srv_q[r,y,sf,sim] * TrueSrvIdx[r,y,sf,sim] * exp(stats::rnorm(1, 0, ObsSrvIdx_SE[r,y,sf])) # Observed survey index w/ lognormal deviations

          #### Survey Compositions -----------------------------------------------------
          # Sample survey ages
          ObsSrvAgeComps <- simulate_comps(r = r,
                                           y = y,
                                           f = sf,
                                           sim = sim,
                                           Exp = SrvIAA,
                                           ISS = ISS_SrvAgeComps,
                                           AgeingError = AgeingError,
                                           comp_like = comp_srvage_like,
                                           ln_theta = ln_SrvAge_theta,
                                           ln_theta_agg = ln_SrvAge_theta_agg,
                                           corr_pars = SrvAge_corr_pars,
                                           corr_pars_agg = SrvAge_corr_pars_agg,
                                           comp_type = SrvAgeComps_Type,
                                           n_sexes = n_sexes,
                                           n_regions = n_regions,
                                           n_cat = n_ages,
                                           Obs = ObsSrvAgeComps,
                                           age_or_len = 0)

          # Sample survey lengths
          if(exists("SizeAgeTrans") && !is.null(SizeAgeTrans)) {
            ObsSrvLenComps <- simulate_comps(r = r,
                                             y = y,
                                             f = sf,
                                             sim = sim,
                                             Exp = SrvIAL,
                                             ISS = ISS_SrvLenComps,
                                             AgeingError = NULL,
                                             comp_like = comp_srvlen_like,
                                             ln_theta = ln_SrvLen_theta,
                                             ln_theta_agg = ln_SrvLen_theta_agg,
                                             corr_pars = SrvLen_corr_pars,
                                             corr_pars_agg = SrvLen_corr_pars_agg,
                                             comp_type = SrvLenComps_Type,
                                             n_sexes = n_sexes,
                                             n_regions = n_regions,
                                             n_cat = n_lens,
                                             Obs = ObsSrvLenComps,
                                             age_or_len = 1)
          } # end if size age transition if availiable

        } # end sf loop
      } # end r loop

      if(UseTagging == 1) {
        for(r in 1:n_regions) {
          ### Release Tags ------------------------------------------------
          # Get indices for tag cohorts in the current year and region
          tag_rel <- which(tag_release_indicator[,1] == r & tag_release_indicator[,2] == y) # Get tag cohort (release event)

          # Release Tags if any events
          if(length(tag_rel) != 0) {
            # Tag Indexing
            tag_rel_region <- tag_release_indicator[tag_rel,1] # tag release region
            tag_rel_yr <- tag_release_indicator[tag_rel,2] # tag release year

            # Release Tagged Fish
            if(!exists("n_tags_rel_input")) n_tags_rel <- round(ObsSrvIdx[tag_rel_region,tag_rel_yr,1,sim] / sum(ObsSrvIdx[,tag_rel_yr,1,sim]) * n_tags) # Number of tags apportioned across regions
            else n_tags_rel <- n_tags_rel_input[tag_rel] # input number of tag releases
            tmp_SrvAgeComps_Prob <- as.vector(SrvIAA[tag_rel_region, tag_rel_yr,,,1,sim]) # Tagged Fish Proportions at age and sex (from first survey)
            tagged_fish <- round((tmp_SrvAgeComps_Prob / sum(tmp_SrvAgeComps_Prob)) * n_tags_rel) # Tagged fish
            Tagged_Fish[tag_rel,,,sim] <- array(tagged_fish, dim = c(n_ages, n_sexes)) # Input into Tagged_Fish array and apply initial tag induced mortality
          }
        } # if using tagging data / simulating
      } # end r loop

      if(UseTagging == 1) {
        #### Generate Tag Recaptures ------------------------------------------------
        for(tag_rel in 1:n_tag_rel_events) {

          # get indexing
          tag_rel_region <- tag_release_indicator[tag_rel,1] # tag release region
          tag_rel_yr <- tag_release_indicator[tag_rel,2] # tag release year

          # Skipping stuff if hasn't occurred yet, or if max liberty
          if(tag_rel_yr > y) next # skip if tagging hasn't occurred
          recap_yr <- y - tag_rel_yr + 1 # get tag liberty
          if(recap_yr > max_liberty) next # skip if max liberty

          # Input tagged fish into available tags for recapture and adjust initial number of tagged fish for tag induced mortality (exponential mortality process)
          if(recap_yr == 1) Tag_Avail[1,tag_rel,tag_rel_region,,,sim] <- Tagged_Fish[tag_rel,,,sim] * exp(-exp(ln_Init_Tag_Mort))

          # Get mortality estimates
          tmp_F <- Get_Tagging_Mortality(tag_selex = tag_selex, # Tag Selex Type
                                         tag_natmort = tag_natmort, # Tag natural mortality
                                         # Reformat fishing mortality
                                         Fmort = array(Fmort[,,,sim,drop = FALSE], dim = c(n_regions, n_yrs, n_fish_fleets)),
                                         # Reformat natural mortality
                                         natmort = array(natmort[,,,,sim,drop = FALSE], dim = c(n_regions, n_yrs, n_ages, n_sexes)),
                                         Tag_Shed = exp(ln_Tag_Shed),
                                         # Reformat fishery selectivity
                                         fish_sel = array(fish_sel[,,,,,sim,drop = FALSE], dim = c(n_regions, n_yrs, n_ages, n_sexes, n_fish_fleets)),
                                         n_regions = n_regions,
                                         n_ages = n_ages,
                                         n_sexes = n_sexes,
                                         n_fish_fleets = n_fish_fleets,
                                         y = y,
                                         what = "F"
          )

          # Get total mortality
          tmp_ZAA <- Get_Tagging_Mortality(tag_selex = tag_selex, # tag selex
                                           tag_natmort = tag_natmort, # tag natmort
                                           # Reformat fishing mortality
                                           Fmort = array(Fmort[,,,sim,drop = FALSE], dim = c(n_regions, n_yrs, n_fish_fleets)),
                                           # Reformat natural mortality
                                           natmort = array(natmort[,,,,sim,drop = FALSE], dim = c(n_regions, n_yrs, n_ages, n_sexes)),
                                           Tag_Shed = exp(ln_Tag_Shed),
                                           # Reformat fishery selectivity
                                           fish_sel = array(fish_sel[,,,,,sim,drop = FALSE], dim = c(n_regions, n_yrs, n_ages, n_sexes, n_fish_fleets)),
                                           n_regions = n_regions,
                                           n_ages = n_ages,
                                           n_sexes = n_sexes,
                                           n_fish_fleets = n_fish_fleets,
                                           y = y,
                                           what = "Z"
          )

          # Discount with tagging time (t_tagging) if it doesn't happen at the start of the year
          if(recap_yr == 1 && t_tagging != 1) tmp_ZAA <- tmp_ZAA * t_tagging

          # Move tagged fish around
          if(t_tagging != 1 && recap_yr == 1) { # Movement does not occur if tagging does not happen at start of year
          } else{
            for(a in 1:n_ages) {
              for(s in 1:n_sexes) {
                Tag_Avail[recap_yr,tag_rel,,a,s,sim] <- Tag_Avail[recap_yr,tag_rel,,a,s,sim] %*% Movement[,,y,a,s,sim]
              } # end s loop
            } # end a loop
          }

          # Apply mortality and ageing to tagged fish
          for(a in 1:n_ages) {
            for(s in 1:n_sexes) {
              if(a < n_ages) { # If not in plus group
                Tag_Avail[recap_yr+1,tag_rel,,a+1,s,sim] <- Tag_Avail[recap_yr,tag_rel,,a,s,sim] * exp(-tmp_ZAA[,1,a,s])
              } else{ # Accumulate plus group here
                Tag_Avail[recap_yr+1,tag_rel,,n_ages,s,sim] <- Tag_Avail[recap_yr+1,tag_rel,,n_ages,s,sim] +
                  Tag_Avail[recap_yr,tag_rel,,n_ages,s,sim]  * exp(-tmp_ZAA[,1,n_ages,s])
              } # end else for in plus group
            } # end s loop
          } # end a loop

          # Apply Baranov's to get predicted recaptures
          Pred_Tag_Recap[recap_yr,tag_rel,,,,sim] <- Tag_Reporting[,y,sim] * (tmp_F[,1,,] / tmp_ZAA[,1,,]) *
                                                     Tag_Avail[recap_yr,tag_rel,,,,sim] * (1 - exp(-tmp_ZAA[,1,,]))

          # Simulate Tag Recoveries
          if(tag_like %in% c(0,1)) {
            for(r in 1:n_regions) {
              for(a in 1:n_ages) {
                for(s in 1:n_sexes) {
                  # Poisson Tag Recovery
                  if(tag_like == 0){
                    Obs_Tag_Recap[recap_yr,tag_rel,r,a,s,sim] <- rpois(n = 1, lambda = Pred_Tag_Recap[recap_yr,tag_rel,r,a,s,sim])
                  }
                  # Negative Binomial Tag Recovery
                  if(tag_like == 1) {
                    Obs_Tag_Recap[recap_yr,tag_rel,r,a,s,sim] <- rnbinom(n = 1, mu = Pred_Tag_Recap[recap_yr,tag_rel,r,a,s,sim], size = exp(ln_tag_theta))
                  }
                } # end s loop
              } # end a loop
            } # end r loop
          } # end if

          # Multinomial or Dirichlet-Multinomial tag recovery (release conditioned)
          if(tag_like %in% c(2,4)) {

            # Get number of initial tags released
            tmp_n_tags_rel <- round(sum(Tagged_Fish[tag_rel,,,sim]))
            # get recapture probabilities ordered by ages, sexes, regions
            tmp_recap <- aperm(Pred_Tag_Recap[recap_yr,tag_rel,,,,sim, drop = FALSE] / tmp_n_tags_rel, c(4,5,3,1,2,6))
            # concatenate recapture and non recapture probabilities
            tmp_probs <- c(tmp_recap, 1 - sum(tmp_recap))

            # simulate
            if(tag_like == 2) tmp_sim_recap <- stats::rmultinom(1, tmp_n_tags_rel, tmp_probs) # multinomial draws
            if(tag_like == 4) tmp_sim_recap <- rdirM(n = 1, N = tmp_n_tags_rel, exp(ln_tag_theta) * tmp_n_tags_rel * tmp_probs) # dirichlet-multinomial draws

            # remove last group (not recaptured) and then reshape into correct format
            tmp_sim_recap <- aperm(array(tmp_sim_recap[-length(tmp_sim_recap)], dim(tmp_recap)), c(4,5,3,1,2,6))
            # input recaptures from multinomial into observed array
            Obs_Tag_Recap[recap_yr,tag_rel,,,,sim] <- tmp_sim_recap

          } # end if for multinomial likelihood (release conditioned)

          # Multinomial or Dirichlet-Mutlinomial tag recovery (recovery conditioned)
          if(tag_like %in% c(3,5)) {

            # Get number of tags to simulate
            tmp_n_tags_recap <- sum(Pred_Tag_Recap[recap_yr,tag_rel,,,,sim]) # Get number of tags to simulate
            # get recapture probabilities ordered by ages, sexes, regions
            tmp_probs <- aperm(Pred_Tag_Recap[recap_yr,tag_rel,,,,sim, drop = FALSE] / tmp_n_tags_recap, c(4,5,3,1,2,6))

            # simulate
            if(tag_like == 3) tmp_sim_recap <- stats::rmultinom(1, tmp_n_tags_recap, tmp_probs) # mutlinomial draws
            if(tag_like == 5) tmp_sim_recap <- rdirM(n = 1, N = tmp_n_tags_recap, exp(ln_tag_theta) * tmp_n_tags_recap * tmp_probs) # dirichlet-multinomial draws

            # reshape into correct format
            tmp_sim_recap <- aperm(array(tmp_sim_recap, dim(tmp_probs)), c(4,5,3,1,2,6))
            # input recaptures from multinomial into observed array
            Obs_Tag_Recap[recap_yr,tag_rel,,,,sim] <- tmp_sim_recap

          } # end if for Multinomial likelihood (recovery conditioned)
        } # end tag_rel loop
      } # end if using tagging data

    ### Compute Recruitment for Next Year ---------------------------------------
    tmp_ln_RecDevs_next <- NULL # Initialize container vector for next year's recruitment deviations
    if(y < n_yrs) { # Get Recruitment Deviations for next year
      for(r in 1:n_regions) {
        if (exists("Rec_input")) {
          if ((y + 1 <= dim(Rec_input)[2]) || (run_feedback == FALSE) || (run_feedback == TRUE && y < feedback_start_yr)) {
            for (s in 1:n_sexes) NAA[r, y+1, 1, s, sim] <- Rec_input[r, y+1, sim] * sexratio[r, y+1, s, sim]
          }
        } else {
          # Global Recruitment Deviations
          if(rec_dd == 0 && is.null(tmp_ln_RecDevs_next)) {
            tmp_ln_RecDevs_next <- stats::rnorm(1, 0, exp(ln_sigmaR[2]))
          }

          # Local Recruitment Deviations
          if(rec_dd == 1) {
            tmp_ln_RecDevs_next <- ln_RecDevs[r,y+1,sim] <- stats::rnorm(1, 0, exp(ln_sigmaR[2]))
          }

          ln_RecDevs[r,y+1,sim] = tmp_ln_RecDevs_next # Input recruitment deviations into vector

          # Get Deterministic Recruitment for next year
          tmp_det_rec_next <- Get_Det_Recruitment(recruitment_model = recruitment_opt,
                                                  recruitment_dd = rec_dd,
                                                  y = y+1,
                                                  rec_lag = rec_lag,
                                                  R0 = sum(R0[,y,sim]), # sum to get global R0
                                                  Rec_Prop = R0[,y,sim] / sum(R0[,y,sim]), # get R0 proportion
                                                  h = h[,y,sim],
                                                  n_regions = n_regions,
                                                  n_ages = n_ages,
                                                  WAA = array(WAA[,y,,1,sim], dim = c(n_regions, n_ages)),
                                                  MatAA = array(MatAA[,y,,1,sim], dim = c(n_regions, n_ages)),
                                                  natmort = array(natmort[,y,,1,sim], dim = c(n_regions, n_ages)),
                                                  SSB_vals = array(SSB[,,sim], dim = c(n_regions, n_yrs))
          )

          # Store next year's recruitment (will be used when next year starts)
          for(s in 1:n_sexes)  NAA[r,y+1,1,s,sim] <- tmp_det_rec_next[r] * exp(ln_RecDevs[r,y+1,sim] - exp(ln_sigmaR[2])^2/2) * sexratio[r,y+1,s,sim]
        }

        Rec[r,y+1,sim] <- sum(NAA[r,y+1,1,,sim]) # Save recruitment estimates
        NAA0[r,y+1,1,,sim] <- NAA[r,y+1,1,,sim] # populate unfished NAA

      } # end r loop
    } # end if

  }) # end simulation environment

  return(invisible(NULL))

}

#' Simulates a static spatial, sex, and age-structured population (no feedback loop)
#'
#' @param output_path path to output simulation objects
#' @param sim_list Simulation list objects
#'
#' @returns a list object with a bunch of simulated values and outputs
#' @export Simulate_Pop_Static
#' @family Simulation Setup
Simulate_Pop_Static <- function(sim_list,
                                output_path = NULL) {

  # Setup simulation environment
  sim_env <- Setup_sim_env(sim_list)

  # Start Simulation
  for (sim in 1:sim_env$n_sims) {
    for (y in 1:sim_env$n_yrs) {
      # Run annual cycle here
      run_annual_cycle(y, sim, sim_env)
    } # end y loop
  } # end sim loop

  # Output simulation outputs as a list
  sim_out <- list(init_F = sim_env$init_F,
                  Fmort = sim_env$Fmort,
                  ln_sigmaC = sim_env$ln_sigmaC,
                  fish_sel = sim_env$fish_sel,
                  fish_q = sim_env$fish_q,
                  ln_RecDevs = sim_env$ln_RecDevs,
                  ln_InitDevs = sim_env$ln_InitDevs,
                  natmort = sim_env$natmort,
                  ZAA = sim_env$ZAA,
                  sexratio = sim_env$sexratio,
                  R0 = sim_env$R0,
                  Rec = sim_env$Rec,
                  WAA = sim_env$WAA,
                  WAA_fish = sim_env$WAA_fish,
                  WAA_srv = sim_env$WAA_srv,
                  MatAA = sim_env$MatAA,
                  ln_sigmaR = sim_env$ln_sigmaR,
                  Movement = sim_env$Movement,
                  Init_NAA = sim_env$Init_NAA,
                  NAA = sim_env$NAA,
                  NAA0 = sim_env$NAA0,
                  Dynamic_SSB0 = sim_env$Dynamic_SSB0,
                  SSB = sim_env$SSB,
                  Total_Biom = sim_env$Total_Biom,
                  TrueCatch = sim_env$TrueCatch,
                  ObsCatch = sim_env$ObsCatch,
                  ObsFishIdx = sim_env$ObsFishIdx,
                  TrueFishIdx = sim_env$TrueFishIdx,
                  ObsFishIdx_SE = sim_env$ObsFishIdx_SE,
                  CAA = sim_env$CAA,
                  CAL = sim_env$CAL,
                  ObsFishAgeComps = sim_env$ObsFishAgeComps,
                  ObsFishLenComps = sim_env$ObsFishLenComps,
                  ObsSrvIdx = sim_env$ObsSrvIdx,
                  TrueSrvIdx = sim_env$TrueSrvIdx,
                  ObsSrvIdx_SE = sim_env$ObsSrvIdx_SE,
                  SrvIAA = sim_env$SrvIAA,
                  SrvIAL = sim_env$SrvIAL,
                  srv_sel = sim_env$srv_sel,
                  srv_q = sim_env$srv_q,
                  ObsSrvAgeComps = sim_env$ObsSrvAgeComps,
                  ObsSrvLenComps = sim_env$ObsSrvLenComps,
                  tag_release_indicator = as.matrix(sim_env$tag_release_indicator),
                  Tag_Reporting = sim_env$Tag_Reporting,
                  Tagged_Fish = sim_env$Tagged_Fish,
                  ln_Init_Tag_Mort = sim_env$ln_Init_Tag_Mort,
                  ln_Tag_Shed = sim_env$ln_Tag_Shed,
                  Tag_Avail = sim_env$Tag_Avail,
                  UseTagging = sim_env$UseTagging,
                  Pred_Tag_Recap = sim_env$Pred_Tag_Recap,
                  Obs_Tag_Recap = sim_env$Obs_Tag_Recap,
                  SizeAgeTrans = if(!is.null(sim_env$SizeAgeTrans)) sim_env$SizeAgeTrans else NULL,
                  AgeingError = sim_env$AgeingError,
                  ISS_FishAgeComps = sim_env$ISS_FishAgeComps,
                  ISS_FishLenComps = sim_env$ISS_FishLenComps,
                  ISS_SrvAgeComps = sim_env$ISS_SrvAgeComps,
                  ISS_SrvLenComps = sim_env$ISS_SrvLenComps,
                  n_sims = sim_env$n_sims,
                  n_regions = sim_env$n_regions,
                  n_years = sim_env$n_yrs,
                  n_ages = sim_env$n_ages,
                  n_lens = if(!is.null(sim_env$n_lens)) sim_env$n_lens else NULL,
                  n_sexes = 1:sim_env$n_sexes,
                  n_fish_fleets = sim_env$n_fish_fleets,
                  n_srv_fleets = sim_env$n_srv_fleets
                  )

  # save RDS file
  if(!is.null(output_path)) saveRDS(sim_out, file = output_path)

  return(sim_out)

} # end function


#' Conduct a Simulation Self Test
#'
#' This function runs a self test of the fitted RTMB model by simulating new
#' datasets under the fitted parameters, refitting the model, and comparing
#' estimated outputs to the true values used for simulation. It can be run
#' sequentially or in parallel.
#'
#' @param data A list containing model data from an RTMB object.
#' @param parameters A list of fitted parameter values from an RTMB object.
#' @param mapping A list specifying parameter mappings from an RTMB object.
#' @param random Character vector specifying random effects.
#' @param rep A list of report values from an RTMB object (`$rep`).
#' @param sd_rep An `sdreport` object from RTMB summarizing parameter uncertainty.
#' @param n_sims Integer. Number of simulation replicates to run.
#' @param newton_loops Integer. Number of Newton loops used in model fitting (default: `3`).
#' @param do_sdrep Logical. If `TRUE`, compute `sdreport` for each fitted replicate (default: `FALSE`).
#' @param do_par Logical. If `TRUE`, run simulations in parallel (default: `FALSE`).
#' @param n_cores Integer. Number of cores to use for parallelization (default: `NULL` = detect automatically).
#' @param output_path Optional file path. If provided, the simulated datasets are written to this location.
#' @param what Character vector. Names of report elements in `rep` to extract and store for each replicate.
#' @family Simulation Setup
#' @return
#' A list with elements corresponding to the requested `what` values, each containing
#' an array of simulation results across replicates. If `do_sdrep = TRUE`, an additional
#' element `"sd_rep"` is included with the list of `sdreport` objects (or `NA` if a replicate fails).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Run a simple self test with 10 simulations, extracting SSB
#' res <- simulation_self_test(
#'   data = model$data,
#'   parameters = model$parameters,
#'   mapping = model$mapping,
#'   random = model$random,
#'   rep = model$rep,
#'   sd_rep = model$sd_rep,
#'   n_sims = 10,
#'   what = "SSB"
#' )
#'
#' str(res$SSB) # look at simulated SSB arrays
#' }
simulation_self_test <- function(data,
                                 parameters,
                                 mapping,
                                 random,
                                 rep,
                                 sd_rep,
                                 n_sims,
                                 newton_loops = 3,
                                 do_sdrep = FALSE,
                                 do_par = FALSE,
                                 n_cores = NULL,
                                 output_path = NULL,
                                 what = c('SSB', 'Rec')
                                 ) {

  missing_names <- setdiff(what, names(rep))
  if(length(missing_names) > 0)  stop(paste("The following elements in 'what' are not found in rep:",  paste(missing_names, collapse = ", ")))
  optim_parameters_list <- get_optim_param_list(parameters, mapping, sd_rep, random) # get optimized parameters in original list format

  # Modify any data weights that are NA to 0
  if(any(is.na(data$Wt_Catch))) data$Wt_Catch[is.na(data$Wt_Catch)] <- 0
  if(any(is.na(data$Wt_FishAgeComps))) data$Wt_FishAgeComps[is.na(data$Wt_FishAgeComps)] <- 0
  if(any(is.na(data$Wt_FishLenComps))) data$Wt_FishLenComps[is.na(data$Wt_FishLenComps)] <- 0
  if(any(is.na(data$Wt_FishIdx))) data$Wt_FishIdx[is.na(data$Wt_FishIdx)] <- 0
  if(any(is.na(data$Wt_SrvAgeComps))) data$Wt_SrvAgeComps[is.na(data$Wt_SrvAgeComps)] <- 0
  if(any(is.na(data$Wt_SrvLenComps))) data$Wt_SrvLenComps[is.na(data$Wt_SrvLenComps)] <- 0
  if(any(is.na(data$Wt_SrvIdx))) data$Wt_SrvIdx[is.na(data$Wt_SrvIdx)] <- 0
  if(any(is.na(data$Wt_Tagging))) data$Wt_Tagging[is.na(data$Wt_Tagging)] <- 0

  # Setup Model Dimensions --------------------------------------------------
  sim_list <- Setup_Sim_Dim(n_sims = n_sims, # number of simulations
                            n_yrs = length(data$years), # number of years
                            n_regions = data$n_regions,  # number of regions
                            n_ages = length(data$ages), # number of ages
                            n_obs_ages = dim(data$ObsFishAgeComps)[3], # number of observed ages
                            n_lens = length(data$lens), # number of lengths
                            n_sexes = data$n_sexes, # number of sexes
                            n_fish_fleets = data$n_fish_fleets, # number of fishery fleets
                            n_srv_fleets = data$n_srv_fleets # number of survey fleets
  )

  # Setup Simulation Containers ---------------------------------------------
  sim_list <- Setup_Sim_Containers(sim_list)

  # Setup Fishing Processes -------------------------------------------------
  ln_sigmaC <- array(NA, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_fish_fleets)) # setup sigmaC container
  # Loop through to populate ln_sigmaC with associated weights
  for(r in 1:sim_list$n_regions) for(f in 1:sim_list$n_fish_fleets) {
    if(!is.vector(data$Wt_Catch)) ln_sigmaC[r,,f] <- log(exp(optim_parameters_list$ln_sigmaC[r,,f]) / sqrt(data$Wt_Catch[r,,f]))
    else ln_sigmaC[r,,f] <- log(exp(optim_parameters_list$ln_sigmaC[r,,f]) / sqrt(data$Wt_Catch))
  }

  # setup fishery simulation processes
  sim_list <- Setup_Sim_Fishing(sim_list = sim_list, # update simulate list
                                ln_sigmaC = ln_sigmaC, # sigmaC
                                Fmort_input = replicate(n = sim_list$n_sims, rep$Fmort), # use fishing mortality from report
                                fish_sel_input = replicate(n = sim_list$n_sims, rep$fish_sel), # use fishery selectivity from report
                                fish_q_input = replicate(n = sim_list$n_sims, rep$fish_q), # use fishery catchability from report
                                ObsFishIdx_SE = data$ObsFishIdx_SE / sqrt(data$Wt_FishIdx), # fishery index uncertainty
                                fish_idx_type = data$fish_idx_type, # fishery index type
                                init_F_val = rep$init_F,

                                # fishery age composition specifications
                                comp_fishage_like = data$FishAgeComps_LikeType, # age comps likelihood
                                FishAgeComps_Type = data$FishAgeComps_Type, # age comps structure
                                ISS_FishAgeComps = replicate(sim_list$n_sims, data$ISS_FishAgeComps[,,,,drop = F] * data$Wt_FishAgeComps), # input sample size
                                ln_FishAge_theta = optim_parameters_list$ln_FishAge_theta[,,,drop = F] , # overdispersion
                                ln_FishAge_theta_agg = optim_parameters_list$ln_FishAge_theta_agg, # aggregated overdispersion for fishery age comps
                                FishAge_corr_pars_agg = optim_parameters_list$FishAge_corr_pars_agg, # correaltion parameters for aggregated fishery age comps
                                FishAge_corr_pars = optim_parameters_list$FishAge_corr_pars[,,,,drop = F], # correlation parameters for fishery age comps

                                # fishery length composition specifications
                                comp_fishlen_like = data$FishLenComps_LikeType, # length comp likelihood
                                FishLenComps_Type = data$FishLenComps_Type, # length comps structure
                                ISS_FishLenComps = replicate(sim_list$n_sims, data$ISS_FishLenComps[,,,,drop = F] * data$Wt_FishLenComps), # input sample size
                                ln_FishLen_theta = optim_parameters_list$ln_FishLen_theta[,,,drop = F], # overdispersion
                                ln_FishLen_theta_agg = optim_parameters_list$ln_FishLen_theta_agg, # aggregated overdispersion
                                FishLen_corr_pars_agg = optim_parameters_list$FishLen_corr_pars_agg, # correaltion parameters for aggregated fish len comps
                                FishLen_corr_pars = optim_parameters_list$FishLen_corr_pars[,,,,drop = F] # correlation parameters for fish len comps
  )


  # Setup Survey Processes --------------------------------------------------
  sim_list <- Setup_Sim_Survey(
    sim_list = sim_list,
    srv_sel_input = replicate(n = sim_list$n_sims, rep$srv_sel),
    srv_q_input = replicate(n = sim_list$n_sims, rep$srv_q),
    ObsSrvIdx_SE = data$ObsSrvIdx_SE / sqrt(data$Wt_SrvIdx), # survey observation error
    srv_idx_type = data$srv_idx_type,
    t_srv = data$t_srv,

    # survey age composition specifications
    comp_srvage_like = data$SrvAgeComps_LikeType, # age comps likelihood
    SrvAgeComps_Type = data$SrvAgeComps_Type, # age comps structure
    ISS_SrvAgeComps = replicate(sim_list$n_sims, data$ISS_SrvAgeComps[,,,,drop = F] * data$Wt_SrvAgeComps), # input sample size
    ln_SrvAge_theta = optim_parameters_list$ln_SrvAge_theta[,,,drop = F] , # overdispersion
    ln_SrvAge_theta_agg = optim_parameters_list$ln_SrvAge_theta_agg, # aggregated overdispersion for survey age comps
    SrvAge_corr_pars_agg = optim_parameters_list$SrvAge_corr_pars_agg, # correaltion parameters for aggregated survey age comps
    SrvAge_corr_pars = optim_parameters_list$SrvAge_corr_pars[,,,,drop = F], # correlation parameters for survey age comps

    # survey length composition specifications
    comp_srvlen_like = data$SrvLenComps_LikeType, # length comp likelihood
    SrvLenComps_Type = data$SrvLenComps_Type, # length comps structure
    ISS_SrvLenComps = replicate(sim_list$n_sims, data$ISS_SrvLenComps[,,,,drop = F] * data$Wt_SrvLenComps), # input sample size
    ln_SrvLen_theta = optim_parameters_list$ln_SrvLen_theta[,,,drop = F], # overdispersion
    ln_SrvLen_theta_agg = optim_parameters_list$ln_SrvLen_theta_agg, # aggregated overdispersion
    SrvLen_corr_pars_agg = optim_parameters_list$SrvLen_corr_pars_agg, # correaltion parameters for aggregated srv len comps
    SrvLen_corr_pars = optim_parameters_list$SrvLen_corr_pars[,,,,drop = F] # correlation parameters for srv len comps
  )

  # Setup Biological Dynamics -----------------------------------------------
  sim_list <- Setup_Sim_Biologicals(
    sim_list = sim_list, # simualtion list
    natmort_input = replicate(n = sim_list$n_sims, rep$natmort), # natuyral mortality
    WAA_input = replicate(n = sim_list$n_sims, data$WAA), # weight at age
    WAA_fish_input = replicate(n = sim_list$n_sims, data$WAA_fish), # fishery weight at age
    WAA_srv_input = replicate(n = sim_list$n_sims, data$WAA_srv), # survey weight at age
    MatAA_input = replicate(n = sim_list$n_sims, data$MatAA), # maturity at age
    AgeingError_input = replicate(n = sim_list$n_sims, data$AgeingError), # ageing error
    SizeAgeTrans_input = replicate(n = sim_list$n_sims, data$SizeAgeTrans) # size age transition matrix
  )

  # Movement
  sim_list$Movement <- replicate(n = sim_list$n_sims, rep$Movement)

  # Setup Recruitment Processes ---------------------------------------------
  sim_list <- Setup_Sim_Rec(
    sim_list = sim_list,
    do_recruits_move = data$do_recruits_move, # whether recruits move
    t_spawn = data$t_spawn, # spawn timing
    init_age_strc = data$init_age_strc, # initilaizing age structure
    h_input = replicate(n = sim_list$n_sims, array(rep$h_trans, dim = c(sim_list$n_regions, sim_list$n_yrs))), # steepness
    R0_input = replicate(n = sim_list$n_sims, expr = array(rep$R0 * rep$Rec_trans_prop, dim = c(sim_list$n_regions, sim_list$n_yrs))), # R0
    sexratio_input = replicate(n = sim_list$n_sims, expr = rep$sexratio), # sex ratio
    ln_sigmaR = optim_parameters_list$ln_sigmaR, # ln_sigmaR
    Rec_input = replicate(n = sim_list$n_sims, expr = rep$Rec), # recruitment time series
    ln_InitDevs_input = replicate(sim_list$n_sims, optim_parameters_list$ln_InitDevs) # init devs
  )

  # Setup Tagging -----------------------------------------------------------
  if(!is.na(sum(data$Tagged_Fish))) n_tags_rel_input <- apply(data$Tagged_Fish, 1, sum) else n_tags_rel_input <- NA
  if(exists("tag_release_indicator", data)) tag_release_indicator <- data$tag_release_indicator  else tag_release_indicator <- NA
  Tag_Reporting_input <- if(!is.null(rep$Tag_Reporting)) replicate(n = sim_list$n_sims, rep$Tag_Reporting) else NULL

  sim_list <- Setup_Sim_Tagging(
    sim_list = sim_list, # simulation list
    max_liberty = data$max_tag_liberty, # maximum tag liberty
    t_tagging = data$t_tagging, # time of tagging
    n_tags_rel_input = n_tags_rel_input * data$Wt_Tagging,  # number of tags to release per event
    tag_release_indicator = tag_release_indicator,  # tag release indicator
    ln_Init_Tag_Mort = optim_parameters_list$ln_Init_Tag_Mort,  # inital tagging mortality
    ln_Tag_Shed = optim_parameters_list$ln_Tag_Shed, # chronic tag shedding
    Tag_Reporting_input = Tag_Reporting_input, # tag reporting rates
    UseTagging = data$UseTagging, # whether or not tagging is used / simulated
    tag_selex = data$tag_selex, # tag selectivity type
    tag_natmort = data$tag_natmort, # tag natural mortality type
    tag_like = data$Tag_LikeType, # tag likelihood
    ln_tag_theta = parameters$ln_tag_theta # tag overdispersion
  )


  # Run Simulation ----------------------------------------------------------

  # storage
  store_res_list <- vector("list", length(what) + 1) # get list
  names(store_res_list) <- c(what, "sd_rep") # name list
  for(j in 1:length(what)) store_res_list[[j]] <- vector("list", n_sims) # stick in n_sims lists into storage

  sim_obj <- Simulate_Pop_Static(sim_list = sim_list, output_path = output_path) # get simulated datasets

  if(do_par == FALSE) {

    for(i in 1:n_sims) {

      tryCatch({

        # set up data stuff
        tmp_data <- data # set up temporary data list
        tmp_pars <- parameters # set up temporary parameter list
        tmp_data$ObsFishIdx <- array(sim_obj$ObsFishIdx[,,,i], dim = dim(tmp_data$ObsFishIdx)) # new fish index
        tmp_data$ObsSrvIdx <- array(sim_obj$ObsSrvIdx[,,,i], dim = dim(tmp_data$ObsSrvIdx)) # new srv index
        tmp_data$ObsCatch <- array(sim_obj$ObsCatch[,,,i], dim = dim(tmp_data$ObsCatch)) # new catch
        tmp_data$ObsFishAgeComps <- array(sim_obj$ObsFishAgeComps[,,,,,i], dim = dim(tmp_data$ObsFishAgeComps)) # new fishery ages
        tmp_data$ObsSrvAgeComps  <- array(sim_obj$ObsSrvAgeComps[,,,,,i], dim = dim(tmp_data$ObsSrvAgeComps)) # new srv ages
        tmp_data$ObsFishLenComps <- array(sim_obj$ObsFishLenComps[,,,,,i], dim = dim(tmp_data$ObsFishLenComps)) # new fishery lens
        tmp_data$ObsSrvLenComps  <- array(sim_obj$ObsSrvLenComps[,,,,,i], dim = dim(tmp_data$ObsSrvLenComps)) # new survey lens

        # setup tagging data stuff if tagging is done
        if(tmp_data$UseTagging == 1) {
          tmp_data$Wt_Tagging <- 1 # set to 1 because ess (original n_tags * Wt_Tagging) already baked into simulated data (tagged fish); avoids double-weighting
          tmp_data$Tagged_Fish <- array(sim_obj$Tagged_Fish[,,,i], dim = dim(tmp_data$Tagged_Fish)) # new tagged fish
          tmp_data$Obs_Tag_Recap <- array(sim_obj$Obs_Tag_Recap[,,,,,i], dim = dim(tmp_data$Obs_Tag_Recap)) # new tag recaps
          tmp_data$tag_release_indicator <- sim_obj$tag_release_indicator # release indicator
        }

        # Fit model
        obj <- fit_model(
          data = tmp_data,
          parameters = parameters,
          mapping = mapping,
          random = random,
          newton_loops = newton_loops,
          silent = TRUE
        )

        # Populate results into store list
        for(j in 1:length(what)) store_res_list[[j]][[i]] <- obj$rep[[what[j]]]

        if(do_sdrep == TRUE) {
          tryCatch({
            obj$sd_rep <- RTMB::sdreport(obj)
            store_res_list[[length(what) + 1]][[i]] <- obj$sd_rep # input sd report
          }, error = function(e) {
            store_res_list[[length(what) + 1]][[i]] <- NA
          })
        }

      }, error = function(e) {
        # Skip failed simulations
        for(j in 1:length(what)) store_res_list[[j]][[i]] <<- NA
        if(do_sdrep == TRUE) store_res_list[[length(what) + 1]][[i]] <<- NA
      })

    } # end i loop

    # Convert result lists to array
    for(j in 1:length(what)) store_res_list[[j]] <- simplify2array(store_res_list[[j]])

  } # not doing parallelization

  if(do_par == TRUE) {

    options(future.globals.maxSize = 5e3 * 1024^2)  # increase parrlalel size
    future::plan(future::multisession, workers = n_cores) # set up cores
    progressr::with_progress({
      p <- progressr::progressor(along = 1:n_sims) # progress bar

      sim_results <- future.apply::future_lapply(1:n_sims, function(i) {

        tryCatch({
          tmp_data <- data # set up temporary data list
          tmp_data$ObsFishIdx <- array(sim_obj$ObsFishIdx[,,,i], dim = dim(tmp_data$ObsFishIdx)) # new fish index
          tmp_data$ObsSrvIdx <- array(sim_obj$ObsSrvIdx[,,,i], dim = dim(tmp_data$ObsSrvIdx)) # new srv index
          tmp_data$ObsCatch <- array(sim_obj$ObsCatch[,,,i], dim = dim(tmp_data$ObsCatch)) # new catch
          tmp_data$ObsFishAgeComps <- array(sim_obj$ObsFishAgeComps[,,,,,i], dim = dim(tmp_data$ObsFishAgeComps)) # new fishery ages
          tmp_data$ObsSrvAgeComps  <- array(sim_obj$ObsSrvAgeComps[,,,,,i], dim = dim(tmp_data$ObsSrvAgeComps)) # new srv ages
          tmp_data$ObsFishLenComps <- array(sim_obj$ObsFishLenComps[,,,,,i], dim = dim(tmp_data$ObsFishLenComps)) # new fishery lens
          tmp_data$ObsSrvLenComps  <- array(sim_obj$ObsSrvLenComps[,,,,,i], dim = dim(tmp_data$ObsSrvLenComps)) # new survey lens

          # setup tagging data stuff if tagging is done
          if(tmp_data$UseTagging == 1) {
            tmp_data$Tagged_Fish <- array(sim_obj$Tagged_Fish[,,,i], dim = dim(tmp_data$Tagged_Fish)) # new tagged fish
            tmp_data$Obs_Tag_Recap <- array(sim_obj$Obs_Tag_Recap[,,,,,i], dim = dim(tmp_data$Obs_Tag_Recap)) # new tag recaps
            tmp_data$tag_release_indicator <- sim_obj$tag_release_indicator # release indicator
          }

          # Fit model
          obj <- fit_model(
            data = tmp_data,
            parameters = parameters,
            mapping = mapping,
            random = random,
            newton_loops = newton_loops,
            silent = TRUE
          )

          # Extract what we need and return
          result <- list()
          for(j in 1:length(what)) result[[what[j]]] <- obj$rep[[what[j]]]

          if(do_sdrep == TRUE) {
            tryCatch({
              obj$sd_rep <- RTMB::sdreport(obj) # get sdreport
              result[[length(what) + 1]] <- obj$sd_rep # input sd report
            }, error = function(e) {
              result[[length(what) + 1]] <- NA
            })
          }

          p() # update progress
          return(result)

        }, error = function(e) {
          # Skip failed simulations
          result <- list()
          for(j in 1:length(what)) result[[what[j]]] <- NA
          if(do_sdrep == TRUE) result[[length(what) + 1]] <- NA

          p() # update progress
          return(result)
        })

      }, future.seed = TRUE)

      future::plan(future::sequential)  # Reset
    })

    # Populate results from parallel run
    for(i in 1:n_sims) for(j in 1:length(what)) store_res_list[[j]][[i]] <- sim_results[[i]][[what[j]]]
    if(do_sdrep == TRUE) for(i in 1:n_sims) store_res_list[[length(what) + 1]][[i]] <- sim_results[[i]][[length(what) + 1]]
    for(j in 1:length(what)) store_res_list[[j]] <- simplify2array(store_res_list[[j]])  # Convert lists to array
  }

  return(store_res_list)

}

#' Extract simulation data into SPoRC format
#'
#' This function subsets and reshapes biological, tagging, fishery, and survey
#' data from a simulation environment for use in SPoRC analyses.
#'
#' @param sim_env A simulation environment / object (list or environment) containing
#'   arrays of biological quantities, tagging information, fishery data,
#'   and survey data.
#' @param y Integer. Number of years to retain (subset from `1:y`).
#' @param sim Integer. Simulation replicate index to extract.
#' @family Simulation Setup
#' @return A named list with the following elements:
#' \describe{
#'   \item{WAA}{Weight-at-age array [region × year × age × sex].}
#'   \item{MatAA}{Maturity-at-age array [region × year × age × sex].}
#'   \item{SizeAgeTrans}{Size–age transition array [region × year × length × age × sex].}
#'   \item{AgeingError}{Ageing error matrix [year × age × error × sim].}
#'   \item{tag_release_indicator}{Tag release indicators (or `NULL` if tagging not used).}
#'   \item{Obs_Tag_Recap}{Observed tag recapture array (or `NULL`).}
#'   \item{Tagged_Fish}{Tagged fish counts (or `NULL`).}
#'   \item{n_tag_cohorts}{Number of tag release cohorts (or `NULL`).}
#'   \item{ObsCatch}{Observed fishery catch array [region × year × fleet].}
#'   \item{ln_sigmaC}{Log Fishery Catch SD [region × year × fleet].}
#'   \item{UseCatch}{Binary indicator array for catch data availability.}
#'   \item{ObsFishIdx}{Observed fishery index array [region × year × fleet].}
#'   \item{ObsFishIdx_SE}{Standard error for fishery index array.}
#'   \item{UseFishIdx}{Binary indicator array for fishery indices.}
#'   \item{ObsFishAgeComps}{Observed fishery age composition array.}
#'   \item{ObsFishLenComps}{Observed fishery length composition array.}
#'   \item{ISS_FishAgeComps}{Implied sample sizes for fishery age compositions.}
#'   \item{ISS_FishLenComps}{Implied sample sizes for fishery length compositions.}
#'   \item{UseFishAgeComps}{Binary indicator array for fishery age comps.}
#'   \item{UseFishLenComps}{Binary indicator array for fishery length comps.}
#'   \item{ObsSrvIdx}{Observed survey index array [region × year × fleet].}
#'   \item{ObsSrvIdx_SE}{Standard error for survey index array.}
#'   \item{UseSrvIdx}{Binary indicator array for survey indices.}
#'   \item{ObsSrvAgeComps}{Observed survey age composition array.}
#'   \item{ObsSrvLenComps}{Observed survey length composition array.}
#'   \item{ISS_SrvAgeComps}{Implied sample sizes for survey age compositions.}
#'   \item{ISS_SrvLenComps}{Implied sample sizes for survey length compositions.}
#'   \item{UseSrvAgeComps}{Binary indicator array for survey age comps.}
#'   \item{UseSrvLenComps}{Binary indicator array for survey length comps.}
#' }
#'
#' @export simulation_data_to_SPoRC
simulation_data_to_SPoRC <- function(sim_env,
                                     y,
                                     sim) {

    # Biologicals
    WAA <- array(sim_env$WAA[,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_ages, sim_env$n_sexes))
    WAA_fish <- array(sim_env$WAA_fish[,1:y,,,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_ages, sim_env$n_sexes, sim_list$n_fish_fleets))
    WAA_srv <- array(sim_env$WAA_srv[,1:y,,,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_ages, sim_env$n_sexes, sim_list$n_srv_fleets))
    MatAA <- array(sim_env$MatAA[,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_ages, sim_env$n_sexes))
    SizeAgeTrans <- if(!is.null(sim_env$SizeAgeTrans)) {
      array(sim_env$SizeAgeTrans[,1:y,,,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_lens, sim_env$n_ages, sim_env$n_sexes))
      } else NULL
    AgeingError <- sim_env$AgeingError[1:y,,,sim]

    # Tagging
    if(sim_env$UseTagging == 1) {
      keep_tag_cohorts <- which(sim_env$tag_release_indicator[,2] %in% 1:y)
      tag_release_indicator <- sim_env$tag_release_indicator[keep_tag_cohorts,,drop = FALSE]
      Obs_Tag_Recap <- sim_env$Obs_Tag_Recap[,keep_tag_cohorts,,,sim, drop = FALSE]
      Tagged_Fish <- sim_env$Tagged_Fish[,keep_tag_cohorts,sim, drop = FALSE]
      n_tag_cohorts <- nrow(tag_release_indicator)
    } else {
      tag_release_indicator = Obs_Tag_Recap = Tagged_Fish = n_tag_cohorts = NULL
    }

    # Fishery Catches
    ObsCatch <- array(sim_env$ObsCatch[,1:y,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_fish_fleets))
    ln_sigmaC <- array(sim_env$ln_sigmaC[,1:y,, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_fish_fleets))
    UseCatch <- array(0, dim = dim(ObsCatch))
    UseCatch[!is.na(ObsCatch) & ObsCatch > 0] <- 1

    # Fishery Indices
    ObsFishIdx <- array(sim_env$ObsFishIdx[, 1:y,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_fish_fleets))
    ObsFishIdx_SE <- array(sim_env$ObsFishIdx_SE[, 1:y,, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_fish_fleets))
    UseFishIdx <- array(0, dim = dim(ObsFishIdx))
    UseFishIdx[!is.na(ObsFishIdx) & ObsFishIdx > 0] <- 1

    # Fishery Compositions
    ObsFishAgeComps <- array(sim_env$ObsFishAgeComps[, 1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), dim(sim_env$AgeingError)[3], sim_env$n_sexes, sim_env$n_fish_fleets))
    ObsFishLenComps <- if(!is.null(sim_env$n_lens)) {
      array(sim_env$ObsFishLenComps[, 1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_lens, sim_env$n_sexes, sim_env$n_fish_fleets))
    } else NULL
    ISS_FishAgeComps <- array(sim_env$ISS_FishAgeComps[, 1:y,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_sexes, sim_env$n_fish_fleets))
    ISS_FishLenComps <- array(sim_env$ISS_FishLenComps[, 1:y,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_sexes, sim_env$n_fish_fleets))
    UseFishAgeComps <- apply(ObsFishAgeComps, c(1,2,5), sum)
    UseFishAgeComps[!is.na(UseFishAgeComps) & UseFishAgeComps > 0] <- 1
    if(!is.null(sim_env$n_lens)) {
      UseFishLenComps <- apply(ObsFishLenComps, c(1,2,5), sum)
      UseFishLenComps[!is.na(UseFishLenComps) & UseFishLenComps > 0] <- 1
    } else UseFishLenComps <- array(0, dim = c(sim_env$n_regions, length(1:y), sim_env$n_fish_fleets))

    # Survey Indices
    ObsSrvIdx <- array(sim_env$ObsSrvIdx[, 1:y,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_srv_fleets))
    ObsSrvIdx_SE <- array(sim_env$ObsSrvIdx_SE[, 1:y,, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_srv_fleets))
    UseSrvIdx <- array(0, dim = dim(ObsSrvIdx))
    UseSrvIdx[!is.na(ObsSrvIdx) & ObsSrvIdx > 0] <- 1

    # Survey Compositions
    ObsSrvAgeComps <- array(sim_env$ObsSrvAgeComps[, 1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), dim(sim_env$AgeingError)[3], sim_env$n_sexes, sim_env$n_srv_fleets))
    ObsSrvLenComps <- if(!is.null(sim_env$n_lens)) {
      array(sim_env$ObsSrvLenComps[, 1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_lens, sim_env$n_sexes, sim_env$n_srv_fleets))
    } else NULL
    ISS_SrvAgeComps <- array(sim_env$ISS_SrvAgeComps[, 1:y,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_sexes, sim_env$n_srv_fleets))
    ISS_SrvLenComps <- array(sim_env$ISS_SrvLenComps[, 1:y,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_sexes, sim_env$n_srv_fleets))
    UseSrvAgeComps <- apply(ObsSrvAgeComps, c(1,2,5), sum)
    UseSrvAgeComps[!is.na(UseSrvAgeComps) & UseSrvAgeComps > 0] <- 1
    if(!is.null(sim_env$n_lens)) {
      UseSrvLenComps <- apply(ObsSrvLenComps, c(1,2,5), sum)
      UseSrvLenComps[!is.na(UseSrvLenComps) & UseSrvLenComps > 0] <- 1
    } else UseSrvLenComps <- array(0, dim = c(sim_env$n_regions, length(1:y), sim_env$n_fish_fleets))

    return(list(WAA = WAA,
                WAA_fish = WAA_fish,
                WAA_srv = WAA_srv,
                MatAA = MatAA,
                SizeAgeTrans = SizeAgeTrans,
                AgeingError = AgeingError,
                tag_release_indicator = tag_release_indicator,
                Obs_Tag_Recap = Obs_Tag_Recap,
                Tagged_Fish = Tagged_Fish,
                n_tag_cohorts = n_tag_cohorts,
                ObsCatch = ObsCatch,
                ln_sigmaC = ln_sigmaC,
                UseCatch = UseCatch,
                ObsFishIdx = ObsFishIdx,
                ObsFishIdx_SE = ObsFishIdx_SE,
                UseFishIdx = UseFishIdx,
                ObsFishAgeComps = ObsFishAgeComps,
                ObsFishLenComps = ObsFishLenComps,
                ISS_FishAgeComps = ISS_FishAgeComps,
                ISS_FishLenComps = ISS_FishLenComps,
                UseFishAgeComps = UseFishAgeComps,
                UseFishLenComps = UseFishLenComps,
                ObsSrvIdx = ObsSrvIdx,
                ObsSrvIdx_SE = ObsSrvIdx_SE,
                UseSrvIdx = UseSrvIdx,
                ObsSrvAgeComps = ObsSrvAgeComps,
                ObsSrvLenComps = ObsSrvLenComps,
                ISS_SrvAgeComps = ISS_SrvAgeComps,
                ISS_SrvLenComps = ISS_SrvLenComps,
                UseSrvAgeComps = UseSrvAgeComps,
                UseSrvLenComps = UseSrvLenComps
    ))

}

#' Predict ISS fishery compositions under fishing mortality
#'
#' Uses historical ISS fishery compositions and fishing mortality rates
#' to estimate ISS compositions in the projection year. Compositions are
#' scaled relative to the historical maximum fishing mortality with
#' linear interpolation between the minimum and maximum observed ISS values.
#' If historical values are not available, defaults to the mean or zero.
#'
#' @param ISS_FishComps Array of ISS fishery compositions with dimensions
#'   `[region, year, sex, fleet, sim]`.
#' @param Fmort Array of fishing mortality rates with dimensions
#'   `[region, year, fleet, sim]`.
#' @param y Integer, projection year index for prediction.
#' @param sim Integer, simulation index.
#'
#' @returns Array with predicted ISS values for year `y`.
#' @keywords internal
#'
predict_sim_fish_iss_fmort <- function(ISS_FishComps,
                                       Fmort,
                                       y,
                                       sim
                                       ) {

  # dimensions
  dims <- dim(ISS_FishComps)
  n_regions <- dims[1]
  n_sexes <- dims[3]
  n_fish_fleets <- dims[4]

  # extract temp vars
  tmp_iss <- ISS_FishComps[, 1:(y-1), , , sim, drop = FALSE]
  tmp_fmort <- Fmort[, 1:(y-1), , sim, drop = FALSE]

  # container
  iss_container <- array(0, dim = c(n_regions, length(1:y), n_sexes, n_fish_fleets))
  iss_container[, 1:(y-1), , ] <- ISS_FishComps[, 1:(y-1), , , sim] # fill in values back

  for(r in 1:n_regions) {
    for(s in 1:n_sexes) {
      for(f in 1:n_fish_fleets) {
        # get ISS and Fmort
        iss_vec <- tmp_iss[r, , s, f, ]
        fmort_vec <- tmp_fmort[r, , f, ]
        # remove zeros/NAs
        valid_idx <- which(iss_vec > 0 & !is.na(iss_vec) & !is.na(fmort_vec))
        if(length(valid_idx) > 0) {
          iss_valid <- iss_vec[valid_idx]
          fmort_valid <- fmort_vec[valid_idx]
          # min and max ISS from conditioning period
          min_iss <- min(iss_valid)
          max_iss <- max(iss_valid)
          # max Fmort from conditioning period
          max_fmort_hist <- max(fmort_valid)
          # new Fmort
          fmort_new <- Fmort[r, y, f, sim]
          # scale ISS proportionally to Fmort relative to historical max
          if(max_fmort_hist > 0 && fmort_new >= 0) {
            # linear scaling between min and max ISS
            scaling_factor <- min(fmort_new / max_fmort_hist, 1)  # cap scaling at 1
            new_iss <- min_iss + scaling_factor * (max_iss - min_iss) # linear scaling
          } else {
            # use mean ISS if conditions not met ...
            new_iss <- mean(iss_valid)
          }
          iss_container[r, y, s, f] <- new_iss
        } else {
          iss_container[r, y, s, f] <- 0
        }
      } # end f loop
    } # end s loop
  } # end r loop

  return(iss_container)
}
