# Req = array(0, dim = n_regions)  # equilibrium recruitment by region
# Req_o = R0 * Rec_Prop  # starting guess

# for(iter in 1:1e3) {
#
#   Req_new = array(0, dim = n_regions)
#
#   for(d in 1:n_regions) {  # recruitment destination region
#
#     # ========== Numerator ==========
#     # 4*h^d*χ^d*R0 * Σ_o S^{o,d}(F)*R^o(F)
#     numer_sum = 0
#     for(o in 1:n_regions) {
#       numer_sum = numer_sum + SB_origin_dest[2, o, d] * Req_o[o]
#     }
#
#     numerator = 4 * h[d] * Rec_Prop[d] * R0 * numer_sum
#
#     # ========== Denominator part 1: (1-h^d)*Σ_o S^{o,d}(0)*χ^o*R0 ==========
#     denom1_sum = 0
#     for(o in 1:n_regions) {
#       denom1_sum = denom1_sum + SB_origin_dest[1, o, d] * Rec_Prop[o] * R0
#     }
#     denom_1 = (1 - h[d]) * denom1_sum
#
#     # ========== Denominator part 2: (5*h^d-1)*Σ_o S^{o,d}(F)*R^o(F) ==========
#     denom2_sum = 0
#     for(o in 1:n_regions) {
#       denom2_sum = denom2_sum + SB_origin_dest[2, o, d] * Req_o[o]
#     }
#     denom_2 = (5 * h[d] - 1) * denom2_sum
#
#     # ========== Solve for R^d ==========
#     denominator = denom_1 + denom_2
#     Req_new[d] = numerator / denominator
#
#   } # end destination loop
#
#   # Check convergence
#   conv = sum(abs(Req_new - Req_o))
#   Req_o = Req_new
#
# } # end iteration loop

#' #' Title Get Local FMSY from a Beverton-Holt (Spatial)
#' #'
#' #' @param pars Parameter List
#' #' @param data Data List
#' #' @keywords internal
#' #' @import RTMB
#' local_BH_Fmsy <- function(pars,
#'                           data) {
#'
#'   "c" <- RTMB::ADoverload("c")
#'   "[<-" <- RTMB::ADoverload("[<-")
#'
#'   RTMB::getAll(pars, data) # get parameters and data
#'
#'   n_regions = dim(fish_sel)[1] # number of regions
#'   n_model_ages = dim(fish_sel)[2] # number of model ages
#'   n_ages = n_model_ages * 100 # get number of ages to iterate through for plus group
#'
#'   # set up containers
#'   SB_age = Nspr = array(0, dim = c(2, n_regions, n_regions, n_ages)) # 2 slots in rows, for unfished, and fished at Fmsy
#'   CAA = array(0, c(n_regions, n_regions, n_ages)) # catch at age
#'   Yield_r = array(0, dim = n_regions) # yield by region
#'   SB_unfished_mat = matrix(0, n_regions, n_regions)  # unfished spawning biomass per recruit
#'   SB_fished_mat = matrix(0, n_regions, n_regions) # fished spawning biomass per recruit
#'   Bmsy_r = array(0, dim = n_regions) # BMSY
#'   B0_r = array(0, dim = n_regions) # unfished B0
#'   SPR_r = array(0, dim = n_regions) # spawning potential ratio
#'
#'   # Extend WAA by repeating the last age
#'   WAA_ext <- cbind(WAA, matrix(WAA[, n_model_ages, drop = FALSE], nrow = n_regions, ncol = n_ages - n_model_ages))
#'
#'   # exponentitate reference points to "solve"
#'   Fmsy = exp(log_Fmsy)
#'
#'   # Set up the initial recruits (1 recruit per area)
#'   for(o in 1:n_regions) {
#'     for(d in 1:n_regions) {
#'       if(o == d) Nspr[1,o,d,1] = Nspr[2,o,d,1]  = 1
#'       else Nspr[1,o,d,1] = Nspr[2,o,d,1] = 0
#'     } # end d loop
#'   } # end o loop
#'
#'   # # Loop through, apply movement first, then decrement recruit
#'   # for(j in 2:n_ages) {
#'   #
#'   #   # Get age index to use for demographics
#'   #   if(j <= n_model_ages) age_idx = j - 1
#'   #   else age_idx = n_model_ages
#'   #
#'   #   # move individuals from origin region and move them around
#'   #   for(o in 1:n_regions) {
#'   #
#'   #     # Get temporary values from origin region
#'   #     tmp_unfished = Nspr[1,o,,j-1]
#'   #     tmp_fished = Nspr[2,o,,j-1]
#'   #
#'   #     # Apply movement
#'   #     if(do_recruits_move == 1 || (do_recruits_move == 0 && j > 2)) {
#'   #       tmp_unfished = t(tmp_unfished) %*% Movement[,,age_idx]
#'   #       tmp_fished = t(tmp_fished) %*% Movement[,,age_idx]
#'   #     }
#'   #
#'   #     tmp_F = apply(F_fract_flt * Fmsy * fish_sel[, age_idx, , drop = FALSE], 1, sum) # get F
#'   #     tmp_Z = tmp_F + natmort[, age_idx] # get Z
#'   #
#'   #     # decrement fish and project forward
#'   #     Nspr[1,o,,j] = tmp_unfished * exp(-1 * natmort[, age_idx])
#'   #     Nspr[2,o,,j] = tmp_fished * exp(-1 * tmp_Z)
#'   #
#'   #   } # end o loop
#'   # } # end j loop
#'   #
#'   #
#'   # # Derive spawning biomass per recruit and yield per recruit quantities by origin and destination region
#'   # for(j in 1:n_ages) {
#'   #
#'   #   # Get age index to use for demographics
#'   #   if(j <= n_model_ages) age_idx = j
#'   #   else age_idx = n_model_ages
#'   #
#'   #   tmp_F = apply(F_fract_flt * Fmsy * fish_sel[,age_idx,,drop = FALSE], 1, sum) # temporary fishing mortality at age
#'   #   tmp_Z = tmp_F + natmort[,age_idx] # temporary total mortality at age
#'   #
#'   #   for(o in 1:n_regions) {
#'   #     for(d in 1:n_regions) {
#'   #
#'   #       SB_age[1,o,d,j] = Nspr[1,o,d,j] * WAA[d,age_idx] * MatAA[d,age_idx] * exp(-t_spwn * natmort[d,age_idx]) # unfished
#'   #       SB_age[2,o,d,j] = Nspr[2,o,d,j] * WAA[d,age_idx] * MatAA[d,age_idx] * exp(-t_spwn * tmp_Z[d]) # fished
#'   #
#'   #       # Get catch at age to derive yield per recruit
#'   #       CAA[o,d,j] = Nspr[2,o,d,j] * (tmp_F[d] / tmp_Z[d]) * (1 - exp(-tmp_Z[d])) # Baranov's
#'   #
#'   #     } # end d loop
#'   #   } # end o loop
#'   #
#'   # } # end j loop
#'
#'   # Loop through, apply movement first, then decrement recruit
#'   for(j in 2:n_ages) {
#'     # Get age index to use for demographics
#'     if(j <= n_model_ages) age_idx = j - 1
#'     else age_idx = n_model_ages
#'
#'     # move individuals from origin region and move them around
#'     for(o in 1:n_regions) {
#'       # Get temporary values from origin region
#'       tmp_unfished = Nspr[1,o,,j-1]
#'       tmp_fished = Nspr[2,o,,j-1]
#'
#'       # Apply movement
#'       if(do_recruits_move == 1 || (do_recruits_move == 0 && j > 2)) {
#'         tmp_unfished = t(tmp_unfished) %*% Movement[,,age_idx]
#'         tmp_fished = t(tmp_fished) %*% Movement[,,age_idx]
#'       }
#'
#'       # Get demographics index for spawning biomass (use j-1 for "before mortality")
#'       if(j-1 <= n_model_ages) age_idx_demo = j - 1
#'       else age_idx_demo = n_model_ages
#'
#'       # Calculate F and Z
#'       tmp_F = apply(F_fract_flt * Fmsy * fish_sel[, age_idx, , drop = FALSE], 1, sum)
#'       tmp_Z = tmp_F + natmort[, age_idx]
#'
#'       # compute SBPR here
#'       for(d in 1:n_regions) {
#'         SB_age[1,o,d,j-1] = tmp_unfished[d] * WAA[d,age_idx_demo] * MatAA[d,age_idx_demo] * exp(-t_spwn * natmort[d,age_idx_demo])
#'         SB_age[2,o,d,j-1] = tmp_fished[d] * WAA[d,age_idx_demo] * MatAA[d,age_idx_demo] * exp(-t_spwn * tmp_Z[d])
#'
#'         # Get catch at age for yield per recruit
#'         CAA[o,d,j-1] = tmp_fished[d] * (tmp_F[d] / tmp_Z[d]) * (1 - exp(-tmp_Z[d]))
#'       }
#'       # decrement fish and project forward
#'       Nspr[1,o,,j] = tmp_unfished * exp(-1 * natmort[, age_idx])
#'       Nspr[2,o,,j] = tmp_fished * exp(-1 * tmp_Z)
#'     } # end o loop
#'   } # end j loop
#'
#'   # Handle the last age (n_ages) separately since loop ends
#'   for(o in 1:n_regions) {
#'     tmp_unfished = Nspr[1,o,,n_ages]
#'     tmp_fished = Nspr[2,o,,n_ages]
#'
#'     age_idx_demo = n_model_ages
#'     tmp_F = apply(F_fract_flt * Fmsy * fish_sel[, age_idx_demo, , drop = FALSE], 1, sum)
#'     tmp_Z = tmp_F + natmort[, age_idx_demo]
#'
#'     for(d in 1:n_regions) {
#'       SB_age[1,o,d,n_ages] = tmp_unfished[d] * WAA[d,age_idx_demo] * MatAA[d,age_idx_demo] * exp(-t_spwn * natmort[d,age_idx_demo])
#'       SB_age[2,o,d,n_ages] = tmp_fished[d] * WAA[d,age_idx_demo] * MatAA[d,age_idx_demo] * exp(-t_spwn * tmp_Z[d])
#'       CAA[o,d,n_ages] = tmp_fished[d] * (tmp_F[d] / tmp_Z[d]) * (1 - exp(-tmp_Z[d]))
#'     }
#'   }
#'
#'   # Remove the old "Derive spawning biomass per recruit" loop entirely
#'
#'   # Determine equilibrium recruitment for destination region
#'   # parse out and compute unfished and fished spawning biomass per recruit
#'   for(o in 1:n_regions) {
#'     for(d in 1:n_regions) {
#'       SB_unfished_mat[o, d] = sum(SB_age[1, o, d, ])  # unfished
#'       SB_fished_mat[o, d] = sum(SB_age[2, o, d, ])  # fished at Fmsy
#'     } # end o loop
#'   } # end d loop
#'
#'   A = 4 * h * Rec_Prop * R0 # define first part of the numerator of BH recruitment
#'   B = rep(0, n_regions) # define first part of the denominator of BH recruitment
#'   for(d in 1:n_regions) B[d] = (1 - h[d]) * sum(SB_unfished_mat[,d] * Rec_Prop * R0)
#'   C = 5 * h - 1 # define second part of the denominator for BH recruitment
#'
#'   # define initial guess to solve for equilibrium recruitment from origin region
#'   Req_o = R0 * Rec_Prop
#'
#'   for(nit in 1:newton_steps) {
#'     # compute equilibrium spawning biomass (SSBR * Req) in destination region
#'     x_vec = as.numeric(t(SB_fished_mat) %*% Req_o)  # function of equilibrium recruitment in origin region (what we are solving for)
#'     numer_vec = A * x_vec # compute numerator of BH
#'     denom_vec = B + (C * x_vec) # compute denominator of BH
#'     g_vec = numer_vec / denom_vec # equilibrium recruitment in destination region
#'
#'     # define root and fine Jacobian
#'     iter_vec = Req_o - g_vec # find values of origin recruitment that are consisitent w/ destination recruitment such that pop'n is in equilibrium
#'
#'     # construct Jacobian for root
#'     # we need J = df (iter_vec)/dReq = dReq/dReq (or I) - dg/dReq
#'     # we basically want to know dg / dReq (how does destination equil rec change as orign equil rec change)
#'     # dg / dReq = (dg / dxk) * (dxk / dReq)
#'     # to get (dg / dxk), use quotient rule of (BH recruitment)
#'     # note that dxk / dReq S_2mat * Req = S_2mat
#'     dg_dxk = (A * B) / (denom_vec^2)
#'     dg_dReq = matrix(0, n_regions, n_regions)
#'     for(d in 1:n_regions) dg_dReq[d, ] = dg_dxk[d] * SB_fished_mat[, d]# now compute to see how destination equilibrium rec changes, as origin equil rec changes
#'
#'     # compute jacobian
#'     J = diag(1, n_regions) - dg_dReq
#'     delta = solve(J, iter_vec) # get step to move towards solution
#'     Req_o = Req_o - delta # newton raphson update
#'   }
#'
#'   # get destination reigon yield
#'   for(d in 1:n_regions) {
#'     tmp = 0 # define temp variable
#'     for(o in 1:n_regions) tmp = tmp + sum(CAA[o, d, ] * WAA_ext[d, 1:n_ages]) * Req_o[o] # get yield to destination
#'     Yield_r[d] = tmp
#'   } # end d loop
#'
#'   # get other derived quantities
#'   for(d in 1:n_regions) {
#'     Bmsy_r[d] = sum(SB_fished_mat[,d] * Req_o)
#'     B0_r[d] = sum(SB_unfished_mat[,d] * R0 * Rec_Prop)
#'     SPR_r[d] = Bmsy_r[d] / B0_r[d]
#'   }
#'
#'   # maximize total yield
#'   Yield_total = sum(Yield_r)
#'   obj_fun = -Yield_total
#'
#'   sum_SB_unfished_mat = sum(SB_unfished_mat)
#'
#'   # RTMB::REPORT(eqrec_prop)
#'   RTMB::REPORT(Fmsy)
#'   RTMB::REPORT(Req_o)
#'   RTMB::REPORT(Bmsy_r)
#'   RTMB::REPORT(SPR_r)
#'   RTMB::REPORT(Yield_r)
#'   RTMB::REPORT(Yield_total)
#'   RTMB::REPORT(dg_dReq)
#'   RTMB::REPORT(B0_r)
#'   RTMB::REPORT(SPR_r)
#'   RTMB::REPORT(iter_vec)
#'   RTMB::REPORT(SB_fished_mat)
#'   RTMB::REPORT(SB_unfished_mat)
#'   RTMB::REPORT(sum_SB_unfished_mat)
#'   RTMB::REPORT(Nspr)
#'   RTMB::REPORT(SB_age)
#'
#'   return(obj_fun)
#' }

#' Title Get Local FMSY from a Beverton-Holt (Spatial)
#'
#' @param pars Parameter List
#' @param data Data List
#' @keywords internal
#' @import RTMB
local_BH_Fmsy <- function(pars,
                          data) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data) # get parameters and data

  n_regions = dim(fish_sel)[1] # number of regions
  n_model_ages = dim(fish_sel)[2] # number of model ages
  n_ages = n_model_ages  # get number of ages to iterate through for plus group

  # set up containers
  SB_age = Nspr = array(0, dim = c(2, n_regions, n_regions, n_ages)) # 2 slots in rows, for unfished, and fished at Fmsy
  CAA = array(0, c(n_regions, n_regions, n_ages)) # catch at age
  Yield_r = array(0, dim = n_regions) # yield by region
  SB_unfished_mat = matrix(0, n_regions, n_regions)  # unfished spawning biomass per recruit
  SB_fished_mat = matrix(0, n_regions, n_regions) # fished spawning biomass per recruit
  Bmsy_r = array(0, dim = n_regions) # BMSY
  B0_r = array(0, dim = n_regions) # unfished B0
  SPR_r = array(0, dim = n_regions) # spawning potential ratio

  # Extend WAA by repeating the last age
  WAA_ext <- cbind(WAA, matrix(WAA[, n_model_ages, drop = FALSE], nrow = n_regions, ncol = n_ages - n_model_ages))

  # exponentitate reference points to "solve"
  Fmsy = exp(log_Fmsy)

  # Set up the initial recruits (1 recruit per area)
  for(o in 1:n_regions) {
    for(d in 1:n_regions) {
      if(o == d) Nspr[1,o,d,1] = Nspr[2,o,d,1]  = 1
      else Nspr[1,o,d,1] = Nspr[2,o,d,1] = 0
    } # end d loop
  } # end o loop

  # Loop through, apply movement first, then decrement recruit
  for(j in 2:n_ages) {

    # Get age index to use for demographics
    if(j <= n_model_ages) age_idx = j - 1
    else age_idx = n_model_ages

    # move individuals from origin region and move them around
    for(o in 1:n_regions) {
      # Get temporary values from origin region
      tmp_unfished = Nspr[1,o,,j-1]
      tmp_fished = Nspr[2,o,,j-1]

      # Apply movement
      if(do_recruits_move == 1 || (do_recruits_move == 0 && j > 2)) {
        tmp_unfished = t(tmp_unfished) %*% Movement[,,age_idx]
        tmp_fished = t(tmp_fished) %*% Movement[,,age_idx]
      }

      # Get age index to use for demographics
      if(j <= n_model_ages) age_idx_demo = j - 1
      else age_idx_demo = n_model_ages


      # Calculate F and Z
      tmp_F = apply(F_fract_flt * Fmsy * fish_sel[, age_idx, , drop = FALSE], 1, sum)
      tmp_Z = tmp_F + natmort[, age_idx]

      # compute SBPR here
      for(d in 1:n_regions) {
        SB_age[1,o,d,j-1] = tmp_unfished[d] * WAA[d,age_idx_demo] * MatAA[d,age_idx_demo] * exp(-t_spwn * natmort[d,age_idx_demo])
        SB_age[2,o,d,j-1] = tmp_fished[d] * WAA[d,age_idx_demo] * MatAA[d,age_idx_demo] * exp(-t_spwn * tmp_Z[d])

        # Get catch at age for yield per recruit
        CAA[o,d,j-1] = tmp_fished[d] * (tmp_F[d] / tmp_Z[d]) * (1 - exp(-tmp_Z[d]))
      }
      # decrement fish and project forward
      Nspr[1,o,,j] = tmp_unfished * exp(-1 * natmort[, age_idx])
      Nspr[2,o,,j] = tmp_fished * exp(-1 * tmp_Z)
    } # end o loop
  } # end j loop

  # Get movement matrices and survival for plus group
  M_penult <- Movement[,, n_model_ages - 1]  # movement for age n_ages-1
  M_plus <- Movement[,, n_model_ages]  # movement for plus group

  age_idx_demo = n_model_ages
  tmp_F_plus = apply(F_fract_flt * Fmsy * fish_sel[, age_idx_demo, , drop = FALSE], 1, sum)
  tmp_Z_plus = tmp_F_plus + natmort[, age_idx_demo]
  tmp_F_penult = apply(F_fract_flt * Fmsy * fish_sel[, n_model_ages - 1, , drop = FALSE], 1, sum)
  tmp_Z_penult = tmp_F_penult + natmort[, n_model_ages - 1]

  s_penult_unfished = exp(-natmort[, n_model_ages - 1])  # survival of age n_ages-1 (unfished)
  s_plus_unfished = exp(-natmort[, n_model_ages])  # survival in plus group (unfished)
  s_penult_fished = exp(-tmp_Z_penult)  # survival of age n_ages-1 (fished)
  s_plus_fished = exp(-tmp_Z_plus)  # survival in plus group (fished)

  I_mat = diag(n_regions)

  # Loop over origin regions
  for(o in 1:n_regions) {
    # UNFISHED: Get age n_ages-1 abundance
    N_penult_unfished = Nspr[1, o, , n_ages - 1]

    # Apply movement to penultimate age, then survival
    N_penult_moved_unfished = as.numeric(t(M_penult) %*% N_penult_unfished)
    source_unfished = N_penult_moved_unfished * s_penult_unfished

    # Transition matrix for plus group (already moved fish that survive)
    T_mat_unfished = diag(s_plus_unfished, n_regions) %*% t(M_plus)

    # Solve for equilibrium plus group abundance
    N_plus_equil_unfished = solve(I_mat - T_mat_unfished, source_unfished)
    Nspr[1, o, , n_ages] = N_plus_equil_unfished

    # FISHED: Same process for fished population
    N_penult_fished = Nspr[2, o, , n_ages - 1]
    N_penult_moved_fished = as.numeric(t(M_penult) %*% N_penult_fished)
    source_fished = N_penult_moved_fished * s_penult_fished

    T_mat_fished = diag(s_plus_fished, n_regions) %*% t(M_plus)
    N_plus_equil_fished = solve(I_mat - T_mat_fished, source_fished)
    Nspr[2, o, , n_ages] = N_plus_equil_fished

    # Calculate spawning biomass and catch for plus group
    for(d in 1:n_regions) {
      SB_age[1, o, d, n_ages] = N_plus_equil_unfished[d] * WAA[d, age_idx_demo] *
        MatAA[d, age_idx_demo] * exp(-t_spwn * natmort[d, age_idx_demo])

      SB_age[2, o, d, n_ages] = N_plus_equil_fished[d] * WAA[d, age_idx_demo] *
        MatAA[d, age_idx_demo] * exp(-t_spwn * tmp_Z_plus[d])

      CAA[o, d, n_ages] = N_plus_equil_fished[d] * (tmp_F_plus[d] / tmp_Z_plus[d]) * (1 - exp(-tmp_Z_plus[d]))
    }
  } # end o loop

  # Determine equilibrium recruitment for destination region
  # parse out and compute unfished and fished spawning biomass per recruit
  for(o in 1:n_regions) {
    for(d in 1:n_regions) {
      SB_unfished_mat[o, d] = sum(SB_age[1, o, d, ])  # unfished
      SB_fished_mat[o, d] = sum(SB_age[2, o, d, ])  # fished at Fmsy
    } # end o loop
  } # end d loop

  A = 4 * h * Rec_Prop * R0 # define first part of the numerator of BH recruitment
  B = rep(0, n_regions) # define first part of the denominator of BH recruitment
  for(d in 1:n_regions) B[d] = (1 - h[d]) * sum(SB_unfished_mat[,d] * Rec_Prop * R0)
  C = 5 * h - 1 # define second part of the denominator for BH recruitment

  # define initial guess to solve for equilibrium recruitment from origin region
  Req_o = R0 * Rec_Prop

  for(nit in 1:newton_steps) {
    # compute equilibrium spawning biomass (SSBR * Req) in destination region
    x_vec = as.numeric(t(SB_fished_mat) %*% Req_o)  # function of equilibrium recruitment in origin region (what we are solving for)
    numer_vec = A * x_vec # compute numerator of BH
    denom_vec = B + (C * x_vec) # compute denominator of BH
    g_vec = numer_vec / denom_vec # equilibrium recruitment in destination region

    # define root and fine Jacobian
    iter_vec = Req_o - g_vec # find values of origin recruitment that are consisitent w/ destination recruitment such that pop'n is in equilibrium

    # construct Jacobian for root
    # we need J = df (iter_vec)/dReq = dReq/dReq (or I) - dg/dReq
    # we basically want to know dg / dReq (how does destination equil rec change as orign equil rec change)
    # dg / dReq = (dg / dxk) * (dxk / dReq)
    # to get (dg / dxk), use quotient rule of (BH recruitment)
    # note that dxk / dReq S_2mat * Req = S_2mat
    dg_dxk = (A * B) / (denom_vec^2)
    dg_dReq = matrix(0, n_regions, n_regions)
    for(d in 1:n_regions) dg_dReq[d, ] = dg_dxk[d] * SB_fished_mat[, d]# now compute to see how destination equilibrium rec changes, as origin equil rec changes

    # compute jacobian
    J = diag(1, n_regions) - dg_dReq
    delta = solve(J, iter_vec) # get step to move towards solution
    Req_o = Req_o - delta # newton raphson update
  }

  # get destination reigon yield
  for(d in 1:n_regions) {
    tmp = 0 # define temp variable
    for(o in 1:n_regions) tmp = tmp + sum(CAA[o, d, ] * WAA_ext[d, 1:n_ages]) * Req_o[o] # get yield to destination
    Yield_r[d] = tmp
  } # end d loop

  # get other derived quantities
  for(d in 1:n_regions) {
    Bmsy_r[d] = sum(SB_fished_mat[,d] * Req_o)
    B0_r[d] = sum(SB_unfished_mat[,d] * R0 * Rec_Prop)
    SPR_r[d] = Bmsy_r[d] / B0_r[d]
  }

  # maximize total yield
  Yield_total = sum(Yield_r)
  obj_fun = -Yield_total

  sum_SB_unfished_mat = sum(SB_unfished_mat)

  # RTMB::REPORT(eqrec_prop)
  RTMB::REPORT(Fmsy)
  RTMB::REPORT(Req_o)
  RTMB::REPORT(Bmsy_r)
  RTMB::REPORT(SPR_r)
  RTMB::REPORT(Yield_r)
  RTMB::REPORT(Yield_total)
  RTMB::REPORT(dg_dReq)
  RTMB::REPORT(B0_r)
  RTMB::REPORT(SPR_r)
  RTMB::REPORT(iter_vec)
  RTMB::REPORT(SB_fished_mat)
  RTMB::REPORT(SB_unfished_mat)
  RTMB::REPORT(sum_SB_unfished_mat)
  RTMB::REPORT(Nspr)
  RTMB::REPORT(SB_age)

  return(obj_fun)
}

library(SPoRC)
devtools::load_all(here("R"))
data <- three_rg_sable_data
rep <- readRDS("/Users/matthewcheng/Desktop/PostDoc/Spatial Assessments and Sablefish/SPoRC/dev/dev_output/3_Region_Model_Sablefish/rep.RDS")
# data <- mlt_rg_sable_data
# rep <- mlt_rg_sable_rep


data_list <- list() # set up data list
# determine years to average over demogrphaics
n_yrs <- length(data$years)
n_avg_yrs = 1
avg_yrs <- (n_yrs - n_avg_yrs + 1):n_yrs
t_spwn = 0

f_ref_pt <- vector() # set up storage
b_ref_pt <- vector() # set up storage
virgin_b_ref_pt <- vector() # set up storage

# determine years to average over demogrphaics
n_yrs <- length(data$years)
avg_yrs <- (n_yrs - n_avg_yrs + 1):n_yrs


# Extract out relevant elements for a given region
n_ages <- length(data$ages) # number of ages to iterate through
n_years <- length(data$years) # number of years
n_regions <- 3 # number of regions
data_list$t_spwn <- t_spwn # specified mortality time up until spawning

# Fleet fraction F
data_list$F_fract_flt <- rep$Fmort[,n_years,,drop = FALSE] / apply(rep$Fmort[,n_years,,drop = FALSE], 1, sum) # get fleet F fraction to derive population level selectivity

# fishery selectivity
fish_sel_avg <- apply(rep$fish_sel[,avg_yrs,,1,,drop = FALSE], c(1, 3, 4, 5), mean)
data_list$fish_sel <- array(fish_sel_avg, dim = c(n_regions, n_ages, data$n_fish_fleets)) # get female selectivity for all fleets

# natural mortality
natmort_avg <- apply(rep$natmort[,avg_yrs,,1,drop = FALSE], c(1, 3, 4), mean)
data_list$natmort <- array(natmort_avg, dim = c(n_regions, n_ages)) # get female natural mortality

# weight at age
WAA_avg <- apply(data$WAA[,avg_yrs,,1,drop = FALSE], c(1, 3, 4), mean)
data_list$WAA <- array(WAA_avg, dim = c(n_regions, n_ages)) # weight at age for females
data_list$WAA[1,] <- data_list$WAA[1,]

# maturity at age
MatAA_avg <- apply(data$MatAA[,avg_yrs,,1,drop = FALSE], c(1, 3, 4), mean)
data_list$MatAA <- array(MatAA_avg, dim = c(n_regions, n_ages)) # maturity at age for females

# Movement
Movement_avg <- apply(rep$Movement[,,avg_yrs,,1,drop = FALSE], c(1,2,4,5), mean)
Movement_avg <- sweep(Movement_avg, c(1, 3, 4), apply(Movement_avg, c(1, 3, 4), sum), "/")
data_list$Movement <- array(Movement_avg, dim = c(n_regions, n_regions, n_ages)) # Movement
# data_list$Movement <- array(replicate(30, diag(1,2)), dim = c(n_regions, n_regions, n_ages)) # Movement
# data_list$Movement <- array(replicate(30, diag(1/3,2)), dim = c(n_regions, n_regions, n_ages)) # Movement

# Recruitment options
data_list$do_recruits_move <- data$do_recruits_move # whether recruits move
data_list$Rec_Prop <- rep$Rec_trans_prop[] / sum(rep$Rec_trans_prop[]) # recruitment proportions

data_list$R0 <- rep$R0  # unfished recruitment
data_list$h <- rep$h_trans
data_list$newton_steps <- 10

par_list <- list() # set up parameter list
par_list$log_Fmsy <- log(rep(0.1, 3)) # Fmsy starting value


# Make adfun object
tmp_obj <- RTMB::MakeADFun(cmb(local_BH_Fmsy, data_list), parameters = par_list, map = NULL, silent = TRUE)
tmp_obj$optim <- stats::nlminb(tmp_obj$par, tmp_obj$fn, tmp_obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
tmp_obj$rep <- tmp_obj$report(tmp_obj$env$last.par.best) # get report
tmp_obj$rep$Fmsy

plot(a)
plot(tmp_obj$rep$SB_age[1,1,1,])

sdreport(tmp_obj)
tmp_obj$rep$Fmsy
tmp_obj$rep$Bmsy_r
tmp_obj$rep$Bmsy_r

# newton_steps <- 6
# log_Fmsy_vals <- rep(log(0.1), 2) # starting
# RTMB::getAll(par_list, data_list) # get parameters and data
# for(i in 1:5) {
#   grad <- numDeriv::jacobian(local_BH_Fmsy_grad, log_Fmsy_vals) # get gradient
#   hess <- numDeriv::hessian(local_BH_Fmsy_grad, log_Fmsy_vals) # get hessian
#   log_Fmsy_vals <- log_Fmsy_vals - as.vector(solve(hess, as.vector(grad)))
#   print(i)
# }

log_Fmsy_vals <- rep(log(0.1), 2) # starting
tmp_obj1 <- RTMB::MakeADFun(cmb(local_BH_Fmsy, data_list), parameters = par_list, map = NULL, silent = TRUE)
for(i in 1:4) {
  grad <- tmp_obj1$gr(x = log_Fmsy_vals)
  hess <- tmp_obj1$he(x = log_Fmsy_vals)
  log_Fmsy_vals <- log_Fmsy_vals - as.vector(solve(hess, as.vector(grad)))
}


f <- expand.grid(data.frame(
  r1 = seq(0.0,1,.05),
  r2 = seq(0.0,1,.05),
  r3 = seq(0,1,.05)
  # r4 = seq(0,1,.1),
  # r5 = seq(0,1,.1)
))

f$y <- NA


library(pbapply)

f$y <- pbapply(f[,1:3], 1, function(row) {
  abs(obj$fn(log(row)))
}, cl = 5)


ggplot(f, aes(x = r1, y = r2, fill = y, z = y)) +
  geom_contour_filled(aes(fill = after_stat(level)), bins = 30, alpha = 1) +
  # geom_point(aes(x = tmp_obj$rep$Fmsy[1],
  #                y = tmp_obj$rep$Fmsy[2]),
  #            size = 4, shape = 21, fill = "white")  +
  coord_cartesian(xlim = c(0.0,0.25), ylim = c(0.0,0.25)) +
  scale_fill_viridis_d(option = 'magma') +
  theme_sablefish()
  # labs(x = 'region 1', y = "region 5")


