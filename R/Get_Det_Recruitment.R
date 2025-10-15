#' Get Deterministic Recruitment
#'
#' Computes deterministic recruitment for each region based on either
#' a mean recruitment model or a Beverton–Holt stock–recruitment relationship.
#'
#' @param recruitment_model Integer flag specifying the recruitment model:
#'   \itemize{
#'     \item `0` = mean recruitment
#'     \item `1` = Beverton–Holt recruitment with steepness
#'   }
#' @param recruitment_dd Integer flag specifying the scale of density dependence:
#'   \itemize{
#'     \item `0` = local (region-specific)
#'     \item `1` = global (shared across regions)
#'   }
#' @param y Current model year (used for SSB lag indexing)
#' @param rec_lag Recruitment lag (number of years between spawning and recruitment)
#' @param R0 Virgin or mean recruitment (global scalar)
#' @param Rec_Prop Vector of recruitment proportions by region (used to allocate global `R0` under local density dependence)
#' @param h Vector of Beverton–Holt steepness values by region
#' @param n_regions Number of spatial regions
#' @param n_ages Number of modeled age classes
#' @param WAA Matrix of weight-at-age by region and age
#' @param MatAA Matrix of maturity-at-age by region and age
#' @param natmort Matrix or vector of natural mortality by region and age
#' @param SSB_vals Matrix of spawning stock biomass (SSB) by region and year
#' @param Movement 3D array of movement probabilities between regions by age (`[origin, destination, age]` or alternatively, `[n_regions, n_regions, age]`)
#' @param do_recruits_move Logical or integer flag (0/1) indicating whether recruits move during their first year
#' @param t_spawn Fraction of the year at which spawning occurs (used for survival to spawning)
#' @param sex_ratio_f Vector of female proportions by region (used to scale initial recruits)
#'
#' @details
#' The function returns region-specific deterministic recruitment estimates
#' based on the chosen recruitment model and density dependence structure.
#'
#' When `recruitment_model = 0`, recruitment is fixed at mean values (`R0 * Rec_Prop`).
#' When `recruitment_model = 1`, Beverton–Holt recruitment is applied using:
#' \deqn{R = \frac{4hR_0SSB}{(1 - h)S_0 + (5h - 1)SSB}}
#' where `S_0` is unfished spawning biomass per recruit, computed separately for each
#' region (local) or summed across all regions (global).
#'
#' @return A numeric vector of length `n_regions` containing deterministic
#' recruitment values for each region.
#'
#' @keywords internal
Get_Det_Recruitment <- function(recruitment_model,
                                recruitment_dd,
                                y,
                                rec_lag,
                                R0,
                                Rec_Prop,
                                h,
                                n_regions,
                                n_ages,
                                WAA,
                                MatAA,
                                natmort,
                                SSB_vals,
                                Movement,
                                do_recruits_move,
                                t_spawn,
                                sex_ratio_f
                                ) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  if(recruitment_model == 0) rec = R0 * Rec_Prop # mean recruitment apportioned across n_regions

  # Beverton-Holt
  if(recruitment_model == 1) {

    # Storage for recruitment and S0
    rec = S0 = rep(0, n_regions)
    # Restructure some stuff
    tmp_natmort <- array(natmort, dim = c(n_regions, n_ages))
    tmp_WAA <- array(WAA, dim = c(n_regions, n_ages))
    tmp_MatAA <- array(MatAA, dim = c(n_regions, n_ages))

    # local density dependence
    if(recruitment_dd == 0) {
      # Calculate unexploited naa per recruit by origin area and destination area
      SB_age = Nspr = array(0, dim = c(n_regions, n_regions, n_ages_ext))
      SB_unfished_mat = array(0, c(n_regions, n_regions))

      # Set up the initial recruits (1 recruit per area)
      for(o in 1:n_regions) {
        for(d in 1:n_regions) {
          if(o == d) Nspr[o,d,1] = 1 * sex_ratio_f[o]
          else Nspr[o,d,1] = 0
        } # end d loop
      } # end o loop

      # Loop through, apply movement first, then decrement recruit
      for(j in 2:n_ages) {
        # move individuals from origin region and move them around
        for(o in 1:n_regions) {
          # Get temporary values from origin region
          tmp_unfished = Nspr[o,,j-1]
          # Apply movement
          if(do_recruits_move == 1 || (do_recruits_move == 0 && j > 2)) {
            tmp_unfished = t(tmp_unfished) %*% Movement[,,j - 1]
          }
          # compute pr quantities here (before mortality)
          for(d in 1:n_regions) SB_age[o,d,j-1] = tmp_unfished[d] * tmp_WAA[d,j - 1] * tmp_MatAA[d,j - 1] * exp(-t_spawn * tmp_natmort[d,j - 1])
          # decrement fish and project forward
          Nspr[o,,j] = tmp_unfished * exp(-1 * tmp_natmort[, j - 1])
        } # end o loop
      } # end j loop

      # Get movement for penultimate and plus group
      M_penult = Movement[,, n_ages - 1]  # movement for age n_ages-1
      M_plus = Movement[,, n_ages]  # movement for plus group

      # compute survival
      s_penult_unfished = exp(-tmp_natmort[, n_ages - 1])  # survival of age n_ages-1 (unfished)
      s_plus_unfished = exp(-tmp_natmort[, n_ages])  # survival in plus group (unfished)
      I_mat = diag(n_regions) # identity matrix to solve

      # Loop over origin regions
      for(o in 1:n_regions) {
        # unfished pr
        N_penult_unfished = Nspr[ o, , n_ages - 1] # get penultimate unfished
        # Apply movement to penultimate age, then survival
        source_unfished = as.numeric(t(M_penult) %*% N_penult_unfished) * s_penult_unfished
        T_mat_unfished = diag(s_plus_unfished, n_regions) %*% t(M_plus) # unfished transition matrix
        # Solve for equilibrium plus group abundance
        N_plus_equil_unfished = solve(I_mat - T_mat_unfished, source_unfished) # (I-T)^-1 %*% source_unfished
        Nspr[o,,n_ages] = N_plus_equil_unfished

        # Calculate spawning biomass and catch for plus group
        for(d in 1:n_regions) {
          SB_age[ o, d, n_ages] = N_plus_equil_unfished[d] * tmp_WAA[d, n_ages] *
            tmp_MatAA[d, n_ages] * exp(-t_spawn * tmp_natmort[d, n_ages])
        }
      } # end o loop

      # Remove the old spawning biomass calculation loop entirely
      # parse out and compute unfished spawning biomass per recruit
      for(o in 1:n_regions) for(d in 1:n_regions) SB_unfished_mat[o, d] = sum(SB_age[o, d, ])
      for(d in 1:n_regions) S0[d] = sum(SB_unfished_mat[,d] * Rec_Prop * R0)
    }

    # global density dependence
    if(recruitment_dd == 1) {
      n_ages <- n_ages
      SB_age = Nspr = array(0, dim = c(n_regions, n_ages))

      # Set up the initial recruits
      Nspr[,1] = Rec_Prop * sex_ratio_f

      # Loop through, apply movement first, then decrement recruit
      for(j in 2:n_ages) {
        # Get temporary values
        tmp_unfished = Nspr[,j-1]
        # Apply movement
        if(do_recruits_move == 1 || (do_recruits_move == 0 && j > 2)) {
          tmp_unfished = t(tmp_unfished) %*% Movement[,,j - 1]
        }
        # compute pr quantities here (before mortality)
        for(d in 1:n_regions) {
          SB_age[d,j-1] = tmp_unfished[d] * WAA[d,j - 1] * MatAA[d,j - 1] * exp(-t_spawn * natmort[d,j - 1])
        }
        # decrement fish and project forward
        Nspr[,j] = tmp_unfished * exp(-1 * natmort[, j- 1])
      } # end j loop

      # Get movement for penultimate and plus group
      M_penult = Movement[,, n_ages - 1]  # movement for age n_ages-1
      M_plus = Movement[,, n_ages]  # movement for plus group
      # compute survival
      s_penult_unfished = exp(-natmort[, n_ages - 1])  # survival of age n_ages-1 (unfished)
      s_plus_unfished = exp(-natmort[, n_ages])  # survival in plus group (unfished)
      I_mat = diag(n_regions) # identity matrix to solve

      # Loop over origin regions
      # unfished pr
      N_penult_unfished = Nspr[, n_ages - 1] # get penultimate unfished
      # Apply movement to penultimate age, then survival
      source_unfished = as.numeric(t(M_penult) %*% N_penult_unfished) * s_penult_unfished
      T_mat_unfished = diag(s_plus_unfished, n_regions) %*% t(M_plus) # unfished transition matrix
      # Solve for equilibrium plus group abundance
      N_plus_equil_unfished = solve(I_mat - T_mat_unfished, source_unfished) # (I-T)^-1 %*% source_unfished
      Nspr[,n_ages] = N_plus_equil_unfished

      # Calculate spawning biomass and catch for plus group
      for(d in 1:n_regions) {
        SB_age[d, n_ages] = N_plus_equil_unfished[d] * WAA[d, n_ages] * MatAA[d, n_ages] * exp(-t_spawn * natmort[d, n_ages])
      }
      S0[] = sum(SB_age) * R0 * Rec_Prop # get virgin biomass
    }

    # get SSB to use to predict recruitment
    if(y <= rec_lag) SSB = S0 else SSB = SSB_vals[,y-rec_lag]

    # Get recruitment based on SSB and R0
    for(r in 1:n_regions) {
      # Local Density Dependence
      if(recruitment_dd == 0) {
        local_R0 = R0 * Rec_Prop[r] # get local R0 based on recruitment proportions
        rec[r] = (4*h[r]*local_R0*SSB[r]) / ((1-h[r])*S0[r] + (5*h[r]-1)*SSB[r]) # get local beverton holt
      }
      # Global Density Dependence
      if(recruitment_dd == 1) {
        rec[r] = (4*h[r]*R0*sum(SSB)) / ((1-h[r])*sum(S0) + (5*h[r]-1)*sum(SSB)) * Rec_Prop[r] # get global beverton holt and then apportion to different regions
      }

    } # end r loop
  } # end Beverton-Holt

    return(rec)
}
