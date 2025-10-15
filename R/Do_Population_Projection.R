#' Do Population Projections
#'
#' @param n_proj_yrs Number of projection years
#' @param n_regions Number of regions
#' @param n_ages Number of ages
#' @param n_sexes Number of sexes
#' @param sexratio Array of recruitment sex ratio (n_regions, n_proj_yrs, n_sexes)
#' @param n_fish_fleets Number of fishery fleets
#' @param do_recruits_move Whether recruits move (0 == don't move, 1 == move)
#' @param recruitment Recruitment matrix dimensioned by n_regions, and n_yrs that we want to summarize across, or condition our projection on
#' @param terminal_NAA Terminal Numbers at Age dimensioned by n_regions, n_ages, n_sexes
#' @param terminal_NAA0 Terminal Unfished Numbers at Age dimensioned by n_regions, n_ages, n_sexes
#' @param terminal_F Terminal fishing mortality rate, dimensioned by n_regions, n_fish_fleets
#' @param natmort Natural mortality, dimensioned by n_regions, n_proj_yrs, n_ages, n_sexes
#' @param WAA Weight at age, dimensioned by n_regions, n_proj_yrs, n_ages, n_sexes
#' @param WAA_fish Weight at age for the fishery, dimensioned by n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets
#' @param MatAA Maturity at age, dimensioned by n_regions, n_proj_yrs, n_ages, n_sexes
#' @param fish_sel Fishery selectivity, dimensioned by n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets
#' @param Movement Movement, dimensioned by n_regions, n_regions, n_proj_yrs, n_ages, n_sexes
#' @param f_ref_pt Fishing mortality reference point dimensioned by n_regions and n_proj_yrs
#' @param b_ref_pt Biological reference point dimensioned by n_regions and n_proj_yrs
#' @param HCR_function Function describing a harvest control rule. The function should always have the following arguments: x, which represents SSB, frp, which takes inputs of fishery reference points, and brp, which takes inputs of biological reference points. Any additional arguments should be specified with defaults or hard coded / fixed within the function.
#' @param recruitment_opt Recruitment simulation option, where options are "inv_gauss", which simulates future recruitment based on the the recruitment values supplied using an inverse gaussian distribution, "mean_rec", which takes the mean of the recruitment values supplied for a given region, and "zero", which assumes that future recruitment does not occur
#' @param fmort_opt Fishing Mortality option, which includes "HCR", which modifies the F reference point using a user supplied HCR_function, or "Input", which uses projected F values supplied by the user.
#' @param t_spawn Fraction time of spawning used to compute projected SSB
#' @param bh_rec_opt A list object containing the following arguments:
#' \describe{
#' \item{recruitment_dd}{A value (0 or 1) indicating global (1) or local density dependence (0). In the case of a single region model, either local or global will give the same results}
#' \item{rec_lag}{A value indicating the number of years lagged that a given year's SSB produces recruits}
#' \item{R0}{The virgin recruitment parameter}
#' \item{Rec_Prop}{Recruitment apportionment values. In a single region model, this should be set at a value of 1. Dimensioned by n_regions}
#' \item{h}{Steepness values for the stock recruitment curve. Dimensioned by n_regions}
#' \item{WAA}{A weight-at-age array dimensioned by n_regions, n_ages, and n_sexes, where the reference year should utilize values from the first year}
#' \item{MatAA}{A maturity at age array dimensioned by n_regions, n_ages, and n_sexes, where the reference year should utilize values from the first year}
#' \item{natmort}{A natural mortality at age array dimensioned by n_regions, n_ages, and n_sexes, where the reference year should utilize values from the first year}
#' \item{SSB}{All SSB values estimated from a given model, dimensioned by n_regions and n_yrs}
#' }
#'
#' @returns A list containing projected F, catch, SSB (and dynamic unfished), and Numbers at Age (and dynamic unfished). (Objects are generally dimensioned in the following order: n_regions, n_yrs, n_ages, n_sexes, n_fleets)
#' @export Do_Population_Projection
#' @family Reference Points and Projections
Do_Population_Projection <- function(n_proj_yrs = 2,
                                     n_regions,
                                     n_ages,
                                     n_sexes,
                                     sexratio,
                                     n_fish_fleets,
                                     do_recruits_move = 0,
                                     recruitment,
                                     terminal_NAA,
                                     terminal_NAA0,
                                     terminal_F,
                                     natmort,
                                     WAA,
                                     WAA_fish,
                                     MatAA,
                                     fish_sel,
                                     Movement,
                                     f_ref_pt = NULL,
                                     b_ref_pt = NULL,
                                     HCR_function = NULL,
                                     recruitment_opt = "inv_gauss",
                                     fmort_opt = 'HCR',
                                     t_spawn,
                                     bh_rec_opt = NULL
                                     ) {


# Error Checking ----------------------------------------------------------

  if(!recruitment_opt %in% c("inv_gauss", "mean_rec", "zero", "bh_rec")) stop("Recruitment options are not specified correctly! Should be inv_gauss, mean_rec, zero, or bh_rec")
  if(!fmort_opt %in% c("HCR", "Input")) stop("Fishing Mortality options are not specified correctly! Should be HCR or Input")
  if(recruitment_opt == "bh_rec") {
    required_fields <- c("recruitment_dd", "rec_lag", "R0", "h", "Rec_Prop", "WAA", "MatAA", "natmort", "SSB")
    diff <- setdiff(required_fields, names(bh_rec_opt)) # find difference
    if(length(diff) > 0) stop(paste("bh_rec_opt is missing the following required fields:", paste(diff)))
  }

# Define Containers -------------------------------------------------------
  fratio <- terminal_F / apply(terminal_F, 1, sum) # get fishing mortality ratio among fleets
  proj_NAA <- array(0, dim = c(n_regions, n_proj_yrs + 1, n_ages, n_sexes))
  proj_NAA0 <- array(0, dim = c(n_regions, n_proj_yrs + 1, n_ages, n_sexes))
  proj_ZAA <- array(0, dim = c(n_regions, n_proj_yrs + 1, n_ages, n_sexes))
  proj_FAA <- array(0, dim = c(n_regions, n_proj_yrs + 1, n_ages, n_sexes, n_fish_fleets))
  proj_CAA <- array(0, dim = c(n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets))
  proj_Catch <- array(0, dim = c(n_regions, n_proj_yrs, n_fish_fleets))
  proj_SSB <- array(0, dim = c(n_regions, n_proj_yrs))
  proj_Dynamic_SSB0 <- array(0, dim = c(n_regions, n_proj_yrs))
  proj_F <- array(0, dim = c(n_regions, n_proj_yrs + 1))

# Start Projection --------------------------------------------------------
  # Input terminal year assessment at age
  proj_NAA[,1,,] <- terminal_NAA
  proj_NAA0[,1,,] <- terminal_NAA0

  for(y in 1:n_proj_yrs) {

# Construct Mortality Processes -------------------------------------------

    # use terminal F in the first year (subsequent years use F derived from reference points and HCR)
    if(y == 1) proj_F[,y] <- rowSums(terminal_F)

    for(a in 1:n_ages) {
      for(s in 1:n_sexes) {
        for(f in 1:n_fish_fleets) {
        # get fishing mortality at age
        proj_FAA[,y,a,s,f] <- proj_F[,y] * fratio[,f] * fish_sel[,y,a,s,f]
      } # end f loop

        # Get Total Mortality at Age
        proj_ZAA[,y,a,s] <- natmort[,y,a,s] + apply(proj_FAA[,y,a,s,,drop = FALSE], c(1:4), sum) # M and sum F across fleets

    } # end s loop
  } # end a loop

# Recruitment Processes ---------------------------------------------------
    if(y > 1) {

      # Get recruitment
      tmp_rec <- switch(recruitment_opt,

                        "inv_gauss" = { # if inverse gaussian
                          sapply(1:n_regions, function(r) rinvgauss_rec(1, recruitment[r,]))
                        },

                        "mean_rec" = { # if mean recruitment
                          sapply(1:n_regions, function(r) mean(recruitment[r,]))
                        },

                        "zero" = { # if zero recruitment
                          rep(0, n_regions)
                        },

                        "bh_rec" = { # if beverton holt recruitment
                          Get_Det_Recruitment(recruitment_model = 1,
                                              recruitment_dd = bh_rec_opt$recruitment_dd,
                                              y = y + dim(bh_rec_opt$SSB)[2],
                                              rec_lag = bh_rec_opt$rec_lag,
                                              R0 = bh_rec_opt$R0,
                                              Rec_Prop = bh_rec_opt$Rec_Prop,
                                              h = bh_rec_opt$h,
                                              n_regions = n_regions,
                                              n_ages = n_ages,
                                              WAA = bh_rec_opt$WAA,
                                              MatAA = bh_rec_opt$MatAA,
                                              natmort = bh_rec_opt$natmort,
                                              SSB_vals = cbind(bh_rec_opt$SSB, proj_SSB),
                                              Movement = bh_rec_opt$Movement,
                                              do_recruits_move = bh_rec_opt$do_recruits_move,
                                              t_spawn = bh_rec_opt$t_spawn,
                                              sex_ratio_f = bh_rec_opt$sex_ratio_f
                                              )
                        }
                        )

      # Apply recruitment to projected NAA
      for(r in 1:n_regions) {
        if(n_sexes == 2) tmp <- tmp_rec[r] * sexratio[r,y,]
        if(n_sexes == 1) tmp <- tmp_rec[r]
        proj_NAA[r,y,1,] <- proj_NAA0[r,y,1,]  <- tmp
      } # end r loop

    }

# Movement Processes ------------------------------------------------------
    # Only apply movement if more than 1 region, or if y > 1 (because terminal NAA already has movement applied)
    if(n_regions > 1 && y > 1) {
      # Recruits don't move
      if(do_recruits_move == 0) {
        # Apply movement after ageing processes - start movement at age 2
        for(a in 2:n_ages) for(s in 1:n_sexes) proj_NAA[,y,a,s] = t(proj_NAA[,y,a,s]) %*% Movement[,,y,a,s] # fished
        for(a in 2:n_ages) for(s in 1:n_sexes) proj_NAA0[,y,a,s] = t(proj_NAA0[,y,a,s]) %*% Movement[,,y,a,s] # unfished
      } # end if recruits don't move
      # Recruits move here
      if(do_recruits_move == 1) {
        for(a in 1:n_ages) for(s in 1:n_sexes) proj_NAA[,y,a,s] = t(proj_NAA[,y,a,s]) %*% Movement[,,y,a,s] # fished
        for(a in 1:n_ages) for(s in 1:n_sexes) proj_NAA0[,y,a,s] = t(proj_NAA0[,y,a,s]) %*% Movement[,,y,a,s] # unfished
      }
    } # only compute if spatial

# Mortality and Ageing ----------------------------------------------------
    proj_NAA[,y+1,2:n_ages,] = proj_NAA[,y,1:(n_ages-1),] * exp(-proj_ZAA[,y,1:(n_ages-1),]) # Exponential mortality for individuals not in plus group (fished)
    proj_NAA[,y+1,n_ages,] = proj_NAA[,y+1,n_ages,] + proj_NAA[,y,n_ages,] * exp(-proj_ZAA[,y,n_ages,]) # Accumulate plus group (fished)
    proj_NAA0[,y+1,2:n_ages,] = proj_NAA0[,y,1:(n_ages-1),] * exp(-natmort[,y,1:(n_ages-1),]) # Exponential mortality for individuals not in plus group (unfished)
    proj_NAA0[,y+1,n_ages,] = proj_NAA0[,y+1,n_ages,] + proj_NAA0[,y,n_ages,] * exp(-natmort[,y,n_ages,]) # Accumulate plus group (unfished)

# Derive Biomass ----------------------------------------------------------
    proj_SSB[,y] = apply(proj_NAA[,y,,1,drop = FALSE] * exp(-proj_ZAA[,y,,1,drop = FALSE] * t_spawn) * WAA[,y,,1,drop = FALSE] * MatAA[,y,,1,drop = FALSE], 1, sum) # Spawning Stock Biomass
    proj_Dynamic_SSB0[,y] = apply(proj_NAA0[,y,,1,drop = FALSE] * exp(-natmort[,y,,1,drop = FALSE] * t_spawn) * WAA[,y,,1,drop = FALSE] * MatAA[,y,,1,drop = FALSE], 1, sum) # Spawning Stock Biomass
    if(n_sexes == 1) {  # If single sex model, multiply SSB calculations by 0.5
      proj_SSB[,y] = proj_SSB[,y] * 0.5
      proj_Dynamic_SSB0[,y] = proj_Dynamic_SSB0[,y] * 0.5
    }

# Derive Catches ----------------------------------------------------------
    for(r in 1:n_regions) {
      for(f in 1:n_fish_fleets) {
        for(a in 1:n_ages) {
          for(s in 1:n_sexes) {
            # Get catch at age with Baranov's
            proj_CAA[r,y,a,s,f] <- (proj_FAA[r,y,a,s,f] / proj_ZAA[r,y,a,s]) * proj_NAA[r,y,a,s] * (1 - exp(-proj_ZAA[r,y,a,s]))
          } # end s loop
        } # end a loop

        # Get total catch
        proj_Catch[r,y,f] <- sum(proj_CAA[r,y,,,f] * WAA_fish[r,y,,,f])

# Project F using HCR and reference points -----------------------------------------------------
        if(fmort_opt == 'HCR') {
          proj_F[r,y+1] <- HCR_function(x = proj_SSB[r,y],
                                        frp = f_ref_pt[r,y],
                                        brp = b_ref_pt[r,y])
        }

# Project F using User Inputs ---------------------------------------------
        if(fmort_opt == 'Input') proj_F[r,y+1] <- f_ref_pt[r,y]

      } # end f loop
    } # end r loop

  } # end y loop

  return(list(proj_F = proj_F,
              proj_Catch = proj_Catch,
              proj_SSB = proj_SSB,
              proj_Dynamic_SSB0 = proj_Dynamic_SSB0,
              proj_NAA = proj_NAA,
              proj_NAA0 = proj_NAA0,
              proj_ZAA = proj_ZAA)
  )

} # end function

