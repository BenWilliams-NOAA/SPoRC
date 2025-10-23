#' Initialize Numbers-at-Age (NAA) for a Population Model
#'
#' This function generates initial numbers-at-age (NAA) for a structured population
#' model across regions, sexes, and age classes. It supports multiple initialization
#' methods, including iterative solution, scalar geometric series, and matrix
#' geometric series, optionally accounting for movement and fishing mortality. Initial
#' age deviations can also be applied.
#'
#' @param init_age_strc Integer specifying the initialization method for the age structure:
#'   - 0: Iterative solution to equilibrium
#'   - 1: Scalar geometric series solution (no movement in plus group)
#'   - 2: Matrix geometric series solution (generalizes scalar solution with movement)
#' @param init_iter Integer; number of iterations to run when `init_age_strc = 0`.
#' @param n_regions Integer; number of spatial regions.
#' @param n_sexes Integer; number of sexes.
#' @param n_ages Integer; number of age classes.
#' @param natmort Array of natural mortality rates with dimensions `[regions, ages, sexes]`.
#' @param init_F Numeric; initial fishing mortality applied (0 for unfished population).
#' @param fish_sel Array of fishery selectivity with dimensions `[regions, ages, sexes, fleets]`.
#' @param R0_r Numeric vector of recruitment values for each region.
#' @param sexratio Array `[regions, sexes]` giving the proportion of each sex in recruitment.
#' @param Movement Array `[regions, regions, ages, sexes]` defining movement probabilities.
#' @param do_recruits_move Integer; 0 = recruits do not move, 1 = recruits move according to `Movement`.
#' @param ln_InitDevs Array `[regions, ages-1]` of log-scale deviations for initial numbers-at-age.
#'
#' @return Array of initial numbers-at-age with dimensions `[regions, ages, sexes]`.
#'
#' @keywords internal
Get_Init_NAA <- function(init_age_strc,
                         init_iter,
                         n_regions,
                         n_sexes,
                         n_ages,
                         natmort,
                         init_F,
                         fish_sel,
                         R0_r,
                         sexratio,
                         Movement,
                         do_recruits_move,
                         NAA,
                         ln_InitDevs) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  # create containers
  Init_NAA = array(0, dim = c(n_regions, n_ages, n_sexes))
  Init_NAA_next_year = array(0, dim = c(n_regions, n_ages, n_sexes))
  NAA = array(0, dim = c(n_regions, n_ages, n_sexes))

  # Iterative Solution
  if(init_age_strc == 0) {
    # initialize age structure (starting point)
    for(r in 1:n_regions) {
      for(s in 1:n_sexes) {
        tmp_cumsum_Z = cumsum(natmort[r,1:(n_ages-1),s] + init_F * fish_sel[r,1:(n_ages-1),s,1])
        Init_NAA[r,,s] = c(R0_r[r] * sexratio[r,s], R0_r[r] * sexratio[r,s] * exp(-tmp_cumsum_Z))
      } # end s loop
    } # end r loop

    # Apply annual cycle and iterate to equilibrium
    for(i in 1:init_iter) {
      for(s in 1:n_sexes) {
        Init_NAA_next_year[,1,s] = R0_r * sexratio[,s] # recruitment
        # movement
        if(do_recruits_move == 0) for(a in 2:n_ages) Init_NAA[,a,s] = t(Init_NAA[,a,s]) %*% Movement[,,a,s] # recruits don't move
        if(do_recruits_move == 1) for(a in 1:n_ages) Init_NAA[,a,s] = t(Init_NAA[,a,s]) %*% Movement[,,a,s] # recruits move
        for(r in 1:n_regions) {
          # ageing and mortality
          Init_NAA_next_year[r,2:n_ages,s] = Init_NAA[r,1:(n_ages-1),s] * exp(-(natmort[r,1:(n_ages-1),s] + (init_F * fish_sel[r,1:(n_ages-1),s,1])))
          # accumulate plus group
          Init_NAA_next_year[r,n_ages,s] = (Init_NAA_next_year[r,n_ages,s]) + (Init_NAA[r,n_ages,s] * exp(-(natmort[r,n_ages,s] + (init_F * fish_sel[r,n_ages,s,1]))))
        } # end r loop
        Init_NAA = Init_NAA_next_year # iterate to next cycle
      } # end s loop
    } # end i loop
    # save result
    NAA[,,] = Init_NAA
  } # end if iterative solution

  # Scalar Geometric Series Solution (no movement in plus group)
  if(init_age_strc == 1) {
    for(r in 1:n_regions) {
      for(s in 1:n_sexes) {
        NAA[r,1,s] = R0_r[r] * sexratio[r,s]  # initialize population
        for(a in 2:(n_ages-1))  NAA[r,a,s] = NAA[r,a-1,s] * exp(-(natmort[r,a-1,s] + (init_F * fish_sel[r,a-1,s,1]))) # iterate
        # Plus group - scalar geometric series
        Z_penult = natmort[r,n_ages-1,s] + (init_F * fish_sel[r,n_ages-1,s,1])
        Z_plus = natmort[r,n_ages,s] + (init_F * fish_sel[r,n_ages,s,1])
        NAA[r,n_ages,s] = NAA[r,n_ages-1,s] * exp(-Z_penult) / (1 - exp(-Z_plus))
      } # end s loop
    } # end r loop
  } # end if

  # Matrix Geometric Series Solution (genearlizes to scalar w/o movement)
  if(init_age_strc == 2) {

    # projection initial abundance forward
    for(i in 1:n_ages) {
      for(s in 1:n_sexes) {
        Init_NAA[,1,s] = R0_r * sexratio[,s] # initialize recruitment
        # movement
        if(do_recruits_move == 0) for(a in 2:n_ages) Init_NAA[,a,s] = t(Init_NAA[,a,s]) %*% Movement[,,a,s] # recruits don't move
        if(do_recruits_move == 1) for(a in 1:n_ages) Init_NAA[,a,s] = t(Init_NAA[,a,s]) %*% Movement[,,a,s] # recruits move
        for(r in 1:n_regions) {
          # ageing and mortality
          Init_NAA[r,2:n_ages,s] = Init_NAA[r,1:(n_ages-1),s] * exp(-(natmort[r,1:(n_ages-1),s] + (init_F * fish_sel[r,1:(n_ages-1),s,1])))
          # accumulate plus group
          Init_NAA[r,n_ages,s] = (Init_NAA[r,n_ages,s]) + (Init_NAA[r,n_ages,s] * exp(-(natmort[r,n_ages,s] + (init_F * fish_sel[r,n_ages,s,1]))))
        } # end r loop
      } # end s loop
    } # end a loop

    # Set up analytical solution for plus group
    for(s in 1:n_sexes) {
      Move_penult = Movement[,,n_ages-1,s] # get movement matrices
      Move_plus = Movement[,,n_ages,s] # get movement matrices
      s_penult = exp(-(natmort[,n_ages-1,s] + (init_F * fish_sel[,n_ages-1,s,1]))) # survival of penultimate age
      s_plus = exp(-(natmort[,n_ages,s] + (init_F * fish_sel[,n_ages,s,1]))) # survival of plus group age
      init_penult <- Init_NAA[,n_ages-1,s] # get initial penultimate age
      source = (t(Move_penult) %*% init_penult) * s_penult # compute forward projection of penultimate age to plus group
      T_mat = diag(s_plus, n_regions) %*% t(Move_plus) # transition matrix
      I_mat = diag(n_regions) # get identity matrix
      plus = solve(I_mat - T_mat, source) # solve to get plus group (I-T)^-1 %*% source
      Init_NAA[,n_ages,s] = plus # input plus group here
    }
    # save result
    NAA = Init_NAA
  }

  # Apply initial age deviations
  for(r in 1:n_regions) {
    for(s in 1:n_sexes) {
      NAA[r,2:n_ages,s] = NAA[r,2:n_ages,s] * exp(ln_InitDevs[r,])
    } # end s loop
  } # end r loop

  return(NAA)

}
