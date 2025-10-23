get_movement_dp_design_matrix <- function(data,
                                          preference_formula,
                                          diffusion_formula
) {
  X_zk = model.matrix(diffusion_formula, data) # diffussion design matrix
  W_zk = model.matrix(preference_formula, data) # preference design matrix
  return(list(
    n_theta = ncol(X_zk),
    n_gamma = ncol(W_zk),
    X_zk = X_zk,
    W_zk = W_zk
  ))
}

Get_Movement <- function(move_type,
                         do_recruits_move,
                         n_regions,
                         n_yrs,
                         n_proj_yrs_devs,
                         n_ages,
                         n_sexes,
                         move_pars,
                         move_devs,
                         use_fixed_movement,
                         Fixed_Movement = NULL,
                         ctmc_move_dat = NULL,
                         preference_formula = NULL,
                         diffusion_formula = NULL,
                         log_move_diffusion_pars,
                         move_preference_pars,
                         area_r,
                         adjacency_mat
                         ) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  move_pen = 0 # initialize movement penalty if used

  # use fixed movement matrix
  if(use_fixed_movement == 1) {
    Movement = Fixed_Movement
  } else if(move_type == 0) { # Unstructured markov movement
    Movement = array(data = 0, dim = c(n_regions, n_regions, n_yrs + n_proj_yrs_devs, n_ages, n_sexes)) # movement "matrix"
    ref_region = 1 # Set up reference region (always set at 0)
    for(r in 1:n_regions) {
      for(y in 1:(n_yrs + n_proj_yrs_devs)) {
        for(a in 1:n_ages) {
          for(s in 1:n_sexes) {

            move_tmp = rep(0, n_regions) # temporary movement vector to store values
            counter = 1  # counter

            for(rr in 1:n_regions) {
              if(rr != ref_region) {
                # extract movement parameters
                if(y <= n_yrs) tmp_move_pars = move_pars[r,counter,y,a,s]
                else tmp_move_pars = move_pars[r,counter,n_yrs,a,s]
                move_tmp[rr] = tmp_move_pars + move_devs[r,counter,y,a,s]
                counter = counter + 1
              } # end if not reference region
            } # end rr loop

            if(use_fixed_movement == 0) Movement[r,,y,a,s] = exp(move_tmp) / sum(exp(move_tmp)) # multinomial logit transform estimated movement

          } # end s loop
        } # end a loop
      } # end y loop
    } # end r loop
  } else if(move_type == 1) { # continuous markov chain movement, NOTE: no projections supported at the moment

    # set up dimensions of movement matrix
    dims = list(from = 1:n_regions, to = 1:n_regions, years = 1:n_yrs, ages = 1:n_ages, sexes = 1:n_sexes)
    Movement = Taxis = Diffusion = array(0, dim = sapply(dims, length),  dimnames = dims)
    loop = expand.grid(dims[-(1:2)]) # get year, age, and sexes to loop through
    if(do_recruits_move == 0) loop = loop[-which(loop$ages == 1),] # remove age 1, if recruits move

    # setup design matrix
    design_mat = get_movement_dp_design_matrix(ctmc_move_dat, preference_formula, diffusion_formula)
    X_zk = design_mat$X_zk # diffussion
    W_zk = design_mat$W_zk # preference

    # diffusion rate from each region
    theta_k = exp(log_move_diffusion_pars) # get difussion parameter
    theta_z = (X_zk %*% theta_k)[,1] # multiply difussion parameter by design matrix
    theta_z = theta_z/area_r[ctmc_move_dat[,'regions']]  # scale difussion matrix by area

    # preference for each region
    gamma_k = move_preference_pars # get preference parameters
    gamma_z = (W_zk %*% gamma_k)[,1] # multiply preference parameters by design matrix

    # Make instantaneous diffusion rate matrix
    for( index in seq_len(nrow(loop)) ){

      # get year, age, and sex specific indices for a given stratum combination
      which_rows = expand.grid( 1:n_regions, loop[index,"years"], loop[index,"ages"], loop[index,"sexes"] )
      which_rows$index = NA
      colnames(which_rows) = c( "regions", names(loop), "index" )

      # match up indices to provided data
      for( i2 in seq_len(nrow(which_rows)) ){
        which_rows$index[i2] = which((which_rows[i2,'regions'] == ctmc_move_dat[,'regions']) &
                                       which_rows[i2, "years"] == ctmc_move_dat[,"years"] &
                                       (which_rows[i2,'ages'] == ctmc_move_dat[,'ages']) &
                                       (which_rows[i2,'sexes'] == ctmc_move_dat[,'sexes']) )
      } # end i2 loop

      # create difussion matrix for strat, year, age, sex combinations
      D_ss = adjacency_mat %*% diag(theta_z[which_rows$index], n_regions) # get corresponding thetas
      diag(D_ss) = -1 * Matrix::colSums(D_ss) # diag to enforce sum to 1
      D_ss = as(D_ss, "sparseMatrix") # force sparse

      # preference for each strata, year, age, sex combination
      gamma_s = gamma_z[which_rows$index] # get corresponding gammas
      pref_s = exp(gamma_s) / sum(exp(gamma_s)) # softmax to keep estimation constrainted
      Z_ss = adjacency_mat * outer( pref_s, pref_s, FUN = "-" ) # h(i) - h(j)
      diag(Z_ss) = -1 * Matrix::colSums(Z_ss) # diag to enforce sum to 1

      # turn continuous rates to movement fractions
      M_ss = Matrix::expm( D_ss + Z_ss  )

      # populate matrices
      Movement[,,loop$years[index],loop$ages[index],loop$sexes[index]] = t(as.matrix(M_ss))
      Taxis[,,loop$years[index],loop$ages[index],loop$sexes[index]] = t(as.matrix(Z_ss))
      Diffusion[,,loop$years[index],loop$ages[index],loop$sexes[index]] = t(as.matrix(D_ss))

      # return penalty (Lagrange multiplier)
      move_pen = move_pen + sum(gamma_s)^2

    } # end index loop
  }

  # recruits don't move
  if(do_recruits_move == 0) Movement[,,,1,] = diag(1, n_regions)

  return(list(Movement = Movement, move_pen = move_pen))
}

