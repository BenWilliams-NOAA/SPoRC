
####################
# Notes
#  * I didn't include the area adjustment ... smaller areas should have higher bulk-diffusion than large areas (higher surface area to area ratio) ... easy to add
#  * I only added the matrix exponential option ... backwards or forwards Euler might be faster in some cases, but prolly not necessary yet
#  * I didn't show the extension to variation over ages/sexes/years/seasons, but that could be added using the formula interface
#  * RTMB is finicky about sparse matrices, so follow up if you're having problems with advector being dropped
#  * I used packages explicitly ... not sure what you're format is
#  * I used a likelihood penalty on preference sum-to-zero constraint (rather than an explicit multivariate-logit) ... its a bit like a Lagrange multiplier and shouldn't affect the MLE, while perhaps being easier numerically ... happy to discuss
####################

# make movement design matrix and get dimensions
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

# movement matrix
make_M <-
function( A_ss,
          parlist,
          area_s = rep(1,nrow(A_ss)), # Equal areas for simplicity of presentation,
          data,
          preference_formula = ~ 0,
          diffusion_formula = ~ 1,
          method = 0  # 0:matexp;  >0:forwards_Euler; <0: backwards_Euler
          ){

  # get design matrices
  design_mat <- get_movement_dp_design_matrix(data, preference_formula, diffusion_formula)
  X_zk <- design_mat$X_zk # diffussion
  W_zk <- design_mat$W_zk # preference

  # diffusion rate from each stratum
  theta_k = exp(parlist$log_theta) # get difussion parameter
  theta_z = (X_zk %*% theta_k)[,1] # multiply difussion parameter by design matrix
  theta_z = theta_z/area_s[data[,'stratum']]  # scale difussion matrix by area

  # preference for each stratum
  gamma_k = parlist$gamma # get preference parameters
  gamma_z = (W_zk %*% gamma_k)[,1] # multiply preference parameters by design matrix

  # set up dimensions of movement matrix
  dims = list( from = strata, to = strata, years = years, age = ages, sexes = sexes)
  M = Z = D = array(0, dim = sapply(dims, length), dimnames = dims )
  loop = expand.grid( dims[-(1:2)] ) # get year, age, and sexes to loop through
  log_penalty = 0 # penalty to enforce sum to zero

  # Make instantaneous diffusion rate matrix
  for( index in seq_len(nrow(loop)) ){

    # get year, age, and sex specific indices for a given stratum combination
    which_rows = expand.grid( strata, loop[index,"years"], loop[index,"age"], loop[index,"sexes"] )
    which_rows$index = NA
    colnames(which_rows) = c( "stratum", names(loop), "index" )

    # match up indices to provided data
    for( i2 in seq_len(nrow(which_rows)) ){
      which_rows$index[i2] = which( (which_rows[i2,'stratum'] == data[,'stratum']) &
                                      which_rows[i2, "years"] == data[,"years"] &
                                      (which_rows[i2,'age'] == data[,'age']) &
                                      (which_rows[i2,'sexes'] == data[,'sexes']) )
    }

    # create difussion matrix for strat, year, age, sex combinations
    D_ss = A_ss %*% Matrix::Diagonal( n = n_s, x = theta_z[which_rows$index] )
    diag(D_ss) = -1 * Matrix::colSums(D_ss) # diag to enforce sum to 1

    # preference for each strata, year, age, sex combination
    gamma_s = gamma_z[which_rows$index]
    pref_s = exp(gamma_s) / sum(exp(gamma_s))
    Z_ss = A_ss * outer( pref_s, pref_s, FUN = "-" )
    diag(Z_ss) = -1 * Matrix::colSums(Z_ss) # diag to enforce sum to 1

    # Make movement fraction matrix
    if(method == 0){
      M_ss = Matrix::expm( D_ss + Z_ss ) # turn continuous rates to movement fractions
    } else{
      stop("Not added yet")
      # Could add other options if it's too slow
    }

    # populate matrices
    M[,,loop$years[index],loop$age[index],loop$sexes[index]] = t(as.matrix(M_ss))
    Z[,,loop$years[index],loop$age[index],loop$sexes[index]] = t(as.matrix(Z_ss))
    D[,,loop$years[index],loop$age[index],loop$sexes[index]] = t(as.matrix(D_ss))

    # return penalty (Lagrange multiplier)
    log_penalty = log_penalty + sum(gamma_s)^2
  }

  return(list(M = M, log_penalty = log_penalty))
}


##################
# EXAMPLE
##################

library(igraph)
library(splines2)

# Geographic domain and age-set
Domain = make_graph(
  ~ 1 - 2,
  2 - 3,
  3 - 4,
  4 - 5
)
ages = 1:30
years = 1:40
sexes = 1:2

# Convenience
strata = names(V(Domain))
n_s = length(strata)
n_a = length(ages)
n_y = length(years)
n_sx = length(sexes)

# make adjacency matrix
A_ss = as_adjacency_matrix(Domain)

# Generate data frame
data = expand.grid(
  stratum = strata,
  years = years,
  age = ages,
  sexes = sexes
)

# Associate with covariates .. these are both static across ages
# data$Temp_C = rnorm(n_s)[match(strata, data$stratum)]
# data$Bathy_m = exp(rnorm(n_s) + log(100))[match(strata,data$stratum)]

# Needs some sanity checks on identifiability, or user understanding theory
diffusion_formula = ~ 1
# preference_formula = ~ 0 + factor(stratum):bSpline(age, df=4, intercept = TRUE):bSpline(years, df = 4, intercept = TRUE):factor(sexes)
preference_formula = ~ 0 + factor(stratum):bSpline(age, df=4, intercept = TRUE)

# Get design matrices to get size of parameters
designs <- get_movement_dp_design_matrix(data, preference_formula, diffusion_formula)

set.seed(123)
parlist = list(
  log_theta = log(0.1) * 1:designs$n_theta,
  gamma = 0.5 * rnorm(designs$n_gamma)
)

out = make_M(
  A_ss = A_ss,
  area_s = c(1,2,1,1,1),
  parlist = parlist,
  data = data,
  preference_formula = preference_formula,
  diffusion_formula = diffusion_formula
)

M = out$M

plot(M[1,1,1,,1])

apply( M, MARGIN=2:5, FUN = sum ) # conserves abundance within stratum combination

# Dominant eigenvector for t(M_ss) is proportional to the stationary distribution
apply(
  M,
  MARGIN = 3:5,
  FUN = function(mat){
    eigen(mat)$vectors[,1] / sum(eigen(mat)$vectors[,1])
  }
)






# GAM Scratch -------------------------------------------------------------

# out_move = make_M_gam(
#   A_ss = A_ss,
#   log_move_theta = log_move_theta,
#   move_gamma = move_gamma,
#   log_lambda = log_lambda,
#   area_s = area_s,
#   move_dim = move_dim,
#   preference_formula = preference_formula,
#   diffusion_formula = diffusion_formula,
#   method = 0
# )

#
# get_movement_gam_design_matrix <- function(data, preference_formula, diffusion_formula) {
#
#   # set up diffusion formula
#   gam_X_zk = mgcv::gam(
#     update.formula(old = diffusion_formula, new = "fake ~ ."),
#     data = cbind("fake" = 0, data),
#     fit = FALSE
#   )
#
#   # get diffusion
#   X_zk = gam_X_zk$X
#
#   # set up preference formula
#   gam_W_zk = mgcv::gam(
#     update.formula(old = preference_formula, new = "fake ~ ."),
#     data = cbind("fake" = 0, data),
#     fit = FALSE
#   )
#
#   # get preference
#   W_zk = gam_W_zk$X
#
#   # return preference function penalties
#   S_list = lapply( seq_along(gam_W_zk$smooth), function(x) gam_W_zk$smooth[[x]]$S )
#   S_kk = Matrix::.bdiag( lapply(S_list, Matrix::.bdiag) )
#   Sdims = unlist( lapply(S_list, FUN = function(List){nrow(List[[1]])}) )
#   Sblock = unlist( lapply(S_list, length) )
#   S = gam_W_zk$S
#
#   return(
#     list(X_zk = X_zk,
#          W_zk = W_zk,
#          S = S,
#          S_kk = S_kk,
#          Sdims = Sdims,
#          Sblock = Sblock)
#   )
# }
#
# get_gam_penalty <- function(gamma_k, Sdims, Sblock, S_kk, log_lambda) {
#
#   # ADoverload to make sure ad class doesn't get dropped
#   "c" <- RTMB::ADoverload("c")
#   "[<-" <- RTMB::ADoverload("[<-")
#
#   nll = 0.0
#   start_gamma = 1   # R is 1-indexed
#   start_k = 1       # R is 1-indexed
#   index_lambda = 1  # R is 1-indexed
#
#   for(z in seq_along(Sdims)) {
#     ngamma_z = Sdims[z]
#
#     # Extract gamma segment
#     gamma_segment = gamma_k[start_gamma:(start_gamma + ngamma_z - 1)]
#     start_gamma = start_gamma + ngamma_z
#
#     # s() splines (single penalty matrix)
#     if(Sblock[z] == 1) {
#       # Extract penalty matrix block
#       S_block = S_kk[start_k:(start_k + ngamma_z - 1),
#                      start_k:(start_k + ngamma_z - 1)]
#
#       # Quadratic form: gamma' * S * gamma
#       quad_form = as.numeric(t(gamma_segment) %*% S_block %*% gamma_segment)
#
#       # Penalty contribution
#       nll = nll - 0.5 * ngamma_z * log_lambda[index_lambda] +
#         0.5 * exp(log_lambda[index_lambda]) * quad_form
#
#       start_k = start_k + ngamma_z
#       index_lambda = index_lambda + 1
#     }
#   }
#
#   return(nll)
# }


# movement matrix
make_M_gam <-
  function( A_ss,
            log_move_theta,
            move_gamma,
            log_lambda,
            area_s = rep(1,nrow(A_ss)), # Equal areas for simplicity of presentation,
            move_dim,
            preference_formula = ~ 0,
            diffusion_formula = ~ 1,
            method = 0  # 0:matexp;  >0:forwards_Euler; <0: backwards_Euler
  ){

    # ADoverload to make sure ad class doesn't get dropped
    "c" <- RTMB::ADoverload("c")
    "[<-" <- RTMB::ADoverload("[<-")

    # get design matrices
    design_mat = get_movement_gam_design_matrix(move_dim, preference_formula, diffusion_formula)
    X_zk = design_mat$X_zk # diffussion
    W_zk = design_mat$W_zk # preference

    # diffusion rate from each stratum
    theta_k = exp(log_move_theta) # get difussion parameter
    theta_z = (X_zk %*% theta_k)[,1] # multiply difussion parameter by design matrix
    theta_z = theta_z/area_s[move_dim[,'stratum']]  # scale difussion matrix by area

    # preference for each stratum
    gamma_k = move_gamma # get preference parameters
    gamma_z = (W_zk %*% gamma_k)[,1] # multiply preference parameters by design matrix

    # set up dimensions of movement matrix
    dims = list( from = strata, to = strata, years = years, age = ages, sexes = sexes)
    M = Z = D = array(0, dim = sapply(dims, length), dimnames = dims )
    loop = expand.grid( dims[-(1:2)] ) # get year, age, and sexes to loop through
    pen = 0 # penalty to enforce sum to zero

    # Make instantaneous diffusion rate matrix
    for( index in seq_len(nrow(loop)) ){

      # get year, age, and sex specific indices for a given stratum combination
      which_rows = expand.grid( strata, loop[index,"years"], loop[index,"age"], loop[index,"sexes"] )
      which_rows$index = NA
      colnames(which_rows) = c( "stratum", names(loop), "index" )

      # match up indices to provided data
      for( i2 in seq_len(nrow(which_rows)) ){
        which_rows$index[i2] = which( (which_rows[i2,'stratum'] == move_dim[,'stratum']) &
                                        which_rows[i2, "years"] == move_dim[,"years"] &
                                        (which_rows[i2,'age'] == move_dim[,'age']) &
                                        (which_rows[i2,'sexes'] == move_dim[,'sexes']) )
      }

      # create difussion matrix for strat, year, age, sex combinations
      D_ss = A_ss %*% diag(theta_z[which_rows$index], n_s)
      diag(D_ss) = -1 * Matrix::colSums(D_ss) # diag to enforce sum to 1
      D_ss = as(D_ss, "sparseMatrix") # force sparse

      # preference for each strata, year, age, sex combination
      gamma_s = gamma_z[which_rows$index]
      pref_s = exp(gamma_s) / sum(exp(gamma_s))
      Z_ss = A_ss * outer( pref_s, pref_s, FUN = "-" )
      diag(Z_ss) = -1 * Matrix::colSums(Z_ss) # diag to enforce sum to 1

      # Make movement fraction matrix
      if(method == 0){
        M_ss = Matrix::expm( D_ss + Z_ss ) # turn continuous rates to movement fractions
      } else{
        stop("Not added yet")
        # Could add other options if it's too slow
      }

      # populate matrices
      M[,,loop$years[index],loop$age[index],loop$sexes[index]] = t(as.matrix(M_ss))
      # M[,,loop$years[index],loop$age[index],loop$sexes[index]] = sweep(as.matrix(M_ss), 1, rowSums(M_ss), "/")
      # Z[,,loop$years[index],loop$age[index],loop$sexes[index]] = as.matrix(Z_ss)
      # D[,,loop$years[index],loop$age[index],loop$sexes[index]] = as.matrix(D_ss)

      # return penalty (Lagrange multiplier)
      pen = pen + sum(gamma_s)^2
    }

    # Return gam penalty
    gam_pen = get_gam_penalty(move_gamma, Sdims =  design_mat$Sdims,
                              Sblock = design_mat$Sblock,
                              S_kk = design_mat$S_kk,
                              log_lambda = log_lambda)

    pen = gam_pen + pen

    return(list(Movement = M, move_pen = pen))
  }

# GAM setup ---------------------------------------------------------------
# data$diffusion_formula = ~ 1
# data$preference_formula = ~ 0 + s(age, by = factor(stratum), bs = 'tp', k = 4)
#
# designs = get_movement_gam_design_matrix(data$move_dim, preference_formula = data$preference_formula,
#                                          diffusion_formula = data$diffusion_formula)
#
#
# # movement and gam parameters
# parameters$log_move_theta = log(0.1)
# parameters$move_gamma = rnorm(ncol(designs$W_zk), 0, 0.1)
# parameters$log_lambda = rep(3, length(designs$S))
# # mapping$log_lambda <- factor(rep(NA, 3))
#
# out = make_M_gam(
#   A_ss = data$A_ss,
#   log_move_theta = parameters$log_move_theta,
#   move_gamma = parameters$move_gamma,
#   log_lambda = parameters$log_lambda,
#   area_s = data$area_s,
#   move_dim = data$move_dim,
#   preference_formula = data$preference_formula,
#   diffusion_formula = data$diffusion_formula,
#   method = 0
# )

# rowSums(out$Z[,,1,1,1] + out$D[,,1,1,1])
# out$M[,,1,1,1]
# rowSums(out$M[,,1,1,1])
#
# Q_ss = out$Z[,,1,1,1] + out$D[,,1,1,1]
# print("Off-diagonal entries:")
# print(Q_ss[row(Q_ss) != col(Q_ss)])
# print("Min off-diagonal:")
# print(min(Q_ss[row(Q_ss) != col(Q_ss)]))

# Movement rates
# reshape2::melt(out$Movement) %>%
#   filter(sexes == 1) %>%
#   mutate(
#     from = paste("from", c(1:3)[from]),
#     to = paste("to", c(1:3)[to])) %>%
#   ggplot(aes(x = age, y = value, group = years)) +
#   geom_line(lwd = 1) +
#   ggh4x::facet_grid2(to~from, scales = 'free', independent = "all")
#
#
# penalty_value <- get_gam_penalty(
#   parameters$move_gamma,
#   Sdims = designs$Sdims,
#   Sblock = designs$Sblock,
#   S_kk = designs$S_kk,
#   log_lambda = parameters$log_lambda
# )
