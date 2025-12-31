# ar1 correlation matrix
get_AR1_CorrMat <- function(n, rho) {
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")
  corrMatrix <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      # Calculate the correlation based on the lag distance
      corrMatrix[i, j] <- rho^(abs(i - j))
    } # end j
  } # end i
  return(corrMatrix)
}

# constant correlation matrix
get_Constant_CorrMat <- function(n, rho) {
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")
  corrMatrix <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      if(i != j) corrMatrix[i, j] <- rho
      else corrMatrix[i, j] <- 1
    } # end j
  } # end i
  return(corrMatrix)
}

# combine function
cmb <- function (f, d)
{
  function(p) f(p, d)
}

# constrain stuff between -1 and 1
rho_trans <- function (x)
{
  2/(1 + exp(-2 * x)) - 1
}

inv_rho_trans <- function(y) {
  log((y + 1)/(1 - y))/2
}

# Simulate ----------------------------------------------------------------
# set.seed(123)
# dimensions
A <- 5 # ages
S <- 2 # sexes
B <- A * S # categories

# construct 2dar1
sd <- 0.3
rho_s <- 0
rho_a <- 0.5
ar1_mat <- get_AR1_CorrMat(A, rho_a)
sex_mat <- get_Constant_CorrMat(S, rho_s)
Sigma <- kronecker(sex_mat, ar1_mat) * sd^2 / (1 - rho_s^2) / (1 - rho_a^2)
Sigma <- Sigma[-nrow(Sigma), -ncol(Sigma)]
image(Sigma)

# generate expected values
props <- runif(B, 0, 1)
props <- props / sum(props)
mu <- log(props[-length(props)])
mu <- mu - log(props[length(props)])

# simulate data
n_samples <- 100
p_matrix <- matrix(NA, nrow = n_samples, ncol = B)
x_matrix <- matrix(NA, nrow = n_samples, ncol = B - 1)
for (i in 1:n_samples) {
  x <- MASS::mvrnorm(1, mu, Sigma)
  x_matrix[i, ] <- x
  p <- exp(x) / (1 + sum(exp(x)))
  p_matrix[i, ] <- c(p, 1 - sum(p))
}


# Model -------------------------------------------------------------------

nll <- function(pars, data) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data)

  nll = 0

  # Build L matrix
  L = matrix(0, nrow = B - 1, ncol = n_factors)
  dmat = matrix(0, B - 1, B - 1) # diagonal matrix

  # build covariance w/ cholesky
  if(covar_type == 'cholesky') {
    idx = 1
    for (j in 1:n_factors) {
      for (i in j:(B-1)) { # fill only lower diag
        L[i, j] = Lvec[idx]
        idx = idx + 1
      } # end i
    } # end j

    diag(dmat) = exp(dvec)^2
    Sigma = (L %*% t(L)) + dmat # lower triangular factorization
  }

  # build covariance w/ kronecker
  if(covar_type == 'kronecker') {
    sd = exp(dvec)^2
    rho_s = rho_trans(logit_rho_s)
    rho_a = rho_trans(logit_rho_a)
    ar1_mat = get_AR1_CorrMat(A, rho_a)
    sex_mat = get_Constant_CorrMat(S, rho_s)
    Sigma = Matrix::kronecker(sex_mat, ar1_mat) * sd
    Sigma = Sigma[-nrow(Sigma), -ncol(Sigma)]
  }

  # multivariate normal
  mu = log(props[-length(props)]) - log(props[length(props)])
  for(i in 1:nrow(p)) {
    tmp_Obs = log(p[i, -length(p[i,])]) - log(p[i, length(p[i,])])
    nll = nll - RTMB::dmvnorm(x = tmp_Obs, mu = mu, Sigma = Sigma, log = TRUE)
  }

  nll = nll + sum(Lvec)^2
  RTMB::REPORT(Sigma)

  return(nll)
}

# Fit Model ---------------------------------------------------------------
n_factors <- 3

# Calculate number of lower pars
n_L_params <- n_factors * (B - 1) - n_factors * (n_factors - 1) / 2

# build data list
data <- list(
  props = props,
  p = p_matrix,
  n_factors = n_factors,
  B = B,
  Sigma_true = Sigma,
  covar_type = 'kronecker'
)

pars <- list(
  dvec = log(0.1),
  Lvec = rnorm(n_L_params, 0, 0.1),
  logit_rho_s = 0,
  logit_rho_a = 0
)

map <- NULL
map <- list(
  # logit_rho_s = factor(NA),
  # logit_rho_a = factor(NA)
  Lvec = factor(rep(NA, length(pars$Lvec)))
)

obj <- RTMB::MakeADFun(
  cmb(nll, data),
  parameters = pars,
  map = map,
  random = NULL,
  silent = FALSE
)


# rhos <- inv_rho_trans(seq(-0.95, 0.95, 0.01))
# nll_save <- vector()
# for(i in 1:length(rhos)) {
#   pars$logit_rho_a <- rhos[i]
#   nll_save[i] <- nll(pars, data)
# }
#
# plot(rho_trans(rhos), nll_save, type = 'l')

# Now, optimize the function
optim <- stats::nlminb(obj$par, obj$fn, obj$gr,
                       control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))

# profile <- TMB::tmbprofile(obj, 'dvec')
# plot(exp(profile$dvec), profile$value, type = 'l')
# abline(v = sd)

# newton steps
# try_improve <- tryCatch(expr =
#                           for(i in 1:5) {
#                             g = as.numeric(obj$gr(optim$par))
#                             h = optimHess(optim$par, fn = obj$fn, gr = obj$gr)
#                             optim$par = optim$par - solve(h,g)
#                             optim$objective = obj$fn(optim$par)
#                           }
#                         , error = function(e){e}, warning = function(w){w})


rep <- obj$report(obj$env$last.par.best)
image(rep$Sigma)
# get sdreport
sdrep <- RTMB::sdreport(obj)
sdrep

max(rep$Sigma)
max(Sigma)

plot(rep$Sigma, Sigma)

# image(Sigma)
# image(rep$Sigma)
# image(Sigma)
# qr(rep$Sigma)$rank

# image(rep$L %*% t(rep$L) + rep$dmat)
# plot(as.vector(rep$Sigma))
# lines(as.vector(Sigma))

# L_est <- rep$L
# rotation <- varimax(L_est)
# rotation$loadings

# rho = seq(-0.8,0.8, 0.01)
# plot(rho, sd^2 / (1 - rho^2) / (1 - rho^2), type = 'l', ylab = 'marginal variance')
