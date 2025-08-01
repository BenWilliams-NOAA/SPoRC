#' ggplot theme for sablefish
#'
#' @return ggplot theme
#' @export theme_sablefish
#' @import ggplot2
theme_sablefish <- function() {
   theme_bw() +
    theme(legend.position = "top",
          strip.text = element_text(size = 17),
          title = element_text(size = 21, color = 'black'),
          axis.text = element_text(size = 15, color = "black"),
          axis.title = element_text(size = 17, color = 'black'),
          legend.text = element_text(size = 15, color = "black"),
          legend.title = element_text(size = 17, color = 'black'))
}


#' Function to fill in an n x n correlation AR(1) matrix
#'
#' @param n Number of bins
#' @param rho correaltion parameter
#'
#' @return correlation matrix for an ar1 process
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' mat <- get_AR1_CorrMat(10, 0.5)
#' }
get_AR1_CorrMat <- function(n, rho) {
  corrMatrix <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      # Calculate the correlation based on the lag distance
      corrMatrix[i, j] <- rho^(abs(i - j))
    } # end i
  } # end j
  return(corrMatrix)
}

#' Constant correlation matrix
#'
#' @param n Number of bins
#' @param rho correaltion parameter
#'
#' @return constant correlation matrix
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' mat <- get_Constant_CorrMat(10, 0.5)
#' }
get_Constant_CorrMat <- function(n, rho) {
  corrMatrix <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      if(i != j) corrMatrix[i, j] <- rho
      else corrMatrix[i, j] <- 1
    } # end i
  } # end j
  return(corrMatrix)
}

#' For combining a parameter and data list in RTMB so a data object can be explicitly defined
#'
#' @param f Parameter list
#' @param d Data list
#' @keywords internal
#' @examples
#' \dontrun{
#'   obj <- RTMB::MakeADFun(cmb(sabie_RTMB, data), parameters = parameters, map = mapping, random = random, silent = TRUE)
#' }
cmb <- function(f, d) {
  function(p) f(p, d)
}

#' Helper function to collect messages
#'
#' @param ... character vector of messages
#' @keywords internal
collect_message <- function(...) {
  messages_list <<- c(messages_list, paste(..., sep = ""))
}

#' Helper function to check package availbility
#'
#' @param pkg package name character
#' @keywords internal
is_package_available <- function(pkg) {
  nzchar(system.file(package = pkg))
}

#' Go from TAC to Fishing Mortality using bisection
#'
#' @param f_guess Initial guess of F
#' @param catch Provided catch values
#' @param NAA Numbers, dimensioned by ages, and sexes
#' @param WAA Weight, dimensioned by ages and sexes
#' @param natmort Natural mortality dimensioned by ages and sex
#' @param fish_sel Fishery selectivity, dimesnioned by ages and sex
#' @param n.iter Number of iterations for bisection
#' @param lb Lower bound of F
#' @param ub Upper bound of F
#'
#' @returns Fishing mortality values
#' @export bisection_F
bisection_F <- function(f_guess,
                        catch,
                        NAA,
                        WAA,
                        natmort,
                        fish_sel,
                        n.iter = 20,
                        lb = 0,
                        ub = 2) {

  range <- vector(length=2) # F range
  range[1] <- lb # Lower bound
  range[2] <- ub # Upper bound

  for(i in 1:n.iter) {

    # Get midpoint of range
    midpoint <- mean(range)

    # Caclulate baranov's
    FAA <- (midpoint * fish_sel)
    ZAA <- FAA + natmort
    pred_catch <- sum((FAA / ZAA * NAA * (1 - exp(-ZAA))) * WAA)

    if(pred_catch < catch) {
      range[1] <- midpoint
      range[2] <- range[2]
    }else {
      range[1] <- range[1]
      range[2] <- midpoint
    }

  } # end i loop

  return(midpoint)
}

#' Post Optimization Model Convergence Checks
#'
#' @param sd_rep sd report list from a `SPoRC` model
#' @param rep report list from a `SPoRC` model
#' @param gradient_tol Value for maximum gradient tolerance to use
#' @param se_tol Value for maximum standard error tolerance to use
#' @param corr_tol Value for maximum correlation tolerance to use
#'
#' @export post_optim_sanity_checks
post_optim_sanity_checks <- function(sd_rep,
                                     rep,
                                     gradient_tol = 1e-3,
                                     se_tol = 100,
                                     corr_tol = 0.99
                                     ) {

  passed_post_sanity_checks <- TRUE

  # check likelihoods are all finite and not NA
  if(!all(is.finite(rep$jnLL))) {
    message("Found Inf in joint log-likelihood, model is not converged!")
    passed_post_sanity_checks <- F
  }

  # check maximum absolute gradients
  max_abs_grad_ndx <- which.max(abs(sd_rep$gradient.fixed))
  max_abs_grad <- abs(sd_rep$gradient.fixed)[max_abs_grad_ndx]
  if(gradient_tol < max_abs_grad) {
    message("Parameter: ", names(sd_rep$par.fixed)[max_abs_grad_ndx], " had absolute gradient = ", max_abs_grad,
            " which was greater than tolerance ", gradient_tol,". This indicates potential non-convergence according to the tolerance.\n")
    passed_post_sanity_checks <- F
  }

  # check hessian
  if(!sd_rep$pdHess) {
    message("Hessian is not positive definite, model is not converged!")
    passed_post_sanity_checks <- F
  }

  # check if standard errors are finite (if finite, then check other stuff)
  if(!all(is.finite(sqrt(diag(sd_rep$cov.fixed))))) {
    message("Found non finite elements in standard errors of parameters, model is not converged!")
    passed_post_sanity_checks <- F
  } else {
    # check if standard errors are big
    if(max(sqrt(diag(sd_rep$cov.fixed))) > se_tol) {
      message("Parameter: ", names(diag(sd_rep$cov.fixed))[which.max(sqrt(diag(sd_rep$cov.fixed)))], " has a standard error = ",
              max(sqrt(diag(sd_rep$cov.fixed))), " which was greated than tolerance ", se_tol, ". This indicates potential non-convergence according to the tolerance. \n")
      passed_post_sanity_checks <- F
    }

    # check if correlations are big
    corr_mat <- cov2cor(sd_rep$cov.fixed)
    diag(corr_mat) <- "Same" # set diagonal to "Same" to remove from max calculations

    # reshape to dataframe
    corr_df <- reshape2::melt(corr_mat) %>%
      dplyr::filter(value != 'Same') %>%
      dplyr::mutate(value = as.numeric(value))

    if(max(abs(corr_df$value)) > corr_tol) {
      message("Parameter pairs: ", corr_df$Var1[which.max(abs(corr_df$value))], " and ", corr_df$Var2[which.max(abs(corr_df$value))], " have a correlation of ", max(abs(corr_df$value)), ". This indicates potential non-convergence according to the tolerance.")
      passed_post_sanity_checks <- F
    }
  }

  if(passed_post_sanity_checks) {
    message("Successfully passed post-optim-sanity checks\n")
  }

  return(passed_post_sanity_checks)

}

#' Helper function for extracting elements from TMB report
#'
#' @param obj TMB report object
#' @param name Name of object to be extracted
#' @keywords internal
safe_extract <- function(obj, name) {
  if (name %in% names(obj) && !is.null(obj[[name]])) {
    return(obj[[name]])
  } else {
    return(0)
  }
}

#' Helper function for extracting parameter information and names from TMB
#'
#' @param parameters Parameter list from setting up TMB object
#' @param mapping Mapping list from setting up TMB object
#' @param sd_rep SD Report from TMB obj
#'
#' @returns A list of dataframes for estimated and non-estimated parameter values.
#' @export get_par_est_info
get_par_est_info <- function(parameters, mapping, sd_rep) {

  # get parameter names
  par_names <- reshape2::melt(parameters) %>%
    dplyr::rename_with(~str_replace(., "^Var(\\d+)$", "Dim\\1")) %>%
    dplyr::rename(Init_Val = value, Par = L1) %>%
    dplyr::group_by(Par) %>%
    # unique parameters based on dimensions of parameter list
    dplyr::mutate(Par_Num = paste(Par, row_number(), sep = "_"))

  # get mapping names
  map_names <- reshape2::melt(mapping) %>%
    dplyr::rename(map = value, Par = L1) %>%
    dplyr::group_by(Par) %>%
    dplyr::mutate(Par_Num = paste(Par, row_number(), sep = "_"),
                  map = as.character(map), # make character
                  map = ifelse(is.na(map), 'NE', map) # denote NAs as NE (for not estimated instead)
    )

  # join parameter names and mapping names
  par_map_names <- par_names %>%
    dplyr::left_join(map_names, by = c("Par", "Par_Num")) %>%
    dplyr::group_by(Par) %>%
    # make sure to turn NAs (they are estimated, but just were not in the mapping list)
    dplyr::mutate(map = ifelse(is.na(map), as.character(row_number()), map))

  # Make sure mapping numbers are sequential to match up with the sdreport
  par_map_names <- par_map_names %>%
    # filter(str_detect(Par, "ln_F_mean")) %>%
    dplyr::group_by(Par) %>%
    dplyr::mutate(map = ifelse(map != 'NE', as.character(cumsum(map != 'NE')), 'NE'),
                  Par_Num_map = paste(Par, map, sep = "_")) # now, make a variable that is consistent with numbering in sdreport

  # now, get estimater parameter names and values
  est_names <- data.frame(Par = names(sd_rep$par.fixed),
                          Est_Val = sd_rep$par.fixed,
                          SE_Val = sqrt(diag(sd_rep$cov.fixed)),
                          Abs_Grad_Val = abs(as.vector(sd_rep$gradient.fixed))) %>%
    dplyr::group_by(Par) %>%
    dplyr::mutate(Par_Num_map = paste(Par, row_number(), sep = '_'))

  # join to get estimated parameters, along with initial starting values
  estimated_pars <- est_names %>%
    dplyr::left_join(par_map_names, by = c("Par", "Par_Num_map")) %>%
    dplyr::select(-c(Par_Num))

  # also get non-estimated parameters
  non_estimated_pars <- par_map_names %>%
    dplyr::filter(map == 'NE')

  return(
    list(est_pars = estimated_pars,
         non_est_pars = non_estimated_pars)
  )
}
