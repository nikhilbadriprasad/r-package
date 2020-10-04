#' Fit Gaussian Process and creates predictions on climate data using optim
#'
#' @param obj An object of class \code{climr} from \code{\link{load_clim}}
#' @param method The type of optimisation method required, either Nelder-Mead (\code{Nelder-Mead}), BFGS(\code{BFGS}),SANN(\code{SANN}) or Brent (\code{Brent})
#'
#' @return Returns a list of class \code{climr_gp_fit} which includes the predicted values as well as the data set and method type used
#' @seealso \code{\link{load_clim}}, \code{\link{plot.climr_gp_fit}},\code{\link{plot.climr_fit}}, \code{\link{fit}}
#' @export
#' @importFrom magrittr "%$%"
#' @importFrom mvtnorm "dmvnorm"
#' @importFrom stats "sd" "optim"
#' @importFrom tibble "tibble"
#'
#' @examples
#' data1 = load_clim('NH')
#' data2 = gp_fit(data1, 'BFGS')

#' @export
gp_fit <- function(obj, method = c('Nelder-Mead', 'BFGS', 'SANN', 'Brent')) {
  UseMethod('gp_fit')
}

# Defining criterion to be minimised in Gaussian process
gp_criterion = function(p,x,y) {
  sig_sq = exp(p[1])
  rho_sq = exp(p[2])
  tau_sq = exp(p[3])
  Mu = rep(0, length(x))
  Sigma = sig_sq * exp( - rho_sq * outer(x, x, '-')^2 ) + tau_sq * diag(length(x))
  ll = dmvnorm(y, Mu, Sigma, log = TRUE)
  return(-ll)
}

#' @export
gp_fit.climr <- function(obj, method = c('Nelder-Mead', 'BFGS', 'SANN', 'Brent')) {

  # Create global variables to avoid annoying CRAN notes
  DJF = Dec = `J-D` = Jan = SON = Year = month = pred = quarter = temp = x = year = NULL

  # Scaling the temperature variable
  y <- scale(obj$clim_year$temp)
  y <- y[,]
  x <- obj$clim_year$year
  x_g <- pretty(x, n = 100)

  # Find what type of fitting method
  fit_arg <- match.arg(method)

  # Fit models and optimising the hyperparameters of the Gaussian process using optim
  if(fit_arg == 'Nelder-Mead') {
    mod <- obj$clim_year %$% optim(rep(0, 3), gp_criterion, x = x, y = y, method = fit_arg)
  } else if(fit_arg == 'BFGS') {
    mod <- obj$clim_year %$% optim(rep(0, 3), gp_criterion, x = x, y = y, method = fit_arg)
  } else if(fit_arg == 'SANN') {
    mod <- obj$clim_year %$% optim(rep(0, 3), gp_criterion, x = x, y = y, method = fit_arg)
  } else if(fit_arg == 'Brent') {
    mod <- obj$clim_year %$% optim(rep(0, 3), gp_criterion, x = x, y = y, method = fit_arg)
  }

  # Extract the results
  sig_sq = exp(mod$par[1])
  rho_sq = exp(mod$par[2])
  tau_sq = exp(mod$par[3])

  # Create covariance matrices
  C = sig_sq * exp( - rho_sq * outer(x_g, x, '-')^2 )
  Sigma = sig_sq * exp( - rho_sq * outer(x, x, '-')^2 ) + tau_sq * diag(length(x))

  # Create predictions
  pred_gp = C %*% solve(Sigma, y)

  # Unscaling the variables
  temp1 <- sd(obj$clim_year$temp)
  mean <- mean(obj$clim_year$temp)
  pred_gp <- pred_gp * temp1
  pred_gp <- pred_gp + mean

  # Output of prdeicted values
  out <- list(model = pred_gp,
              data = obj$clim_year,
              method = fit_arg)
  class(out) <- 'climr_gp_fit'

  invisible(out)
}
