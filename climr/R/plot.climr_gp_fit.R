#' Plot climr_gp_fit output
#'
#' @param x Output from the \code{\link{gp_fit}} function
#' @param time_grid  An optional time grid over which to produce fitted values of the model
#' @param ... Other arguments to plot (not currently implemented)
#'
#' @return A nice plot to illustrate the output of gp_fit as a scatterplot, and the smoothed Gaussian process regression line
#' @seealso \code{\link{load_clim}}, \code{\link{gp_fit}}, \code{\link{plot.climr_fit}}, \code{\link{fit}}
#' @export
#' @import ggplot2
#' @importFrom stats "sd" "optim"
#' @importFrom tibble "tibble"
#' @importFrom viridis "scale_color_viridis"
#'
#' @examples
#' data1 = load_clim('NH')
#' data2 = gp_fit(data1, 'BFGS')
#' plot(data2)
#' @export
plot.climr_gp_fit <- function(x, time_grid = pretty(x$data$year, n = 100), ...) {

  # Create global variables to avoid annoying CRAN notes
  DJF = Dec = `J-D` = Jan = SON = Year = month = pred = quarter = temp = year = NULL

  # Create a nice plot from the output of gp_fit

  # Get the data set out
  df <- x$data

  # For predicted values based on the time_grid
  if(x$method == 'Nelder-Mead') {
    fits <- tibble(time_grid, x$model)
  } else if(x$method == 'BFGS') {
    fits <- tibble(time_grid, x$model)
  } else if(x$method == 'SANN') {
    fits <- tibble(time_grid, x$model)
  } else if(x$method == 'Brent') {
    fits <- tibble(time_grid, x$model)
  }

  # Finally create the plot
  ggplot(df, aes(year, temp)) +
    geom_point(aes(colour = temp)) +
    theme_bw() +
    xlab('Year') +
    ylab('Temperature anomaly') +
    geom_line(data = fits, aes(x = time_grid, y = x$model, colour = x$data$temp)) +
    theme(legend.position = 'None') +
    scale_color_viridis()
}
