#' Plot the GLD based in its lambda values.
#'
#' This function read all the csv that are stored in a directory and create an NxMxS array, where N and M are sparial dimensions
#' and S is the number of simulations on each spatial point.
#' @author Noel Moreno Lemus
#' @param n the directory where the .csv files are stored.
#' @param L a pattern to read the .csv (e.g. pattern="^[h]").
#' @param param (default = "fmkl")
#' @examples
#' suq2.plot.gld(L = c(0, 2, 0.14, 0.14))
#'
#' @export
#'

suq2.plot.gld <- function(n = 1000, L = c(0, 2, 4, 4), param = "fmkl"){
  # Generate an uniform sampling and sort it.
  Y = sort(runif(n))
  # Compute x = Q(Y)
  if (identical(param, "rs")){
    x <- L[1] + (Y^L[3] - (1-Y)^L[4])/L[2]
    # Compute f(x) = PDF(x)
    y <- L[2]/(L[3]*Y^(L[3]-1) + L[4]*(1-Y)^(L[4]-1))
  }
  if(identical(param, "fmkl")) {
    A = ((Y^L[3])-1)/L[3]
    B = (((1-Y)^L[4])-1)/L[4]
    x <- L[1] + (1/L[2])*(A - B)
    # Compute f(x) = PDF(x)
    y <- L[2]/(Y^(L[3]-1) + (1-Y)^(L[4]-1))
  }
  df = data.frame(x, y)
  gp = ggplot(df, aes(df$x, df$y)) + geom_line() + theme_bw() +
    labs(y="Density", x="x")
  gp
}
