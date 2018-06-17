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

suq2.plot.multigld <-function(L){

  s = dim(L)
  temp = gldToPlot(L[1,])
  df = data.frame(x = temp$x, y = temp$y)
  gp <- ggplot(data = df, aes(x = df$x, y = df$y)) + theme_bw() + geom_line() +
    labs(y="Density", x="x")

  for(i in 2:s[1]) {
    temp <- gldToPlot(L[i,])
    df <- data.frame(x = temp$x, y = temp$y)
    gp <- gp + geom_line(data = df, aes(x = x, y = y))
  }
  gp
}

suq2.plot.by.cluster.1D <- function(cl, cluster_number, n){
  lista = array(0, dim = c(n, 4))
  count = 1;
  i = 1;

  dim_clusters = length(clusters)

  for(i in 1:dim_clusters){
    #for(j in 1:dim_clusters[2]){
      if (clusters[i]==cluster_number && count <= n){
        lista[count, ] = c(0, lambda2[i], lambda3[i], lambda4[i])
        count = count + 1;
      }
    #}
  }

  gld.plot.multigld(lista)
}

suq2.plot.l3l4 <- function(cl, colors = cl$clusters$cluster) {
  df = data.frame(cl$x)
  gp = ggplot(df, aes(x = df$X1, y = df$X2)) + geom_point(col = colors, size = 3) + theme_bw() +
    # geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
    labs(title=TeX('Clusters in $\\lambda_{3}-\\lambda_{4}'), y=TeX('$\\lambda_{4}'), x=TeX('$\\lambda_{3}'))
  gp
}

suq2.plot.clustersByIndex <- function(cl){
  x = sequence(length(cl$clusters$cluster))
  df = data.frame(x, y = cl$clusters$cluster)
  ggplot(df, aes(x = df$x, y = df$y)) + geom_point(col = cl$clusters$cluster, size = 3) +
    theme_bw() + labs(title="Cluster number by dataset index", y="Cluster Number", x="Index")
  # plot(cl$clusters$cluster, col = cl$clusters$cluster, pch = 16, xlab = "Index", ylab = "Cluster Number", main = "Cluster number by dataset index")
}
