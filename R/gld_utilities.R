#' Plot the GLD based in the Lambda values
#'
#' This function allows you to show image with scale.
#' @param L Lambda values of the GLD.
#' @keywords lmoments
#' @examples
#' L = c(0, 2, 0.25, 1.5)
#' suq2.plot.gldPlot(L)
#'
#' @export
#'

suq2.plot.gldPlot <- function(L, param = "fmkl"){
  # Create a Y vector between 0 and 1.
  # Y <- c(0.0001, 0.001*(1:9), 0.01*(1:99), .99+0.001*(1:9), .9999)
  Y = sort(runif(1000))
  # Compute x = Q(Y)
  x <- L[1] + (Y^L[3] - (1-Y)^L[4])/L[2]
  # Compute f(x) = PDF(x)
  y <- L[2]/(L[3]*Y^(L[3]-1) + L[4]*(1-Y)^(L[4]-1))

  if(identical(param, "fmkl")) {
    A = ((Y^L[3])-1)/L[3]
    B = (((1-Y)^L[4])-1)/L[4]
    x <- L[1] + (1/L[2])*(A - B)
    # Compute f(x) = PDF(x)
    y <- L[2]/(Y^(L[3]-1) + (1-Y)^(L[4]-1))
  }
  plot(x, y, type="l")
}

#' Return (x, y) values of the GLD based in the Lambda values
#'
#' This function allows you to return the (x,y) values of the GLD, based in the Lambda values.
#' @param L Lambda values of the GLD.
#' @examples
#' L = c(0, 2, 0.25, 1.5)
#' suq2.plot.gldToPlot(L)
#'
#' @export
#'

suq2.plot.gldToPlot <- function(L, param = "fmkl"){
  # Create a Y vector between 0 and 1.
  # Y <- c(0.0001, 0.001*(1:9), 0.01*(1:99), .99+0.001*(1:9), .9999)
  Y = sort(runif(1000))
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

  result <- list("x" = x, "y" = y)
  return (result)
}

#' Compare if two GLDs belongs to the same distribution
#'
#' This function return D distance from the ks-test that compare if both GLDs belongs to the same distribution.
#' @author Noel Moreno Lemus
#' @param L1 Lambda values of the first GLD.
#' @param L2 Lambda values of the second GLD.
#' @examples
#' L1 = c(0, 2, 0.25, 1.5)
#' L2 = c(0, 2, 0.3, 1.75)
#' D <- suq2.utils.gldComparison(L1, L2)
#'
#' @export
#'

suq2.utils.gldComparison <- function(L1, L2, param = "fmkl", no.test = 1000, len = floor(0.9*no.test), alpha = 0.05) {
    sample1 <- rgl(len * no.test, L1[1], L1[2], L1[3], L1[4], param)
    sample2 <- rgl(len * no.test, L2[1], L2[2], L2[3], L2[4], param)

    # test and sample.fitted to use the same names used in fun.diag.ks.g
    test <- split(sample1, 1:no.test)
    sample.fitted <- split(sample2, 1:no.test)

    result.o <- sum(sapply(1:no.test, function(i, test, sample.fitted)
    ks.gof(test[[i]], sample.fitted[[i]])$p.value, test, sample.fitted) > alpha)
    return(result.o)
}

#' Compare if two GLDs belongs to the same distribution based in its KL.dist
#'
#' This function return the KL.dist of two GLDs.
#' @author Noel Moreno Lemus
#' @param L1 Lambda values of the first GLD.
#' @param L2 Lambda values of the second GLD.
#' @examples
#' L1 = c(0, 2, 0.25, 1.5)
#' L2 = c(0, 2, 0.3, 1.75)
#' D <- suq2.utils.gldComparisonKL(L1, L2)
#'
#' @export
#'

suq2.utils.gldComparisonKL <- function(L1, L2, param = "fmkl", no.test = 1000, len = floor(0.9*no.test), alpha = 0.05) {
  sample1 <- rgl(len * no.test, L1[1], L1[2], L1[3], L1[4], param)
  sample2 <- rgl(len * no.test, L2[1], L2[2], L2[3], L2[4], param)

  # test and sample.fitted to use the same names used in fun.diag.ks.g
  # test <- split(sample1, 1:no.test)
  # sample.fitted <- split(sample2, 1:no.test)

  D = KL.dist(sample1, sample2, k=5)
  a = 0
  count = 0
  for (i in 1:5) {
    if (D[i] != Inf && !is.nan(D[i])){
      a = a + D[i]
      count = count + 1
    }
  }
  a = a/count

  return(a)
}

#' Function to be used as a distance function in clustering algoritms
#'
#' This function return the distances between all the centroids and the elements of the dataset, using a KS-test as a measure of the distance.
#' @author Noel Moreno Lemus
#' @param x dataset.
#' @param centers the centroids to be analized.
#' @examples
#' TODO
#'
#' @export
#'

suq2.utils.distGLDComparison <- function(x, centers) {
  if (ncol(x) != ncol(centers))
    stop(sQuote("x")," and ", sQuote("centers")," must have the same number of columns")
  z <- matrix(0, nrow = nrow(x), ncol = nrow(centers))

  alpha = 0.05
  no.test = 1000
  len = floor(0.9*no.test)

  for (k in 1:nrow(centers)) {
    sample2 <- rgl(len * no.test, 0, 2, centers[k, 1], centers[k, 2], "fmkl")
    sample.fitted <- split(sample2, 1:no.test)
    for (j in 1:nrow(x)) {
      sample1 <- rgl(len * no.test, 0, 2, x[j, 1], x[j, 2], "fmkl")
      test <- split(sample1, 1:no.test)
      result.o <- sum(sapply(1:no.test, function(i, test, sample.fitted)
        ks.gof(test[[i]], sample.fitted[[i]])$p.value, test, sample.fitted) > alpha)

      z[j, k] <- (1000 - result.o)
    }
  }
  z
}

#' Function to be used as a distance function in clustering algoritms
#'
#' This function return the distances between all the centroids and the elements of the dataset, using KL-divergence as a measure of the distance.
#' @author Noel Moreno Lemus
#' @param x dataset.
#' @param centers the centroids to be analized.
#' @examples
#' TODO
#'
#' @export
#'
suq2.utils.distGLDComparisonKL <- function(x, centers) {
  if (ncol(x) != ncol(centers))
    stop(sQuote("x")," and ", sQuote("centers")," must have the same number of columns")
  z <- matrix(0, nrow = nrow(x), ncol = nrow(centers))

  alpha = 0.05
  no.test = 1000
  len = floor(0.9*no.test)

  for (k in 1:nrow(centers)) {
    sample2 <- rgl(no.test, 0, 2, centers[k, 1], centers[k, 2], "fmkl")
    sample.fitted <- split(sample2, 1:no.test)
    for (j in 1:nrow(x)) {
      sample1 <- rgl(no.test, 0, 2, x[j, 1], x[j, 2], "fmkl")
      test <- split(sample1, 1:no.test)
      D = mean(KL.dist(sample1, sample2, k=2))
      z[j, k] <- D
    }
  }
  z
}

suq2.utils.distEuclideanM <- function (x, centers){
  if (ncol(x) != ncol(centers))
    stop(sQuote("x"), " and ", sQuote("centers"), " must have the same number of columns")
  z <- matrix(0, nrow = nrow(x), ncol = nrow(centers))
  for (k in 1:nrow(centers)) {
    z[, k] <- sqrt(colSums((t(x) - centers[k, ])^2))
  }
  z
}

#' Evaluate if the centroid of a cluster is a good representative of the other members of this cluster
#'
#' This function return a list of D distances from the ks-test that compare the centroid of a cluster wiht n members of this cluster.
#' @author Noel Moreno Lemus
#' @param cluster_number The cluster ID we are interesting in.
#' @param n Number of elements of the cluster to analize.
#' @param centroid The lambda values of the centroid of the cluster.
#' @examples
#'
#' Ds <- suq2.clustering.gldClusterComparison(3, 60, c(0, 2, 1.5, 1.3))
#'
#' @export
#'

suq2.clustering.gldClusterComparison <- function(cluster_number, n, centroid) {
  lista = list()
  count = 1;
  i = 1;
  dim_clusters = dim(clusters)
  for(i in 1:dim_clusters[1]){
    for(j in 1:dim_clusters[2]){
      if (clusters[i,j]==cluster_number && count <= n){
        D = gldComparison(centroid, c(0, lambda2[i,j], lambda3[i,j], lambda4[i,j]))
        lista[[length(lista)+1]] = D;
        count = count + 1;
      }
    }
  }
  as.numeric(lista)
}

suq2.clustering.gldClusterComparison1D <- function(cluster_number, n, centroid) {
  lista = list()
  count = 1;
  i = 1;
  dim_clusters = length(clusters)
  for(i in 1:dim_clusters){
      if (clusters[i]==cluster_number && count <= n){
        D = gldComparison(c(0, centroid), c(0, lambda2[i], lambda3[i], lambda4[i]))
        lista[[length(lista)+1]] = D;
        count = count + 1;
      }
  }
  as.numeric(lista)
}

suq2.clustering.gldClusterComparisonKL <- function(cluster_number, n, centroid) {
  library("FNN")
  lista = list()
  count = 1;
  i = 1;
  dim_clusters = dim(clusters)
  for(i in 1:dim_clusters[1]){
    for(j in 1:dim_clusters[2]){
      if (clusters[i,j]==cluster_number && count <= n){
        sample1 = rgl(10000, centroid)
        sample2 = rgl(10000, c(0, lambda2[i,j], lambda3[i,j], lambda4[i,j]))
        D = mean(KL.dist(sample1, sample2))
        #D = gldComparison(centroid, c(0, lambda2[i,j], lambda3[i,j], lambda4[i,j]))
        lista[[length(lista)+1]] = D;
        count = count + 1;
      }
    }
  }
  as.numeric(lista)
}

suq2.clustering.gldClusterComparisonKL1D <- function(cluster_number, n, centroid) {
  library("FNN")
  lista = list()
  count = 1;
  i = 1;
  dim_clusters = length(clusters)
  for(i in 1:dim_clusters){

      if (clusters[i]==cluster_number && count <= n){
        sample1 = rgl(10000, c(0, lambda2[i], centroid))
        sample2 = rgl(10000, c(0, lambda2[i], lambda3[i], lambda4[i]))
        D = mean(KL.dist(sample1, sample2))
        #D = gldComparison(centroid, c(0, lambda2[i,j], lambda3[i,j], lambda4[i,j]))
        lista[[length(lista)+1]] = D;
        count = count + 1;
      }

  }
  as.numeric(lista)
}

suq2.clustering.klDivergenceInClusters <- function(centers, no_clusters){
  kld = array(0, dim = c(no_clusters))
  for (i in 1:no_clusters) {
    D = gldClusterComparisonKL(i, 100, centers[i,])
    kld[i] = mean(D[!is.inf(D)])
  }
  plot(kld, pch = 16, col = c(seq(no_clusters)), main = "KL-divergence inside each Cluster", xlab = "Cluster Number", ylab = "KL-divergence")
  lines(kld)
}

suq2.clustering.DInClusters <- function(centers, no_clusters){
 DD = array(0, dim = c(no_clusters))
  for (i in 1:16) {
    D = gldClusterComparison1D(i, 34, centers[i,])
    DD[i] = mean(D)
  }
  plot(DD, pch = 16, col = c(seq(no_clusters)), main = "KS D distance inside each Cluster", xlab = "Cluster Number", ylab = "D-distance")
  lines(DD)
}


#' Clustering of the GLDs in function of its l2, l3 and l4 values
#'
#' TODO.
#' @author Noel Moreno Lemus
#' @param lambdas A matrix of n x m x 4 of all the lambda values.
#' @param no_clusters Number of clusters.
#'
#' @examples
#'
#' suq2.clustering.gldClustering(lambdas, 10, TRUE)
#'
#' @export
#'
suq2.clustering.gldClustering <- function(lambdas, no_clusters, l234 = TRUE){
  lambda1 = lambdas[,,1]
  lambda2 = lambdas[,,2]
  lambda3 = lambdas[,,3]
  lambda4 = lambdas[,,4]

  dimension = dim(lambdas)
  dimTotal = dimension[1]*dimension[2]

  lambdas34 = array(0, dim = c(dimTotal, 6))
  if (l234){
    x = array(0, dim = c(dimTotal, 3))
  } else {
    x = array(0, dim = c(dimTotal, 2))
  }


  count = 0;
  for(i in 1:dimension[1]){
    for(j in 1:dimension[2]){
      count = count + 1;
      print(count)
      lambdas34[count,1] = i;
      lambdas34[count,2] = j;
      lambdas34[count,3] = lambda1[i, j];
      lambdas34[count,4] = lambda2[i, j];
      lambdas34[count,5] = lambda3[i, j];
      lambdas34[count,6] = lambda4[i, j];
      #x[count,1] = lambda1[i, j];
      if (l234){
        x[count,1] = lambda2[i, j];
        x[count,2] = lambda3[i, j];
        x[count,3] = lambda4[i, j];
      } else {
        x[count,1] = lambda3[i, j];
        x[count,2] = lambda4[i, j];
      }

    }
  }

  cl <- kmeans(x, no_clusters)

  #cl <- dbscan(x, eps = 0.5)
  clusters = array(0, dim = c(dimension[1], dimension[2]))
  for(k in 1:dimTotal){
    clusters[lambdas34[k,1], lambdas34[k,2]] = cl$cluster[k]
  }

  image_display(clusters)

  results = list("clusters" = clusters, "x" = x, "cl" = cl)
  return(results)
}

suq2.clustering.gldClustering1D <- function(lambdas, no_clusters, l234 = TRUE){
  lambda1 = lambdas[,,1]
  lambda2 = lambdas[,,2]
  lambda3 = lambdas[,,3]
  lambda4 = lambdas[,,4]

  dimension = length(lambda1)

  dimTotal = dimension


  if (l234){
    x = array(0, dim = c(dimTotal, 3))
  } else {
    x = array(0, dim = c(dimTotal, 2))
  }


  count = 0;
  for(i in 1:dimension){

      count = count + 1;
      #print(count)
      #x[count,1] = lambda1[i, j];
      if (l234){
        x[count,1] = lambda2[i];
        x[count,2] = lambda3[i];
        x[count,3] = lambda4[i];
      } else {
        x[count,1] = lambda3[i];
        x[count,2] = lambda4[i];
      }


  }

  cl <- kmeans(x, no_clusters)

  #cl <- dbscan(x, eps = 0.5)

  results = list("clusters" = cl, "x" = x)
  return(results)
}

#' Plot the clusters in l3-l4 space.
#'
#' TODO.
#' @author Noel Moreno Lemus
#' @param clusters An n x m matrix wiht the clusters by positions.
#' @param x Array used to create the clusters in function @method gldClustering.
#'
#' @examples
#'
#' suq2.plot.gldClustersL3L4(clusters, x)
#'
#' @export
#'
suq2.plot.gldClustersL3L4 <- function (clusters, x) {
  library(latex2exp)
  plot(lambda3, lambda4, type = "n", xlab = TeX('$\\lambda_{3}'), ylab = TeX('$\\lambda_{4}'), main = TeX('Clusters in $\\lambda_{3}-\\lambda_{4}'))  # setting up coord. system

  no_clusters = max(clusters)

  #legend_list = list()

  for (i in 1:no_clusters) {
    plot.clusters.l3l4(i, col = i)
    #legend_list[[length(legend_list)+1]] = paste("Cluster ", i)
  }
  #legend("right", 95, legend_list, col=c(1:no_clusters), lty=1, cex=0.8)
}

suq2.plot.gldClustersL2L3L4 <- function(clusters, x){
  library("plotly")
  a = data.frame(x[,1], x[,2], x[,3], clusters)

  p <- plot_ly(a, x = a$x...1., y = a$x...2., z = a$x...3., color = a$clusters, marker = list(size = 2, color = clusters, width = 1), mode = 'markers') %>%
    layout(yaxis = list(title = 'lambda_3'),
           xaxis = list(title = 'lambda_2'))
  p
}
