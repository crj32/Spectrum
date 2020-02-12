## kernel functions for spectrum

#' CNN_kernel: fast adaptive density-aware kernel
#'
#' @param mat Matrix: matrix should have samples as columns and rows as features
#' @param NN Numerical value: the number of nearest neighbours to use when calculating local sigma
#' @param NN2 Numerical value: the number of nearest neighbours to use when calculating common nearest neighbours
#'
#' @return A kernel matrix
#' @export
#'
#' @examples
#' CNN_kern <- CNN_kernel(blobs[,1:50])
CNN_kernel <- function(mat, NN = 3, NN2 = 7) {
  n <- ncol(mat)
  ## need N nearest neighbour distance per sample (kn), 
  ## and names of NN2 nearest neigbours (nbs)
  nbs <- list()
  dm <- Rfast::Dist(t(mat))
  dimnames(dm) <- list(colnames(mat), colnames(mat))
  kn <- c()
  ## find kth nearest neighbour for each sample 1...N
  for (i in seq_len(n)) {
    # sort the vector to retrieve the N nearest neighbour and the names of the NN2 nearest neighbours
    sortedvec <- sort.int(dm[i, ])
    # append the NNth nearest neighbour distance
    kn <- c(kn, sortedvec[NN + 1])
    # append the names of the NN2 nearest neighbours
    nbs[[i]] <- names(sortedvec[2:(NN2+1)])
    names(nbs)[[i]] <- names(sortedvec)[1]
  }
  ## make the symmetrical matrix of kth nearest neighbours distances
  sigmamatrix <- kn %o% kn
  ## calculate the kernel using the local statistics (scale and density) of each sample pair
  out <- matrix(nrow = n, ncol = n)  # don't overwrite functions
  # calculate the numerator beforehand
  upper <- -dm^2
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      # shared nearest neighbours between ith and jth sample
      cnns <- length(intersect(nbs[[i]], nbs[[j]]))
      upperval <- upper[i, j]
      # retrieve sigma
      localsigma <- sigmamatrix[i, j]
      # calculate local affinity
      out[i, j] <- exp(upperval / (localsigma * (cnns + 1)))
    }
  }
  # reflect, diag = 1
  out <- pmax(out, t(out), na.rm = TRUE)
  diag(out) <- 1
  ## add colnames back on
  colnames(out) <- colnames(mat)
  row.names(out) <- colnames(mat)
  ## return kernel
  return(out)
}

#' rbfkernel_b: fast self-tuning kernel
#'
#' @param mat Matrix: matrix should have samples as columns and rows as features
#' @param K Numerical value: the number of nearest neighbours to use when calculating local sigma
#' @param sigma Numerical value: a global sigma, usually left to 1 which has no effect
#'
#' @return A kernel matrix
#' @export
#'
#' @examples
#' stsc_kern <- rbfkernel_b(blobs[,1:50])
rbfkernel_b <- function (mat, K = 3, sigma = 1) { # calculate gaussian kernel with local sigma
  n <- ncol(mat)
  NN <- K # nearest neighbours (2-3)
  dm <- Rfast::Dist(t(mat))
  kn <- c() # find kth nearest neighbour for each sample 1...N
  for (i in seq_len(n)) {
    sortedvec <- as.numeric(sort.int(dm[i, ]))
    sortedvec <- sortedvec[!sortedvec == 0]
    kn <- c(kn, sortedvec[NN])
  }
  sigmamatrix <- kn %o% kn # make the symmetrical matrix of kth nearest neighbours distances
  upper <- -dm^2 # calculate the numerator beforehand
  out <- matrix(nrow = n, ncol = n)
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      lowerval <- sigmamatrix[i, j] # retrieve sigma
      upperval <- upper[i, j]
      out[i, j] <- exp(upperval / (lowerval * sigma)) # calculate local affinity
    }
  }
  # reflect, diag = 1
  out <- pmax(out, t(out), na.rm = TRUE)
  diag(out) <- 1
  # add names back on
  colnames(out) <- colnames(mat)
  row.names(out) <- colnames(mat)
  ## return kernel
  return(out)
}

#' ng_kernel: Kernel from the Ng spectral clustering algorithm
#' 
#' This is the kernel from the Ng spectral clustering algorithm. It takes a global
#' sigma which requires tuning for new datasets in most cases. It is possible to use
#' the sigma_finder function to find a sigma for a dataset. Sigma is assumed to be
#' squared already.
#' 
#' @param data Data frame or matrix: with points as columns, features as rows
#' @param sigma Numerical value: a global sigma that controls the drop off in affinity
#' 
#' @references Ng, Andrew Y., Michael I. Jordan, and Yair Weiss. "On spectral clustering: 
#' Analysis and an algorithm." Advances in neural information processing systems. 2002.
#'
#' @return A similarity matrix of the input data
#' @export
#'
#' @examples
#' ng_similarity <- ng_kernel(brain[[1]])
ng_kernel <- function(data,sigma=0.1){
  dm <- Rfast::Dist(t(data))
  K <- exp(-dm^2/sigma)
  #
  colnames(K) <- colnames(data)
  row.names(K) <- colnames(data)
  #
  return(K)
}

#' sigma_finder: heuristic to find sigma for the Ng kernel
#' 
#' This is a heuristic to find the sigma for the kernel from the Ng spectral clustering algorithm. 
#' It returns a global sigma. It uses the mean K nearest neighbour distances of all samples to
#' determine sigma.
#' 
#' @param mat Data frame or matrix: with points as columns, features as rows
#' @param NN Numerical value: the number of nearest neighbours to use (default=3)
#' 
#' @return A global sigma
#' @export
#'
#' @examples
#' sig <- sigma_finder(blobs)
sigma_finder <- function(mat,NN=3){
  # get a global sigma for the Ng kernel using the KNN distances of every sample
  n <- ncol(mat)
  dm <- Rfast::Dist(t(mat))
  dimnames(dm) <- list(colnames(mat), colnames(mat))
  kn <- c()
  ## find kth nearest neighbour for each sample 1...N
  for (i in seq_len(n)) {
    # sort the vector to retrieve the N nearest neighbour distances
    sortedvec <- sort.int(dm[i, ])
    # append the NNth nearest neighbour distance
    kn <- c(kn, sortedvec[NN + 1])
  }
  # get a global sigma
  globalsigma <- mean(kn)*mean(kn) # sigma^2
  #
  message(paste('estimated sigma: ',globalsigma))
  return(globalsigma)
}
