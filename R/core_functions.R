### core functions for Spectrum

#' integrate_similarity_matrices: integrate similarity matrices using a tensor product graph
#' linear combination and diffusion technique
#' 
#' Given a list of similarity matrices this function will integrate them running
#' the Shu algorithm, also can reduce noise if the input is a list consisting of
#' a single matrix.
#' 
#' @param kernellist A list of similarity matrices: those to be integrated
#' @param KNNs_p Numerical value: number of nearest neighbours for KNN graph (default=10, suggested=10-20)
#' @param diffusion_iters Numerical value: number of iterations for graph diffusion (default=4, suggested=2-6)
#' @param method Character: either TPG (see reference below) or mean (default=TPG)
#' 
#' @references Shu, Le, and Longin Jan Latecki. "Integration of single-view graphs with 
#' diffusion of tensor product graphs for multi-view spectral clustering." Asian Conference 
#' on Machine Learning. 2016.
#' 
#' @return An integrated similarity matrix
#' @export
#'
#' @examples
#' i_test <- integrate_similarity_matrices(misslfilled,method='mean')
integrate_similarity_matrices <- function(kernellist,KNNs_p=10,diffusion_iters=4,
                                          method='TPG'){
  A <- Reduce('+', kernellist) # construct A by linear combination
  ### diffusion on KNN graph
  if (method=='TPG'){
    ## get KNN graph
    for (col in seq(1,ncol(A))){
      KNNs <- head(rev(sort(A[,col])),(KNNs_p+1)) # find the KNNs (10 default)
      tokeep <- names(KNNs)
      A[!(names(A[,col])%in%tokeep),col] <- 0
    }
    A <- A/rowSums(A) # row normalise A
    ## diffusion iterations
    Qt <- A
    im <- matrix(ncol=ncol(A),nrow=ncol(A))
    im[is.na(im)] <- 0
    diag(im) <- 1
    for (t in seq(1,diffusion_iters)){ # diffusion_iterations (4 default)
      Qt <- A%*%Qt%*%t(A)+im
    }
    A2 <- t(Qt)
  }else if (method=='mean'){
    # if we are not doing TPG method use simple mean A
    A2 <- A/length(kernellist)
  }
  #
  return(A2)
}

#' estimate_k: estimate K using the eigengap or multimodality gap heuristics
#' 
#' This function will try to estimate K given a similarity matrix. Generally the
#' maximum eigengap is preferred, but on some data examining the distribution
#' of the eigenvectors as in the multimodality gap heuristic may be beneficial.
#' 
#' @param A2 Data frame or matrix: a similarity matrix
#' @param maxk Numerical value: maximum number of K to be considered
#' @param showplots Character value: whether to show the plot on the screen
#' 
#' @return A data frame containing the eigenvalues and dip-test statistics of the 
#' eigenvectors of the graph Laplacian
#' @export
#'
#' @examples
#' k_test <- estimate_k(missl[[1]])
estimate_k <- function(A2, maxk=10, showplots=TRUE){
  ### calculate graph Laplacian of A2 (both)
  dv <- 1/sqrt(rowSums(A2))
  l <- dv * A2 %*% diag(dv)
  decomp <- eigen(l)
  ### heuristics
  # mgap
  xi <- decomp$vectors[,1:(maxk+1)]
  res <- EM_finder(xi,silent=TRUE)
  zi <- as.numeric(res[,2])
  d1 <- data.frame('K'=seq(1,maxk+1),'Z'=res[1:(maxk+1),2])
  if (showplots){
    plot_multigap(d1,maxk=maxk,dotsize=3,fontsize=18)
  }
  # egap
  evals <- as.numeric(decomp$values)
  diffs <- diff(evals)
  diffs <- diffs[-1]
  optkb <- which.max(abs(diffs[1:maxk-1]))+1 # put a cap on the maximum number of clusters
  message(paste('egap optimal K:',optkb))
  nn <- maxk+1
  d <- data.frame(K=seq(1,maxk+1),evals=evals[1:nn],z=d1$Z)
  if (showplots){
    plot_egap(d,maxk=maxk,dotsize=3,fontsize=18)
  }
  #
  return(d)
}

#' cluster_similarity: cluster a similarity matrix using the Ng method
#' 
#' This function performs clustering of a similarity matrix following the method
#' of Ng or of Melia. We recommend using the Ng method with GMM to cluster the 
#' eigenvectors instead of k-means.
#' 
#' @param A2 Data frame or matrix: a similarity matrix
#' @param k Numerical value: the number of clusters
#' @param clusteralg Character value: GMM or km clustering algorithm (suggested=GMM)
#' @param specalg Character value: Ng or Melia variant of spectral clustering (default=Ng)
#' 
#' @references Ng, Andrew Y., Michael I. Jordan, and Yair Weiss. "On spectral clustering: 
#' Analysis and an algorithm." Advances in neural information processing systems. 2002.
#' @references Meila, Marina, et al. "Spectral Clustering: a Tutorial for the 2010’s." Handbook 
#' of Cluster Analysis. CRC Press, 2016. 1-23.
#'
#' @return A numeric vector of cluster assignments
#' @export
#'
#' @examples
#' ng_similarity <- cluster_similarity(missl[[1]],k=8)
cluster_similarity <- function(A2, k=k, clusteralg='GMM', specalg='Ng'){
  ## get the eigenvectors of graph
  if (specalg == 'Ng'){
    message('Ng spectral clustering.')
    # get norm graph Laplacian and eigensystem
    dv <- 1/sqrt(rowSums(A2))
    l <- dv * A2 %*% diag(dv)
    decomp <- eigen(l)
    # get k eigenvectors
    xi <- decomp$vectors[,1:k]
    yi <- xi/sqrt(rowSums(xi^2))
    yi[which(!is.finite(yi))] <- 0
  }else if (specalg == 'Melia'){
    message('Melia spectral clustering.')
    # get transition matrix and eigensystem
    tm <- A2/rowSums(A2)
    decomp <- eigen(tm)
    yi <- decomp$vectors[,1:k]
    # added row scaling from Ng
    yi <- yi/sqrt(rowSums(yi^2))
  }
  ## cluster the eigenvectors
  if (clusteralg == 'GMM'){
    message('Gaussian Mixture Modelling clustering.')
    ### GMM
    gmm <- ClusterR::GMM(yi, k, verbose = F, seed_mode = "random_spread") # use random spread          
    pr <- ClusterR::predict_GMM(yi, gmm$centroids, gmm$covariance_matrices, gmm$weights)
    names(pr)[3] <- 'cluster'
    if (0 %in% pr$cluster){
      pr$cluster <- pr$cluster+1
    }
    pr <- pr$cluster
  }else if (clusteralg == 'km'){
    message('K-means clustering.')
    ### k means
    pr <- kmeans(yi, k)
    pr <- pr$cluster
  }
  return(pr)
}

### search multimodality diffs to find best k
findk <- function(res,maxk=maxk,thresh=4,frac=2){
  v <- diff(res[,2])
  v <- v[1:maxk] 
  # parameters for search used in paper
  # thresh <- 4
  # frac <- 2
  for (e in seq(2,length(v))){
    if (e == 2){ # store element
      saved <- v[e]
      index <- 2
    }else{
      currentindex <- e
      xx <- currentindex-index # how far are we ahead of minima
      if (xx >= thresh){
        # if we have a max followed by 4 negs, then stop search
        break
      }
      if (v[e] < saved | v[e] < saved/frac){
        # we have found a new max diff, so replace cache and index
        saved <- v[e]
        index <- e
      }
    }
  }
  optk <- index
  return(optk)
}

### tuner for locally adaptive density aware kernel
kernfinder_mine <- function(data,maxk=10,fontsize=fontsize,silent=silent,
                            showres=showres,dotsize=dotsize){ 
  if (silent == FALSE){
    message('finding optimal NN kernel parameter by examining eigenvector distributions')
  }
  rr <- c()
  for (param in seq(1,10)){
    if (silent == FALSE){
      message(paste('tuning kernel NN parameter:',param))
    }
    kern <- CNN_kernel(data,NN=param,NN2=7)
    kern[which(!is.finite(kern))] <- 0 # deal with possible NaNs
    ## calculate difference between most multimodal eigenvectors and background
    dv <- 1/sqrt(rowSums(kern)) # D = diag(1/sqrt(rowSums(A)))
    l <- dv * kern %*% diag(dv) # L = D%*%A%*%D
    xi <- eigen(l)$vectors
    res <- matrix(nrow=ncol(xi),ncol=2)
    for (ii in seq(1,ncol(xi))){
      r <- diptest::dip.test(xi[,ii], simulate.p.value = FALSE, B = 2000)
      res[ii,1] <- r$p.value
      res[ii,2] <- r$statistic
    }
    ## chose kernel via maximum difference
    diffs <- diff(res[,2])
    diffs <- diffs[-1] # remove difference corresponding to first eigenvector
    tophit <- diffs[1:(maxk+1)][which.min((diffs[1:(maxk+1)]))]
    rr <- c(rr, tophit)
  }
  ## find index that yields lowest diff
  optimalparam <- which.min(rr)
  if (silent == FALSE){
    message(paste('optimal NN:'),optimalparam)
  }
  ## print tuning res
  d <- data.frame(x=seq(1,10),y=rr)
  py <- ggplot2::ggplot(data=d,aes(x=x,y=rr)) + ggplot2::geom_point(colour = 'black', size = dotsize) +
    ggplot2::theme_bw() +
    ggplot2::geom_line() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = fontsize, colour = 'black'),
          axis.text.x = ggplot2::element_text(size = fontsize, colour = 'black'),
          axis.title.x = ggplot2::element_text(size = fontsize),
          axis.title.y = ggplot2::element_text(size = fontsize),
          legend.text = ggplot2::element_text(size = fontsize),
          legend.title = ggplot2::element_text(size = fontsize),
          plot.title = ggplot2::element_text(size = fontsize, colour = 'black', hjust = 0.5),
          panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::ylab('D') +
    ggplot2::xlab('NN') +
    ggplot2::scale_x_continuous(limits=c(1,10),breaks=c(seq(1,10,by=1)))
  if (showres == TRUE){
    print(py)
  }
  ##
  return(optimalparam)
}

### tuner for self tuning kernel
kernfinder_local <- function(data,maxk=maxk,silent=silent,dotsize=dotsize,
                             fontsize=fontsize,showres=showres){ 
  if (silent == FALSE){
    message('finding optimal kernel NN parameter by examining eigenvectors')
  }
  rr <- c()
  for (param in seq(1,10)){
    if (silent == FALSE){
      message(paste('tuning NN parameter:',param))
    }
    kern <- rbfkernel_b(data,K=param,sigma=1) # sigma (0.75)
    ## calculate difference between most multimodal eigenvectors and background
    dv <- 1/sqrt(rowSums(kern)) # D = diag(1/sqrt(rowSums(A)))
    l <- dv * kern %*% diag(dv) # L = D%*%A%*%D
    xi <- eigen(l)$vectors
    res <- matrix(nrow=ncol(xi),ncol=2)
    for (ii in seq(1,ncol(xi))){
      r <- diptest::dip.test(xi[,ii], simulate.p.value = FALSE, B = 2000)
      res[ii,1] <- r$p.value
      res[ii,2] <- r$statistic
    }
    ## chose kernel via maximum difference
    diffs <- diff(res[,2])
    diffs <- diffs[-1] # remove difference corresponding to first eigenvector
    tophit <- diffs[1:(maxk+1)][which.min((diffs[1:(maxk+1)]))]
    rr <- c(rr, tophit)
  }
  ## find index that yields lowest diff
  optimalparam <- which.min(rr)
  if (silent == FALSE){
    message(paste('optimal NN:'),optimalparam)
  }
  ## print tuning res
  d <- data.frame(x=seq(1,10),y=rr)
  py <- ggplot2::ggplot(data=d,aes(x=x,y=rr)) + ggplot2::geom_point(colour = 'black', size = dotsize) +
    ggplot2::theme_bw() +
    ggplot2::geom_line() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = fontsize, colour = 'black'),
                   axis.text.x = ggplot2::element_text(size = fontsize, colour = 'black'),
                   axis.title.x = ggplot2::element_text(size = fontsize),
                   axis.title.y = ggplot2::element_text(size = fontsize),
                   legend.text = ggplot2::element_text(size = fontsize),
                   legend.title = ggplot2::element_text(size = fontsize),
                   plot.title = ggplot2::element_text(size = fontsize, colour = 'black', hjust = 0.5),
                   panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::ylab('D') +
    ggplot2::xlab('NN') +
    ggplot2::scale_x_continuous(limits=c(1,10),breaks=c(seq(1,10,by=1)))
  if (showres == TRUE){
    print(py)
  }
  ## this is the maximum method, better to look at the distribution
  optimalparam <- which.min(rr)
  if (silent == FALSE){
    message(paste('optimal NN:'),optimalparam)
  }
  return(optimalparam)
}

### evaluate multimodality of eigenvectors to get max gap
EM_finder <- function(xi,silent=silent){ # accepts the eigenvector decomposition of L
  if (silent == FALSE){
    message('finding informative eigenvectors...')
  }
  res <- matrix(nrow=ncol(xi),ncol=2)
  for (ii in seq(1,ncol(xi))){
    r <- diptest::dip.test(xi[,ii], simulate.p.value = FALSE, B = 2000)
    res[ii,1] <- r$p.value
    res[ii,2] <- r$statistic
  }
  if (silent == FALSE){
    message('done.')
  }
  return(res)
}

#' harmonise_ids: works on a list of similarity matrices to add entries of NA where
#' there are missing observations between views
#' 
#' Simply adds a column and row of NA with the missing ID for data imputation. The
#' similarity matrix requires row and column IDs present for this to work.
#' 
#' @param l A list of similarity matrices: those to be harmonised.
#' 
#' @return A list of harmonised similarity matrices.
#' @export
#'
#' @examples
#' h_test <- harmonise_ids(missl) 
harmonise_ids <- function(l){
  allids <- c()
  for (view in seq(1,length(l))){
    allids <- c(allids,colnames(l[[view]]))
  }
  allids <- unique(allids)
  l2 <- list()
  for (view in seq(1,length(l))){
    l2[[view]] <- setdiff(allids,colnames(l[[view]]))
  }
  ## add ids to each platform, get in correct order
  for (view in seq(1,length(l))){
    namevector <- l2[[view]]
    if (length(namevector) == 0){
      namevector <- c()
    }
    # NXN so require column + row to be added
    ln <- matrix(ncol=length(namevector),nrow=nrow(l[[view]]))
    row.names(ln) <- row.names(l[[view]])
    colnames(ln) <- namevector
    l[[view]] <- cbind(l[[view]],ln)
    #
    ln2 <- matrix(nrow=length(namevector),ncol=ncol(l[[view]]))
    row.names(ln2) <- namevector
    colnames(ln2) <- c(row.names(ln),namevector)
    l[[view]] <- rbind(l[[view]],ln2)
    # sort both
    l[[view]] <- l[[view]][allids,allids]
  }
  return(l)
}

#' mean_imputation: mean imputation function for multi-view spectral clustering
#' with missing data
#' 
#' Works on a list of similarity matrices to impute missing values using the mean 
#' from the other views.
#' 
#' @param l A list of data frames: all those to be included in the imputation. 
#' 
#' @return A list of completed data frames.
#' @export
#'
#' @examples
#' m_test <- mean_imputation(misslfilled)
mean_imputation <- function(l){
  N <- length(l)
  ## get missing indices of matrix
  missing <- c()
  for (i in seq(1,length(l))){
    if (any(is.na(l[[i]]))){
      missing <- c(missing,i)
    }
  }
  ## for N
  l <- lapply(l, function(x) {if(any(class(x)=="data.frame")) as.matrix(x) else x})
  mn <- seq(1,N)
  for (mm in missing){
    # make a list of everything but missing matrix and take mean matrix
    X <- l[-mm]
    Y <- do.call(cbind, X)
    Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
    Z <- apply(Y, c(1, 2), mean, na.rm = TRUE)
    # replace NA with mean
    l[[mm]][is.na(l[[mm]])] <- Z[is.na(l[[mm]])]
  }
  if (is.null(missing)){
    message('sorry, nothing to impute.')
  }else{
    message('imputed.')
  }
  return(l)
}