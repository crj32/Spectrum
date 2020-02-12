#' Spectrum: Fast Adaptive Spectral Clustering for Single and Multi-view Data
#'
#' Spectrum is a self-tuning spectral clustering method for single or multi-view data. Spectrum uses a new type of adaptive
#' density-aware kernel that strengthens connections between points that share common nearest neighbours in the graph. 
#' For integrating multi-view data and reducing noise a tensor product graph data integration and diffusion procedure is used. 
#' Spectrum analyses eigenvector variance or distribution to determine the number of clusters. Spectrum is well suited for a wide 
#' range of data, including both Gaussian and non-Gaussian structures.
#'
#' @param data Data frame or list of data frames: contains the data with points to cluster as columns and rows as features. For multi-view data a list of dataframes is to be supplied with the samples in the same order.
#' @param method Numerical value: 1 = default eigengap method (Gaussian clusters), 2 = multimodality gap method (Gaussian/ non-Gaussian clusters), 3 = no automatic method (see fixk param)
#' @param maxk Numerical value: the maximum number of expected clusters (default=10). This is data dependent, do not set excessively high.
#' @param fixk Numerical value: if we are just performing spectral clustering without automatic selection of K, set this parameter and method to 3
#' @param silent Logical flag: whether to turn off messages
#' @param showres Logical flag: whether to show the results on the screen
#' @param diffusion Logical flag: whether to perform graph diffusion to reduce noise (default=TRUE)
#' @param kerneltype Character string: 'density' (default) = adaptive density aware kernel, 'stsc' = Zelnik-Manor self-tuning kernel
#' @param NN Numerical value: kernel param, the number of nearest neighbours to use sigma parameters (default=3)
#' @param NN2 Numerical value: kernel param, the number of nearest neighbours to use for the common nearest neigbours (default = 7)
#' @param showpca Logical flag: whether to show pca when running on one view
#' @param frac Numerical value: optk search param, fraction to find the last substantial drop (multimodality gap method param)
#' @param thresh Numerical value: optk search param, how many points ahead to keep searching (multimodality gap method param)
#' @param fontsize Numerical value: controls font size of the ggplot2 plots
#' @param dotsize Numerical value: controls the dot size of the ggplot2 plots
#' @param tunekernel Logical flag: whether to tune the kernel, only applies for method 2 (default=FALSE)
#' @param clusteralg Character string: clustering algorithm for eigenvector matrix (GMM or km)
#' @param FASP Logical flag: whether to use Fast Approximate Spectral Clustering (for v. high sample numbers)
#' @param FASPk Numerical value: the number of centroids to compute when doing FASP
#' @param krangemax Numerical value: the maximum K value to iterate towards when running a range of K
#' @param runrange Logical flag: whether to run a range of K or not (default=FALSE), puts Kth results into Kth element of list
#' @param diffusion_iters Numerical value: number of diffusion iterations for the graph (default=4)
#' @param KNNs_p Numerical value: number of KNNs when making KNN graph (default=10, suggested=10-20)
#' @param missing Logical flag: whether to impute missing data in multi-view analysis (default=FALSE)
#'
#' @return A list, containing: 
#' 1) cluster assignments, in the same order as input data columns 
#' 2) eigenvector analysis results (either eigenvalues or dip test statistics)
#' 3) optimal K
#' 4) final similarity matrix
#' 5) eigenvectors and eigenvalues of graph Laplacian
#' @export
#'
#' @examples
#' res <- Spectrum(brain[[1]][,1:50])

Spectrum <- function(data,method=1,silent=FALSE,showres=TRUE,diffusion=TRUE,
                     kerneltype=c('density','stsc'),maxk=10,NN=3,NN2=7,
                     showpca=FALSE,frac=2,thresh=7,
                     fontsize=18,dotsize=3,tunekernel=FALSE,clusteralg='GMM',
                     FASP=FALSE,FASPk=NULL,fixk=NULL,krangemax=10,
                     runrange=FALSE,diffusion_iters=4,KNNs_p=10,missing=FALSE){
  
  ###
  kerneltype <- match.arg(kerneltype)
  
  ### error handling
  if (!inherits(data,'list')){
    datalist <- list(data) # just a single data source
  }else{
    datalist <- data
  }
  if ( length(datalist) > 1 & FASP == TRUE ) {
    stop("Error: FASP method works for only a single view")
  }
  if (is.null(FASPk) == TRUE & FASP == TRUE){
    stop("Error: FASP method requires a number of centroids to compute")
  }
  if (runrange == TRUE & method == 3){
    stop("Error: cannot run a range of K whilst method=3")
  }
  if (is.null(fixk) == TRUE & method == 3){
    stop("Error: need to set the value of K using the fixk parameter for method 3")
  }
  
  ### initial messages
  if (silent == FALSE){
    message('***Spectrum***')
    message(paste('detected views:',length(datalist)))
    message(paste('method:',method))
    message(paste('kernel:',kerneltype))
  }
  
  ### if running FASP calculate the centroids
  if (FASP){
    message('running with FASP data compression')
    cs <- kmeans(t(datalist[[1]]),centers=FASPk)
    csx <- cs$centers # these are the centroids to cluster
    cas <- cs$cluster # now samples are assigned to centroids
    datalist[[1]] <- data.frame(t(csx))
  }
  
  ### calculate individual kernels
  kernellist <- list()
  for (platform in seq(1,length(datalist))){
    ### calculate kernel
    if (silent == FALSE){
      message(paste('calculating similarity matrix',platform))
    }
    if (kerneltype == 'stsc'){
      if (method == 2){
        if (tunekernel){
          NN <- kernfinder_local(datalist[[platform]],maxk=maxk,silent=silent,fontsize=fontsize,
                                 dotsize=dotsize,showres=showres)
        }
      }
      kerneli <- rbfkernel_b(datalist[[platform]],K=NN,sigma=1)
    }else if (kerneltype == 'density'){ 
      if (method == 2){
        if (tunekernel){
          NN <- kernfinder_mine(datalist[[platform]],maxk=maxk,silent=silent,
                                showres=showres,fontsize=fontsize,dotsize=dotsize)
        }
      }
      kerneli <- CNN_kernel(datalist[[platform]],NN=NN,NN2=7)
    }
    if (silent == FALSE){
      message('done.')
    }
    ### save kernel in list
    #colnames(kerneli) <- colnames(datalist[[platform]])
    #row.names(kerneli) <- colnames(datalist[[platform]])
    kernellist[[platform]] <- kerneli
  }
  
  ### impute missing data
  if (missing){
    message('imputing missing data...')
    kernellist <- harmonise_ids(kernellist)
    kernellist <- mean_imputation(kernellist)
    message('done.')
  }
  
  ### fuse and truncate/ make KNN graph (both)
  if (silent == FALSE){
    message('combining similarity matrices if > 1 and making kNN graph...')
  }
  A <- Reduce('+', kernellist) # construct A by linear combination
  
  ### diffusion on TPG using the Shu algorithm
  if (diffusion == TRUE){
    ## get KNN graph
    for (col in seq(1,ncol(A))){
      KNNs <- head(rev(sort(A[,col])),(KNNs_p+1)) # find the KNNs (10 default)
      tokeep <- names(KNNs)
      A[!(names(A[,col])%in%tokeep),col] <- 0
    }
    A <- A/rowSums(A) # row normalise A
    if (silent == FALSE){
      message('done.')
    }
    ## diffusion iterations
    if (silent == FALSE){
      message('diffusing on tensor product graph...')
    }
    Qt <- A
    im <- matrix(ncol=ncol(A),nrow=ncol(A))
    im[is.na(im)] <- 0
    diag(im) <- 1
    for (t in seq(1,diffusion_iters)){ # diffusion_iterations (4 default)
      Qt <- A%*%Qt%*%t(A)+im
    }
    A2 <- t(Qt)
    if (silent == FALSE){
      message('done.')
    }
    #diag(A2) <- 0 # removing self similarity makes no difference
  }else if (diffusion == FALSE){
    # if we are not doing TPG method use simple mean A
    A2 <- A/length(datalist)
  }
  
  ### calculate graph Laplacian of A2 (both)
  if (silent == FALSE){
    message('calculating graph laplacian (L)...')
  }
  dv <- 1/sqrt(rowSums(A2))
  l <- dv * A2 %*% diag(dv)
  
  ### eigengap heuristic
  if (method == 1){
    if (silent == FALSE){
      message('getting eigensystem of L...')
    }
    decomp <- eigen(l)
    if (silent == FALSE){
      message('done.')
      message('examining eigenvalues to select K...')
    }
    evals <- as.numeric(decomp$values)
    diffs <- diff(evals)
    diffs <- diffs[-1]
    optk <- which.max(abs(diffs[1:maxk-1]))+1 # put a cap on the maximum number of clusters
    if (silent == FALSE){
      message(paste('optimal K:',optk))
    }
    nn <- maxk+1
    d <- data.frame(K=seq(1,maxk+1),evals=evals[1:nn])
    if (showres == TRUE){
      plot_egap(d,maxk=maxk,dotsize=dotsize,
                fontsize=fontsize)
    }
  }else if (method == 2){ # multimodality gap heuristic
    if (silent == FALSE){
      message('getting eigensystem of L...')
    }
    decomp <- eigen(l)
    if (silent == FALSE){
      message('done.')
      message('examining eigenvector distributions to select K...')
    }
    xi <- decomp$vectors[,1:(maxk+1)]
    res <- EM_finder(xi,silent=silent)
    d <- data.frame('K'=seq(1,maxk+1),'Z'=res[1:(maxk+1),2])
    if (showres == TRUE){
      plot_multigap(d,maxk=maxk,dotsize=dotsize,
                    fontsize=fontsize)
    }
    optk <- findk(res,maxk=maxk,frac=frac,thresh=thresh)
    if (silent == FALSE){
      message(paste('optimal K:',optk))
    }
  }else if (method == 3){ # no automatic K selection method
    decomp <- eigen(l)
    optk <- fixk
  } 
  
  ## cluster either a range of K or just optimal K
  results <- list()
  if (runrange){
    for (tk in seq(2,krangemax)){
      # get range of k results
      ### select optimal eigenvectors
      xi <- decomp$vectors[,1:tk]
      # normalise rows
      yi <- xi/sqrt(rowSums(xi^2))
      # replace NA values (row = 0) with zeros
      yi[which(!is.finite(yi))] <- 0
      if (clusteralg == 'GMM'){
        ### GMM
        gmm <- ClusterR::GMM(yi, tk, verbose = F, seed_mode = "random_spread") # use random spread          
        pr <- ClusterR::predict_GMM(yi, gmm$centroids, gmm$covariance_matrices, gmm$weights)
        names(pr)[3] <- 'cluster'
        if (0 %in% pr$cluster){
          pr$cluster <- pr$cluster+1
        }
        if (silent == FALSE){
          message('clustered.')
        }
      }else if (clusteralg == 'km'){
        ### k means
        pr <- kmeans(yi, tk)
        if (silent == FALSE){
          message('clustered.')
        }
      }
      #
      if (method != 3){
        if (FASP){
          ### if running FASP assign original samples to centroid clusters
          casn <- cas
          casn <- casn[seq_along(casn)] <- pr$cluster[as.numeric(casn[seq_along(casn)])]
          names(casn) <- names(cas) 
          ### return results
          results[[tk]] <- list('allsample_assignments'=casn,'centroid_assignments'=pr$cluster,
                          'eigenvector_analysis'=d,'K'=tk,'similarity_matrix'=A2,
                          'eigensystem'=decomp)
        }else{
          ### return results
          results[[tk]] <- list('assignments'=pr$cluster,'eigenvector_analysis'=d,
                          'K'=tk,'similarity_matrix'=A2,'eigensystem'=decomp)
        }
      }
    }
  }else{
    # just optimal k/ forced k
    ### select optimal eigenvectors
    xi <- decomp$vectors[,1:optk]
    # normalise rows
    yi <- xi/sqrt(rowSums(xi^2))
    # replace NA values (row = 0) with zeros
    yi[which(!is.finite(yi))] <- 0
    
    if (clusteralg == 'GMM'){
      ### GMM
      if (silent == FALSE){
        message('doing GMM clustering...')
      }
      gmm <- ClusterR::GMM(yi, optk, verbose = F, seed_mode = "random_spread") # use random spread          
      pr <- ClusterR::predict_GMM(yi, gmm$centroids, gmm$covariance_matrices, gmm$weights)
      names(pr)[3] <- 'cluster'
      if (0 %in% pr$cluster){
        pr$cluster <- pr$cluster+1
      }
      if (silent == FALSE){
        message('done.')
      }
    }else if (clusteralg == 'km'){
      ### k means
      if (silent == FALSE){
        message('doing k-means clustering...')
      }
      pr <- kmeans(yi, optk)
      if (silent == FALSE){
        message('done.')
      }
    }
    
    ### display clusters
    if (length(datalist) == 1 && showres == TRUE){ # for one data source only
      if (showpca == TRUE){
        pca(datalist[[1]],labels=as.factor(pr$cluster),axistextsize=fontsize,legendtextsize=fontsize,dotsize=dotsize)
      }
    }
    
    if (method != 3){
      if (FASP){
        ### if running FASP assign original samples to centroid clusters
        casn <- cas
        casn <- casn[seq_along(casn)] <- pr$cluster[as.numeric(casn[seq_along(casn)])]
        names(casn) <- names(cas) 
        ### return results
        results <- list('allsample_assignments'=casn,'centroid_assignments'=pr$cluster,
                        'eigenvector_analysis'=d,'K'=optk,'similarity_matrix'=A2,
                        'eigensystem'=decomp)
      }else{
        ### return results
        results <- list('assignments'=pr$cluster,'eigenvector_analysis'=d,
                        'K'=optk,'similarity_matrix'=A2,'eigensystem'=decomp)
      }
    }else if (method == 3){
      if (FASP){
        ### if running FASP assign original samples to centroid clusters
        casn <- cas
        casn <- casn[seq_along(casn)] <- pr$cluster[as.numeric(casn[seq_along(casn)])]
        names(casn) <- names(cas) 
        ### return results
        results <- list('allsample_assignments'=casn,'centroid_assignments'=pr$cluster,
                        'K'=optk,'similarity_matrix'=A2,
                        'eigensystem'=decomp)
      }else{
        ### return results
        results <- list('assignments'=pr$cluster,
                        'K'=optk,'similarity_matrix'=A2,'eigensystem'=decomp)
      }
    }
  }
  
  if (silent == FALSE){
    message('finished.')
  }
  
  return(results)
}

