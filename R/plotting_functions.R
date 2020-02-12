## plotting functions for spectrum

if(getRversion() >= "2.15.1")  utils::globalVariables(c('K','PC1','PC2','X1','X2','Z','evals','x'))

### plot eigengap
plot_egap <- function(d,fontsize=22,width=26,maxk=maxk,dotsize=dotsize){
  py <- ggplot2::ggplot(data=d,aes(x=K,y=evals)) + ggplot2::geom_point(colour = 'grey', size = dotsize) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = fontsize, colour = 'black'),
                   axis.text.x = ggplot2::element_text(size = fontsize, colour = 'black'),
                   axis.title.x = ggplot2::element_text(size = fontsize),
                   axis.title.y = ggplot2::element_text(size = fontsize),
                   legend.text = ggplot2::element_text(size = fontsize),
                   legend.title = ggplot2::element_text(size = fontsize),
                   plot.title = ggplot2::element_text(size = fontsize, colour = 'black', hjust = 0.5),
                   panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::ylab('Eigenvalue') +
    ggplot2::xlab('Eigenvector') +
    ggplot2::scale_x_continuous(limits=c(1,maxk+1),breaks=c(seq(1,maxk+1,by=1)))
  print(py)
}

#' pca: A pca function
#'
#' @param mydata Data frame or matrix: matrix or data frame with samples as columns, features as rows
#' @param labels Factor: to label the plot with colours
#' @param axistextsize Numerical value: axis text size
#' @param legendtextsize Numerical value: legend text size
#' @param dotsize Numerical value: dot size
#' 
#' @return A pca plot object
#' @export
#'
#' @examples
#' ex_pca <- pca(blobs[,1:50])
pca <- function(mydata, labels=FALSE, dotsize = 3, axistextsize = 18, legendtextsize = 18){
  pca1 = prcomp(t(mydata))
  scores <- data.frame(pca1$x) # PC score matrix
  p <- ggplot2::ggplot(data = scores, aes(x = PC1, y = PC2) ) + ggplot2::geom_point(aes(colour = factor(labels)), size = dotsize) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = axistextsize, colour = 'black'),
                   axis.text.x = ggplot2::element_text(size = axistextsize, colour = 'black'),
                   axis.title.x = ggplot2::element_text(size = axistextsize),
                   axis.title.y = ggplot2::element_text(size = axistextsize),
                   legend.title = ggplot2::element_text(size = legendtextsize),
                   legend.text = ggplot2::element_text(size = legendtextsize)) 
  print(p)
}

### plot eigengap
plot_multigap <- function(d,fontsize=fontsize,width=26,maxk=maxk,
                          dotsize=dotsize){
  py <- ggplot2::ggplot(data=d,aes(x=K,y=Z)) + ggplot2::geom_point(colour = 'grey', size=dotsize) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = fontsize, colour = 'black'),
                   axis.text.x = ggplot2::element_text(size = fontsize, colour = 'black'),
                   axis.title.x = ggplot2::element_text(size = fontsize),
                   axis.title.y = ggplot2::element_text(size = fontsize),
                   legend.text = ggplot2::element_text(size = fontsize),
                   legend.title = ggplot2::element_text(size = fontsize),
                   plot.title = ggplot2::element_text(size = fontsize, colour = 'black', hjust = 0.5),
                   panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::ylab('Z') +
    ggplot2::xlab('Eigenvector') +
    ggplot2::scale_x_continuous(limits=c(1,maxk+1),breaks=c(seq(1,maxk+1,by=1)))
  print(py)
}

#' kernel_pca: A kernel pca function
#'
#' @param datam Dataframe or matrix: a data frame with samples as columns, rows as features, or a kernel matrix
#' @param labels Factor: to label the plot with colours
#' @param axistextsize Numerical value: axis text size
#' @param legendtextsize Numerical value: legend text size
#' @param dotsize Numerical value: dot size
#' @param similarity Logical flag: whether the input is a similarity matrix or not
#' 
#' @return A kernel PCA plot
#' @export
#'
#' @examples
#' ex_kernel_pca <- kernel_pca(blobs[,1:50], similarity=FALSE)
kernel_pca <- function(datam, labels = FALSE, axistextsize = 18, legendtextsize = 18, dotsize = 3,
                       similarity = TRUE){
  if (similarity == FALSE){
    km <- CNN_kernel(datam)
  }else{
    km <- datam
  }
  # input a NXN kernel
  m <- nrow(km)
  # center kernel matrix
  kc <- t(t(km - colSums(km)/m) - rowSums(km)/m) + sum(km)/m^2
  # compute eigenvectors
  res <- eigen(kc/m,symmetric=TRUE)
  # multiply each eigenvector by the square root of its eigenvalue
  features <- m
  ret <- suppressWarnings(t(t(res$vectors[,1:features])/sqrt(res$values[1:features])))
  # do plot
  scores <- data.frame(ret)
  colnames(scores)[1:2] <- c('PC1','PC2')
  p1 <- ggplot2::ggplot(data = scores, aes(x = PC1, y = PC2) ) + ggplot2::geom_point(aes(colour = factor(labels)),size=dotsize) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(size = axistextsize, colour = 'black'),
          axis.text.x = ggplot2::element_text(size = axistextsize, colour = 'black'),
          axis.title.x = ggplot2::element_text(size = axistextsize),
          axis.title.y = ggplot2::element_text(size = axistextsize),
          legend.title = ggplot2::element_text(size = legendtextsize),
          legend.text = ggplot2::element_text(size = legendtextsize)) 
  return(p1)
}
