#' 8 blob like structures
#'
#' A simulated dataset of 8 Gaussian blobs. Simulated using 
#' the 'clusterlab' CRAN package.
#'
#' @format A data frame with 10 rows and 800 variables
"blobs"

#' A brain cancer dataset
#'
#' A dataset containing The Cancer Genome Atlas expression
#' data. From this publication https://tcga-data.nci.nih.gov/docs/publications/lgggbm_2016/.
#' The first data frame is a 5133X150 RNA-seq data matrix, the second is a 262X150
#' miRNA-seq data matrix, the third is 45X150 protein array data matrix. The data was
#' all pre-normalised then subject to log transform.
#'
#' @format A list of data frames
#' @source \url{https://gdac.broadinstitute.org/}
"brain"

#' Two spirals wrapped around one another
#'
#' Simulated data using the 'mlbench' CRAN package.
#'
#' @format A data frame with 2 rows and 180 variables
"spirals"

#' Three concentric circles
#'
#' Simulated data using the 'clusterSim' CRAN package.
#'
#' @format A data frame with 2 rows and 540 variables
"circles"

#' A list of the blob data as similarity matrices with
#' a missing entry in one
#'
#' Two copies of a simulated dataset of 8 Gaussian blobs in a 
#' list converted to a similarity matrix, but one has a 
#' missing observation.
#'
#' @format A list of two data frames
"missl"

#' A list of the blob data as similarity matrices with
#' a missing entry in one filled with NAs
#'
#' Two copies of a simulated dataset of 8 Gaussian blobs in a 
#' list converted to a similarity matrix, but one has a 
#' missing observation filled with NAs.
#'
#' @format A list of two data frames
"misslfilled"


