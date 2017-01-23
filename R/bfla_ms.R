#' List of known chemicals in the chemical signatures of Bombus flavifrons
#'
#' @description
#' This dataset is comprised of substances and their corresponding retention times that have been identified in a chemical dataset from cephalic labial gland secretion of the bumblebee *Bombus flavifrons* (Dellicour and Lecocq 2013). Samples are the same as in \code{\link{bfla}}.
#'
#' @details
#' Data was extracted from Table 4 of the supporting information from Dellicour and Lecocq (2013). We removed all rows that contain retention times of unknown substances to allow the data beeing used for a conservative estimation of alignment errors of the dataset \code{\link{bfla}}.
#'
#' @format
#' A table containing data on known substances in the checmical signature of 11 samples of *Bombus flavifrons*. Missing values indicate the absence of a certain substance in a sample.
#'
#' @source
#' http://onlinelibrary.wiley.com/store/10.1002/jssc.201300388/asset/supinfo/jssc3437-sup-0001-TableS1.zip?v=1&s=57d5d58273d1d4207e70c72cecd5bba4b1fe95a1
#'
#' @references
#' Dellicour, S. and Lecocq, T. (2013), GCALIGNER 1.0: An alignment program to compute a multiple sample comparison data matrix from large eco-chemical datasets obtained by GC. J. Sep. Science, 36: 3206–3209. doi:10.1002/jssc.201300388
#'
#' @keywords datasets
#' @name bfla_ms
#'
#' @docType data
NULL
