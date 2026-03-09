#' eClosure: Multiple testing methods based on the e-Closure principle
#'
#' Provides functions for applying the closed e-Benjamini-Hochberg (eBH)
#' and closed Benjamini-Yekutieli (BY) procedures within the e-Closure
#' framework for simultaneous FDR control in multiple hypothesis testing.
#'
#' @section Main functions:
#' The two core functions are:
#' \itemize{
#'   \item \code{\link{closedeBH}} — applies the closed eBH procedure
#'   \item \code{\link{closedBY}} — applies the closed BY procedure
#' }
#'
#' @author Jelle Goeman
#'
#' @docType package
#' @name eClosure
#' @references
#' Xu, Z., Solari, A., Fischer, L., de Heide, R., Ramdas, A., & Goeman, J. (2025).
#' Bringing closure to false discovery rate control: A general principle for
#' multiple testing. arXiv preprint arXiv:2509.02517.
"_PACKAGE"