#' @useDynLib eClosure, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Closed eBH procedure for simultaneous FDR control
#'
#' @description
#' Applies the closed testing version of the e-BH (e-values Benjamini-Hochberg)
#' procedure. The standard eBH procedure controls the false discovery rate (FDR)
#' at level \eqn{\alpha} but only provides a single set of rejections. The closed
#' eBH procedure provides **simultaneous FDR control**: for every set of
#' hypotheses, it determines whether that set can be reported as discoveries
#' while maintaining FDR control at level \eqn{\alpha}, regardless of which
#' other sets were inspected.
#'
#' @details
#' The closed eBH procedure is based on the concept of **mean consistency**. A
#' set \eqn{R} of hypotheses is mean consistent — and therefore a valid
#' simultaneous rejection — if and only if:
#' \deqn{
#'   \frac{1}{|S|} \sum_{i \in S} e_i \;\geq\; \frac{|S \cap R|}{|R| \cdot \alpha}
#' }
#' is satisfied jointly for all \eqn{S \subseteq [m]}, when \eqn{m} hypotheses are tested. 
#' In practice, this condition guarantees
#' that \eqn{R} is a valid closed-testing rejection, providing post-hoc FDR
#' control: you may report any mean-consistent set as your discovery set without
#' inflating the FDR above \eqn{\alpha}, even if the choice of set was
#' data-driven.
#'
#' The function has two modes:
#'
#' - **Set-checking mode** (when `set` is supplied): Returns `TRUE` if the
#'   specified set is mean consistent (i.e., can be reported as a valid
#'   simultaneous rejection at level \eqn{\alpha}), and `FALSE` otherwise.
#'
#' - **Discovery mode** (when `set = NULL`): Returns the size \eqn{r} of the
#'   largest mean-consistent set. The \eqn{r} hypotheses with the largest
#'   e-values always form one such set. This gives the maximum number of
#'   hypotheses that can be reported while maintaining simultaneous FDR control.
#'
#' Note that mean consistency is not a monotone property: a set of size \eqn{r}
#' being mean consistent does not imply that all smaller sets are as well.
#' The exact algorithm therefore checks all set sizes, while the approximate
#' algorithm (`approximate = TRUE`) uses a faster bisection strategy that may
#' occasionally underestimate the largest consistent set.
#'
#' @param e Numeric vector of e-values, one per hypothesis. E-values must be
#'   non-negative; each is interpreted as the evidence against its null
#'   hypothesis. The e-values should have expectation at most 1 under the null 
#'   hypothesis.
#' @param set Optional subsetting vector for `e` (logical, index or negative index), 
#'   indicating
#'   which hypotheses belong to the set to be checked for mean consistency.
#'   If `NULL` (the default), the function instead returns the size of the
#'   largest mean-consistent set.
#' @param alpha Numeric scalar in \eqn{[0, 1]}. The target FDR level. Defaults
#'   to `0.05`.
#' @param approximate Logical. If `FALSE` (the default), uses an exact algorithm
#'   that is guaranteed to find the largest mean-consistent set. If `TRUE`, uses
#'   a faster approximate algorithm that may occasionally return a smaller set.
#'   The approximate method is recommended for exploratory analyses or large
#'   inputs where computation time is a concern.
#'
#' @return
#' - If `set` is supplied: a single logical value. `TRUE` indicates that the
#'   specified set is mean consistent and can be reported as a simultaneous
#'   rejection at FDR level \eqn{\alpha}. `FALSE` indicates it cannot.
#'
#' - If `set = NULL`: a single non-negative integer \eqn{r}. The \eqn{r}
#'   hypotheses with the largest e-values form a valid simultaneous rejection
#'   set. A return value of `0` means no non-empty set can be rejected.
#'
#' @examples
#' set.seed(42)
#' # 20 null hypotheses (e ~ Exp(1)) and 10 non-nulls (e ~ Exp(0.1), larger on average)
#' e <- c(rexp(20, rate = 1), rexp(10, rate = 0.1))
#'
#' # --- Discovery mode ---
#' # Find the maximum number of simultaneous rejections at FDR level 5%
#' r <- closedeBH(e, alpha = 0.05)
#' cat("Largest simultaneous rejection set:", r, "\n")
#'
#' # The r hypotheses with the largest e-values form a valid discovery set
#' discovery_set <- e >= sort(e, decreasing = TRUE)[r]
#' cat("E-values in discovery set:", round(sort(e[discovery_set], decreasing = TRUE), 2), "\n")
#'
#' # --- Set-checking mode ---
#' # Check whether a researcher-defined set is a valid simultaneous rejection
#' candidate_set <- e > 3
#' closedeBH(e, set = candidate_set, alpha = 0.05)
#'
#' # --- Exact vs. approximate ---
#' r_exact  <- closedeBH(e, alpha = 0.05, approximate = FALSE)
#' r_approx <- closedeBH(e, alpha = 0.05, approximate = TRUE)
#' cat("Exact:", r_exact, "  Approximate:", r_approx, "\n")
#'
#' @references
#' Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate:
#' a practical and powerful approach to multiple testing.
#' \emph{Journal of the Royal Statistical Society: Series B}, 57(1), 289–300.
#'
#' Wang, R., & Ramdas, A. (2022). False discovery rate control with
#' e-values. \emph{Journal of the Royal Statistical Society: Series B},
#' 84(3), 822–852.
#' 
#' Xu, Z., Solari, A., Fischer, L., de Heide, R., Ramdas, A., & Goeman, J. (2025). 
#' Bringing closure to false discovery rate control: A general principle for multiple testing. 
#' arXiv preprint arXiv:2509.02517.
#'
#' @seealso
#' [p.adjust()] for non-simultaneous multiple testing corrections.
#'
#' @export
closedeBH <- function(e, set = NULL, alpha = 0.05, approximate = FALSE) {
  m <- length(e)
  
  if (!is.null(set)) {
    # Check if a specific set is mean consistent
    set_indices <- (1:m)[set]
    r <- length(set_indices)
    
    # Empty set is always mean consistent
    if (r == 0) return(TRUE)
    
    # Sort e-values: first those outside the set, then those inside
    e_in <- sort(e[set_indices])
    e_out <- sort(e[-set_indices])
    
    # Compute cumulative sums
    cum_e <- cumsum(c(0, e_out, e_in))
    
    # Call C++ function
    result <- meanconsistent_wrapper_cpp(cum_e, r, alpha)
    return(result)
    
  } else {
    # Find the largest mean consistent set
    e_sorted <- sort(e)
    
    if (approximate) {
      result <- largestmeanconsistent_approximate_cpp(e_sorted, alpha)
    } else {
      result <- largestmeanconsistent_cpp(e_sorted, alpha)
    }
    
    return(result)
  }
}