#' @useDynLib eClosure, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

# Compute ell_alpha = -W_{-1}(-alpha/e), the Lambert W correction factor.
# su_cpp handles this internally; closedSu must apply it before calling
# find_largest_r_cpp and cSu_check_cpp (which work with a plain threshold).
.ell_alpha <- function(alpha) {
  # ell_alpha(0) = +Inf by the limit W_{-1}(0-) -> -Inf; guard explicitly.
  if (alpha <= 0) return(Inf)
  # W_{-1}(x) via Newton iterations
  x  <- -alpha / exp(1)
  lx <- log(-x)
  w  <- lx - log(-lx)
  for (i in seq_len(8)) {
    ew    <- exp(w)
    denom <- ew * (w + 1)
    if (abs(denom) < 1e-300) break
    w_new <- w - (w * ew - x) / denom
    if (abs(w_new - w) < 1e-15 * abs(w_new)) break
    w <- w_new
  }
  -w   # ell_alpha = -W_{-1}(-alpha/e)
}

#' Closed Su procedure for simultaneous FDR control
#'
#' @description
#' Applies the closed testing improvement of the Su (2018) procedure. The
#' standard Su procedure controls the false discovery rate (FDR) at level
#' \eqn{\alpha} under the PRDN assumption but only provides a single set of
#' rejections. The closed Su procedure provides simultaneous FDR control:
#' for every set of hypotheses, it determines whether that set can be reported
#' as discoveries while maintaining FDR control at level \eqn{\alpha},
#' regardless of which other sets are inspected.
#'
#' @details
#' A set \eqn{R} of hypotheses is closed-Su-significant — and therefore a
#' valid simultaneous rejection — if and only if, for every
#' \eqn{1 \le u \le |R|} and \eqn{0 \le v \le m - |R|}, the combined sorted
#' vector \eqn{t_1 \ge \cdots \ge t_{u+v}} of the \eqn{u} largest p-values
#' in \eqn{R} and the \eqn{v} largest p-values outside \eqn{R} satisfies
#' \deqn{
#'   \min_i \left\{ t_i \cdot (u+v) - c \cdot (u+v-i+1) \right\} \le 0,
#'   \quad c = \frac{|R|\alpha}{u \cdot \ell_\alpha},
#' }
#' where \eqn{\ell_\alpha = -W_{-1}(-\alpha/e)} is the Lambert W correction
#' factor (Xu et al., 2025, Section 6.2).
#'
#' The function has two modes:
#'
#' - **Set-checking mode** (when `set` is supplied): Returns `TRUE` if the
#'   specified set is closed-Su-significant, and `FALSE` otherwise.
#'
#' - **Discovery mode** (when `set = NULL`): Returns the size \eqn{r} of the
#'   largest closed-Su-significant set consisting of the \eqn{r} smallest
#'   p-values.
#'
#' @param p Numeric vector of p-values, one per hypothesis.
#' @param set Optional subsetting vector for `p` (logical, positive index, or
#'   negative index). If `NULL` (the default), returns the size of the largest
#'   closed-Su-significant set.
#' @param alpha Numeric scalar in \eqn{[0, 1]}. FDR level. Default `0.05`.
#' @param approximate Logical. If `FALSE` (the default), uses an exact algorithm
#'   that checks every set size from largest to smallest and is guaranteed to
#'   find the largest closed-Su-significant set. If `TRUE`, uses a faster
#'   bisection strategy seeded by the Su lower bound (BH at
#'   \eqn{\alpha/\ell_\alpha}). The approximate method may occasionally
#'   underestimate the largest significant set, but is recommended for
#'   exploratory analyses or large inputs where computation time is a concern.
#'
#' @return
#' - If `set` is supplied: a single logical (`TRUE`/`FALSE`).
#' - If `set = NULL`: a single non-negative integer \eqn{r} (0 = no rejection).
#'
#' @examples
#' set.seed(42)
#' p <- c(runif(20), rbeta(10, 0.1, 1))
#'
#' # Discovery mode
#' r <- closedSu(p, alpha = 0.05)
#'
#' # Set-checking mode
#' closedSu(p, set = p < 0.01, alpha = 0.05)
#'
#' @references
#' Su, W. J. (2018). The FDR-Linking Theorem. arXiv:1812.08965.
#'
#' Xu, Z., Solari, A., Fischer, L., de Heide, R., Ramdas, A., & Goeman, J.
#' (2025). Bringing closure to false discovery rate control.
#' arXiv:2509.02517.
#'
#' @seealso [closedBY()], [closedeBH()]
#'
#' @export
closedSu <- function(p, set = NULL, alpha = 0.05, approximate = FALSE) {
  
  m        <- length(p)
  alpha_su <- alpha / .ell_alpha(alpha)
  
  # ------------------------------------------------------------------
  # SET MODE  (approximate has no effect here)
  # ------------------------------------------------------------------
  if (!is.null(set)) {
    set_indices <- seq_len(m)[set]
    r           <- length(set_indices)
    if (r == 0) return(TRUE)
    
    p_combined <- c(sort(p[-set_indices], decreasing = TRUE),
                    sort(p[ set_indices], decreasing = TRUE))
    return(cSu_check_cpp(p_combined, r, alpha_su))
  }
  
  # ------------------------------------------------------------------
  # DISCOVERY MODE
  # ------------------------------------------------------------------
  p_sorted <- sort(p, decreasing = TRUE)
  
  if (approximate) {
    find_largest_r_approximate_cpp(p_sorted, alpha_su)
  } else {
    find_largest_r_cpp(p_sorted, alpha_su)
  }
}