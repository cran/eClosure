# R wrapper for the C++ closedBH implementation.
# The _cpp functions (cBH_cpp, cBH_adjust_cpp) are compiled via Rcpp and called
# here. The next.r function (r.closed) and its helper (r.variant) are built into
# the C++ source; these wrappers only forward the user-facing arguments.

#' Closed Benjamini-Hochberg procedure for FDR control
#'
#' @description
#' Applies the e-Closure version of the Benjamini-Hochberg (BH) procedure.
#' `closedBH()` returns the number of rejections at a given level, while
#' `closedBH.adjust()` returns adjusted p-values, one per hypothesis.
#'
#' @details
#' `closedBH()` returns the size \eqn{r} of the closedBH-significant set. The
#' \eqn{r} hypotheses with the smallest p-values form one such set. This gives
#' the maximum number of hypotheses that can be reported while maintaining FDR
#' control.
#'
#' `closedBH.adjust()` returns the adjusted p-value of each hypothesis: the
#' smallest level \eqn{\alpha} at which that hypothesis is among the rejections
#' of `closedBH()`. Because the rejection set is always the hypotheses with the
#' smallest p-values, a hypothesis of sorted rank \eqn{k} is rejected at level
#' \eqn{\alpha} exactly when `closedBH()` returns at least \eqn{k}; the adjusted
#' p-value is the smallest such \eqn{\alpha}. Consequently
#' `sum(closedBH.adjust(p, cap = FALSE) <= alpha)` equals
#' `closedBH(p, alpha = alpha)`.
#'
#' All of the heavy lifting (the breadth-first traversal over set sizes and the
#' built-in `next.r` recursion) is performed in C++; these wrappers only sort
#' out the user-facing arguments and forward them. `closedBH()` runs in
#' \eqn{O(m + r^2)} time; `closedBH.adjust()` runs in \eqn{O(m^2)} time and
#' \eqn{O(m)} memory.
#'
#' @param p Numeric vector of p-values, one per hypothesis. Values must lie
#'   in \eqn{[0, 1]}.
#' @param alpha Numeric scalar in \eqn{[0, 1]}. The target FDR level. Defaults
#'   to `0.05`.
#' @param cap Logical. If `TRUE` (the default), each adjusted p-value is raised
#'   to be at least its own raw p-value, so that adjustment never makes a
#'   p-value smaller. If `FALSE`, the raw adjusted values are returned.
#'
#' @return
#' `closedBH()` returns a single non-negative integer \eqn{r}: the \eqn{r}
#' hypotheses with the smallest p-values form a valid rejection set, and a
#' value of `0` means no non-empty set can be rejected.
#'
#' `closedBH.adjust()` returns a numeric vector of adjusted p-values, in the
#' same order as the input `p`, each lying in \eqn{[0, 1]}.
#'
#' @examples
#' p <- c(
#'   0.0001, 0.0004, 0.0019, 0.0095, 0.0201,
#'   0.0278, 0.0298, 0.0344, 0.0459, 0.3240,
#'   0.4262, 0.5719, 0.6528, 0.7590, 1.0000
#' )
#' closedBH(p, alpha = 0.05)
#' closedBH(p, alpha = 0.10)
#'
#' padj <- closedBH.adjust(p)
#' sum(padj <= 0.05)             # matches closedBH(p, alpha = 0.05)
#'
#' @references
#' Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate:
#' a practical and powerful approach to multiple testing.
#' \emph{Journal of the Royal Statistical Society: Series B}, 57(1), 289–300.
#'
#' Xu, Z., Solari, A., Fischer, L., de Heide, R., Ramdas, A., & Goeman, J. (2025).
#' Bringing closure to false discovery rate control: A general principle for
#' multiple testing. arXiv preprint arXiv:2509.02517.
#'
#' @seealso
#' [closedBY()] for the analogous Benjamini-Yekutieli procedure.
#' [closedSu()] for the analogous Su procedure.
#' [closedeBH()] for the analogous procedure based on e-values.
#' [p.adjust()] for standard non-simultaneous multiple testing corrections.
#'
#' @export
closedBH = function(p, alpha = 0.05) {
  cBH_cpp(as.numeric(p), as.numeric(alpha))
}

#' @rdname closedBH
#' @export
closedBH.adjust = function(p, cap = TRUE) {
  cBH_adjust_cpp(as.numeric(p), as.logical(cap))
}