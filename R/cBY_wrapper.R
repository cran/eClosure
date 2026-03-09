#' @useDynLib eClosure, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

# R wrappers for the C++ closedBY implementation.
# The two _cpp functions (BY_cpp, cBY_check_cpp) are compiled via Rcpp
# and called here.  All R-level logic (harmonic caching, sorting, bisection)
# lives in this file; only the inner loops are in C++.

# ---------------------------------------------------------------------------
# Harmonic-number cache
# Stored in a private environment so the vector survives across calls but is
# not visible in the global workspace.
# ---------------------------------------------------------------------------
eClosure.env = new.env(parent = emptyenv())
.onLoad = function(libname, pkgname) {
  assign("harmonic", value = 0, envir = eClosure.env)
}

# Retrieve harmonic numbers H(0), H(1), ..., H(m) as a 0-indexed numeric
# vector (harmonic[i+1] in R == harmonic[i] in C++, i.e. H(i)).
# Extends the cached vector on demand so each value is computed at most once.
get_harmonic = function(m) {
  harmonic = get("harmonic", envir = eClosure.env)
  oldm = length(harmonic) - 1          # current maximum n for which H(n) is stored
  if (oldm < m) {
    harmonic[(oldm+2):(m+1)] = harmonic[oldm+1] + cumsum(1/((oldm+1):m))
    assign("harmonic", harmonic, envir = eClosure.env)
  }
  harmonic  # 0-indexed: harmonic[1]=H(0)=0, harmonic[2]=H(1)=1, ...
            # but in C++ we pass this directly and index as harmonic[s] = H(s)
}

#' Closed Benjamini-Yekutieli procedure for simultaneous FDR control
#'
#' @description
#' Applies the closed testing version of the Benjamini-Yekutieli (BY)
#' procedure. The standard BY procedure controls the false discovery rate (FDR)
#' at level \eqn{\alpha} under arbitrary dependence but only provides a single
#' set of rejections. The closed BY procedure provides **simultaneous FDR
#' control**: for every set of hypotheses, it determines whether that set can
#' be reported as discoveries while maintaining FDR control at level
#' \eqn{\alpha}, regardless of which other sets are inspected.
#'
#' @details
#' The closed BY procedure is based on a local e-value for every
#' intersection hypothesis. A set \eqn{R} of hypotheses is closedBY-significant —
#' and therefore a valid simultaneous rejection — if and only if, for every
#' subset \eqn{S \subseteq [m]}, the local e-value exceeds
#' \eqn{|S \cap R|/|R|\alpha}.
#' This guarantees post-hoc FDR control: you may report any closedBY-significant
#' set as your discovery set without inflating the FDR above \eqn{\alpha},
#' even if the choice of set was data-driven.
#'
#' The function has two modes:
#'
#' - **Set-checking mode** (when `set` is supplied): Returns `TRUE` if the
#'   specified set is closedBY-significant (i.e., can be reported as a valid
#'   simultaneous rejection at level \eqn{\alpha}), and `FALSE` otherwise.
#'
#' - **Discovery mode** (when `set = NULL`): Returns the size \eqn{r} of the
#'   largest closedBY-significant set. The \eqn{r} hypotheses with the smallest
#'   p-values always form one such set. This gives the maximum number of
#'   hypotheses that can be reported while maintaining simultaneous FDR control. 
#'   In particular, the set \eqn{R} consisting of the \eqn{r} smallest p-values is
#'   closedBY-significant.
#'
#' Note that closedBY significance is not a monotone property: a set of size \eqn{r}
#' being closedBY-significant does not imply that all smaller sets are as well.
#' The exact algorithm therefore checks all set sizes, while the approximate
#' algorithm (`approximate = TRUE`) uses a faster bisection strategy that may
#' occasionally underestimate the largest significant set.
#'
#' @param p Numeric vector of p-values, one per hypothesis. Values must lie
#'   in \eqn{[0, 1]}; exact zeros are replaced internally by a small positive
#'   constant to avoid numerical issues.
#' @param set Optional subsetting vector for `p` (logical, index, or negative
#'   index), indicating which hypotheses belong to the set to be checked for
#'   closedBY significance. If `NULL` (the default), the function instead returns
#'   the size of the largest closedBY-significant set.
#' @param alpha Numeric scalar in \eqn{[0, 1]}. The target FDR level. Defaults
#'   to `0.05`.
#' @param approximate Logical. If `FALSE` (the default), uses an exact algorithm
#'   that is guaranteed to find the largest closedBY-significant set. If `TRUE`,
#'   uses a faster approximate algorithm based on bisection that may
#'   occasionally return a smaller set. The approximate method is recommended
#'   for exploratory analyses or large inputs where computation time is a
#'   concern.
#'
#' @return
#' - If `set` is supplied: a single logical value. `TRUE` indicates that the
#'   specified set is closedBY-significant and can be reported as a simultaneous
#'   rejection at FDR level \eqn{\alpha}. `FALSE` indicates it cannot.
#'
#' - If `set = NULL`: a single non-negative integer \eqn{r}. The \eqn{r}
#'   hypotheses with the smallest p-values form a valid simultaneous rejection
#'   set. A return value of `0` means no non-empty set can be rejected.
#'
#' @examples
#' set.seed(42)
#' # 20 null hypotheses (p ~ Uniform(0,1)) and 10 non-nulls (p ~ Beta(0.1, 1), smaller on average)
#' p <- c(runif(20), rbeta(10, 0.1, 1))
#'
#' # --- Discovery mode ---
#' # Find the maximum number of simultaneous rejections at FDR level 5%
#' r <- closedBY(p, alpha = 0.05)
#' cat("Largest simultaneous rejection set:", r, "\n")
#'
#' # The r hypotheses with the smallest p-values form a valid discovery set
#' discovery_set <- p <= sort(p)[r]
#' cat("P-values in discovery set:", round(sort(p[discovery_set]), 4), "\n")
#'
#' # --- Set-checking mode ---
#' # Check whether a researcher-defined set is a valid simultaneous rejection
#' candidate_set <- p < 0.01
#' closedBY(p, set = candidate_set, alpha = 0.05)
#'
#' # --- Exact vs. approximate ---
#' r_exact  <- closedBY(p, alpha = 0.05, approximate = FALSE)
#' r_approx <- closedBY(p, alpha = 0.05, approximate = TRUE)
#' cat("Exact:", r_exact, "  Approximate:", r_approx, "\n")
#'
#' @references
#' Benjamini, Y., & Yekutieli, D. (2001). The control of the false discovery
#' rate in multiple testing under dependency.
#' \emph{The Annals of Statistics}, 29(4), 1165–1188.
#'
#' Goeman, J. J., & Solari, A. (2011). Multiple testing for exploratory
#' research. \emph{Statistical Science}, 26(4), 584–597.
#'
#' Xu, Z., Solari, A., Fischer, L., de Heide, R., Ramdas, A., & Goeman, J. (2025).
#' Bringing closure to false discovery rate control: A general principle for
#' multiple testing. arXiv preprint arXiv:2509.02517.
#'
#' @seealso
#' [p.adjust()] for standard p-value-based non-simultaneous multiple testing corrections,
#' including the BY procedure (`method = "BY"`).
#' [closedeBH()] for the analogous procedure based on e-values.
#'
#' @export
closedBY = function(p, set = NULL, alpha = 0.05, approximate = FALSE) {
  # ---------------------------------------------------------------------------
  # closedBY
  #
  # Main user-facing function.  Two modes of operation:
  #
  #   set mode    (set != NULL): returns TRUE/FALSE — is `set` closedBY-significant?
  #   no-set mode (set == NULL): returns an integer — largest closedBY-significant
  #                              set size (0 if none).
  #
  # Arguments:
  #   p           - vector of m p-values (need not be sorted on input)
  #   set         - integer/logical index vector identifying the set of interest,
  #                 or NULL for no-set mode
  #   alpha       - significance level (default 0.05)
  #   approximate - if TRUE, use bisection over r rather than exhaustive search
  #                 (ignored in set mode); may underestimate the true maximum
  # ---------------------------------------------------------------------------
  
  m = length(p)
  harmonic = get_harmonic(m)   # length m+1; harmonic[s+1] = H(s) in R,
                               # equivalently harmonic[s] = H(s) in 0-indexed C++

  # Replace exact zeros to avoid division issues in the C++ ceiling computation
  p[p == 0] = 1e-12

  # ------------------------------------------------------------------
  # SET MODE: check whether a given set is closedBY-significant
  # ------------------------------------------------------------------
  if (!is.null(set)) {
    set = seq_len(m)[set]        # coerce logical/integer to integer indices
    r   = length(set)

    # Arrange p so that p[1:(m-r)] are the outside p-values (sorted decreasing)
    # and p[(m-r+1):m] are the inside p-values (sorted decreasing).
    # This matches the layout expected by cBY_check_cpp.
    p = c(sort(p[-set], decreasing = TRUE),
          sort(p[ set], decreasing = TRUE))

    return(cBY_check_cpp(p, r, harmonic, 0L, alpha)$res)
  }

  # ------------------------------------------------------------------
  # NO-SET MODE: find the largest r such that the top-r set is significant
  # ------------------------------------------------------------------
  # Sort p decreasing once; the "set" tested is always {1, ..., r} in this
  # sorted order, so the p_out block is p[1:(m-r)] and p_in is p[(m-r+1):m].
  p = sort(p, decreasing = TRUE)

  # ---- Non-approximate: exhaustive search from r = m down to 0 ----
  if (!approximate) {
    warm_s = 0L
    for (r in m:1) {
      check = cBY_check_cpp(p       = p,
                            r       = r,
                            harmonic = harmonic,
                            warm_s  = warm_s,
                            alpha   = alpha)
      if (check$res) return(r)     
      warm_s = check$at_s
    }
    return(0L) # check$res is always TRUE for r == 0
  }

  # ---- Approximate: bisection over r ----
  # Start from the BY lower bound (guaranteed achievable), then try to
  # improve via binary search.  Risk: may miss the true maximum if it lies
  # in a non-monotone region (hence "approximate").
  r_BY  = BY_cpp(p = p, harmonic = harmonic, alpha = alpha)
  best_r = r_BY

  interval = c(r_BY + 1L, m)
  warm_s   = 0L

  while (diff(interval) >= 0) {
    r     = sum(interval) %/% 2L
    check = cBY_check_cpp(p        = p,
                          r        = r,
                          harmonic = harmonic,
                          warm_s   = warm_s,
                          alpha    = alpha)
    if (check$res) {              # better r found; search upper half
      best_r   = r
      interval = c(r + 1L, interval[2])
      # Keep warm_s from last *failure* for continuity of the hint
    } else {                      # r failed; search lower half
      interval = c(interval[1], r - 1L)
      warm_s   = check$at_s
    }
  }

  return(best_r)
}
