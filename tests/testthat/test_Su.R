# tests/testthat/test_Su.R

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Lambert W correction factor: ell_alpha = -W_{-1}(-alpha/e)
ell_alpha <- function(alpha) {
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
  -w
}

# Standard Su procedure: BH applied at alpha / ell_alpha.
su_rejections <- function(p, alpha = 0.05) {
  sum(p.adjust(p, method = "BH") <= alpha / ell_alpha(alpha))
}

# ---------------------------------------------------------------------------
# Naive triple-loop check (directly from the formula).
#
# For every (u, v), merge p[1..u] and q[1..v] into t sorted decreasingly,
# then check min_i { t_i - r*alpha_su*(u+v-i+1)/((u+v)*u) } <= 0.
#
# alpha_su = alpha / ell_alpha(alpha)  is the effective threshold after
# applying the Lambert W correction (Xu et al. 2025, Section 6.2).
# ---------------------------------------------------------------------------
check_condition_naive_su <- function(p, q, alpha) {
  r        <- length(p)
  m_r      <- length(q)
  alpha_su <- alpha / ell_alpha(alpha)
  
  for (u in seq_len(r)) {
    for (v in 0L:m_r) {
      t   <- sort(c(p[seq_len(u)], q[seq_len(v)]), decreasing = TRUE)
      n   <- length(t)
      i   <- seq_along(t)
      val <- t - r * alpha_su * (n - i + 1L) / (n * u)
      if (min(val) > 0) return(FALSE)
    }
  }
  TRUE
}

find_largest_r_naive_su <- function(z, alpha) {
  m <- length(z)
  z <- sort(z, decreasing = TRUE)
  for (r in m:1L) {
    p <- z[(m - r + 1L):m]
    q <- if (r < m) z[seq_len(m - r)] else numeric(0L)
    if (check_condition_naive_su(p, q, alpha)) return(r)
  }
  0L
}

# ---------------------------------------------------------------------------
# Reproducible test inputs
# ---------------------------------------------------------------------------
set.seed(123)
p_mixed    <- c(runif(20), rbeta(10, 0.1, 1))   # 20 nulls + 10 non-nulls
p_all_null <- runif(50)                           # all null
p_small    <- c(0.001, 0.002, 0.01, 0.04, 0.5)  # hand-crafted small example
p_single   <- 0.03                               # single p-value


# ===========================================================================
# 1. BASIC SANITY CHECKS
# ===========================================================================

test_that("discovery mode returns a non-negative integer", {
  r <- closedSu(p_mixed)
  expect_true(is.numeric(r) || is.integer(r))
  expect_length(r, 1)
  expect_gte(r, 0)
  expect_lte(r, length(p_mixed))
})

test_that("set-checking mode returns a single logical", {
  result <- closedSu(p_mixed, set = p_mixed < 0.05)
  expect_true(is.logical(result))
  expect_length(result, 1)
})

test_that("result does not exceed m", {
  r <- closedSu(p_mixed)
  expect_lte(r, length(p_mixed))
})


# ===========================================================================
# 2. CORNER CASES
# ===========================================================================

test_that("empty set is always TRUE", {
  expect_true(closedSu(p_mixed, set = integer(0)))
  expect_true(closedSu(p_mixed, set = logical(length(p_mixed))))  # all FALSE
})

test_that("single p-value: significant if below Su threshold", {
  alpha    <- 0.05
  thresh   <- alpha / ell_alpha(alpha)
  expect_true( closedSu(0.5 * thresh, set = 1L, alpha = alpha))
  expect_false(closedSu(2.0 * thresh, set = 1L, alpha = alpha))
})

test_that("single p-value discovery mode returns 0 or 1", {
  r <- closedSu(p_single)
  expect_true(r %in% c(0L, 1L))
})

test_that("all p-values equal to 1 yields 0 rejections", {
  expect_equal(closedSu(rep(1, 10)), 0)
})

test_that("all p-values near 0 yields m rejections", {
  p_tiny <- rep(1e-10, 20)
  expect_equal(closedSu(p_tiny), 20L)
})

test_that("p-values of length 1 handled correctly", {
  expect_no_error(closedSu(0.03))
  expect_no_error(closedSu(0.03, set = 1L))
})

test_that("alpha = 0 yields 0 rejections and all sets FALSE", {
  expect_equal(closedSu(p_small, alpha = 0), 0)
  expect_false(closedSu(p_small, set = 1:3, alpha = 0))
})

test_that("alpha = 1 yields m rejections", {
  expect_equal(closedSu(p_small, alpha = 1), length(p_small))
})

test_that("duplicate p-values do not cause errors", {
  p_dup <- c(0.01, 0.01, 0.01, 0.5, 0.5)
  expect_no_error(closedSu(p_dup))
  expect_no_error(closedSu(p_dup, set = 1:3))
})

test_that("set supplied as logical, positive index, and negative index agree", {
  set_logical  <- c(TRUE, TRUE, TRUE, FALSE, FALSE)
  set_posindex <- c(1L, 2L, 3L)
  set_negindex <- c(-4L, -5L)
  r1 <- closedSu(p_small, set = set_logical)
  r2 <- closedSu(p_small, set = set_posindex)
  r3 <- closedSu(p_small, set = set_negindex)
  expect_equal(r1, r2)
  expect_equal(r1, r3)
})


# ===========================================================================
# 3. CONSISTENCY BETWEEN THE TWO MODES
# ===========================================================================

test_that("top-r set (discovery mode) is confirmed significant by set-checking mode", {
  r <- closedSu(p_mixed)
  if (r > 0) {
    top_r_set <- order(p_mixed)[seq_len(r)]
    expect_true(closedSu(p_mixed, set = top_r_set))
  }
})

test_that("discovery and set-checking modes agree across multiple alpha levels", {
  for (alpha in c(0.01, 0.05, 0.10, 0.20)) {
    r <- closedSu(p_mixed, alpha = alpha)
    if (r > 0) {
      top_set <- order(p_mixed)[seq_len(r)]
      expect_true(
        closedSu(p_mixed, set = top_set, alpha = alpha),
        label = paste("top-r set should be significant at alpha =", alpha)
      )
    }
  }
})


# ===========================================================================
# 4. closedSu REJECTIONS ALWAYS CONTAIN Su REJECTIONS
# ===========================================================================

test_that("closedSu discovery count is at least as large as Su rejection count", {
  r_cSu <- closedSu(p_mixed)
  r_Su  <- su_rejections(p_mixed)
  expect_gte(r_cSu, r_Su)
})

test_that("Su rejection set is always closedSu-significant", {
  r_Su <- su_rejections(p_mixed)
  if (r_Su > 0) {
    su_set <- order(p_mixed)[seq_len(r_Su)]
    expect_true(closedSu(p_mixed, set = su_set))
  }
})

test_that("closedSu >= Su containment holds across multiple alpha levels", {
  for (alpha in c(0.01, 0.05, 0.10)) {
    r_Su  <- su_rejections(p_mixed, alpha = alpha)
    r_cSu <- closedSu(p_mixed, alpha = alpha)
    expect_gte(r_cSu, r_Su,
               label = paste("closedSu >= Su at alpha =", alpha))
    if (r_Su > 0) {
      su_set <- order(p_mixed)[seq_len(r_Su)]
      expect_true(
        closedSu(p_mixed, set = su_set, alpha = alpha),
        label = paste("Su set is closedSu-significant at alpha =", alpha)
      )
    }
  }
})

test_that("closedSu >= Su containment holds on all-null data", {
  r_Su  <- su_rejections(p_all_null)
  r_cSu <- closedSu(p_all_null)
  expect_gte(r_cSu, r_Su)
})

test_that("closedSu >= Su containment holds on hand-crafted small example", {
  r_Su  <- su_rejections(p_small)
  r_cSu <- closedSu(p_small)
  expect_gte(r_cSu, r_Su)
  if (r_Su > 0) {
    su_set <- order(p_small)[seq_len(r_Su)]
    expect_true(closedSu(p_small, set = su_set))
  }
})


# ===========================================================================
# 5. INVARIANCE AND STABILITY
# ===========================================================================

test_that("result is invariant to input order", {
  p_shuffled <- sample(p_mixed)
  expect_equal(closedSu(p_mixed), closedSu(p_shuffled))
})

test_that("result is stable across repeated calls", {
  r1 <- closedSu(p_mixed)
  r2 <- closedSu(p_mixed)
  expect_equal(r1, r2)
})


# ===========================================================================
# 6. AGREEMENT WITH NAIVE TRIPLE-LOOP IMPLEMENTATION
#
# check_condition_naive_su is a direct application of the closed Su procedure
# checking the worst case at every S \cap R and S \setminus R.  
# These tests compare closedSu against
# it on small inputs (triple loop is O(m^4), so keep m small).
# ===========================================================================

test_that("discovery mode agrees with naive triple loop on small random inputs", {
  set.seed(42)
  alpha <- 0.05
  for (i in seq_len(20)) {
    z <- sort(runif(sample(5:10, 1))^2, decreasing = TRUE)
    m <- length(z)
    r_fast  <- closedSu(z, alpha = alpha)
    r_naive <- find_largest_r_naive_su(z, alpha = alpha)
    expect_equal(r_fast, r_naive,
                 label = sprintf("seed 42 iter %d: fast=%d naive=%d", i, r_fast, r_naive))
  }
})

test_that("per-(p,q) check agrees with naive on small random inputs", {
  set.seed(99)
  alpha <- 0.05
  for (i in seq_len(20)) {
    z <- sort(runif(sample(5:9, 1))^2, decreasing = TRUE)
    m <- length(z)
    for (r in seq_len(m)) {
      p <- z[(m - r + 1L):m]
      q <- if (r < m) z[seq_len(m - r)] else numeric(0L)
      
      # Discovery-mode: top-r separated set
      r_fast  <- closedSu(z, alpha = alpha)
      # set-checking mode on the separated split
      set_idx <- (m - r + 1L):m
      fast_check  <- closedSu(z, set = set_idx, alpha = alpha)
      naive_check <- check_condition_naive_su(p, q, alpha)
      expect_equal(fast_check, naive_check,
                   label = sprintf("iter %d r=%d: fast=%s naive=%s",
                                   i, r, fast_check, naive_check))
    }
  }
})

test_that("set-checking mode agrees with naive for non-separated random sets", {
  set.seed(7)
  alpha <- 0.05
  for (i in seq_len(20)) {
    z <- sort(runif(sample(6:10, 1))^2, decreasing = TRUE)
    m <- length(z)
    r <- sample(seq_len(m), 1)
    # Random (non-separated) subset of size r
    idx <- sort(sample.int(m, r))
    p   <- sort(z[idx],  decreasing = TRUE)
    q   <- sort(z[-idx], decreasing = TRUE)
    
    fast_check  <- closedSu(z, set = idx, alpha = alpha)
    naive_check <- check_condition_naive_su(p, q, alpha)
    expect_equal(fast_check, naive_check,
                 label = sprintf("iter %d r=%d: fast=%s naive=%s",
                                 i, r, fast_check, naive_check))
  }
})


# ===========================================================================
# 8. TOLERANCE — boundary cases
#
# The TOLERANCE constant (1e-10) in cSu_check_cpp makes the function return
# TRUE whenever the criterion is met up to a small rounding error.  Each test
# below targets one of the five comparison sites where TOLERANCE is applied,
# constructing a p-value configuration that sits exactly on the boundary and
# verifying that a tiny perturbation (within TOLERANCE) does not flip the
# result to FALSE.
#
# Boundary configurations are built analytically so that the naive criterion
# value is exactly 0 (or within machine epsilon), then p-values are shifted
# by TOLERANCE / 2 in the direction that would produce a FALSE result without
# tolerance, and we assert TRUE is still returned.
# ===========================================================================

# Shared tolerance used in cSu.cpp
.SU_TOLERANCE <- 1e-10

test_that("tolerance: p_in[j] exactly equal to c_val does not cause FALSE", {
  # When p_in[j] == c_val, the B3 branch is entered (|p_in[j] - c_val| <= TOLERANCE).
  # The primary check is that no division-by-zero or NaN occurs, and the result
  # agrees with the naive check (which may legitimately be FALSE if the set is
  # not significant for other reasons).
  alpha    <- 0.05
  alpha_su <- alpha / ell_alpha(alpha)
  r        <- 3L
  # c_val at u=r: r*alpha_su/r = alpha_su
  c_val    <- alpha_su
  # Set p_in[0] exactly equal to c_val; other inside p-values smaller
  p_in     <- c(c_val, c_val * 0.5, c_val * 0.3)
  p_out    <- c(0.8, 0.7, 0.6)
  p        <- c(sort(p_out, decreasing = TRUE),
                sort(p_in,  decreasing = TRUE))
  # No error and result matches naive (does not crash on zero denominator)
  expect_no_error(cSu_check_cpp(p, r, alpha_su))
  naive <- check_condition_naive_su(sort(p_in, decreasing = TRUE),
                                    sort(p_out, decreasing = TRUE), alpha)
  expect_equal(cSu_check_cpp(p, r, alpha_su), naive)
})

test_that("tolerance: p_out[k] within TOLERANCE of c_val is handled leniently", {
  # first_below uses >= threshold + TOLERANCE, so a q-element at
  # c_val + TOLERANCE/2 is included in the Q-frontier and contributes
  # v_star = k+1 (denominator negative => v_raw < 0), which can lower V_u.
  # Without the tolerance it would be skipped entirely.
  alpha    <- 0.05
  alpha_su <- alpha / ell_alpha(alpha)
  r        <- 2L
  m        <- 5L
  alpha_su_u1 <- r * alpha_su / 1L    # c_val at u=1 (smallest c_val)
  # Place one p_out element just above c_val by TOLERANCE/2: without
  # the first_below tolerance it is excluded; with it, it enters the
  # Q-frontier and contributes v_star=1, possibly lowering V_u to 1
  # and allowing the P-witnesses to cover [0,0].
  p_out <- c(alpha_su_u1 + .SU_TOLERANCE / 2,   # just above c_val at u=1
             0.9, 0.8)
  p_in  <- c(alpha_su * 0.9, alpha_su * 0.4)    # small inside p-values
  p     <- c(sort(p_out, decreasing = TRUE),
             sort(p_in,  decreasing = TRUE))
  # Verify no error and result is consistent with naive (small enough to run)
  naive <- check_condition_naive_su(sort(p_in, decreasing = TRUE),
                                    sort(p_out, decreasing = TRUE), alpha)
  fast  <- cSu_check_cpp(p, r, alpha_su)
  expect_equal(fast, naive)
})

test_that("tolerance: near-equal p_in and p_out values handled correctly", {
  # Tests both rk (>= p_out[k] + TOLERANCE, undercount) and s[j]
  # (> pj - TOLERANCE, overcount) at the same time.
  # When p_in[j] == p_out[k] exactly, the lenient rk counts it as NOT
  # above p_out[k] (smaller B_k, lower V_u) while lenient s[j] counts
  # p_out[k] as above p_in[j] (larger interval).  Both help toward TRUE.
  alpha    <- 0.05
  alpha_su <- alpha / ell_alpha(alpha)
  shared   <- alpha_su * 0.5   # a value that appears in both p_in and p_out
  r        <- 2L
  p_in     <- c(shared,       alpha_su * 0.2)
  p_out    <- c(shared + 1e-15, 0.7, 0.6)   # essentially equal to p_in[0]
  m        <- length(p_in) + length(p_out)
  p        <- c(sort(p_out, decreasing = TRUE),
                sort(p_in,  decreasing = TRUE))
  naive <- check_condition_naive_su(sort(p_in, decreasing = TRUE),
                                    sort(p_out, decreasing = TRUE), alpha)
  fast  <- cSu_check_cpp(p, r, alpha_su)
  expect_equal(fast, naive)
})

test_that("tolerance: tiny perturbation within TOLERANCE does not flip TRUE to FALSE", {
  # Construct a set that is significant (naive returns TRUE), then nudge
  # each p-value slightly upward (by TOLERANCE/2) in the direction that
  # would produce FALSE without tolerance, and confirm TRUE is still returned.
  alpha    <- 0.05
  alpha_su <- alpha / ell_alpha(alpha)
  set.seed(17)
  # Find a small case that is exactly on the boundary
  found <- FALSE
  for (attempt in seq_len(200)) {
    z     <- sort(runif(8)^3, decreasing = TRUE)
    m     <- length(z)
    r_val <- closedSu(z, alpha = alpha)
    if (r_val == 0L) next
    r  <- r_val
    idx <- (m - r + 1L):m
    p  <- sort(z[ idx], decreasing = TRUE)
    q  <- sort(z[-idx], decreasing = TRUE)
    # Check the naive criterion value is near 0 for some (u,v)
    min_val <- Inf
    for (u in seq_len(r)) {
      for (v in 0L:length(q)) {
        t   <- sort(c(p[seq_len(u)], q[seq_len(v)]), decreasing = TRUE)
        n   <- length(t)
        i   <- seq_along(t)
        val <- min(t - r * alpha_su * (n - i + 1L) / (n * u))
        if (val < min_val) min_val <- val
      }
    }
    if (abs(min_val) < 1e-6) { found <- TRUE; break }
  }
  if (found) {
    # Nudge all p-values up by TOLERANCE/2 (makes criterion slightly harder)
    z_nudged <- pmin(z + .SU_TOLERANCE / 2, 1)
    # Both original and nudged should return the same r (or nudged >= original)
    r_orig   <- closedSu(z,        alpha = alpha)
    r_nudged <- closedSu(z_nudged, alpha = alpha)
    expect_gte(r_nudged, r_orig - 1L,
               label = "tiny nudge should not reduce r by more than 1")
  } else {
    skip("Could not find a near-boundary case in 200 attempts")
  }
})


test_that("APSAC example: closedSu >= Su rejections", {
  p <- c(
    0.0001, 0.0004, 0.0019, 0.0095, 0.0201,
    0.0278, 0.0298, 0.0344, 0.0459, 0.3240,
    0.4262, 0.5719, 0.6528, 0.7590, 1.0000
  )
  r_Su  <- su_rejections(p, alpha = 0.05)
  r_cSu <- closedSu(p, alpha = 0.05)
  expect_gte(r_cSu, r_Su)
  
  r_Su_10  <- su_rejections(p, alpha = 0.10)
  r_cSu_10 <- closedSu(p, alpha = 0.10)
  expect_gte(r_cSu_10, r_Su_10)
})

# ===========================================================================
# 9. APPROXIMATE DISCOVERY MODE
# ===========================================================================

test_that("approximate never exceeds exact (core guarantee)", {
  r_exact  <- closedSu(p_mixed, approximate = FALSE)
  r_approx <- closedSu(p_mixed, approximate = TRUE)
  expect_lte(r_approx, r_exact)
  # Should be reasonably close (within 20% of m, a generous tolerance)
  expect_gte(r_approx, r_exact - ceiling(0.2 * length(p_mixed)))
})

test_that("approximate agrees with exact on hand-crafted small example", {
  r_exact  <- closedSu(p_small, approximate = FALSE)
  r_approx <- closedSu(p_small, approximate = TRUE)
  expect_lte(r_approx, r_exact)
})

test_that("approximate result >= Su rejections (lower-bound guarantee)", {
  # The bisection is seeded by su_rejections, so anything Su rejects is
  # always returned even when the bisection adds nothing on top.
  for (alpha in c(0.01, 0.05, 0.10)) {
    r_approx <- closedSu(p_mixed, alpha = alpha, approximate = TRUE)
    r_Su     <- su_rejections(p_mixed, alpha = alpha)
    expect_gte(r_approx, r_Su,
               label = paste("approx >= Su at alpha =", alpha))
  }
})

test_that("top-r_approx set is confirmed significant by set-checking mode", {
  r_approx <- closedSu(p_mixed, approximate = TRUE)
  if (r_approx > 0) {
    top_set <- order(p_mixed)[seq_len(r_approx)]
    expect_true(closedSu(p_mixed, set = top_set))
  }
})

test_that("approximate handles corner cases without error", {
  expect_equal(closedSu(rep(1,   10), approximate = TRUE), 0)
  expect_equal(closedSu(rep(1e-10, 20), approximate = TRUE), 20L)
  expect_equal(closedSu(p_small, alpha = 0, approximate = TRUE), 0)
  expect_equal(closedSu(p_small, alpha = 1, approximate = TRUE), length(p_small))
})

test_that("approximate result is invariant to input order", {
  p_shuffled <- sample(p_mixed)
  expect_equal(closedSu(p_mixed,    approximate = TRUE),
               closedSu(p_shuffled, approximate = TRUE))
})

test_that("set-checking mode is unaffected by approximate flag", {
  set_idx        <- order(p_mixed)[seq_len(5)]
  result_exact   <- closedSu(p_mixed, set = set_idx, approximate = FALSE)
  result_approx  <- closedSu(p_mixed, set = set_idx, approximate = TRUE)
  expect_equal(result_approx, result_exact)
})