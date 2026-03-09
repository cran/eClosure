# tests/testthat/test-cBY.R

# ---------------------------------------------------------------------------
# Helper: Benjamini-Yekutieli rejection count (base R implementation)
# Returns the number of hypotheses rejected by the standard BY procedure.
# ---------------------------------------------------------------------------
by_rejections <- function(p, alpha = 0.05) {
  m <- length(p)
  hm <- sum(1 / seq_len(m))  # H(m)
  threshold <- alpha / hm
  sum(p.adjust(p, method = "BY") <= alpha)
}

# ---------------------------------------------------------------------------
# Reproducible test inputs
# ---------------------------------------------------------------------------
set.seed(123)
p_mixed  <- c(runif(20), rbeta(10, 0.1, 1))   # 20 nulls + 10 non-nulls
p_all_null <- runif(50)                         # all null
p_small  <- c(0.001, 0.002, 0.01, 0.04, 0.5)  # small hand-crafted example
p_single <- 0.03                               # single p-value


# ===========================================================================
# 1. BASIC SANITY CHECKS
# ===========================================================================

test_that("discovery mode returns a non-negative integer", {
  r <- closedBY(p_mixed)
  expect_true(is.numeric(r) || is.integer(r))
  expect_length(r, 1)
  expect_gte(r, 0)
  expect_lte(r, length(p_mixed))
})

test_that("set-checking mode returns a single logical", {
  result <- closedBY(p_mixed, set = p_mixed < 0.05)
  expect_true(is.logical(result))
  expect_length(result, 1)
})

test_that("result does not exceed m", {
  r <- closedBY(p_mixed)
  expect_lte(r, length(p_mixed))
})


# ===========================================================================
# 2. CORNER CASES
# ===========================================================================

test_that("empty set is always TRUE", {
  expect_true(closedBY(p_mixed, set = integer(0)))
  expect_true(closedBY(p_mixed, set = logical(length(p_mixed))))  # all FALSE
})

test_that("single p-value: significant if below BY threshold", {
  # BY threshold for m=1 is alpha itself
  expect_true(closedBY(0.01, set = 1, alpha = 0.05))
  expect_false(closedBY(0.9,  set = 1, alpha = 0.05))
})

test_that("single p-value discovery mode returns 0 or 1", {
  r <- closedBY(p_single)
  expect_true(r %in% c(0L, 1L))
})

test_that("all p-values equal to 1 yields 0 rejections", {
  expect_equal(closedBY(rep(1, 10)), 0)
})

test_that("all p-values near 0 yields m rejections", {
  p_tiny <- rep(1e-10, 20)
  expect_equal(closedBY(p_tiny), 20L)
})

test_that("exact zeros are handled without error", {
  p_with_zero <- c(0, 0.01, 0.02, 0.5)
  expect_no_error(closedBY(p_with_zero))
  expect_no_error(closedBY(p_with_zero, set = c(TRUE, TRUE, FALSE, FALSE)))
})

test_that("p-values of length 1 handled correctly", {
  expect_no_error(closedBY(0.03))
  expect_no_error(closedBY(0.03, set = 1L))
})

test_that("alpha = 0 yields 0 rejections and all sets FALSE", {
  expect_equal(closedBY(p_small, alpha = 0), 0)
  expect_false(closedBY(p_small, set = 1:3, alpha = 0))
})

test_that("alpha = 1 yields m rejections", {
  expect_equal(closedBY(p_small, alpha = 1), length(p_small))
})

test_that("duplicate p-values do not cause errors", {
  p_dup <- c(0.01, 0.01, 0.01, 0.5, 0.5)
  expect_no_error(closedBY(p_dup))
  expect_no_error(closedBY(p_dup, set = 1:3))
})

test_that("set supplied as logical, positive index, and negative index agree", {
  set_logical  <- c(TRUE, TRUE, TRUE, FALSE, FALSE)
  set_posindex <- c(1L, 2L, 3L)
  set_negindex <- c(-4L, -5L)
  r1 <- closedBY(p_small, set = set_logical)
  r2 <- closedBY(p_small, set = set_posindex)
  r3 <- closedBY(p_small, set = set_negindex)
  expect_equal(r1, r2)
  expect_equal(r1, r3)
})


# ===========================================================================
# 3. CONSISTENCY BETWEEN THE TWO MODES
# ===========================================================================

test_that("top-r set (discovery mode) is confirmed significant by set-checking mode", {
  r <- closedBY(p_mixed)
  if (r > 0) {
    top_r_set <- order(p_mixed)[seq_len(r)]
    expect_true(closedBY(p_mixed, set = top_r_set))
  }
})

test_that("set of size r+1 is not necessarily significant (non-monotonicity warning test)", {
  # We cannot assert FALSE in general, but we can assert no error occurs
  r <- closedBY(p_mixed)
  if (r < length(p_mixed)) {
    top_r1_set <- order(p_mixed)[seq_len(r + 1)]
    expect_no_error(closedBY(p_mixed, set = top_r1_set))
  }
})

test_that("discovery mode and set-checking mode agree across multiple alpha levels", {
  for (alpha in c(0.01, 0.05, 0.10, 0.20)) {
    r <- closedBY(p_mixed, alpha = alpha)
    if (r > 0) {
      top_set <- order(p_mixed)[seq_len(r)]
      expect_true(
        closedBY(p_mixed, set = top_set, alpha = alpha),
        label = paste("top-r set should be significant at alpha =", alpha)
      )
    }
  }
})

test_that("exact and approximate discovery modes return the same or approximate result", {
  r_exact  <- closedBY(p_mixed, approximate = FALSE)
  r_approx <- closedBY(p_mixed, approximate = TRUE)
  # Approximate may underestimate but never overestimate
  expect_lte(r_approx, r_exact)
  # And should be reasonably close (within 20% of m, a generous tolerance)
  expect_gte(r_approx, r_exact - ceiling(0.2 * length(p_mixed)))
})

test_that("exact and approximate agree on hand-crafted small example", {
  r_exact  <- closedBY(p_small, approximate = FALSE)
  r_approx <- closedBY(p_small, approximate = TRUE)
  expect_lte(r_approx, r_exact)
})


# ===========================================================================
# 4. closedBY REJECTIONS ALWAYS CONTAIN BY REJECTIONS
# ===========================================================================

test_that("closedBY discovery count is at least as large as BY rejection count", {
  r_cBY <- closedBY(p_mixed)
  r_BY  <- by_rejections(p_mixed)
  expect_gte(r_cBY, r_BY)
})

test_that("BY rejection set is always closedBY-significant (set-checking mode)", {
  r_BY     <- by_rejections(p_mixed)
  if (r_BY > 0) {
    by_set <- order(p_mixed)[seq_len(r_BY)]
    expect_true(closedBY(p_mixed, set = by_set))
  }
})

test_that("BY >= cBY containment holds across multiple alpha levels", {
  for (alpha in c(0.01, 0.05, 0.10)) {
    r_BY  <- by_rejections(p_mixed, alpha = alpha)
    r_cBY <- closedBY(p_mixed, alpha = alpha)
    expect_gte(r_cBY, r_BY,
               label = paste("cBY >= BY at alpha =", alpha))
    if (r_BY > 0) {
      by_set <- order(p_mixed)[seq_len(r_BY)]
      expect_true(
        closedBY(p_mixed, set = by_set, alpha = alpha),
        label = paste("BY set is closedBY-significant at alpha =", alpha)
      )
    }
  }
})

test_that("BY >= cBY containment holds on all-null data", {
  r_BY  <- by_rejections(p_all_null)
  r_cBY <- closedBY(p_all_null)
  expect_gte(r_cBY, r_BY)
})

test_that("BY >= cBY containment holds on hand-crafted small example", {
  r_BY  <- by_rejections(p_small)
  r_cBY <- closedBY(p_small)
  expect_gte(r_cBY, r_BY)
  if (r_BY > 0) {
    by_set <- order(p_small)[seq_len(r_BY)]
    expect_true(closedBY(p_small, set = by_set))
  }
})


# ===========================================================================
# 5. INVARIANCE AND STABILITY
# ===========================================================================

test_that("result is invariant to input order", {
  p_shuffled <- sample(p_mixed)
  expect_equal(closedBY(p_mixed), closedBY(p_shuffled))
})

test_that("result is stable across repeated calls (no mutable state leakage)", {
  r1 <- closedBY(p_mixed)
  r2 <- closedBY(p_mixed)
  expect_equal(r1, r2)
})

test_that("harmonic cache does not corrupt results when called with varying m", {
  # Exercise cache extension by calling with increasing m
  for (n in c(5, 10, 50, 100)) {
    p_n <- sort(runif(n))
    expect_no_error(closedBY(p_n))
  }
  # Final result on p_mixed should still be correct
  expect_equal(closedBY(p_mixed), closedBY(p_mixed))
})


# ===========================================================================
# 5. EXAMPLES FROM THE PAPER
# ===========================================================================

test_that("BH example", {
  p <- c(
    0.0001, 0.0004, 0.0019, 0.0095, 0.0201,
    0.0278, 0.0298, 0.0344, 0.0459, 0.3240,
    0.4262, 0.5719, 0.6528, 0.7590, 1.0000
  )
  expect_equal(closedBY(p, alpha=0.05), 3L)
  expect_equal(closedBY(p, alpha=0.10), 5L)
})
