test_that("empty set is always mean consistent", {
  e <- c(0.5, 1.0, 2.0)
  expect_true(closedeBH(e, set = logical(3), alpha = 0.05))
  expect_true(closedeBH(e, set = integer(0), alpha = 0.05))
})

test_that("empty input returns 0", {
  expect_equal(closedeBH(numeric(0), alpha = 0.05), 0L)
})

test_that("single hypothesis: rejected iff e >= 1/alpha", {
  alpha <- 0.05
  # e = 1/alpha is exactly the boundary (just rejected)
  expect_equal(closedeBH(1 / alpha, alpha = alpha), 1L)
  # e just below boundary: not rejected
  expect_equal(closedeBH(1 / alpha - 1e-6, alpha = alpha), 0L)
  # e well above boundary
  expect_equal(closedeBH(100, alpha = alpha), 1L)
  # e = 0: not rejected
  expect_equal(closedeBH(0, alpha = alpha), 0L)
})

test_that("all e-values equal to 1 (no evidence): no rejections", {
  e <- rep(1, 20)
  expect_equal(closedeBH(e, alpha = 0.05), 0L)
})

test_that("all e-values equal to 1/alpha: all hypotheses rejected", {
  alpha <- 0.05
  e <- rep(1 / alpha, 10)
  expect_equal(closedeBH(e, alpha = alpha), 10L)
})

test_that("all e-values equal to 1/(2*alpha): exactly half evidence, check boundary", {
  # With all e-values equal to c, mean consistency for a set of size r requires
  # c >= (r + k) / (r * alpha) for all k in 0..m-r.
  # Worst case is k = m - r, giving c >= m / (r * alpha).
  # So for a set of size r to be consistent: r >= m / (c * alpha).
  # With c = 1/(2*alpha), we need r >= 2m. So no non-trivial set passes for m > 1.
  e <- rep(1 / (2 * 0.05), 10)
  expect_equal(closedeBH(e, alpha = 0.05), 0L)
})

test_that("alpha = 1: threshold collapses, all non-negative e-values rejected", {
  # With alpha = 1, mean consistency requires sum(e_S) / |S| >= (|S| + k) / |S|
  # i.e. mean(e_S) >= 1 + k/|S| for all k. At k=0, mean(e_S) >= 1.
  # Any set whose mean e-value >= 1 should be consistent.
  e <- rep(2, 5)
  expect_equal(closedeBH(e, alpha = 1), 5L)
  # e-values all below 1: nothing rejected
  e_low <- rep(0.5, 5)
  expect_equal(closedeBH(e_low, alpha = 1), 0L)
})

test_that("discovery mode: return value is a non-negative integer <= m", {
  set.seed(1)
  e <- rexp(50)
  r <- closedeBH(e, alpha = 0.05)
  expect_type(r, "integer")
  expect_gte(r, 0L)
  expect_lte(r, length(e))
})

test_that("set-checking mode: return value is a single logical", {
  set.seed(2)
  e <- rexp(20)
  result <- closedeBH(e, set = e > 1, alpha = 0.05)
  expect_type(result, "logical")
  expect_length(result, 1)
})

# --- Consistency between modes ---

test_that("top-r set from discovery mode is mean consistent in set-checking mode", {
  # 5 strong signals (500 >> m/alpha = 20/0.05 = 400) guarantee rejections
  e <- c(rep(1, 15), rep(500, 5))
  r <- closedeBH(e, alpha = 0.05)
  expect_gt(r, 0)  # explicitly assert rejections exist rather than silently skipping
  top_r_set <- e >= sort(e, decreasing = TRUE)[r]
  expect_true(closedeBH(e, set = top_r_set, alpha = 0.05))
})

test_that("a set larger than r is not necessarily mean consistent", {
  # 2 strong signals, 18 nulls at 1. The full set of 20 fails mean consistency
  # because the nulls drag the mean down. m/alpha = 20/0.05 = 400, so e=500
  # ensures the top-2 set is rejected by eBH, but the full set is not consistent.
  e <- c(rep(1, 18), 500, 500)
  r <- closedeBH(e, alpha = 0.05)
  expect_lt(r, length(e))  # explicitly assert not everything is rejected
  expect_false(closedeBH(e, set = rep(TRUE, length(e)), alpha = 0.05))
})

test_that("set-checking and discovery modes are consistent across alpha values", {
  # Use strong signals that guarantee rejections even at alpha = 0.01
  # Strictest threshold: m/alpha = 25/0.01 = 2500, so e = 3000 ensures rejection
  e <- c(rep(1, 20), rep(3000, 5))
  for (alpha in c(0.01, 0.05, 0.1, 0.2)) {
    r <- closedeBH(e, alpha = alpha)
    expect_gt(r, 0, label = paste("expected rejections at alpha =", alpha))
    top_r <- rank(-e, ties.method = "first") <= r
    expect_true(
      closedeBH(e, set = top_r, alpha = alpha),
      label = paste("top-r set not consistent at alpha =", alpha)
    )
  }
})

# --- Consistency with standard eBH ---

# Helper: standard eBH on e-values
# Returns number of rejections
eBH_reference <- function(e, alpha) {
  m <- length(e)
  e_sorted <- sort(e, decreasing = TRUE)
  for (i in seq_along(e_sorted)) {
    if (e_sorted[i] * i >= m / alpha) {
      return(i)
    }
  }
  return(0L)
}

test_that("closedeBH rejects at least as many as standard eBH", {
  set.seed(5)
  for (trial in 1:20) {
    e <- rexp(30, rate = runif(1, 0.2, 1.5))
    alpha <- sample(c(0.01, 0.05, 0.10), 1)
    r_ceBH <- closedeBH(e, alpha = alpha)
    r_eBH  <- eBH_reference(e, alpha)
    expect_gte(
      r_ceBH, r_eBH,
      label = paste("trial", trial, "alpha =", alpha)
    )
  }
})

test_that("eBH rejections are always mean consistent", {
  set.seed(6)
  for (trial in 1:20) {
    e <- rexp(30, rate = runif(1, 0.2, 1.5))
    alpha <- 0.05
    r_eBH <- eBH_reference(e, alpha)
    if (r_eBH > 0) {
      top_eBH <- e >= sort(e, decreasing = TRUE)[r_eBH]
      expect_true(
        closedeBH(e, set = top_eBH, alpha = alpha),
        label = paste("eBH rejection set not consistent in trial", trial)
      )
    }
  }
})

# --- Exact vs approximate ---

test_that("approximate result never exceeds exact result", {
  set.seed(7)
  for (trial in 1:30) {
    e <- rexp(50, rate = runif(1, 0.1, 2))
    alpha <- 0.05
    r_exact  <- closedeBH(e, alpha = alpha, approximate = FALSE)
    r_approx <- closedeBH(e, alpha = alpha, approximate = TRUE)
    expect_lte(
      r_approx, r_exact,
      label = paste("approximate exceeds exact in trial", trial)
    )
  }
})

test_that("exact and approximate agree on clear signal and clear null cases", {
  # Clear signal: all very large e-values
  e_signal <- rep(1000, 20)
  expect_equal(
    closedeBH(e_signal, alpha = 0.05, approximate = FALSE),
    closedeBH(e_signal, alpha = 0.05, approximate = TRUE)
  )
  # Clear null: all e-values equal to 1
  e_null <- rep(1, 20)
  expect_equal(
    closedeBH(e_null, alpha = 0.05, approximate = FALSE),
    closedeBH(e_null, alpha = 0.05, approximate = TRUE)
  )
})

# --- Monotonicity in alpha ---

test_that("more rejections at larger alpha (weakly monotone)", {
  set.seed(8)
  e <- rexp(40, rate = 0.5)
  alphas <- c(0.01, 0.05, 0.10, 0.20)
  rejections <- sapply(alphas, function(a) closedeBH(e, alpha = a))
  expect_true(
    all(diff(rejections) >= 0),
    label = "rejections should be non-decreasing in alpha"
  )
})

# --- Invariance to input ordering ---

test_that("result does not depend on the order of e-values", {
  set.seed(9)
  e <- rexp(30, rate = 0.4)
  e_shuffled <- sample(e)
  expect_equal(closedeBH(e, alpha = 0.05), closedeBH(e_shuffled, alpha = 0.05))
})

test_that("set-checking result does not depend on order within set or outside set", {
  set.seed(10)
  e <- rexp(20, rate = 0.4)
  set_logical <- e > median(e)
  # Shuffle the e-values while keeping set membership intact
  perm <- sample(length(e))
  expect_equal(
    closedeBH(e, set = set_logical, alpha = 0.05),
    closedeBH(e[perm], set = set_logical[perm], alpha = 0.05)
  )
})

# --- Scaling invariance ---

test_that("scaling all e-values by a constant scales the number of rejections predictably", {
  # If we multiply all e-values by k > 1, we can only reject more (or equal).
  set.seed(11)
  e <- rexp(30, rate = 0.5)
  r1 <- closedeBH(e, alpha = 0.05)
  r2 <- closedeBH(e * 10, alpha = 0.05)
  expect_gte(r2, r1)
})

test_that("scaling all e-values by k is equivalent to scaling alpha by k", {
  # mean consistency condition scales linearly: multiplying e by k is the
  # same as multiplying alpha by k (up to boundary effects from TOLERANCE).
  set.seed(12)
  e <- rexp(20, rate = 0.5)
  k <- 2
  expect_equal(
    closedeBH(e * k, alpha = 0.05),
    closedeBH(e, alpha = 0.05 * k)
  )
})