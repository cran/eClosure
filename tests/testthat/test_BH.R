# tests/testthat/test_BH.R

TOLERANCE <- 1e-10

# ---------------------------------------------------------------------------
# Naive reference implementation: O(m^2) time and memory, very close to the
# math. Builds the full r-, a- and v-matrices over all (k, s) and returns the
# largest k whose cumulative-maximum column values v[k, ] dominate k for every
# set size s. Used as the ground-truth oracle for the fast C++ implementation.
# ---------------------------------------------------------------------------
cBH.naive <- function(p, alpha = 0.05, next.r) {
  if (alpha == 0) return(length(p == 0))
  p <- sort(p)
  m <- length(p)
  
  k <- matrix(1:m, m, m, byrow = FALSE)
  s <- matrix(1:m, m, m, byrow = TRUE)
  
  # make r
  r <- matrix(0, m, m)
  r[k + s <= m + 1] <- k[k + s <= m + 1]
  r[1, m] <- m
  for (run.s in 2:m) {
    for (run.k in (m - run.s + 2):m)
      r[run.k, run.s] <- next.r(r[run.k - 1, run.s], m, run.s, run.k)
  }
  
  # make a
  a <- matrix(0, m, m)
  a[k + s <= m + 1] <- (k * alpha / s)[k + s <= m + 1]
  a[k + s >  m + 1] <- ((k + s - m) * r * alpha / ((r + s - m) * s))[k + s > m + 1]
  
  # make thresholded r
  p <- matrix(p, m, m, byrow = FALSE)
  tr <- r * (p <= a + TOLERANCE)
  # make v
  v <- apply(tr, 2, cummax)
  
  # get r
  rs <- which(apply(v >= k, 1, all))
  r <- ifelse(length(rs) == 0, 0, max(rs))
  return(r)
}

# next.r / r.closed built into the package, replicated here for the oracle.
next.r.ref <- function(r, m, s, k) {
  cap <- if (s < m) trunc(k * s / (m - k) + TOLERANCE) else m
  ms  <- m - s
  ksm <- k - ms
  denom <- ksm * ms - r
  b <- if (denom > 0) trunc(min(m, ms * (ksm - 1) * r / denom) + TOLERANCE) else m
  min(cap, b)
}

naive <- function(p, alpha = 0.05) {
  if (length(p) == 1L) {
    # the matrix construction assumes m >= 2; handle m = 1 directly
    if (alpha == 0) return(0L)
    return(as.integer(sort(p)[1] <= alpha + TOLERANCE))
  }
  cBH.naive(p, alpha, next.r.ref)
}


# ===========================================================================
# 1. EXAMPLES FROM THE PAPER
# ===========================================================================

test_that("BH example", {
  p <- c(
    0.0001, 0.0004, 0.0019, 0.0095, 0.0201,
    0.0278, 0.0298, 0.0344, 0.0459, 0.3240,
    0.4262, 0.5719, 0.6528, 0.7590, 1.0000
  )
  expect_equal(closedBH(p, alpha = 0.05), 4L)
  expect_equal(closedBH(p, alpha = 0.10), 9L)
})


# ===========================================================================
# 2. BASIC SANITY CHECKS
# ===========================================================================

test_that("discovery mode returns a non-negative integer within [0, m]", {
  set.seed(1)
  p <- runif(30)^4
  r <- closedBH(p)
  expect_length(r, 1)
  expect_gte(r, 0)
  expect_lte(r, length(p))
})


# ===========================================================================
# 3. CORNER CASES
# ===========================================================================

test_that("all p-values equal to 1 yields 0 rejections", {
  expect_equal(closedBH(rep(1, 10)), 0)
})

test_that("all p-values near 0 yields m rejections", {
  expect_equal(closedBH(rep(1e-10, 20)), 20L)
})

test_that("alpha = 0 yields 0 rejections", {
  expect_equal(closedBH(runif(8)^4, alpha = 0), 0)
})

test_that("single p-value discovery mode returns 0 or 1", {
  expect_equal(closedBH(0.03, alpha = 0.05), 1L)
  expect_equal(closedBH(0.9,  alpha = 0.05), 0L)
})

test_that("duplicate p-values do not cause errors", {
  expect_no_error(closedBH(c(0.01, 0.01, 0.01, 0.5, 0.5)))
})

test_that("result is invariant to input order", {
  set.seed(7)
  p <- runif(40)^4
  expect_equal(closedBH(p), closedBH(sample(p)))
})


# ===========================================================================
# 4. 'set' ARGUMENT IS IGNORED WITH A WARNING (not yet implemented)
# ===========================================================================

# to add later if set argument implemented

# ===========================================================================
# 5. REGRESSION TESTS: full-set column and dead start indices
#
# These exercise the full-set column s = m, which only becomes binding at
# larger alpha. Two distinct issues lived here:
#   (a) the s = m column must be seeded with r = m, not r = k;
#   (b) a breadth-first pass may start its k-scan at m-s+1 on an index that an
#       earlier pass already removed, so removal must be guarded by liveness to
#       avoid unlinking a dead index twice and corrupting the list.
# ===========================================================================

test_that("regression: r[1, m] = m boundary (m = 2, alpha = 0.5)", {
  p <- c(0.023632, 0.775595)
  expect_equal(closedBH(p, alpha = 0.5), naive(p, 0.5))
  expect_equal(closedBH(p, alpha = 0.5), 2L)
})

test_that("regression: full-set column drives last to 0 (m = 3, alpha = 0.5)", {
  p <- c(0.635327, 0.195443, 0.661134)
  expect_equal(closedBH(p, alpha = 0.5), naive(p, 0.5))
  expect_equal(closedBH(p, alpha = 0.5), 0L)
})

test_that("regression: dead start index in a later s-pass (m = 4, alpha = 0.5)", {
  # The s = m pass starts at k = 1 after smaller-set passes have already removed
  # interior indices; without a liveness guard a dead index is unlinked twice and
  # last settles at a stale survivor (1) instead of 0.
  p <- c(0.1876125, 0.3714447, 0.4229739, 0.6972259)
  expect_equal(closedBH(p, alpha = 0.5), naive(p, 0.5))
  expect_equal(closedBH(p, alpha = 0.5), 0L)
})


# ===========================================================================
# 6. RANDOMIZED AGREEMENT WITH THE NAIVE REFERENCE
#
# Non-trivial inputs: runif(m)^4 (and other powers) push most p-values toward 0
# so that r is a substantial, non-degenerate fraction of m. Swept across a wide
# range of alpha, including large alpha where the s = m column is binding.
# ===========================================================================

test_that("matches naive on runif(m)^4 inputs across alpha", {
  set.seed(2024)
  for (i in 1:150) {
    m     <- sample(2:60, 1)
    p     <- runif(m)^4
    alpha <- sample(c(0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 0.9), 1)
    expect_equal(
      closedBH(p, alpha = alpha),
      naive(p, alpha),
      label = paste0("m=", m, " alpha=", alpha)
    )
  }
})

test_that("matches naive across varied power-law inputs", {
  set.seed(99)
  for (i in 1:150) {
    m     <- sample(2:50, 1)
    pw    <- sample(c(1, 2, 3, 4, 6), 1)
    p     <- runif(m)^pw
    alpha <- sample(c(0.01, 0.05, 0.1, 0.2, 0.5), 1)
    expect_equal(
      closedBH(p, alpha = alpha),
      naive(p, alpha),
      label = paste0("m=", m, " pw=", pw, " alpha=", alpha)
    )
  }
})

test_that("matches naive with a strong-signal / null mixture", {
  set.seed(11)
  for (i in 1:80) {
    m    <- sample(5:60, 1)
    nsig <- sample(0:m, 1)
    p    <- runif(m)
    if (nsig > 0) p[seq_len(nsig)] <- runif(nsig)^8   # strong signals
    alpha <- sample(c(0.01, 0.05, 0.1, 0.2), 1)
    expect_equal(
      closedBH(p, alpha = alpha),
      naive(p, alpha),
      label = paste0("m=", m, " nsig=", nsig, " alpha=", alpha)
    )
  }
})


# ===========================================================================
# 7. closedBH CONTAINS BH
# ===========================================================================

test_that("closedBH discovery count is at least the BH rejection count", {
  set.seed(5)
  for (alpha in c(0.01, 0.05, 0.10)) {
    p <- runif(50)^4
    bh <- sum(p.adjust(p, method = "BH") <= alpha)
    expect_gte(closedBH(p, alpha = alpha), bh)
  }
})


# ===========================================================================
# 8. ADJUSTED P-VALUES
# ===========================================================================

# Naive O(m^2) reference for adjusted p-values, built directly from the r- and
# a-matrices. Ground-truth oracle for closedBH.adjust.
cBH.adjust.naive <- function(p, next.r, cap = TRUE) {
  ord <- order(p)
  p <- p[ord]
  m <- length(p)
  
  k <- matrix(1:m, m, m, byrow = FALSE)
  s <- matrix(1:m, m, m, byrow = TRUE)
  
  r <- matrix(0, m, m)
  r[k + s <= m + 1] <- k[k + s <= m + 1]
  r[1, m] <- m
  if (m >= 2) {
    for (run.s in 2:m)
      for (run.k in (m - run.s + 2):m)
        r[run.k, run.s] <- next.r(r[run.k - 1, run.s], m, run.s, run.k)
  }
  
  a <- matrix(0, m, m)
  a[k + s <= m + 1] <- (k / s)[k + s <= m + 1]
  a[k + s >  m + 1] <- ((k + s - m) * r / ((r + s - m) * s))[k + s > m + 1]
  
  p.adjust <- rep(0, m)
  for (run.s in 1:m) {
    ratio <- rep(Inf, m)
    pos <- a[, run.s] > 0
    ratio[pos] <- p[pos] / a[pos, run.s]
    for (run.k in 1:m) {
      qualifying <- which((1:run.k) <= run.k & r[1:run.k, run.s] >= run.k)
      alpha.needed <- if (length(qualifying) > 0) min(ratio[qualifying]) else Inf
      p.adjust[run.k] <- max(p.adjust[run.k], alpha.needed)
    }
  }
  
  if (m >= 2)
    for (kk in (m - 1):1)
      p.adjust[kk] <- min(p.adjust[kk], p.adjust[kk + 1])
  
  p.adjust <- pmin(p.adjust, 1)
  if (cap) p.adjust <- pmax(p.adjust, p)
  
  p.adjust[ord] <- p.adjust
  p.adjust
}

adjust.naive <- function(p, cap = TRUE) {
  m <- length(p)
  if (m == 1) {
    val <- min(p, 1)
    if (cap) val <- max(val, p)
    return(val)
  }
  cBH.adjust.naive(p, next.r.ref, cap = cap)
}


test_that("adjusted p-values match the naive reference (cap = FALSE)", {
  set.seed(2024)
  for (i in 1:120) {
    m  <- sample(2:40, 1)
    pw <- sample(c(1, 2, 3, 4, 6), 1)
    p  <- runif(m)^pw
    expect_equal(
      closedBH.adjust(p, cap = FALSE),
      adjust.naive(p, cap = FALSE),
      tolerance = 1e-8,
      label = paste0("m=", m, " pw=", pw)
    )
  }
})

test_that("adjusted p-values match the naive reference (cap = TRUE)", {
  set.seed(77)
  for (i in 1:120) {
    m  <- sample(2:40, 1)
    pw <- sample(c(1, 2, 3, 4, 6), 1)
    p  <- runif(m)^pw
    expect_equal(
      closedBH.adjust(p, cap = TRUE),
      adjust.naive(p, cap = TRUE),
      tolerance = 1e-8,
      label = paste0("m=", m, " pw=", pw)
    )
  }
})

test_that("thresholding adjusted p-values reproduces closedBH", {
  # The defining property: a hypothesis is rejected at level alpha exactly when
  # its adjusted p-value (uncapped) is <= alpha; the count must equal closedBH.
  set.seed(13)
  for (i in 1:100) {
    m     <- sample(1:60, 1)
    p     <- runif(m)^4
    padj  <- closedBH.adjust(p, cap = FALSE)
    alpha <- sample(c(0.01, 0.05, 0.1, 0.2, 0.5), 1)
    expect_equal(
      sum(padj <= alpha + 1e-10),
      closedBH(p, alpha = alpha),
      label = paste0("m=", m, " alpha=", alpha)
    )
  }
})

test_that("adjusted p-values are in [0, 1] and respect the cap", {
  set.seed(5)
  p <- runif(50)^4
  padj_t <- closedBH.adjust(p, cap = TRUE)
  padj_f <- closedBH.adjust(p, cap = FALSE)
  expect_true(all(padj_t >= 0 & padj_t <= 1))
  expect_true(all(padj_f >= 0 & padj_f <= 1))
  expect_true(all(padj_t >= p - 1e-10))           # cap = TRUE never below raw p
  expect_length(padj_t, length(p))
})

test_that("adjusted output preserves input order", {
  set.seed(9)
  p <- runif(20)^4
  perm <- sample(seq_along(p))
  expect_equal(
    closedBH.adjust(p)[perm],
    closedBH.adjust(p[perm])
  )
})