#include <Rcpp.h>
using namespace Rcpp;

// ---------------------------------------------------------------------------
// Tolerance for floating-point comparisons (mirrors R's tolerance = 1e-10)
// ---------------------------------------------------------------------------
static const double TOLERANCE = 1e-10;

// ---------------------------------------------------------------------------
// BY_cpp
//
// A basic Benjamini-Yekutieli step-up procedure.
// Returns the number of rejections (0 if none).
//
// Arguments:
//   p        - NumericVector of p-values sorted DECREASING (length m)
//   harmonic - NumericVector of harmonic numbers; harmonic[i] = H(i),
//              so harmonic[0] = 0, harmonic[1] = 1, harmonic[2] = 1.5, ...
//              (length >= m+1, 1-indexed in R becomes 0-indexed here)
//   alpha    - significance level
//
// R indexing note:
//   R: harmonic[m+1]  →  C++: harmonic[m]   (harmonic is 0-indexed here)
//   R: p[i]           →  C++: p[i-1]
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
int BY_cpp(NumericVector p, NumericVector harmonic, double alpha) {
  int m = p.size();
  // harmonic[m] in 0-indexed C++ == harmonic[m+1] in 1-indexed R
  double factor = (double)m * harmonic[m];

  for (int i = 0; i < m; ++i) {
    // In R: p[i] * factor / (m-i+1) <= alpha + tolerance
    // i is 0-indexed; rank from largest = i+1, so denominator = m - i
    if (p[i] * factor / (double)(m - i) <= alpha + TOLERANCE)
      return m - i; // number of rejections
  }
  return 0;
}

// ---------------------------------------------------------------------------
// cBY_check_cpp
//
// Checks whether a specific set of cardinality r is cBY-significant,
// given a vector of p-values p of length m.
//
// Significance condition: for all S ⊆ [m],
//   Σ_{i∈S} e_{i,s}  ≥  u / (r · α)
// where u = |S ∩ set|, s = |S|, and
//   e_{i,s} = 1 / (α · ⌈s · H(s) · p_i / α⌉)   if ·H(s)·p_i ≤ α
//           = 0                                     otherwise
//
// The common factor 1/α appears on both sides and is cancelled throughout.
//
// Input layout of p (both segments sorted DECREASING):
//   p[0 .. m-r-1]   : p-values OUTSIDE the set  ("p_out")
//   p[m-r .. m-1]   : p-values INSIDE the set   ("p_in")
//
// The strategy checks all sizes s = 1..m in breadth-first (bisection) order
// so that violating values of s — which tend to cluster — are found quickly.
//
// For each s, the worst-case S consists of:
//   - the u largest p-values from inside the set  (p[m-r .. m-r+u-1])
//   - the v = s-u largest p-values from outside   (p[0 .. v-1])
// We iterate u from minu to maxu, reusing partial sums.
//
// Returns:
//   res   - true if all checks pass (set IS significant)
//   at_s  - if res=false, the value of s where a violation was found
//           (used as a warm-start for the next call)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List cBY_check_cpp(NumericVector p,   // length m, layout described above
                   int r,             // cardinality of the set
                   NumericVector harmonic, // 0-indexed; harmonic[s] = H(s)
                   int warm_s,        // warm-start s (0 = no warm start)
                   double alpha) {

  // Edge case: empty set is always non-significant (vacuously passes)
  if (r == 0) return List::create(Named("res") = true, Named("at_s") = 0);

  int m = p.size();
  double invalpha = 1.0 / alpha;

  // e_out[i+1] caches the e-value for the (i+1)-th largest p_out entry
  // (1-indexed logic from R preserved: e_out has size m-r+1, index 0 unused)
  // We allocate m-r+1 elements; index 0 is a sentinel (always 0).
  std::vector<double> e_out(m - r + 1, 0.0);

  // ------------------------------------------------------------------
  // Breadth-first bisection queue over s ∈ [1, m]
  // Each element is a pair {lo, hi} representing a sub-interval to process.
  // We always pick the midpoint of the interval, matching R's queue logic.
  // ------------------------------------------------------------------
  // Use a deque for efficient front-pop
  std::deque<std::pair<int,int>> queue;
  queue.push_back({1, m});

  bool use_warm = (warm_s != 0);

  while (!queue.empty()) {
    auto interval = queue.front();
    queue.pop_front();
    int lo = interval.first, hi = interval.second;

    int s;
    if (use_warm) {
      s = warm_s;
      use_warm = false;
    } else {
      s = (lo + hi) / 2; // integer division, same as R's %/%
    }

    // harmonic[s] in 0-indexed C++ == harmonic[s+1] in 1-indexed R
    double Hs     = harmonic[s];          // H(s)
    double factor = (double)s * Hs;       // s * H(s)
    double factoralpha = factor * invalpha; // s * H(s) / α

    // u = |S ∩ set|; v = s - u = |S \ set|
    // Constraints: max(1, s-(m-r)) ≤ u ≤ min(s, r)
    int minu = std::max(1, s - (m - r));
    int maxu = std::min(s, r);
    // maxv = s - minu  (not used explicitly, but v = s - u)

    // ------------------------------------------------------------------
    // Sufficient condition to skip this s (and all smaller s):
    //
    // If p_in[1] · s · H(s) ≤ maxu · α  then e_in[1] ≥ 1/(maxu·α),
    // so Σ e_{i,s} ≥ u/(maxu·α) ≥ u/(r·α) and the condition holds.
    //
    // p_in[1] in R (1-indexed, largest in set) = p[m-r] in C++ (0-indexed).
    //
    // Moreover, if this holds at s it holds for all s' < s:
    //   - s ≤ r : condition is p_in[1]·H(s) ≤ α, easier for smaller H(s)
    //   - s > r : condition is p_in[1]·s·H(s) ≤ r·α, LHS decreasing in s
    // So we only recurse into the upper sub-interval.
    // ------------------------------------------------------------------
    if (p[m - r] * factor <= (double)maxu * alpha) {
      if (hi > s) queue.push_back({s + 1, hi});
      continue;
    }

    // For s=1 and the sufficient condition failed: automatic violation.
    if (s == 1) return List::create(Named("res") = false, Named("at_s") = 1);

    // ------------------------------------------------------------------
    // Compute e_out values for indices v = maxv down to 1.
    // maxv = s - minu.  We need e_out[1..maxv], i.e. p_out[0..maxv-1].
    // We iterate from i = maxv-1 down to 0 (0-indexed p_out).
    //
    // R code:  for (i in rev(seq_len(maxv)))   → i goes maxv, maxv-1, ..., 1
    //          e_out[i+1] stores the e-value for the i-th largest p_out
    //          p[i] in R (1-indexed) = p[i-1] in C++ (0-indexed)
    //
    // Zero-out stale entries for positions that are no longer valid
    // (i.e. where factor * p_out[i] >= alpha).  We stop zeroing when we
    // hit a slot that is already 0 (from previous iteration or init).
    // ------------------------------------------------------------------
    int maxv = s - minu;
    double sum_e_out = 0.0;

    for (int i = maxv - 1; i >= 0; --i) {
      // p_out[i] = p[i] (0-indexed)
      if (Hs * p[i] > alpha + TOLERANCE) {
        // p[i] is too large; zero out stale e_out entries
        // (slot index = i+1 in R convention = i+1 in our 1-indexed e_out)
        int j = i; // j is the 0-indexed p_out position, slot = j+1
        while (j >= 0 && e_out[j + 1] != 0.0) {
          e_out[j + 1] = 0.0;
          --j;
        }
        break;
      }
      double next_e = invalpha / std::ceil(p[i] * factoralpha);
      e_out[i + 1] = next_e;
      sum_e_out += next_e;
    }

    // ------------------------------------------------------------------
    // Pre-accumulate sum_e_in for u = 1 .. minu-1
    // p_in[i] in R (1-indexed) = p[m-r + i - 1] in C++ (0-indexed)
    // ------------------------------------------------------------------
    double sum_e_in = 0.0;
    for (int i = minu - 2; i >= 0; --i) {
      // i is 0-indexed offset into p_in; p_in[i] = p[m-r+i]
      double pi = p[m - r + i];
      if (Hs * pi  > alpha + TOLERANCE) break;
      sum_e_in += invalpha / std::ceil(pi * factoralpha);
    }

    // ------------------------------------------------------------------
    // Main loop over u = minu .. maxu
    // ------------------------------------------------------------------
    for (int u = minu; u <= maxu; ++u) {
      int v = s - u;

      // p_in[u] in R (1-indexed) = p[m-r + u - 1] in C++ (0-indexed)
      double pu_in = p[m - r + u - 1];
      double next_e_in = 0.0;
      if (Hs * pu_in < alpha + TOLERANCE) {
        // e-value is non-zero only if factor*p < alpha (strictly)
        next_e_in = invalpha / std::ceil(pu_in * factoralpha);
      }
      sum_e_in += next_e_in;

      double critical = invalpha * (double)u / (double)r;

      if (sum_e_in > critical - TOLERANCE) {
        // Condition met: because p_in is sorted decreasing, later elements
        // are larger, giving larger e-values, so the condition holds for
        // all remaining u as well.  Break out of u-loop.
        break;
      }

      // Condition not yet met from p_in alone; check with p_out contribution
      if (sum_e_in + sum_e_out < critical - TOLERANCE) {
        return List::create(Named("res") = false, Named("at_s") = s);
      }

      // Remove e_out contribution for slot v+1 (which won't be in S for u+1)
      // e_out[v+1] stores e-value of p_out[v] (0-indexed); v = s-u, so slot = v+1
      sum_e_out -= e_out[v + 1];
    }

    // Recurse into both sub-intervals
    if (lo < s) queue.push_back({lo, s - 1});
    if (hi > s) queue.push_back({s + 1, hi});
  }

  // All checks passed
  return List::create(Named("res") = true, Named("at_s") = 0);
}
