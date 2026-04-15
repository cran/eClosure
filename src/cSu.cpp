// cSu.cpp  — Closed Su procedure for simultaneous FDR control
//
// Exported functions:
//   su_cpp              — standard Su lower bound (BH at alpha / ell_alpha)
//   find_largest_r_cpp  — discovery mode: full r-loop in C++ (no R overhead)
//   cSu_check_sep_cpp   — single-r check, exposed for testing / diagnostics
//   cSu_check_cpp       — set-checking mode: arbitrary set (general case)
//
// Input layout (both modes)
// -------------------------
// p[0 .. m-r-1]  — p-values OUTSIDE the set   (m-r largest, sorted dec)
// p[m-r .. m-1]  — p-values INSIDE  the set   (r smallest,  sorted dec)
//
// Speedups in the separated case
// --------------------------------
// 1. Global P-witness sufficient condition (O(1), checked once per r):
//    r*alpha >= p_in[0]*m  =>  p_in[0] alone witnesses every (u,v) pair.
//
// 2. Per-u P-witness sufficient condition (O(1) per u):
//    r*alpha >= p_in[0]*(m_r+u) => skip Q-frontier entirely for this u.
//
// 3. Fast P-witness lower bound:
//    Check j=u-1 first (smallest p_in, largest Vp_j).
//
// 4. BFS over u — violations found in O(log r) checks per failing r.
//
// 5. Binary search for Q-frontier start (first_below):
//    p_out sorted dec => p_out[k] >= c_val forms a prefix;
//    skip it in O(log m_r) rather than O(k_start).
//
// 6. r-loop entirely in C++ (find_largest_r_cpp):
//    Eliminates R<->C++ copy of warm_cache on every iteration.
//
// 7. Combined check + cache-update (check_and_update_warm):
//    The warm-start Q-frontier check and update_v_star_cache iterate the
//    SAME [k_start, m_r) range computing identical values.  Merging them
//    into one loop halves the dominant cost.  Removing the early-exit
//    also lets the compiler auto-vectorize the division-heavy inner loop.

#include <Rcpp.h>
#include <deque>
#include <utility>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

using namespace Rcpp;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static const double SU_TOL   = 1e-9;   // rounding guard in ifloor/iceil
static const double TOLERANCE = 1e-10;  // floating-point comparison tolerance
// (matches cBY.cpp and ceBH.cpp)

static inline int ifloor_su(double x) {
  if (std::isnan(x) || x > (double)std::numeric_limits<int>::max())
    return std::numeric_limits<int>::max();
  if (x < -(double)std::numeric_limits<int>::max())
    return -std::numeric_limits<int>::max();
  return static_cast<int>(std::floor(x + SU_TOL));
}

static inline int iceil_su(double x) {
  if (std::isnan(x) || x > (double)std::numeric_limits<int>::max())
    return std::numeric_limits<int>::max();
  return static_cast<int>(std::ceil(x - SU_TOL));
}

static double lambert_w_neg(double x) {
  double lx = std::log(-x);
  double w  = lx - std::log(-lx);
  for (int i = 0; i < 8; ++i) {
    double ew    = std::exp(w);
    double denom = ew * (w + 1.0);
    if (std::abs(denom) < 1e-300) break;
    double w_new = w - (w * ew - x) / denom;
    if (std::abs(w_new - w) < 1e-15 * std::abs(w_new)) break;
    w = w_new;
  }
  return w;
}

static inline double compute_ell_alpha(double alpha) {
  return -lambert_w_neg(-alpha / std::exp(1.0));
}

// First k in [0, n) where arr[k] < threshold + TOLERANCE (arr sorted dec).
// Using threshold + TOLERANCE rather than threshold means elements that are
// within TOLERANCE of the boundary are included, matching cBY's convention.
// Returns n if all values are >= threshold + TOLERANCE.
static inline int first_below(const double* arr, int n, double threshold) {
  int lo = 0, hi = n;
  while (lo < hi) {
    int mid = (lo + hi) / 2;
    if (arr[mid] >= threshold + TOLERANCE) lo = mid + 1;
    else                                   hi = mid;
  }
  return lo;
}

// P-witness reach for a given u.
// Returns the maximum Vp_j = floor(c_val*(u-j)/p_in[j] - u) over j in [0, u-1].
// Speedup 3: checks j=u-1 first (smallest p_in => largest Vp_j).
static inline int p_witness_reach(
    int           u,
    const double* p_in,
    double        c_val)
{
  if (p_in[u - 1] == 0.0) return std::numeric_limits<int>::max();
  int reach = ifloor_su(c_val / p_in[u - 1] - (double)u);   // j=u-1, a=1
  
  for (int j = 0; j < u - 1; ++j) {
    if (p_in[j] == 0.0) return std::numeric_limits<int>::max();
    int a    = u - j;
    int Vp_j = ifloor_su(c_val * (double)a / p_in[j] - (double)u);
    if (Vp_j > reach) {
      reach = Vp_j;
      // No early exit here to keep the function simple;
      // callers do the v_max comparison.
    }
  }
  return reach;
}

// ===========================================================================
// Separated case — Q-frontier helpers
// ===========================================================================

// Compute v*(k, u) = max(k+1, ceil(c_val*k/(c_val - p_out[k]) - u)).
// Called for one k; p_out[k] < c_val assumed.
static inline int v_star_k(int k, const double* p_out, double c_val, int u) {
  int    k1    = k + 1;
  double v_raw = (k1 == 1) ? 0.0
  : c_val * (double)k / (c_val - p_out[k]) - (double)u;
  return std::max(k1, iceil_su(v_raw));
}

// ---------------------------------------------------------------------------
// check_one_u_sep
//
// Check one u value in the separated case (no cache update).
// Used by BFS (step 2) where we don't need to update the cache.
// ---------------------------------------------------------------------------
static bool check_one_u_sep(
    int                      u,
    const double*            p_in,
    const double*            p_out,
    double                   c_val,
    double                   r_alpha,
    int                      m_r,
    int&                     out_V_u)
{
  // Speedup 2: per-u P-witness sufficient condition
  if (p_in[0] == 0.0 ||
      r_alpha >= p_in[0] * (double)(m_r + u) - SU_TOL) {
    out_V_u = m_r + 1;
    return true;
  }
  
  int k_start = first_below(p_out, m_r, c_val);
  
  // Q-frontier with early exit (no cache needed for BFS)
  int V_u = m_r + 1;
  for (int k = k_start; k < m_r; ++k) {
    int vs = v_star_k(k, p_out, c_val, u);
    if (vs < V_u) {
      V_u = vs;
      if (V_u <= 0) break;
    }
  }
  
  out_V_u = V_u;
  if (V_u <= 0) return true;
  
  int v_max = std::min(V_u - 1, m_r);
  return p_witness_reach(u, p_in, c_val) >= v_max;
}

// ---------------------------------------------------------------------------
// check_and_update_warm
//
// Combined Q-frontier check AND cache update for warm_failing_u.
//
// Key insight: check_one_u_sep and update_v_star_cache previously made two
// independent passes over [k_start, m_r), each computing identical v*(k,u)
// values.  This function merges them into ONE pass, halving the dominant
// cost.  Removing the early-exit also enables compiler auto-vectorization.
//
// cache is updated in-place:
//   cache[k] = max(cache[k], v*(k, warm_u))   for k in [k_start, m_r)
//   cache[k] = m_r + 1                         for k in [0, k_start)
// cache is also resized to m_r (adding a new element if m_r grew by 1).
//
// Returns true if the condition is satisfied for warm_u.
// ---------------------------------------------------------------------------
static bool check_and_update_warm(
    int               warm_u,
    const double*     p_in,
    const double*     p_out,
    double            c_val,
    double            r_alpha,
    int               m_r,
    std::vector<int>& cache,   // in/out: updated in place
    int&              out_V_u)
{
  // Resize cache to new m_r (grows by 1 each r step)
  // New element defaults to k+1 = m_r, the natural lower bound.
  if ((int)cache.size() < m_r)
    cache.resize(m_r, m_r);
  
  // Speedup 2: per-u P-witness sufficient condition
  if (p_in[0] == 0.0 ||
      r_alpha >= p_in[0] * (double)(m_r + warm_u) - SU_TOL) {
    out_V_u = m_r + 1;
    return true;
  }
  
  int k_start = first_below(p_out, m_r, c_val);
  
  // Fill prefix (p_out[k] >= c_val): these never achieve frontier minimum
  std::fill(cache.begin(), cache.begin() + k_start, m_r + 1);
  
  // Single pass: compute v*(k, warm_u) for k in [k_start, m_r),
  // apply cache lower bound, update cache, track V_u.
  // No early exit => full loop => compiler can auto-vectorize.
  int V_u = m_r + 1;
  for (int k = k_start; k < m_r; ++k) {
    int vs = std::max(cache[k], v_star_k(k, p_out, c_val, warm_u));
    cache[k] = vs;
    if (vs < V_u) V_u = vs;
  }
  
  out_V_u = V_u;
  if (V_u <= 0) return true;
  
  int v_max = std::min(V_u - 1, m_r);
  return p_witness_reach(warm_u, p_in, c_val) >= v_max;
}

// ---------------------------------------------------------------------------
// compute_v_star_sep
//
// Compute v*(k, u) fresh for all k — called when a NEW failing_u is found.
// ---------------------------------------------------------------------------
static std::vector<int> compute_v_star_sep(
    const double* p_out,
    int           m_r,
    double        c_val,
    int           u)
{
  std::vector<int> cache(m_r);
  int k_start = first_below(p_out, m_r, c_val);
  std::fill(cache.begin(), cache.begin() + k_start, m_r + 1);
  for (int k = k_start; k < m_r; ++k)
    cache[k] = v_star_k(k, p_out, c_val, u);
  return cache;
}

// ---------------------------------------------------------------------------
// Result struct: avoids Rcpp List allocation in the inner r-loop
// ---------------------------------------------------------------------------
struct SepResult {
  bool             sat;
  int              failing_u;
  std::vector<int> cache;
};

// ---------------------------------------------------------------------------
// check_r_sep — internal single-r check (no Rcpp types)
// ---------------------------------------------------------------------------
static SepResult check_r_sep(
    const double*    p_in,
    const double*    p_out,
    int              r,
    int              m_r,
    double           r_alpha,
    int              warm_failing_u,
    std::vector<int> warm_cache)    // passed by value; modified in place
{
  // Speedup 1: global P-witness sufficient condition (hardest u = r)
  if (p_in[0] == 0.0 ||
      r_alpha >= p_in[0] * (double)(m_r + r) - SU_TOL)
    return {true, 0, {}};
  
  // ------------------------------------------------------------------
  // Step 1: combined warm-start check + cache update (speedup 7).
  // Merges what were previously two separate O(m_r) passes into one.
  // ------------------------------------------------------------------
  if (warm_failing_u >= 1 && warm_failing_u <= r) {
    double c_val = r_alpha / (double)warm_failing_u;
    int    V_u;
    bool   sat = check_and_update_warm(warm_failing_u, p_in, p_out,
                                       c_val, r_alpha, m_r,
                                       warm_cache, V_u);
    if (!sat)
      // Cache already updated in check_and_update_warm
      return {false, warm_failing_u, std::move(warm_cache)};
  }
  
  // ------------------------------------------------------------------
  // Step 2: BFS over u in [1, r], skipping warm_failing_u.
  // Uses check_one_u_sep (no cache update) since any failure here
  // will get a fresh cache via compute_v_star_sep.
  // ------------------------------------------------------------------
  static const std::vector<int> empty_cache_unused;  // unused — kept for signature
  std::deque<std::pair<int,int>> queue;
  queue.push_back({1, r});
  
  while (!queue.empty()) {
    auto interval = queue.front();
    queue.pop_front();
    int lo = interval.first;
    int hi = interval.second;
    int u  = (lo + hi) / 2;
    
    if (u == warm_failing_u) {
      if (lo < u) queue.push_back({lo, u - 1});
      if (hi > u) queue.push_back({u + 1, hi});
      continue;
    }
    
    double c_val = r_alpha / (double)u;
    int    V_u;
    bool sat = check_one_u_sep(u, p_in, p_out, c_val, r_alpha, m_r, V_u);
    
    if (!sat)
      return {false, u, compute_v_star_sep(p_out, m_r, c_val, u)};
    
    if (lo < u) queue.push_back({lo, u - 1});
    if (hi > u) queue.push_back({u + 1, hi});
  }
  
  return {true, 0, {}};
}

// ===========================================================================
// Exported functions
// ===========================================================================

// [[Rcpp::export]]
int su_cpp(NumericVector p, double alpha) {
  int    m         = p.size();
  double ell_alpha = compute_ell_alpha(alpha);
  double factor    = (double)m * ell_alpha;
  for (int i = 0; i < m; ++i) {
    if (p[i] * factor <= (double)(m - i) * alpha + SU_TOL)
      return m - i;
  }
  return 0;
}

// ---------------------------------------------------------------------------
// find_largest_r_cpp — discovery mode, full r-loop in C++
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
int find_largest_r_cpp(NumericVector p, double alpha) {
  int           m  = p.size();
  const double* pd = p.begin();
  
  int              warm_failing_u = 0;
  std::vector<int> cache;
  
  for (int r = m; r >= 1; --r) {
    int           m_r     = m - r;
    const double* p_in    = pd + m_r;
    const double* p_out   = pd;
    double        r_alpha = (double)r * alpha;
    
    SepResult res = check_r_sep(p_in, p_out, r, m_r, r_alpha,
                                warm_failing_u, std::move(cache));
    if (res.sat) return r;
    
    warm_failing_u = res.failing_u;
    cache          = std::move(res.cache);
  }
  
  return 0;
}

// ---------------------------------------------------------------------------
// cSu_check_sep_cpp — single-r check (Rcpp interface, for diagnostics)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List cSu_check_sep_cpp(NumericVector p,
                       int           r,
                       int           warm_failing_u,
                       IntegerVector warm_cache,
                       double        alpha) {
  int           m     = p.size();
  int           m_r   = m - r;
  const double* p_in  = p.begin() + m_r;
  const double* p_out = p.begin();
  double        r_alpha = (double)r * alpha;
  
  std::vector<int> cache(warm_cache.begin(), warm_cache.end());
  
  SepResult res = check_r_sep(p_in, p_out, r, m_r, r_alpha,
                              warm_failing_u, std::move(cache));
  
  return List::create(
    Named("res")       = res.sat,
    Named("failing_u") = res.failing_u,
    Named("cache")     = IntegerVector(res.cache.begin(), res.cache.end()));
}

// ---------------------------------------------------------------------------
// cSu_check_cpp — set-checking mode, arbitrary set
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
bool cSu_check_cpp(NumericVector p, int r, double alpha) {
  int           m     = p.size();
  int           m_r   = m - r;
  const double* p_out = p.begin();
  const double* p_in  = p.begin() + m_r;
  
  if (r == 0) return true;
  
  // Precompute s[j] = #{k : p_out[k] > p_in[j] - TOLERANCE}
  // Lenient: near-equal outside values are counted as above, giving a
  // larger s[j] and a wider Case A interval.
  std::vector<int> s(r);
  for (int j = 0; j < r; ++j) {
    double pj = p_in[j];
    int lo = 0, hi = m_r;
    while (lo < hi) {
      int mid = (lo + hi + 1) / 2;
      if (p_out[mid - 1] > pj - TOLERANCE) lo = mid; else hi = mid - 1;
    }
    s[j] = lo;
  }
  
  std::vector<int> rk(std::max(m_r, 1), 0);
  std::vector<int> iL(2 * r + 1), iR(2 * r + 1);
  double r_alpha = (double)r * alpha;
  
  for (int u = r; u >= 1; --u) {
    double c_val = r_alpha / (double)u;
    
    if (p_in[0] == 0.0 ||
        r_alpha >= p_in[0] * (double)(m_r + u) - SU_TOL)
      continue;
    
    if (m_r > 0) {
      int ptr = 0;
      for (int k = 0; k < m_r; ++k) {
        while (ptr < u && p_in[ptr] >= p_out[k] + TOLERANCE) ++ptr;
        rk[k] = ptr;
      }
    }
    
    int k_start = first_below(p_out, m_r, c_val);
    int V_u = m_r + 1;
    for (int k = k_start; k < m_r; ++k) {
      int    k1  = k + 1;
      int    B_k = k + rk[k];
      double v_raw = (B_k == 0) ? 0.0
      : c_val * (double)B_k / (c_val - p_out[k]) - (double)u;
      int vs = std::max(k1, iceil_su(v_raw));
      if (vs < V_u) {
        V_u = vs;
        if (V_u <= 0) break;
      }
    }
    if (V_u <= 0) continue;
    
    int v_max = std::min(V_u - 1, m_r);
    
    int ni = 0;
    for (int j = 0; j < u; ++j) {
      int j1 = j + 1;
      int a  = u - j;
      int sj = s[j];
      int Aj = a - sj;
      
      if (sj > 0) {
        int ir_A;
        if (p_in[j] == 0.0) {
          ir_A = std::min(sj - 1, v_max);
        } else {
          ir_A = std::min({sj - 1, v_max,
                          ifloor_su(c_val * (double)a / p_in[j] - (double)u)});
        }
        if (ir_A >= 0) { iL[ni] = 0; iR[ni] = ir_A; ++ni; }
      }
      
      int il_B = sj;
      if (il_B > v_max) continue;
      
      if (p_in[j] < c_val + TOLERANCE) {                // B1 (or near-equal)
        double num = p_in[j] * (double)u - c_val * (double)Aj;
        double den = c_val - p_in[j];
        int il = (den <= 0.0) ? il_B   // near-equal: suffix from il_B
        : std::max(il_B, iceil_su(num / den));
        if (il <= v_max) { iL[ni] = il; iR[ni] = v_max; ++ni; }
      } else if (p_in[j] > c_val + TOLERANCE) {          // B2
        double num2 = c_val * (double)Aj - p_in[j] * (double)u;
        if (num2 >= -TOLERANCE) {
          int ir = std::min(v_max, ifloor_su(num2 / (p_in[j] - c_val)));
          if (il_B <= ir) { iL[ni] = il_B; iR[ni] = ir; ++ni; }
        }
      } else {                                            // B3: |p_in[j] - c_val| <= TOLERANCE
        if (j1 == 1 && sj == 0) { iL[ni] = 0; iR[ni] = v_max; ++ni; }
      }
    }
    
    if (ni == 0) return false;
    
    std::vector<int> ord(ni);
    for (int i = 0; i < ni; ++i) ord[i] = i;
    std::sort(ord.begin(), ord.end(),
              [&](int a, int b) { return iL[a] < iL[b]; });
    
    if (iL[ord[0]] > 0) return false;
    int reach = iR[ord[0]];
    if (reach < v_max) {
      for (int idx = 1; idx < ni; ++idx) {
        int i = ord[idx];
        if (iL[i] > reach + 1) return false;
        if (iR[i] > reach) {
          reach = iR[i];
          if (reach >= v_max) break;
        }
      }
      if (reach < v_max) return false;
    }
  }
  
  return true;
}

// ---------------------------------------------------------------------------
// find_largest_r_approximate_cpp — approximate discovery mode
//
// Uses su_cpp as a lower bound (BH at alpha/ell_alpha), then bisects upward.
// Analogous to largestmeanconsistent_approximate_cpp in ceBH.cpp.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
int find_largest_r_approximate_cpp(NumericVector p, double alpha) {
  int m = p.size();
  
  // Lower bound: su_cpp already applies the ell_alpha correction internally,
  // but here p is sorted decreasing and alpha is already alpha_su (passed from R).
  // So we call su_cpp with the raw alpha to get the plain-BH lower bound at
  // alpha/ell_alpha.  Equivalently: count BH rejections at alpha_su.
  int min_r = su_cpp(p, alpha);
  
  if (min_r == m) return m;
  
  // Compute cumulative sums (not needed here; check_r_sep works on raw p)
  // Bisection over r in [min_r + 1, m]
  int lo = min_r + 1, hi = m;
  int best_r = min_r;
  
  int              warm_failing_u = 0;
  std::vector<int> cache;
  
  while (lo <= hi) {
    int mid_r = (lo + hi) / 2;
    int m_r   = m - mid_r;
    
    const double* pd    = p.begin();
    const double* p_in  = pd + m_r;
    const double* p_out = pd;
    double        r_alpha = (double)mid_r * alpha;
    
    SepResult res = check_r_sep(p_in, p_out, mid_r, m_r, r_alpha,
                                warm_failing_u, std::move(cache));
    
    if (res.sat) {
      best_r = mid_r;
      lo     = mid_r + 1;
      // On success, reset warm state — failure point from a larger r
      // is not a reliable warm start when moving to an even larger r.
      warm_failing_u = 0;
      cache          = {};
    } else {
      hi             = mid_r - 1;
      warm_failing_u = res.failing_u;
      cache          = std::move(res.cache);
    }
  }
  
  return best_r;
}
