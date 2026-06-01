#include <Rcpp.h>
#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// ---------------------------------------------------------------------------
// Tolerance for floating-point comparisons (mirrors R's TOLERANCE = 1e-10)
// ---------------------------------------------------------------------------
static const double TOLERANCE = 1e-10;

// ---------------------------------------------------------------------------
// r.closed (built in; this is the next.r function)
//
// Self-contained: combines the "cap" and the truncated variant bound without
// calling a separate r.variant helper.
//
//   cap   = trunc(k*s/(m-k) + TOLERANCE)              if s < m
//         = m                                          otherwise
//
//   ms    = m - s
//   ksm   = k - ms
//   denom = ksm * ms - r
//   b     = trunc(min(m, ms*(ksm-1)*r / denom) + TOLERANCE)   if denom > 0
//         = m                                                  otherwise
//
//   return min(cap, b)
// ---------------------------------------------------------------------------
static inline double r_closed(double r, int m, int s, int k) {
  double cap = (s < m) ? std::trunc((double)k * s / (double)(m - k) + TOLERANCE)
    : (double)m;
  
  int ms  = m - s;
  int ksm = k - ms;
  double denom = (double)ksm * ms - r;
  double b;
  if (denom > 0.0)
    b = std::trunc(std::min((double)m,
                            (double)ms * (ksm - 1) * r / denom) + TOLERANCE);
  else
    b = (double)m;
  
  return std::min(cap, b);
}

// ---------------------------------------------------------------------------
// cBH_cpp
//
// Closed Benjamini-Hochberg inner loop.  Breadth-first traversal over set
// sizes s (the s.queue bisection), with a doubly linked list (prv/nxt/last)
// tracking the candidate rejections shared across passes.
//
// A given index may be dropped during one s-pass and then re-encountered by a
// later pass that starts its consecutive k-scan at m-s+1 (which may already be
// dead).  Removal is therefore guarded by a liveness flag: only an index still
// in the list is unlinked, preventing a second unlink from corrupting the
// chain (and with it the final value of last).  The v / r / a recursion still
// runs over consecutive k, since the cumulative maximum legitimately depends
// on the r-values of already-removed indices.
//
// Arguments:
//   p     - NumericVector of p-values (need not be sorted; sorted internally)
//   alpha - significance level
//
// Returns:
//   the size of the largest closedBH-significant set (0 if none).
//
// Indexing note: the R reference is 1-indexed.  Here the linked-list arrays
// prv/nxt are kept 1-indexed (index 0 == NULL) to mirror the R logic exactly,
// while p is accessed as p[k-1] for R's p[k].
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
int cBH_cpp(NumericVector p_in, double alpha) {
  int m = p_in.size();
  
  // alpha == 0: no rejections possible.
  if (alpha == 0.0) return 0;
  
  // Sort input p-values and incorporate alpha now to avoid later divisions.
  std::vector<double> p(p_in.begin(), p_in.end());
  std::sort(p.begin(), p.end());
  for (int i = 0; i < m; ++i) p[i] /= alpha;
  
  // ------------------------------------------------------------------
  // Ordered set as a doubly linked list (1-indexed; 0 == NULL).
  //   prv[k] = predecessor of k, nxt[k] = successor of k
  //   alive[k] = whether k is still in the list
  //   last   = largest remaining index (the eventual return value)
  // ------------------------------------------------------------------
  std::vector<int> prv(m + 1), nxt(m + 1);
  std::vector<char> alive(m + 1, 1);
  for (int k = 1; k <= m; ++k) { prv[k] = k - 1; nxt[k] = k + 1; }
  nxt[m] = 0;
  int last = m;
  
  // ------------------------------------------------------------------
  // Breadth-first traversal queue of (from.s, to.s) intervals.
  // A moving head index gives FIFO behaviour (m/2, 3m/4, m/4, ...).
  // ------------------------------------------------------------------
  std::vector<std::pair<int,int> > s_queue;
  s_queue.push_back(std::make_pair(1, m));
  size_t head = 0;
  
  while (head < s_queue.size()) {
    int from_s = s_queue[head].first;
    int to_s   = s_queue[head].second;
    ++head;
    
    int s = (from_s + to_s + 1) / 2;
    from_s = std::max(from_s, m - last + 1);
    if (to_s > s)   s_queue.push_back(std::make_pair(s + 1, to_s));
    if (from_s < s) s_queue.push_back(std::make_pair(from_s, s - 1));
    
    int ms = m - s;
    
    // Starting values.  For the full set (s == m, ms == 0) the bound r is m,
    // not k; for every other column the generic start r = k applies.
    int k = ms + 1;
    double r = (s == m) ? (double)m : (double)k;
    double a = (double)k / s;
    int v = k - 1;
    
    while (k <= last) {
      // Update cumulative maximum.
      if (p[k - 1] <= a + TOLERANCE) v = std::max(v, (int)r);
      
      // If v < k and k is still in the list, drop it: remove(k).
      // The alive guard prevents a dead index (a possible start point of this
      // pass) from being unlinked a second time and corrupting the chain.
      if (v < k && alive[k]) {
        alive[k] = 0;
        if (prv[k]) nxt[prv[k]] = nxt[k];
        if (nxt[k]) {
          prv[nxt[k]] = prv[k];
        } else {
          last = prv[k];
          if (last == 0) return 0;   // early stop
        }
      }
      
      if (v >= last) break;          // v >= last > k for the rest of the loop
      
      // Update.
      k = k + 1;
      r = r_closed(r, m, s, k);
      a = (double)(k - ms) * r / ((r - (double)ms) * (double)s);
    }
  }
  
  return last;
}

// ---------------------------------------------------------------------------
// cBH_adjust_cpp
//
// Adjusted p-values for the closed Benjamini-Hochberg procedure: for each
// hypothesis, the smallest alpha at which it is rejected (i.e. at which cBH
// returns at least that hypothesis's sorted rank).
//
// O(m^2) time, O(m) memory.  Each set size s is processed in a single pass
// over k, maintaining a monotone deque of the windowed minimum ratio
//   p_(k') / a(k', s)   over the qualifying window { k' <= k : r(k', s) >= k }.
// The qualifying window is a contiguous suffix because r(., s) is nondecreasing
// in k', so its left edge only advances as k grows.
//
// Arguments:
//   p_in  - NumericVector of p-values (need not be sorted; sorted internally)
//   cap   - if true, raise each adjusted p-value to at least its raw p-value
//
// Returns:
//   NumericVector of adjusted p-values in the original input order.
//
// Indexing note: the working arrays are 1-indexed (index 0 unused / NULL) to
// mirror the R reference; p is accessed as p[k-1] for R's p[k].
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector cBH_adjust_cpp(NumericVector p_in, bool cap) {
  int m = p_in.size();
  NumericVector out(m);
  if (m == 0) return out;
  
  // Order of the input (ascending), and the sorted p-values.  A stable sort
  // matches R's order(), so tied p-values keep their input order and the
  // result is reproducible and identical to the R reference.
  std::vector<int> ord(m);
  for (int i = 0; i < m; ++i) ord[i] = i;
  std::stable_sort(ord.begin(), ord.end(),
                   [&](int x, int y) { return p_in[x] < p_in[y]; });
  std::vector<double> p(m);
  for (int i = 0; i < m; ++i) p[i] = p_in[ord[i]];
  
  // Adjusted p-values in sorted order; built up via per-(k, s) maxima.
  std::vector<double> padj(m + 1, 0.0);   // 1-indexed; padj[k] for rank k
  
  // Reusable deque (parallel index/ratio arrays) and horizon history.
  // dq holds positions i with strictly increasing ratio, front = current min.
  std::vector<int>    dq_idx(m + 1);
  std::vector<double> dq_rat(m + 1);
  std::vector<double> r_hist(m + 1);      // r(ms+i, s) for i = 1..s
  
  for (int s = m; s >= 1; --s) {
    int ms = m - s;
    
    int dq_front = 1, dq_back = 0;        // empty deque
    int left = 1;                         // left edge of qualifying window
    
    // Starting values for k = ms + 1 (where k + s = m + 1, so a = k / s).
    // The full set (s == m) seeds the horizon at m rather than k.
    double r_cur = (s == m) ? (double)m : (double)(ms + 1);
    int    k = ms + 1;
    double a = (double)k / s;
    
    for (int i = 1; i <= s; ++i) {
      // r_cur and a hold the values for the current k = ms + i.
      r_hist[i] = r_cur;
      // a > 0 and finite are invariants here, so no guard is needed.
      double ratio = p[k - 1] / a;
      
      // Push onto the deque, evicting larger ratios from the back.
      while (dq_back >= dq_front && dq_rat[dq_back] >= ratio) --dq_back;
      ++dq_back;
      dq_idx[dq_back] = i;
      dq_rat[dq_back] = ratio;
      
      // Advance the window's left edge: require r_hist[j] >= k.
      while (left <= i && r_hist[left] < (double)k) ++left;
      
      // Drop deque entries that have fallen out of the window.
      while (dq_front <= dq_back && dq_idx[dq_front] < left) ++dq_front;
      
      if (dq_front > dq_back) {
        // No qualifying k' remains: k cannot survive this s, and since the
        // window only shrinks, neither can any larger k.
        for (int j = k; j <= m; ++j) padj[j] = 1.0;
        break;
      }
      
      double alpha_needed = dq_rat[dq_front];
      if (alpha_needed > padj[k]) padj[k] = alpha_needed;
      
      // Update r and a for the next k = ms + i + 1.
      if (i < s) {
        ++k;
        r_cur = r_closed(r_cur, m, s, k);
        a = (double)(k - ms) * r_cur / ((r_cur - (double)ms) * (double)s);
      }
    }
  }
  
  // cBH returns the largest survivor, so padj[k] = min over k' >= k.
  for (int k = m - 1; k >= 1; --k) {
    if (padj[k + 1] < padj[k]) padj[k] = padj[k + 1];
  }    
  
  // Cap at 1, and (optionally) raise to at least the raw p-value; then
  // scatter back into the original input order.
  for (int k = 1; k <= m; ++k) {
    double val = padj[k];
    if (val > 1.0) val = 1.0;
    if (cap && val < p[k - 1]) val = p[k - 1];
    out[ord[k - 1]] = val;
  }
  
  return out;
}
