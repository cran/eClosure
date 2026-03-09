#include <Rcpp.h>
#include <deque>
#include <utility>
#include <algorithm>

using namespace Rcpp;

// Tolerance threshold to avoid numerical imprecision problems
// Values below this threshold are considered negative
const double TOLERANCE = 1e-12;

// Structure to hold the result of meanconsistency check
// When meanconsistent is false, at_i_in and at_i_out indicate
// where the negative value was found (useful for warm starts)
struct MeanConsistentResult {
  bool meanconsistent;
  int at_i_in;
  int at_i_out;
};

// Helper function to compute the target function for meanconsistency
// This implements a DC (Difference of Convex) function:
//   target(i_in, i_out) = sum_e - penalty
// where:
//   sum_e = sum of i_in smallest e-values in set + i_out smallest e-values outside set
//   penalty = i_in * (i_in + i_out) / (r * alpha)
//
// The function is convex in sum_e (linear, hence convex) and convex in penalty (quadratic)
// Therefore target = convex - convex, which is generally non-convex
// However, for fixed i_in, the function is unimodal in i_out, enabling bisection search
inline double compute_target(const NumericVector& cum_e, 
                             int i_in, 
                             int i_out, 
                             int set_offset, 
                             double div_ra) {
  // sum_e uses cumulative sums for efficiency:
  // cum_e[set_offset + i_in] gives sum of first (m-r) + i_in values
  // cum_e[set_offset] gives sum of first (m-r) values (all outside set)
  // cum_e[i_out] gives sum of first i_out values outside set
  // Therefore: sum of i_in smallest in set = cum_e[set_offset+i_in] - cum_e[set_offset]
  //            sum of i_out smallest outside = cum_e[i_out]
  double sum_e = cum_e[set_offset + i_in] - cum_e[set_offset] + cum_e[i_out];
  
  // Penalty is quadratic: i_in * (i_in + i_out) / (r * alpha)
  double penalty = i_in * (i_in + i_out) * div_ra;
  
  return sum_e - penalty;
}

// Check if a set of size r is mean consistent
// 
// Mean consistency means: target(i_in, i_out) >= 0 for all valid (i_in, i_out)
// where i_in ranges from 0 to r and i_out ranges from 0 to m-r
//
// Algorithm strategy:
// 1. Iterate over all i_in values using breadth-first search (to find violations quickly)
// 2. For each i_in, minimize over i_out using bisection (exploiting unimodality)
// 3. Early termination: return FALSE as soon as any negative value is found
// 4. Warm start: if provided, check the warm start location first (useful when checking
//    multiple r values sequentially, as failure points often persist across similar r)
//
// Parameters:
//   cum_e: cumulative sums of sorted e-values (length m+1, starting with 0)
//   r: size of the set to check
//   alpha: significance level
//   warm_i_in, warm_i_out: optional warm start from previous r value
//
// Returns:
//   MeanConsistentResult with:
//   - meanconsistent: true if all values non-negative, false otherwise
//   - at_i_in, at_i_out: location of negative value (if found)
MeanConsistentResult meanconsistent_cpp(const NumericVector& cum_e, 
                                        int r, 
                                        double alpha,
                                        int warm_i_in = -1, 
                                        int warm_i_out = -1) {
  // Base case: empty set is always mean consistent
  if (r == 0) {
    return {true, -1, -1};
  }
  
  // Pre-calculate constants to avoid repeated computation
  int m = cum_e.size() - 1;  // Total number of hypotheses
  int set_offset = m - r;  // Index from where "in set" e-values start in cum_e
  double div_ra = 1.0 / (r * alpha);  // Reciprocal of r*alpha for penalty calculation
  
  // Breadth-first search over i_in values
  // We use a queue of intervals [start, end] and always check the midpoint first
  // This strategy finds negative values quickly by exploring the middle regions first
  // Note: We start at i_in=1, not i_in=0, because:
  //   At i_in=0: sum_e >= 0 (sum of non-negative e-values) and penalty = 0
  //   Therefore target(0, i_out) = sum_e >= 0 always, no need to check
  
  bool warm = (warm_i_in >= 0);  // Flag to use warm start on first iteration
  std::deque<std::pair<int, int>> queue_i_in;
  queue_i_in.push_back({1, r});
  
  while (!queue_i_in.empty()) {
    // Get the next interval to check
    auto interval_i_in = queue_i_in.front();
    queue_i_in.pop_front();
    
    // Choose midpoint: use warm start if available, otherwise use interval midpoint
    int mid_i_in;
    if (warm) {
      // Warm start: use the i_in from previous r (capped at current r)
      mid_i_in = std::min(warm_i_in, r);
    } else {
      // Normal case: use midpoint of interval
      mid_i_in = (interval_i_in.first + interval_i_in.second) / 2;
    }
    
    // For this fixed i_in value, minimize over i_out using bisection
    // The target function is unimodal in i_out, so we can use bisection
    std::pair<int, int> interval_i_out = {0, m - r};
    
    while (interval_i_out.first <= interval_i_out.second) {
      if (interval_i_out.first == interval_i_out.second) {
        // Interval reduced to single point - check it
        double val = compute_target(cum_e, mid_i_in, interval_i_out.first, 
                                    set_offset, div_ra);
        if (val < -TOLERANCE) {
          // Found a negative value - not mean consistent
          return {false, mid_i_in, interval_i_out.first};
        }
        // This i_out value passes, exit the i_out loop
        break;
      } else {
        // Interval contains at least 2 points - check two adjacent points
        // to determine which direction the minimum lies
        int mid1_i_out;
        if (warm) {
          // Warm start: use i_out from previous r
          // Cap at m-r-1 (not m-r) because we'll check mid2_i_out = mid1_i_out + 1
          // and need mid2_i_out <= m-r
          mid1_i_out = std::min(warm_i_out, m - r - 1);
        } else {
          // Normal case: use midpoint of interval
          mid1_i_out = (interval_i_out.first + interval_i_out.second) / 2;
        }
        
        // Check two adjacent points
        double target1 = compute_target(cum_e, mid_i_in, mid1_i_out, 
                                       set_offset, div_ra);
        if (target1 < -TOLERANCE) {
          return {false, mid_i_in, mid1_i_out};
        }
        
        int mid2_i_out = mid1_i_out + 1;
        double target2 = compute_target(cum_e, mid_i_in, mid2_i_out, 
                                       set_offset, div_ra);
        if (target2 < -TOLERANCE) {
          return {false, mid_i_in, mid2_i_out};
        }
        
        // Determine which half contains the minimum based on the slope
        if (target1 > target2) {
          // Function is decreasing: minimum must be in upper half
          interval_i_out = {mid2_i_out + 1, interval_i_out.second};
        } else {
          // Function is increasing: minimum must be in lower half
          interval_i_out = {interval_i_out.first, mid1_i_out - 1};
        }
      }
      
      // After first iteration, disable warm start
      warm = false;
    }
    
    // This i_in value passed (no negative values found for any i_out)
    // Subdivide the interval and add parts to queue for breadth-first exploration
    
    if (mid_i_in > interval_i_in.first) {
      // Add left part to queue if there are points to the left of mid_i_in
      queue_i_in.push_back({interval_i_in.first, mid_i_in - 1});
    }
    
    if (mid_i_in < interval_i_in.second) {
      // Add right part to queue if there are points to the right of mid_i_in
      queue_i_in.push_back({mid_i_in + 1, interval_i_in.second});
    }
  }
  
  // All i_in values passed - the set is mean consistent
  return {true, -1, -1};
}

// Find the largest mean consistent set
// 
// Since mean consistency is NOT a monotone property (a larger set passing
// does not guarantee a smaller set passes), we must check all set sizes.
// We iterate from r=m down to r=1, returning the first size that passes.
//
// Optimization: We use warm starts by passing the failure location from r
// to r-1, since failure points often persist or shift only slightly between
// consecutive r values.
//
// Parameters:
//   e: sorted vector of e-values
//   alpha: significance level
//
// Returns:
//   The largest r such that the set of r largest e-values is mean consistent
//   (returns 0 if no non-empty set is mean consistent)
//
// [[Rcpp::export]]
int largestmeanconsistent_cpp(NumericVector e, double alpha) {
  int m = e.size();
  
  // Compute cumulative sums for efficiency
  // cum_e[0] = 0, cum_e[k] = sum of first k e-values
  NumericVector cum_e(m + 1);
  cum_e[0] = 0.0;
  for (int i = 0; i < m; i++) {
    cum_e[i + 1] = cum_e[i] + e[i];
  }
  
  // Initialize warm start values
  int i_in = -1;
  int i_out = -1;
  
  // Check each set size from largest to smallest
  for (int r = m; r >= 1; r--) {
    MeanConsistentResult mc = meanconsistent_cpp(cum_e, r, alpha, i_in, i_out);
    
    if (mc.meanconsistent) {
      // Found the largest mean consistent set
      return r;
    }
    
    // Update warm start for next iteration
    i_in = mc.at_i_in;
    i_out = mc.at_i_out;
  }
  
  // No non-empty set is mean consistent
  return 0;
}

// Find the number of rejections by the eBH (e-values Benjamini-Hochberg) procedure
//
// The eBH procedure rejects hypothesis i if:
//   e[i] * (m - i + 1) / m >= 1 / alpha
//
// This provides a lower bound for the largest mean consistent set, since
// any set rejected by eBH is guaranteed to be mean consistent.
//
// Parameters:
//   e: sorted vector of e-values (assumes sorted in increasing order)
//   alpha: significance level
//
// Returns:
//   Number of rejections by eBH (0 if none)
//
// [[Rcpp::export]]
int eBH_cpp(NumericVector e, double alpha) {
  int m = e.size();
  double apm = alpha / m;
  
  // Check each hypothesis in order
  for (int i = 0; i < m; i++) {
    int rev_i = m - i;  // Number of hypotheses at or above current rank
    
    // Check eBH rejection criterion
    if (e[i] * rev_i * apm >= 1.0 + TOLERANCE) {
      return rev_i;
    }
  }
  
  return 0;
}

// Find the largest mean consistent set approximately (faster but less accurate)
//
// This function assumes that mean consistency behaves approximately like a
// monotone property: if a set of size r is NOT mean consistent, then larger
// sets are also unlikely to be mean consistent. While this is not strictly true,
// it often holds in practice.
//
// Strategy:
// 1. Use eBH to find a lower bound (eBH rejections are always mean consistent)
// 2. Use bisection search on r between this lower bound and m
// 3. Use warm starts to speed up each mean consistency check
//
// Parameters:
//   e: sorted vector of e-values
//   alpha: significance level
//
// Returns:
//   Approximate largest r such that the set is mean consistent
//
// [[Rcpp::export]]
int largestmeanconsistent_approximate_cpp(NumericVector e, double alpha) {
  int m = e.size();
  
  // Get lower bound from eBH
  int min_r = eBH_cpp(e, alpha);
  
  // If eBH rejects all hypotheses, no need to search further
  if (min_r == m) {
    return m;
  }
  
  // Setup for bisection search
  std::pair<int, int> interval_r = {min_r + 1, m};
  
  // Compute cumulative sums
  NumericVector cum_e(m + 1);
  cum_e[0] = 0.0;
  for (int i = 0; i < m; i++) {
    cum_e[i + 1] = cum_e[i] + e[i];
  }
  
  int best_r = min_r;
  int i_in = -1;
  int i_out = -1;
  
  // Bisection search on r
  while (interval_r.first <= interval_r.second) {
    int mid_r = (interval_r.first + interval_r.second) / 2;
    
    MeanConsistentResult mc = meanconsistent_cpp(cum_e, mid_r, alpha, i_in, i_out);
    
    if (mc.meanconsistent) {
      // This size passes - try larger sizes
      best_r = mid_r;
      interval_r = {mid_r + 1, interval_r.second};
    } else {
      // This size fails - try smaller sizes
      interval_r = {interval_r.first, mid_r - 1};
      
      // Update warm start (only when we find a failure)
      i_in = mc.at_i_in;
      i_out = mc.at_i_out;
    }
  }
  
  return best_r;
}

// Wrapper function to check mean consistency of a specific set
// This is called by the R wrapper ceBH when a specific set is provided
//
// Parameters:
//   cum_e: cumulative sums of e-values (with out-of-set values first, then in-set)
//   r: size of the set
//   alpha: significance level
//
// Returns:
//   true if the set is mean consistent, false otherwise
//
// [[Rcpp::export]]
bool meanconsistent_wrapper_cpp(NumericVector cum_e, int r, double alpha) {
  MeanConsistentResult result = meanconsistent_cpp(cum_e, r, alpha);
  return result.meanconsistent;
}
