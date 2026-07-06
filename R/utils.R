# =============================================================================
# Utility functions
# =============================================================================
# This file contains small helper functions shared by different fitting methods.
# =============================================================================


# -----------------------------------------------------------------------------
# Truncate probabilities to avoid numerical issues in log or logit calculations.
# -----------------------------------------------------------------------------
# Inputs:
#   x   : Numeric vector of probabilities.
#   eps : Lower and upper truncation level.
#
# Output:
#   Numeric vector with values restricted to [eps, 1 - eps].
# -----------------------------------------------------------------------------
clamp_prob <- function(x, eps = 1e-12) {
  pmin(pmax(x, eps), 1 - eps)
}


# -----------------------------------------------------------------------------
# Stable log-sum-exp calculation.
# -----------------------------------------------------------------------------
# Computes log(sum(exp(x))) in a numerically stable way.
# This is useful when summing terms on the log scale, such as log Bayes factors.
#
# Input:
#   x : Numeric vector on the log scale.
#
# Output:
#   log(sum(exp(x))).
# -----------------------------------------------------------------------------
log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}
