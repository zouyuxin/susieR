#' @export
#'
susie_rss_adjustLD = function (z, R, n, bhat, shat, var_y,
                      z_ld_weight = 0,
                      estimate_residual_variance = FALSE,
                      prior_variance = 50,
                      check_prior = TRUE,...) {

  if (estimate_residual_variance)
    warning_message("For estimate_residual_variance = TRUE, please check ",
                    "that R is the \"in-sample\" LD matrix; that is, the ",
                    "correlation matrix obtained using the exact same data ",
                    "matrix X that was used for the other summary ",
                    "statistics. Also note, when covariates are included in ",
                    "the univariate regressions that produced the summary ",
                    "statistics, also consider removing these effects from ",
                    "X before computing R.",style = "hint")

  # Check input R.
  if (missing(z))
    p <- length(bhat)
  else
    p <- length(z)
  if (nrow(R) != p)
      stop(paste0("The dimension of R (",nrow(R)," x ",ncol(R),") does not ",
                  "agree with expected (",p," x ",p,")"))

  # Check input n.
  if (!missing(n))
    if (n <= 1)
      stop("n must be greater than 1")

  # Check inputs z, bhat and shat. Note that bhat is no longer used
  # after this step.
  if (sum(c(missing(z),missing(bhat) || missing(shat))) != 1)
    stop("Please provide either z or (bhat, shat), but not both")
  if (missing(z)) {
    if (length(shat) == 1)
      shat = rep(shat,length(bhat))
    if (length(bhat) != length(shat))
      stop("The lengths of bhat and shat do not agree")
    if (anyNA(bhat) || anyNA(shat))
      stop("bhat, shat cannot have missing values")
    if (any(shat <= 0))
      stop("shat cannot have zero or negative elements")
    z = bhat/shat
  }
  if (length(z) < 1)
    stop("Input vector z should have at least one element")
  z[is.na(z)] = 0

  # When n is provided, compute the PVE-adjusted z-scores.
  if (!missing(n)) {
    adj = (n-1)/(z^2 + n - 2)
    z   = sqrt(adj) * z
  }

  # Modify R by z_ld_weight; this modification was designed to ensure
  # the column space of R contained z, but susie_suff_stat does not
  # require this, and is no longer recommended.
  if (z_ld_weight > 0) {
    warning_message("As of version 0.11.0, use of non-zero z_ld_weight is no longer ",
            "recommended")
    R = muffled_cov2cor((1-z_ld_weight)*R + z_ld_weight*tcrossprod(z))
    R = (R + t(R))/2
  }

  # Call susie_suff_stat. We call susie_suff_stat in two different
  # ways depending on whether n is provided.
  if (missing(n)) {

    # The sample size (n) is not provided, so use unadjusted z-scores.
    # The choice of n=2, yty=1 is mostly arbitrary except in that it
    # ensures var(y) = yty/(n-1) = 1, and because of this
    # scaled_prior_variance = prior_variance.
    warning_message("Providing the sample size (n), or even a rough estimate of n, ",
            "is highly recommended. Without n, the implicit assumption is ",
            "n is large (Inf) and the effect sizes are small (close to zero).")
    s = susie_suff_stat_adjustXtX(XtX = R,Xty = z,n = 2,yty = 1,
                        scaled_prior_variance = prior_variance,
                        estimate_residual_variance = estimate_residual_variance,
                        standardize = FALSE,check_prior = check_prior, ...)
  } else {

    # The sample size (n) is provided, so use PVE-adjusted z-scores.
    if (!missing(shat) & !missing(var_y)) {

      # var_y, shat (and bhat) are provided, so the effects are on the
      # *original scale*.
      XtXdiag = var_y * adj/(shat^2)
      XtX = t(R * sqrt(XtXdiag)) * sqrt(XtXdiag)
      XtX = (XtX + t(XtX))/2
      Xty = z * sqrt(adj) * var_y / shat
    } else {

      # The effects are on the *standardized* X, y scale.
      XtX = (n-1)*R
      Xty = sqrt(n-1)*z
      var_y = 1
    }
    s = susie_suff_stat_adjustXtX(XtX = XtX,Xty = Xty,n = n,yty = (n-1)*var_y,
                        estimate_residual_variance = estimate_residual_variance,
                        check_prior = check_prior,...)
  }
  s$zR_outliers = sort(s$correct_zR_discrepancy$outlier_index)
  return(s)
}

