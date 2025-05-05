update_each_effect_ss_adjustXtX = function (Xty, s_init,
                                  estimate_prior_variance = FALSE,
                                  estimate_prior_method = "optim",
                                  check_null_threshold = 0, 
                                  check_outlier = TRUE,
                                  dist_thres = 0.1,
                                  LD_thres = 0.7, 
                                  likelihood_thres = 0.5, 
                                  weight_thres = 0.5, region = "local",
                                  outlier_method = "IQR", use_post_mean_l = TRUE) {
  if (!estimate_prior_variance)
    estimate_prior_method = "none"
  
  # Repeat for each effect to update.
  s = s_init
  L = nrow(s$alpha)
  N = attr(s$XtX,"d")[1]
  if (L > 0) {
    for (l in 1:L) {
      # print(paste0("l=", l))
      # Compute residuals.
      s$XtXr = s$XtXr - s$XtX %*% (s$alpha[l,] * s$mu[l,])
      XtR = Xty - s$XtXr
      res = single_effect_regression_ss(as.matrix(XtR), attr(s$XtX,"d"), s$V[l],
              s$sigma2, s$pi, estimate_prior_method, check_null_threshold)

      if (s$correct_zR_discrepancy$to_correct) {
        tmp_s = s
        tmp_s$alpha[l,] = res$alpha
        tmp_s$mu[l,] = res$mu
        tmp_s$mu2[l,] = res$mu2
        tmp_s$V[l] = res$V
        post_mean = susie_get_posterior_mean(tmp_s)
        post_mean_l = res$alpha * res$mu

        if (L == 1) {
          prev_cs = NULL
        } else {
          tmp_s_r = tmp_s
          tmp_s_r$alpha = tmp_s_r$alpha[-l, , drop = FALSE]
          tmp_s_r$V = tmp_s_r$V[-l]
          prev_cs = unique(unname(unlist(susie_get_cs(tmp_s_r, Xcorr = s$XtX/N, coverage = 0.75)$cs)))
          # print("=====PREV CS======")
          # print(prev_cs)
        }

        # For the current effect, get the c_index
        local_c_index = get_non_zero_effects_proxy(matrix(res$alpha, nrow = 1), LD = cov2cor(s$XtX), 
          dist_thres = dist_thres, LD_thres = LD_thres)
        # print("Local cindex:")
        # print("remove already identified outliers:")
        local_c_index = setdiff(local_c_index, s$correct_zR_discrepancy$outlier_index)
        # print(local_c_index)

        # Initialize weight and update pattern: there are 3 update patterns, normal, outlier or adjust
        weight = 0
        update_pattern = "normal"
        outliers = c()
        
        if(region == "local"){
          c_index = local_c_index
        }else{
          c_index = c(local_c_index, prev_cs)
        }

        if(length(local_c_index) > 1){
          # Determine the update pattern
          # original distance when weight = 0
          if (use_post_mean_l) {
            dist_original = distance_calc(weight = 0, c_index = c_index,
              as.matrix(XtR), s$XtX, post_mean_l)
            weight = solve_weight(c_index, XtR, s$XtX, post_mean_l)
            dist_update = distance_calc(weight = weight, c_index = c_index,
              as.matrix(XtR), s$XtX, post_mean_l)
          } else {
            dist_original = distance_calc(weight = 0, c_index = c_index,
              as.matrix(Xty), s$XtX, post_mean)
            weight = solve_weight(c_index, Xty, s$XtX, post_mean)
            dist_update = distance_calc(weight = weight, c_index = c_index,
              as.matrix(Xty), s$XtX, post_mean)
          }
          # print(paste0("dist_original:", dist_original))
          # print(paste0("dist_update:", dist_update))
          dist_ratio = abs(dist_update - dist_original)
          
          # print(paste0("Weight:", weight))
          # print(paste("distance_ratio:", dist_ratio))

          # weight large, change large -- outlier
          if (check_outlier == FALSE & l == 1) {
            # print("first one, no outlier check")
            if(weight > 1e-5 & dist_ratio >= likelihood_thres){
              update_pattern = "adjust"
            } else {
              update_pattern = "normal"
            }
          } else {
            if(weight >= weight_thres & dist_ratio >= likelihood_thres){
              update_pattern = "outlier"
            }else if(weight > 1e-5 & weight < weight_thres & dist_ratio >= likelihood_thres){
              update_pattern = "adjust"
            }else{
              if (use_post_mean_l) {
                outlier_index = test_outlier(XtR, s$XtX, post_mean_l, local_c_index, outlier_method)
              } else {
                outlier_index = test_outlier(Xty, s$XtX, post_mean, local_c_index, outlier_method)
              }
              if(length(outlier_index) > 0){
                update_pattern = "outlier"
                outliers = local_c_index[outlier_index]
              }else{
                update_pattern = "normal"
              }
            }
          }
        }
        # print(paste("update_pattern:", update_pattern))

        if(update_pattern == "adjust"){
          # adjust by weight
          XtX_adj <- XtX_update(s$XtX, c_index, N, weight)
          # print("=============Update this effect====================")
          # print(paste0("Weight:", weight))
          attr(XtX_adj,"d") = diag(XtX_adj)
          attr(XtX_adj,"scaled:scale") = attr(s$XtX,"scaled:scale")
          s$XtX = XtX_adj
          # Keep record of weight
          s$correct_zR_discrepancy$weight[[l]] = c(s$correct_zR_discrepancy$weight[[l]], weight)
        } else if(update_pattern == "outlier") {
          if (length(outliers) == 0) {
            s = outlier_removal(s, outlier = local_c_index)
          } else {
            s = outlier_removal(s, outlier = outliers)
          }
          res = single_effect_regression_ss(as.matrix(XtR), attr(s$XtX,"d"), s$V[l],
            s$sigma2,s$pi,estimate_prior_method,check_null_threshold)
        }
      }

      # Update the variational estimate of the posterior mean.
      s$mu[l,]    = res$mu
      s$alpha[l,] = res$alpha
      s$mu2[l,]   = res$mu2
      s$V[l]      = res$V
      s$lbf[l]    = res$lbf_model
      s$lbf_variable[l,] = res$lbf
      s$KL[l]     = -res$lbf_model +
        SER_posterior_e_loglik_ss(attr(s$XtX,"d"),XtR,s$sigma2,
                                  res$alpha * res$mu,res$alpha * res$mu2)

      s$XtXr = s$XtXr + s$XtX %*% (s$alpha[l,] * s$mu[l,])
    }

    if (s$correct_zR_discrepancy$to_correct) {
      if (identical(s_init$XtX, s$XtX) & setequal(s$correct_zR_discrepancy$outlier_index, s_init$correct_zR_discrepancy$outlier_index)) {
        s$correct_zR_discrepancy$outlier_stable_count = s$correct_zR_discrepancy$outlier_stable_count + 1
        # print(paste0("stb num:", s$correct_zR_discrepancy$outlier_stable_count))
        if (s$correct_zR_discrepancy$outlier_stable_count > s$correct_zR_discrepancy$outlier_stabilize) {
          # Wrap it up
          cat("zR mismatch correction done!\n")
          s$correct_zR_discrepancy$to_correct = FALSE
          s$force_iterate = FALSE
        }
      } else {
        s$correct_zR_discrepancy$outlier_stable_count = 0
      }
    }
  }
  s$XtXr = unname(as.matrix(s$XtXr))
  return(s)
}

test_flat <- function(p0, dist_thres) {
  q <- 1 / length(p0)
  js <- js_divergence_p_q(p0, q)
  return(js >= dist_thres)
}

# Jensen-Shannon Divergence
js_divergence_p_q <- function(p, q) {
  m <- (p + q) / 2
  return((kl_divergence_p_q(p, m) + kl_divergence_p_q(q, m)) / 2)
}

# KL Divergence
kl_divergence_p_q <- function(p, q) {
  p[which(p < .Machine$double.xmin)] <- .Machine$double.xmin
  q[which(q < .Machine$double.xmin)] <- .Machine$double.xmin
  sum(p * (log(p)-log(q)))
}

# Get the LD adjacent variant by max alpha's index and LD_thres
get_LD_adj = function(alpha, LD, LD_thres){
  max_index = which.max(alpha)
  ##print("MAX var")
  ##print(max_index)
  LD_vec = LD[max_index, ]
  index = which(abs(LD_vec) > LD_thres)
  return(index)
}

get_non_zero_effects_proxy = function(alpha, LD, dist_thres, LD_thres) {
  # If there are signal in l effect, then set sig_flag = true, then get the variants in LD.
  sig_flag = test_flat(alpha, dist_thres  = dist_thres)
  if(sig_flag){
    indices = get_LD_adj(alpha, LD, LD_thres = LD_thres)
    return(unique(indices))
  }else{
    return(NULL)
  }
}

# use weight to adjust the region to be adjusted
XtX_update = function(XtX, c_index, N, weight){
    # Construct weight matrix, W1 for original matrix, W2 for identity matrix 
    w1 <- rep(1, nrow(XtX))
    w2 <- rep(0, nrow(XtX))
    w1[c_index] <- sqrt(1 - weight)
    w2[c_index] <- sqrt(N * weight)
    W2 <- diag(w2)
    LD_adj <- t(XtX * w1) * w1 + W2^2
    return(LD_adj)
}


distance_calc <- function(weight, c_index, Xty, XtX, post_mean) {

  N = attr(XtX,"d")[1]
  XtX_sub = XtX[c_index, c_index]
  eigenld = eigen(XtX_sub, symmetric = TRUE)
  eigenld$values[eigenld$values < 1e-8] = 0

  tmp = (1 - weight) * eigenld$values + weight * diag(XtX_sub)
  XtX_adj = XtX_update(XtX, c_index, N, weight)
  include_index = which(tmp > 0)
  ztv = crossprod((Xty - XtX_adj %*% post_mean)[c_index], eigenld$vectors[,include_index]) # 1 x p
  distance = 0.5 * sum(log(tmp[include_index])) + 
    0.5 * tcrossprod(ztv / (tmp[include_index]), ztv)
  return(distance)
}

solve_weight <- function(c_index, Xty, XtX, post_mean) {
  estimate <- optim(par = 0.2, fn = distance_calc, 
    c_index = c_index, Xty = Xty, XtX = XtX, post_mean = post_mean,
    method = "Brent", lower = 0, upper = 1, control = list(reltol = 1e-5))$par
  return(estimate)
}

test_outlier = function(Xty, XtX, post_mean, c_index, method = "IQR"){
    if(method == "IQR"){
        expect_z = (XtX %*% post_mean)[c_index]
        diff_z = Xty[c_index] - expect_z
        Q1 <- quantile(diff_z, 0.25)
        Q3 <- quantile(diff_z, 0.75)
        IQR <- Q3 - Q1
        # Calculate bounds to filter outliers --- not necessary
        lower_bound <- Q1 - 3 * IQR
        upper_bound <- Q3 + 3 * IQR
        outliers <- which((diff_z < min(lower_bound, -2) | diff_z > max(upper_bound, 2)))
        # if (length(outliers) > 0) {
        #   print("outlier diff z summary:")
        #   print(summary(diff_z[outliers]))
        # }
    }else{
        dentist_s = function(z1,z2,r12) {
          (z1 - r12 * z2)^2 / (1 - r12^2)
        }

        is_sign_flip = function(z1,z2,r12) {
          ifelse(sign(z1) * sign(z2) * sign(r12) < 0, TRUE, FALSE)
        }
        r2=0.6
        p=1E-4
        chisq_cutoff = qchisq(1-p, df = 1)
        max_idx = which.max(abs(Xty[c_index]))
        Xty_max = Xty[c_index][max_idx]
        r_max = XtX[c_index, c_index][max_idx,]
        stats_filter = sapply(1:length(c_index), function(i) dentist_s(Xty_max, Xty[c_index][i], r_max[i]))
        stats_filter = (stats_filter > chisq_cutoff)
        r2_filter = sapply(1:length(c_index), function(i) XtX[c_index, c_index][max_idx, i]^2)
        r2_filter = (r2_filter > r2)
        sign_filter = sapply(1:length(c_index), function(i) is_sign_flip(Xty_max, Xty[c_index][i], r_max[i]))
        combined_filter = (stats_filter & r2_filter) | (sign_filter & r2_filter)
        outliers = which(combined_filter)
    }
    return(outliers)
}

outlier_removal = function(s, outlier){
  s$correct_zR_discrepancy$outlier_index = union(s$correct_zR_discrepancy$outlier_index, outlier)  
  # print("=============Outlier=============")
  # print(s$correct_zR_discrepancy$outlier_index)
  s$pi[s$correct_zR_discrepancy$outlier_index] = 0
  s$pi[s$pi>0]= s$pi[s$pi>0] / sum(s$pi)
  return(s)
}

find_index = function(s){
  mu_sum = rowSums(s$mu)
  keep_index = which(mu_sum != 0)
  return(keep_index)
}