#' @export
processOutput = function(result, coef_name, cutoff, infection_ind, plots){
  time = result$time
  vh = result$vh
  se = sqrt(result$var_vh)
  lower = result$ci_lower
  upper = result$ci_upper
  cov_beta = result$Covariance
  beta = result$beta
  p = length(beta)
  n_type = ncol(vh)

  ## covariates
  covariates = matrix(data = NA, nrow = p, ncol = 7)
  colnames(covariates) = c("coef", "se(coef)", "z", "Pr(>|z|)", "exp(coef)", "lower .95", "upper .95")
  covariates[,1] = beta
  covariates[,2] = sqrt(diag(cov_beta))
  covariates[,3] = beta / sqrt(diag(cov_beta))
  covariates[,4] = 2 * (1 - pnorm(abs(covariates[,3])))
  covariates[,5] = exp(beta)
  covariates[,6] = exp(beta - qnorm(0.975) * sqrt(diag(cov_beta)))
  covariates[,7] = exp(beta + qnorm(0.975) * sqrt(diag(cov_beta)))
  rownames(covariates) = coef_name

  ## vaccine
  VE = NULL
  for (i in 1:n_type) {
    tmp_mat = matrix(data = 0, nrow = length(time), ncol = 5)
    colnames(tmp_mat) = c("time", "VE", "se", "lower .95", "upper .95")
    tmp_mat[,1] = time
    tmp_mat[,2] = vh[,i]
    tmp_mat[,3] = se[,i]
    tmp_mat[,4] = lower[,i]
    tmp_mat[,5] = upper[,i]
    if(i >= infection_ind){
      tmp_mat = tmp_mat[which(time>=cutoff),]
    }
    VE[[i]] = tmp_mat
  }

  ## plots
  plot_res = NULL
  if(plots){
    plot_res = VEplot(vaccine = VE, infection_ind = infection_ind)
  }

  output = NULL
  output[["covariates"]] = covariates
  output[["vaccine"]] = VE
  output[["plots"]] = plot_res
  return(output)
}
