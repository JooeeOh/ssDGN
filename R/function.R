generate3_X <- function (n = n, v_min, v_max, mu_sen, mu_mod, mu_res, structure = "random") {
  
  group_levels <- c("sensitive","moderate","resistant")

  Omega_list <- vector("list", num_sub)
  for (net in 1:num_sub) {
    Omega_list[[net]] <- vector("list", length(group_levels))
    for (g_idx in 1:length(group_levels)) {
      gname <- group_levels[g_idx]
      sim_g <- huge::huge.generator(n = n/3, d = p_sub, graph = group_types[gname], verbose = FALSE)
      omega <- as.matrix(sim_g$omega)
      # diag(theta) <- 0
      # theta <- (theta + t(theta)) / 2
      Omega_list[[net]][[g_idx]] <- omega
    }
  }
  
  v_vec <- sort(runif(n, v_min, v_max))
  
  # 앞 절반 서브네트워크
  X1 <- matrix(NA, n, (p/2))
  for (net in 1:(num_sub/2)) {
    X_new <- matrix(NA, n, p_sub)
    for (i in 1:n) {
      g_idx <- match(group_labels[i], group_levels)
      omega <- Omega_list[[net]][[g_idx]]
      v_i   <- v_vec[i]
      Omega_i <- corpcor::make.positive.definite(v_i * omega, tol = 1e-10)
      
      Sigma_i <- tryCatch(solve(Omega_i), error = function(e) chol2inv(chol(Omega_i)))
      Sigma_i <- (Sigma_i + t(Sigma_i))/2
      if (!corpcor::is.positive.definite(Sigma_i)) {
        Sigma_i <- corpcor::make.positive.definite(Sigma_i, tol = 1e-10)
      }
      mu_i <- switch(group_labels[i],
                     "sensitive" = rep(mu_sen, p_sub),
                     "moderate"  = rep(mu_mod, p_sub),
                     "resistant" = rep(mu_res, p_sub))
      X_new[i, ] <- MASS::mvrnorm(1, mu = mu_i, Sigma = Sigma_i)
    }
    X1[, ((net - 1) * p_sub + 1):(net * p_sub)] <- X_new
  }
  
  # 뒤 절반 서브네트워크
  sim <- huge::huge.generator(n = n, d = p_sub, graph = structure, verbose = FALSE)
  Sigma_com <- solve(sim$omega)
  X2 <- do.call(cbind, replicate(num_sub/2,
                                 MASS::mvrnorm(n, mu = rep(mu_mod, p_sub), Sigma = Sigma_com), simplify = FALSE))
  cbind(X1, X2)
}
                          
make_square <- function (mx) {
  p <- nrow(mx)
  res <- matrix(0, nrow = p, ncol = p)
  for (i in 1:p) {
    res[i, ] <- append(mx[i, ], values = 0, after = (i - 1))
  }
  return(res)
}

KL <- function(mu_a, sigma_a, mu_b, sigma_b, eps = 1e-10){
  d <- length(mu_a)
  
  kld <- log(det(sigma_b) / det(sigma_a)) -
    d +
    sum(diag(solve(sigma_b) %*% sigma_a)) +
    t(mu_a -mu_b) %*% solve(sigma_b) %*% (mu_a -mu_b)
  return(kld)
}
