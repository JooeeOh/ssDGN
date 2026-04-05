rm(list = ls())
source("function.R")
library(huge)
library(MASS)
library(corpcor)
library(ggplot2)
library(tictoc)

# setup
T.lambda <- c(2, 2.1, 2.2)
T.alpha  <- seq(from = 0.4, to = 0.6, length = 3)
T.h      <- c(1, 1.1, 1.2)
s0 <- 1e-3
eps_var <- 1e-6

n        <- 30
p        <- 50
num_sub  <- 10
p_sub    <- p / num_sub           
num_pm   <- 100

group_size   <- n / 3              
group_labels <- rep(c("sensitive","moderate","resistant"), each = group_size)
group_types  <- c(sensitive = "band", moderate = "scale-free", resistant = "random")
strt <- "random"
sub_labels   <- rep(paste0("sub", 1:num_sub), each = p_sub)

mu_sen <- 0
mu_mod <- 5
mu_res <- 10

v_min <- 0
v_max <- 3

set.seed(145182)
m <- sort(runif(n, 0, 3))
X <- generate3_X(n = n, v_min = v_min, v_max = v_max, mu_sen = mu_sen, mu_mod = mu_mod, mu_res = mu_res, structure = strt)

parameterMatrix <- NULL
for(i in 1:length(T.lambda)){
  for(j in 1:length(T.alpha)){
    for(l in 1:length(T.h)){
      parameterMatrix <-rbind(parameterMatrix,c(T.lambda[i],T.alpha[j],T.h[l]))
    }
  }
}

opt <- vector("list", num_sub)
for (net in 1:num_sub) {
  opt[[net]] <- vector("list", p_sub)
  for (i in 1:p_sub) {
    opt[[net]][[i]] <- matrix(NA, (n*2/3), 3)
  }
}

kl_vec <- vector("list", num_sub)
kl_ratio <- rep(NA, num_sub)
cov_list_list <- vector("list", num_sub)
mu_list_list <- vector("list", num_sub)

for (net in 1:num_sub) {
  cov_list_list[[net]] <- vector("list", (n*2/3))
  mu_list_list[[net]] <- vector("list", (n*2/3))
  for (i in 1:n) {
    cov_list_list[[net]][[(n*2/3)]] <- matrix(NA, p_sub, p_sub)
    mu_list_list[[net]][[(n*2/3)]] <- rep(NA, p_sub)
  }
}
######################################
#### estimation
rg <- which(group_labels == "sensitive" | group_labels == "resistant")

for (net in 1:num_sub) {
  # 초기화
  beta_list <- vector("list", (n*2/3))
  mu_list   <- vector("list", (n*2/3))
  sig2_list <- vector("list", (n*2/3))
  
  for (u in 1:(n*2/3)) {
    beta_list[[u]] <- matrix(NA, nrow = p_sub, ncol = p_sub - 1)
    mu_list[[u]] <- rep(NA, p_sub)
    sig2_list[[u]] <- rep(NA, p_sub)
  }
  
  # X for subnetwork
  X_sub <- X[, ((net - 1)*p_sub + 1):(net*p_sub)]
  
  # hyperparameter selection
  for (l in 1:p_sub) {
    r <- X_sub[, -l]
    t <- X_sub[, l]
    
    PARA <- matrix(NA, ncol = ncol(parameterMatrix) + (n*2/3), nrow = nrow(parameterMatrix))
    PARA[, 1:3] <- parameterMatrix
    
    for (tp in 1:nrow(PARA)) {
      AIC_vec <- rep(NA, (n*2/3))
      for (rg_i in 1:(n*2/3)) {
        i <- rg[rg_i]
        K <- diag(exp(-1 / PARA[tp,3] * (m - m[i])^2))
        K.r <- sqrt(K) %*% r
        K.t <- sqrt(K) %*% t
        
        fit1 = HDeconometrics::ic.glmnet(K.r, K.t, crit = "aic", alpha = PARA[tp,2], lambda = PARA[tp,1])
        AIC_vec[rg_i] <- fit1$ic[2]
      }
      PARA[, -c(1:3)] <- AIC_vec
    }
    idx <- apply(PARA[, -c(1:3)], 2, which.min)
    opt[[net]][[l]] <- parameterMatrix[idx, ]
    
    for (rg_i in 1:(n*2/3)) {
      i <- rg[rg_i]
      K <- diag(exp(-1 / opt[[net]][[l]][rg_i, 3] * (m - m[i])^2))
      K.r <- sqrt(K) %*% r
      K.t <- sqrt(K) %*% t
      
      fit <- glmnet::glmnet(K.r, K.t, lambda = opt[[net]][[l]][rg_i, 1], alpha = opt[[net]][[l]][rg_i, 2], standardize = FALSE)
      b <- as.numeric(fit$beta)
      
      res <- as.numeric(K.t - K.r %*% b)
      s2 <- sum(res^2) / (n - 1)
      # w  <- exp(-(m - m[i])^2) 
      # s2 <- sum(w * res^2) / sum(w)
      
      beta_list[[rg_i]][l, ] <- b
      sig2_list[[rg_i]][l] <- s2
    }
  }
  
  # Covariance estimation
  cov_list <- vector("list", (n*2/3))
  for (rg_i in 1:(n*2/3)) {
    i <- rg[rg_i]
    Bfull <- make_square(beta_list[[rg_i]])
    Omega <- matrix(NA, p_sub, p_sub)
    for (l in 1:p_sub) {
      if (sig2_list[[rg_i]][l] < eps_var){
        sig2_list[[rg_i]][l] <- eps_var
      }
      Omega[l, l] <- 1 / sig2_list[[rg_i]][l]
      Omega[l, -l] <- - Bfull[l, -l] * Omega[l, l]
    }
    
    Omega <- (Omega + t(Omega)) / 2
    Omega <- Omega + diag(1e-5, p_sub)
    if (!is.positive.definite(Omega)) {
      Omega <- corpcor::make.positive.definite(Omega)
    }
    if (!is.positive.definite(Omega)) {
      Omega <- corpcor::make.positive.definite(Omega)
    }
    Sigma <- tryCatch(solve(Omega), error = function(e) chol2inv(chol(Omega)))
    # Sigma <- solve(Omega)
    Sigma <- (Sigma + t(Sigma)) / 2
    if (!corpcor::is.positive.definite(Sigma)) {
      Sigma <- corpcor::make.positive.definite(Sigma)
    }
    cov_list[[rg_i]] <- Sigma
  }
  
  mu_list <- lapply(rg, function (i) {
    X_sub[i, ]
  })
  
  mu_list_list[[net]] <- mu_list
  cov_list_list[[net]] <- cov_list
  
  ########################################## kl
  
  ref_idx <- length(rg)
  mu_ref <- mu_list[[ref_idx]]
  sigma_ref <- cov_list[[ref_idx]]
  
  kl_vec[[net]] <- sapply(
    1:(length(rg) - 1),
    function (i) KL(mu_ref, sigma_ref,
                    mu_list[[i]], cov_list[[i]])
  )
  
  kl_ratio[net] <- mean(kl_vec[[net]][1:group_size]) / mean(kl_vec[[net]][(group_size + 1):(length(rg) - 1)])
  # cat(sprintf("Estimation for subnetwork %d\n", net))
}

# toc()
par(mfrow = c(2, 5))
for (net in 1:num_sub) {
  plot(1:(length(rg) - 1), kl_vec[[net]],
       main = paste0("subnetwork ", net),
       ylab = "KL divergence", xlab = "sample index")
  abline(v = 40, col = "gray", lty = 2)
}
############################################# kl pm
pm_df <- matrix(NA, ncol = num_sub, nrow = num_pm)
for (net in 1:num_sub) {
  for (pm in 1:num_pm) {
    pm_idx_res <- sample(setdiff(1:length(rg), ref_idx), group_size - 1)
    pm_idx_sen <- setdiff(setdiff(1:length(rg), ref_idx), pm_idx_res)
    
    kl_res_pm <- sapply(pm_idx_res, function (i) {
      KL(mu_ref, sigma_ref,
         mu_list_list[[net]][[i]], cov_list_list[[net]][[i]])
    })
    
    kl_sen_pm <- sapply(pm_idx_sen, function (i) {
      KL(mu_ref, sigma_ref,
         mu_list_list[[net]][[i]], cov_list_list[[net]][[i]])
    })
    pm_df[pm, net] <- mean(kl_sen_pm) / mean(kl_res_pm)
  }
  # cat(sprintf("Permutation completed for subnetwork %d\n", net))
}

# pvalue
pvalue <- rep(NA, num_sub)
for (net in 1:num_sub) {
  pvalue[net] <- length(which(kl_ratio[net] < pm_df[, net]))
  pvalue[net] <- pvalue[net] / 200
}

pvalue



