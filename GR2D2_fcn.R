# Dailin Gan, March 2022

# load necessary packages
library(abind)
library(stats)
library(dplyr)
library(MASS)
library(statmod)
library(GIGrvg)
library(Matrix)
library(Rcpp)
# library(Rfast)
library(pracma)
library(matrixcalc) # Check positive definiteness
library(matlib)

#  Sampling for graphical R2D2
GR2D2_E <- function(S, n, burnin, nmc, eps = 1e-20){
  # Input:
  #   S = Y'*Y: sample covariance matrix * n
  #   n: sample size
  #   burnin, nmc: number of MCMC burnins and saved samples
  
  # Output:
  #   omega_save: (p by p by nmc) matrices of saved posterior samples of precision matrix.
  #   psi_save: (p by p by nmc) vector of saved samples of psi
  #   phi_save: (p by p by nmc) vector of saved samples of phi (local tuning parameter)
  #   w_save: (1 by nmc) vector of saved samples of w (global tuning parameter)
  
  p <- nrow(S)
  print(p)
  omega_save <- array(numeric(), c(p, p, 0))
  psi_save <- array(numeric(), c(p, p, 0))
  phi_save <- array(numeric(), c(p, p, 0))
  w_save <- matrix(0, nrow = 1, ncol = nmc)
  sample_size <- n
  # Set initial values
  #   Initial values 1
  b = 0.5
  r = 1
  C = 1
  p_n = p*(p-1)/2
  alpha_pi = C/((p_n^(b/2))*(sample_size^(r*b/2))*log(sample_size))
  
  a = alpha_pi*p*(p-1)/2
  #   Initial values 2
  Omega <- diag(p)
  Sigma <- diag(p)
  Psi <- matrix(1, nrow = p, ncol = p)
  Phi <- matrix(1, nrow = p, ncol = p)
  w <- 1
  xi <- 1
  
  # Main iterations
  for (iter in 1: (burnin + nmc)) {
    print(iter)
    for (i in 1:p){
      # sample gamma
      s_ii <- S[i, i]
      gamma <- rgamma(1, shape = sample_size/2+1, rate = s_ii/2)
      
      # Value adjustment for Sigma
      Sigma[which(Sigma > 1/eps)] <- 1/eps  # New line------------------------
      Sigma[which(Sigma < -1/eps)] <- -1/eps  # New line------------------------
      Sigma[which(Sigma < eps & Sigma > 0)] <- eps  # New line------------------------
      Sigma[which(Sigma > -eps & Sigma < 0)] <- -eps  # New line------------------------
      
      sigma_vec_ii <- as.matrix(Sigma[-i, i])
      inverse_Omega_ii <- Sigma[-i, -i] - sigma_vec_ii %*% t(sigma_vec_ii)/Sigma[i, i]
      
      # sample beta
      Lambda <- Psi*Phi
      
      Lambda[which(Lambda > 1/eps)] <- 1/eps  # New line------------------------
      Lambda[which(Lambda < -1/eps)] <- -1/eps  # New line------------------------
      Lambda[which(Lambda < eps & Lambda > 0)] <- eps  # New line------------------------
      Lambda[which(Lambda > -eps & Lambda < 0)] <- -eps  # New line------------------------
      
      lambda_vec_ii <- Lambda[-i, i]
      Lambda_star <- diag(lambda_vec_ii)
      D_inv <- s_ii*inverse_Omega_ii + 2*solve(w*Lambda_star, tol = 1e-10000)
      
      # Value adjustment for D_inv
      # D_inv[which(D_inv > 1/eps)] <- 1/eps  # New line------------------------
      # D_inv[which(D_inv < eps)] <- eps  # New line------------------------

      if (!is.positive.definite(D_inv)){ # New line------------------------
        D_inv_PD <- nearPD(D_inv)
        D_inv <- D_inv_PD$mat %>% as.matrix()
      }
      
      D_inv.chol <- chol(D_inv)
      
      D_inv[which(D_inv > 1/eps)] <- 1/eps  # New line------------------------
      D_inv[which(D_inv < eps & D_inv > 0)] <- eps  # New line------------------------
      D_inv[which(D_inv > -eps & D_inv < 0)] <- -eps  # New line------------------------
      D_inv[which(D_inv < -1/eps )] <- -1/eps  # New line------------------------
      
      mu_i <- -mldivide(D_inv, S[-i, i])
      
      # Value adjustment for D_inv
      # mu_i <- -solve(D_inv, tol = 1e-10000)%*%S[-i, i] # New line------------------------
      
      # mu_i <- -Ginv(D_inv, tol = eps)%*%S[-i, i] # Use generalized inverse------------------------
      
      D_inv.chol[which(D_inv.chol > 1/eps)] <- 1/eps  #s New line------------------------
      D_inv.chol[which(D_inv.chol < -1/eps)] <- -1/eps  #s New line------------------------
      D_inv.chol[which(D_inv.chol < eps & D_inv.chol > 0)] <- eps  # New line------------------------
      D_inv.chol[which(D_inv.chol > -eps & D_inv.chol < 0)] <- -eps  # New line------------------------
      
      beta_vec <- mu_i +  mldivide(D_inv.chol, mvrnorm(n = 1, rep(0, p-1), diag(p-1)))
      # beta_vec <- mu_i + Ginv(D_inv.chol, tol = eps)%*%mvrnorm(n = 1, rep(0, p-1), diag(p-1)) # Use generalized inverse------------------------
      
      Omega[i, i] <- gamma + t(beta_vec)%*%inverse_Omega_ii%*%beta_vec
      Omega[-i, i] <- beta_vec
      Omega[i, -i] <- beta_vec
      
      # update Psi
      for (row in 1:p) {
        Psi[row, i] <- rinvgauss(1, mean = sqrt(Phi[row, i]*w/2)/abs(Omega[row, i]), shape = 1)
        Psi[row, i] <- 1/Psi[row, i]
        Psi[i, row] <- 1/Psi[row, i]
      }
      diag(Psi) <- 1
      
      # update Sigma
      Sigma[-i, -i] <- inverse_Omega_ii + (inverse_Omega_ii%*%beta_vec)%*%t(inverse_Omega_ii%*%beta_vec)/gamma
      Sigma[-i, i] <- -(inverse_Omega_ii%*%beta_vec)/gamma
      Sigma[i, -i] <- -(inverse_Omega_ii%*%beta_vec)/gamma
      Sigma[i, i] <- 1/gamma
    }
    
    T_summary <- matrix(0, nrow = p, ncol = p)
    sum_T <- 0
    psi_mn <- 2*xi
    lambda_mn <- alpha_pi - 1/2
    
    tem <- 2*(Omega)^(2)/Psi
    tem[which(tem > 1/eps)] <- 1/eps # New lines----------
    tem[which(tem < -1/eps)] <- -1/eps # New lines----------
    tem[which(tem < eps & tem > 0)] <- eps # New lines----------
    tem[which(tem > -eps & tem < 0)] <- -eps # New lines----------
    
    for (m in 1:(p-1)) {
      for (n in (m+1):p) {
        Chi_mn <- tem[m, n]
        T_mn <- rgig(n = 1, lambda = lambda_mn, chi = Chi_mn, psi = psi_mn)
        T_summary[m, n] <- T_mn
        T_summary[n, m] <- T_mn
        sum_T <- sum_T + T_mn
      }
    }
    Phi <- diag(p) + T_summary/sum_T
    
    # sample w
    Chi_w <- 0
    for (m in 1:(p-1)) {
      for (n in (m+1):p) {
        Chi_w <- Chi_w + 2*(Omega[m, n])^(2)/(Phi[m, n]*Psi[m, n])
      }
    }
    psi_w <- 2*xi
    lambda_w <- a - p*(p-1)/4
    
    if(Chi_w > 1/eps) { # New lines----------
      Chi_w <- 1/eps
    } else {
      Chi_w <- Chi_w
    }
    
    w <- rgig(n = 1, lambda = lambda_w, chi = Chi_w, psi = psi_w)
    
    # sample xi
    xi <- rgamma(1, shape = (a + b), rate = w+1)
    
    # save Omega, psi, phi, w
    if (iter > burnin){
      omega_save <- abind(omega_save, Omega)
      psi_save <- abind(psi_save, Psi)
      phi_save <- abind(phi_save, Phi)
      w_save[iter] <- w
    }
  }
  return(list("omega_save" = omega_save, "psi_save" = psi_save, "phi_save" = phi_save, "w_save" = w_save))
}


# Record observations based on credible intervals
adjusted_fcn <- function(p, estimate, nmc, CI_p = 0.95){
  adjusted <- matrix(0, nrow = p, ncol = p)
  lower_p <- (1-CI_p)/2
  upper_p <- 1-lower_p
  for (i in 1:nrow(adjusted)) {
    for (j in 1:ncol(adjusted)) {
      chosen_entry <- estimate[i, j, 1:nmc]
      q_lower <- quantile(chosen_entry, probs = lower_p)%>% as.numeric()
      q_upper <- quantile(chosen_entry, probs = upper_p)%>% as.numeric()
      if (q_lower < 0 & 0 < q_upper){
        adjusted[i,j] <- 0
      } else{
        adjusted[i,j] <- median(chosen_entry)
      }
    }
  }
  return(adjusted)
}

