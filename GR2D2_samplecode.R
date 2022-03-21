
rm(list = ls())
source("./GR2D2_fcn.R")

set.seed(7)
# Generate different precision matrix structures
Hubs_precision_matrix <- function(dim_p, group_num, fixed_val, num_fix_values){
  group_dim <- dim_p/group_num
  hubs_precision_matrix <- matrix(0, nrow = dim_p, ncol = dim_p)
  
  for (i in 1:group_num) {
    group_i <- diag(group_dim)
    group_inx <- c(1:group_num)
    group_inx <- group_inx[group_inx!=i]
    chosen_rows <- sample(group_inx, num_fix_values, replace = FALSE)
    for (j in chosen_rows) {
      group_i[i, j] <- fixed_val
      group_i[j, i] <- fixed_val
    }
    hubs_precision_matrix[(group_dim*(i-1)+1):(group_dim*i), (group_dim*(i-1)+1):(group_dim*i)] <- group_i
  }
  return(hubs_precision_matrix)
}

Cliques_precision_matrix <- function(dim_p, group_num, num_fixed_val, fixed_val){
  group_dim <- dim_p/group_num
  cliques_precision_matrix <- matrix(0, nrow = dim_p, ncol = dim_p)
  
  for (i in 1:group_num) {
    group_i <- diag(group_dim)
    chosen_rows <- sample(c(1:group_num), num_fixed_val, replace = FALSE)
    for (j in 1:(length(chosen_rows)-1)) {
      for (k in (j+1):length(chosen_rows)) {
        group_i[chosen_rows[j], chosen_rows[k]] <- fixed_val
        group_i[chosen_rows[k], chosen_rows[j]] <- fixed_val
      }
    }
    cliques_precision_matrix[(group_dim*(i-1)+1):(group_dim*i), (group_dim*(i-1)+1):(group_dim*i)] <- group_i
  }
  return(cliques_precision_matrix)
}

Generate_S <- function(precision_matrix, sample_size){
  Sigma <- spdinv(precision_matrix)
  Sigma_chol <- chol(Sigma)
  mu <- c(rep(0, dim_p))
  dim_p <- nrow(precision_matrix)
  sd_normal <- t(mvrnorm(n = sample_size, mu, diag(dim_p)))
  Y <- t(Sigma_chol%*%sd_normal)#n*p
  S <- t(Y)%*%Y
  return(S)
}

burnin <- 5
nmc <- 5
CI <- 0.95

dim_p <- 100
n <- 50

# Generate hubs precision matrix
group_num_hubs <- 10
fixed_val_hubs <- 0.25
num_fix_values <- 4
hubs_precision <- Hubs_precision_matrix(dim_p, group_num_hubs, fixed_val_hubs, num_fix_values)
hubs_S <- Generate_S(hubs_precision, n)

group_num_cliques <- 10
num_fixed_val <- 3

# Generate cliques positive precision matrix
fixed_val_cliques1 <- -0.45
cliques1_precision <- Cliques_precision_matrix(dim_p, group_num_cliques, num_fixed_val, fixed_val_cliques1)
cliques1_S <- Generate_S(cliques1_precision, n)
# Generate cliques negative precision matrix
fixed_val_cliques2 <- 0.75
cliques2_precision <- Cliques_precision_matrix(dim_p, group_num_cliques, num_fixed_val, fixed_val_cliques2)
cliques2_S <- Generate_S(cliques2_precision, n)


# For hubs
GR2D2_hubs <- GR2D2_E(hubs_S, n, burnin, nmc, eps = 1e-20)
omega_save <- GR2D2_hubs[["omega_save"]]
omega_save_mean <- apply(omega_save, c(1,2), mean)
omega_save_adjusted <- adjusted_fcn(dim_p, omega_save, nmc, CI)

# For cliques1
GR2D2_cliques1 <- GR2D2_E(cliques1_S, n, burnin, nmc, eps = 1e-20)
omega_save <- GR2D2_cliques1[["omega_save"]]
omega_save_mean <- apply(omega_save, c(1,2), mean)
omega_save_adjusted <- adjusted_fcn(dim_p, omega_save, nmc, CI)

# For cliques2
GR2D2_cliques2 <- GR2D2_E(cliques2_S, n, burnin, nmc, eps = 1e-20)
omega_save <- GR2D2_cliques2[["omega_save"]]
omega_save_mean <- apply(omega_save, c(1,2), mean)
omega_save_adjusted <- adjusted_fcn(dim_p, omega_save, nmc, CI)

