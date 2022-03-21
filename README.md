# Graphical R2D2 (GR2D2)
## Description
Draw Monte Carlo samples from the posterior distribution based on the graphical R2D2 prior, to estimate the precision matrix for multivariate Gaussian data.

## Input explanations
`S`: Sample covariance matrix * n with dimension p by p

`n`: Sample size

`burnin`: Number of MCMC burnins

`nmc`: Number of MCMC saved samples

## Usage
1. Obtain the estimation 

`GR2D2_hubs <- GR2D2_E(hubs_S, n, burnin, nmc, eps = 1e-20)`

`omega_save <- GR2D2_hubs[["omega_save"]]`

2. Obtain the mean values based on nmc samples

`omega_save_mean <- apply(omega_save, c(1,2), mean)`

3. Use credible intervals to find observations

`omega_save_adjusted <- adjusted_fcn(dim_p, omega_save, nmc, CI)`

## Output of function `GR2D2_E`
`omega_save`: A (p by p by nmc) matrices of saved posterior samples of precision matrix.

`psi_save`: A (p by p by nmc) vector of saved samples of psi

`phi_save`: A (p by p by nmc) vector of saved samples of phi (local tuning parameter)

`w_save` : A (1 by nmc) vector of saved samples of w (global tuning parameter)

## Examples
See GR2D2_samplecode.R, in which we provide three different precision matrix structures.
