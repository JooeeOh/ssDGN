# Workflow
### Step 1. Generate data
- Set the simulation parameters.
- Generate the modulator vector and expression matrix.
- Split variables into subnetworks.

### Step 2. Build the hyperparameter grid
- Create all combinations of:
  - `lambda`
  - `alpha`
  - `h`

### Step 3. Fit sample-specific models
- For each subnetwork and node:
  - fit kernel-weighted elastic net models
  - evaluate all hyperparameter combinations by AIC
  - select the best combination for each target sample
- Refit the model using the selected hyperparameters.

### Step 4. Estimate covariance structures
- Construct the precision matrix from the fitted coefficients and residual variances.
- Convert it to a covariance matrix for each target sample.

### Step 5. Compare distributions
- Compute KL divergence values between the reference sample and each non-reference sample.
- Calculate `kl_ratio` from the computed KL values.

### Step 6. Perform permutation testing
- Generate permutation-based KL ratios.
- Compute empirical p-values for each subnetwork.

---

# `KL()` function

### Arguments
- `mu_a`: estimated mean vector for the reference sample
- `sigma_a`: estimated covariance matrix for the reference sample
- `mu_b`: estimated mean vector for a non-reference sample
- `sigma_b`: estimated covariance matrix for a non-reference sample
- `eps`: stability parameter included in the function signature

### Returns
- a numeric value representing the KL divergence from the reference sample distribution to a non-reference sample distribution

