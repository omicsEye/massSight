data {
  int<lower=1> N;         // Number of data points
  int<lower=1> K;         // Number of clusters
  matrix[N, 2] data;      // Data matrix (MZ and RT data)
}

parameters {
  simplex[K] mixing_proportions;  // Mixing proportions for each cluster
  vector[2] means[K];             // Means for each cluster
  cov_matrix[2] covariances[K];   // Covariances for each cluster
}

model {
  // Priors for the means and covariances
  for (k in 1:K) {
    means[k] ~ normal(0, 100);
    covariances[k] ~ inv_wishart(4, diag_matrix(rep_vector(100, 2)));
  }

  // Likelihood for each data point
  for (n in 1:N) {
    vector[K] log_cluster_probability;
    for (k in 1:K) {
      log_cluster_probability[k] = log(mixing_proportions[k]) + multi_normal_lpdf(data[n] | means[k], covariances[k]);
    }
    target += log_sum_exp(log_cluster_probability);  // Log-sum-exp trick for numerical stability
  }
}
