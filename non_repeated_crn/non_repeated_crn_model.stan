functions {
  vector lower_elements(matrix M, int tri_size){
    int n = rows(M);
    int counter = 1;
    vector[tri_size] out;
    for (i in 2:n){
      for (j in 1:(i - 1)) {
        out[counter] = M[i, j];
        counter += 1;
      }
    }
    return out;
  }
}

data { 
  int<lower=1> N; //n of obs
  int<lower=1> I; //n of ind
  int<lower=1> Y; //n of years
  int<lower=1> D; //n of traits
  array[N] int id; //index of individuals
  array[N] int years; //index of years
  vector[N] climate; //env covariate
  vector[Y] climate_y; //env covariate - correlation submodel
  array[N] real growth; //growth data
  array[N] int productivity; //fecundity data
} 

transformed data {
  int<lower=0> off_diag = (D * (D - 1))/2;
}

parameters { 
  real mu_growth; //intercept - growth submodel
  real mu_productivity; //intercept - fecundity submodel
  real beta_g; //slope climate - growth submodel
  real beta_p; //slope climate - fecundity submodel
  real mu_r; //intercept - correlation submodel
  real beta_r; //slope climate - correlation submodel
  array[Y] corr_matrix[D] R; //year-specific corr matrices
  vector[Y] z_Y; //cluster avg corr
  real<lower=0> sd_Y; //sd of cluster avg
  array[Y] vector<lower=0>[D] sd_W; //sd of ind effects
  array[Y] matrix[I,D] z_W; //year-specific individual random effects
}

transformed parameters {
  //linear predictor of cluster-specific correlations constrained to [0,1]
  vector[Y] eta_Y = inv_logit(mu_r + beta_r * climate_y + sd_Y * z_Y);
  //(0,1) to (-1,1) scale for correlation
  vector[Y] r_Y = 2 * eta_Y - 1;
}

model {
  array[Y] matrix[D,D] P; //cluster cov matrices
  array[Y] matrix[I,D] W_Y; //community-specific individual-level effects
  
  //individual-level effect scaled by year-specific correlations

  for(y in 1:Y) {
    target += normal_lpdf(lower_elements(R[y], off_diag) | r_Y[y], 0.01); 
    //scaled covariance matrix of individual effects (for each year)
    P[y] = diag_matrix(sd_W[y]) * R[y] * diag_matrix(sd_W[y]);
    //scaled individual random effects (for each year)
    W_Y[y] = z_W[y] * diag_pre_multiply(sd_W[y], cholesky_decompose(R[y]))';
    to_vector(z_W[y]) ~ std_normal();
    sd_W[y] ~ exponential(2);
    }
    
  for(n in 1:N){  
    real prob_g = mu_growth + beta_g * climate[n] + col(W_Y[years[n]],1)[id[n]];
    real lambda_p = mu_productivity + beta_p * climate[n] + col(W_Y[years[n]],2)[id[n]];
    target += normal_lpdf(growth[n]| prob_g, 0.1);
    target += poisson_lpmf(productivity[n]| exp(lambda_p));
    }

  //priors
  mu_growth ~ normal(0,1);
  mu_productivity ~ normal(0,1);
  mu_r ~ normal(0,1);
  beta_g ~ normal(0,1);
  beta_p ~ normal(0,1);
  beta_r ~ normal(0,1);
  z_Y ~ std_normal();
  sd_Y ~ exponential(2);
}

