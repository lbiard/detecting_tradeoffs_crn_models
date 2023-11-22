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
  int<lower=1> N; //n of obs offspring - growth submodel
  int<lower=1> M; //n of obs - fecundity submodel
  int<lower=1> I; //n of individuals (number of mothers)
  int<lower=1> Y; //n of years
  int<lower=1> D; //n of traits
  array[N] int id_g; //index of individuals - growth submodel
  array[M] int id_f; //index of individuals - fecundity submodel
  array[N] int years_g; //index of years - growth sunmodel
  array[M] int years_f; //index of years - fecundity submodel
  vector[N] climate_g; //total winter snow - growth submodel
  vector[N] climate_g2; //max June temperature - growth submodel
  vector[N] age_g; //age of mother - growth submodel
  vector[N] mass_mom_g; // mass of mother - growth submodel
  vector[M] climate_f; //total winter snow - fecundity submodel
  vector[M] age_f; //age of mother - fecundity submodel
  vector[M] mass_mom_f; // mass of mother - fecundity submodel
  vector[Y] climate_y; //total winter snow - correlation submodel
  vector[Y] climate_y2; //max June temperature - correlation submodel
  array[N] real growth; //offspring mass data
  array[M] int productivity; //fecundity data
} 

transformed data {
  int<lower=0> off_diag = (D * (D - 1))/2;
}

parameters { 
  real mu_growth; //intercept - growth submodel
  real mu_productivity; //intercept - fecundity submodel
  real beta_g; //slope total winter snow - growth submodel
  real beta_gt; //slope max June temperature - growth submodel
  real beta_g2; //slope age of mother - growth submodel
  real beta_g3; //slope mass of mother - growth submodel
  real beta_f; //slope total winter snow - fecundity submodel
  real beta_f2; //slope age of mother - fecundity submodel
  real beta_f3; //slope mass of mother - fecundity submodel
  real mu_r; //intercept year - correlation submodel
  real beta_r; //slope total winter snow - correlation submodel
  real beta_rt; //slope max June temperature - correlation submodel
  real<lower=0> sigma; //within-litter variance
  array[Y] corr_matrix[D] R; //year-specific corr matrices
  vector[Y] z_Y; //cluster avg corr
  real<lower=0> sd_Y; //sd of cluster avg
  array[Y] vector<lower=0>[D] sd_W; //sd of ind effects
  array[Y] matrix[I,D] z_W; //year-specific individual random effects
}

transformed parameters {
  //linear predictor of cluster-specific correlations constrained to [0,1]
  vector[Y] eta_Y = inv_logit(mu_r + beta_r * climate_y + beta_rt * climate_y2 + sd_Y * z_Y);
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
    
  //offspring mass submodel
    
  for(n in 1:N){  
    real prob_g = mu_growth + beta_g * climate_g[n] + beta_gt * climate_g2[n] + beta_g2 * age_g[n] + beta_g3 * mass_mom_g[n] + col(W_Y[years_g[n]],1)[id_g[n]];
    target += normal_lpdf(growth[n]| prob_g, sigma);
    }

  //fecundity submodel

  for(m in 1:M){ 
    real lambda_p = mu_productivity + beta_f * climate_f[m] + beta_f2 * age_f[m] + beta_f3 * mass_mom_f[m] + col(W_Y[years_f[m]],2)[id_f[m]];
    target += poisson_lpmf(productivity[m]| exp(lambda_p));
  }
  
  //priors
  mu_growth ~ normal(0,1);
  mu_productivity ~ normal(0,1);
  mu_r ~ normal(0,1);
  beta_g ~ normal(0,1);
  beta_gt ~ normal(0,1);
  beta_g2 ~ normal(0,1);
  beta_g3 ~ normal(0,1);
  beta_f ~ normal(0,1);
  beta_f2 ~ normal(0,1);
  beta_f3 ~ normal(0,1);
  beta_r ~ normal(0,1);
  beta_rt ~ normal(0,1);
  z_Y ~ std_normal();
  sd_Y ~ exponential(2);
  sigma ~ exponential(2);
}

