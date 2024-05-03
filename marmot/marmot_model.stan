functions {

//functions are used fom prior work by
//Dan Schrage (https://gitlab.com/dschrage/rcovreg)
  
  real sum_square_x(matrix x, int i, int j) {
    int j_prime;
    real sum_x = 0;
    if(j==1) return(sum_x);

    j_prime = 1;
    while(j_prime < j) {
      sum_x += x[i,j_prime]^2;
      j_prime += 1;
    }
    return(sum_x);
  }
  
  matrix lkj_to_chol_corr(row_vector constrained_reals, int ntrait) {
    int z_counter;
    matrix[ntrait,ntrait] x;

    z_counter = 1;
    x[1,1] = 1;
    for(j in 2:ntrait) {
      x[1,j] = 0;
    }
    for(i in 2:ntrait) {
      for(j in 1:ntrait) {
        if(i==j) {
          x[i,j] = sqrt(1 - sum_square_x(x, i, j));
        } else if(i > j) {
          x[i,j] = constrained_reals[z_counter]*sqrt(1 - sum_square_x(x, i, j));
          z_counter += 1;
        } else {
          x[i,j] = 0;
        }
      }
    }
    return(x);
  }
}

data {
  int<lower=1> N; //total number of observations _ growth
  int<lower=1> M; //total number of observations - fecundity
  int<lower=1> C; //total number of environmental contexts (years)
  int<lower=1> I; //total number of subjects
  int<lower=1> D; //total number of traits/dimensions
  int<lower=0> P_y; //total number of environmental predictors (+ intercept) on correlation
  int<lower=0> P_g; //total number of environmental predictors (+ intercept) on growth
  int<lower=0> P_f; //total number of environmental predictors (+ intercept) on fecundity
  
  array[N] int<lower=0> id_g; //index linking observations to individuals - growth
  array[N] int<lower=0> c_id_g; //index linking observations to contexts - growth
  array[N] int<lower=0> idc_g; //index linking individuals to positions in cmat - growth
  array[N] int<lower=0> id_g_lm; //index linking individuals to positions in cmat - growth
  
  array[M] int<lower=0> id_f; //index linking observations to individuals - fecundity
  array[M] int<lower=0> c_id_f; //index linking observations to contexts - fecundity
  array[M] int<lower=0> idc_f; //index linking individuals to positions in cmat - fecundity
  array[M] int<lower=0> id_f_lm; //index linking individuals to positions in cmat - fecundity
  
  matrix[C,P_y] X; //environmental predictor matrix (+ intercept) on correlation
  matrix[N,P_g] X_g; //environmental predictor matrix (+ intercept) on growth
  matrix[M,P_f] X_f; //environmental predictor matrix (+ intercept) on fecundity
  matrix[I,I] A; //relatedness matrix
  
  int<lower=1> cm; //max number of individuals observed in a context
  array[C, cm] int cmat; //matrix with all individuals observed in each context (row)
  array [C] int<lower=0> cn; //count of individuals observed per context
  int<lower=1> cnt; //total number of individuals across contexts
  
  array[N] real growth; //offspring mass data
  array[M] int productivity; //fecundity data
}

transformed data{
  matrix[I, I] LA = cholesky_decompose(A);
  int ncor = (D*(D-1))/2; //unique cov/cor parameters
  // Compute, thin, and then scale QR decomposition
  matrix[C, P_y] Q = qr_thin_Q(X) * sqrt(C-1);
  matrix[P_y, P_y] R = qr_thin_R(X) / sqrt(C-1);
  matrix[P_y, P_y] R_inv = inverse(R);
  
  matrix[N, P_g] Q_g = qr_thin_Q(X_g) * sqrt(N-1);
  matrix[P_g, P_g] R_g = qr_thin_R(X_g) / sqrt(N-1);
  matrix[P_g, P_g] R_inv_g = inverse(R_g);
  
  matrix[M, P_f] Q_f = qr_thin_Q(X_f) * sqrt(M-1);
  matrix[P_f, P_f] R_f = qr_thin_R(X_f) / sqrt(M-1);
  matrix[P_f, P_f] R_inv_f = inverse(R_f);
}

parameters { 
  //fixed effects
  vector[P_g] B_mq_g; //RN of means - growth
  vector[P_f] B_mq_f; //RN of means - fecundity
  matrix[P_y, ncor] B_cpcq; //RN of canonical partial correlations

  //random effects
  matrix[cnt, D] Z_G; //all context-specific additive values
  real<lower=0> sd_E; //residual standard deviation (within litter variance) - growth
  array[C] vector<lower=0>[D] sd_G; //sd of ind effects
  
  real<lower=0> sd_g; // sd from season random intercept - growth
  real<lower=0> sd_f; // sd from season random intercept - fecundity
  vector[C] z_season_g;
  vector[C] z_season_f;
  
}

model {
  //predicted values from reaction norms
  //growth
  vector[N] mu_growth =  Q_g * B_mq_g;
  
  //fecundity
  vector[M] mu_fecundity =  Q_f * B_mq_f;
                       
  //correlations (expressed as canonical partial correlations)
  matrix[C, ncor] cpc_G = tanh(Q * B_cpcq);
  
  //scale non-focal univariate random effects
  vector[C] re_season_g = z_season_g * sd_g;
  vector[C] re_season_f = z_season_f * sd_f;

 //initialize mean linear predictors
  vector[N] mu_g = mu_growth[id_g_lm] + re_season_g[c_id_g];
  vector[M] mu_f = mu_fecundity[id_f_lm] + re_season_f[c_id_f];

  //scale context-specific multivariate additive genetic effects
  matrix[cnt, D] mat_G;
  int pos = 1; //keep track of position 1:cnt
  for(c in 1:C){
      mat_G[pos:(pos+cn[c]-1)] = 
      LA[cmat[c,1:cn[c]],cmat[c,1:cn[c]]] * Z_G[pos:(pos+cn[c]-1)] * diag_pre_multiply(sd_G[c],lkj_to_chol_corr(cpc_G[c], D))';
      pos = pos + cn[c];   
  }
        
//add context-specific genetic effects to linear predictors
  for(n in 1:N){
  mu_g[n]  += col(mat_G,1)[idc_g[n]];
  }
  
  for(m in 1:M){
  mu_f[m]  += col(mat_G,2)[idc_f[m]]; 
  }        
                  
//likelihood 
  growth ~ normal(mu_g, sd_E);
  productivity ~ poisson(exp(mu_f));

//priors
  to_vector(B_mq_g) ~ normal(0,1);
  to_vector(B_mq_f) ~ normal(0,1);
  B_cpcq[1,1] ~ normal(0,0.5);
  B_cpcq[2,1] ~ normal(0,0.5);
  B_cpcq[3,1] ~ normal(0,0.5);
  B_cpcq[4,1] ~ normal(0,0.5);
  to_vector(Z_G) ~ std_normal();
  
  z_season_g ~ std_normal();
  z_season_f ~ std_normal();
  sd_g ~ exponential(2);
  sd_f ~ exponential(2);

  sd_E ~ exponential(2);
  
  for(c in 1:C){
  sd_G[c] ~ exponential(2);
  }
}

generated quantities{
  vector[P_g] B_m_g; //mean RN parameters for X
  vector[P_f] B_m_f; //mean RN parameters for X
  matrix[P_y,ncor] B_cpc; //partial correlation RN parameters for X

  B_m_g= R_inv_g * B_mq_g;
  B_m_f= R_inv_f * B_mq_f;

  for(d in 1:ncor){
    B_cpc[,d]= R_inv * B_cpcq[,d];
    }

  // Posterior predictive check for fecundity and offpsring mass
  vector[N] mu_growth_bis =  Q_g * B_mq_g;
  vector[M] mu_fecundity_bis =  Q_f * B_mq_f;
  matrix[C, ncor] cpc_G_bis = tanh(Q * B_cpcq);
  vector[C] re_season_g_bis = z_season_g * sd_g;
  vector[C] re_season_f_bis = z_season_f * sd_f;
  vector[N] mu_g_bis = mu_growth_bis[id_g_lm] + re_season_g_bis[c_id_g];
  vector[M] mu_f_bis = mu_fecundity_bis[id_f_lm] + re_season_f_bis[c_id_f];
  matrix[cnt, D] mat_G_bis;
  int pos = 1; //keep track of position 1:cnt
  for(c in 1:C){
      mat_G_bis[pos:(pos+cn[c]-1)] = 
      LA[cmat[c,1:cn[c]],cmat[c,1:cn[c]]] * Z_G[pos:(pos+cn[c]-1)] * diag_pre_multiply(sd_G[c],lkj_to_chol_corr(cpc_G_bis[c], D))';
      pos = pos + cn[c];   
  }
  array[N] real y_rep_g;
  for(n in 1:N){
  mu_g_bis[n]  += col(mat_G_bis,2)[idc_g[n]]; 
  y_rep_g[n] = normal_rng(mu_g_bis[n], sd_E);
  }
  array[M] int y_rep_f;
  for(m in 1:M){
  mu_f_bis[m]  += col(mat_G_bis,2)[idc_f[m]]; 
  y_rep_f[m] = poisson_rng(exp(mu_f_bis[m]));
  }     

}
