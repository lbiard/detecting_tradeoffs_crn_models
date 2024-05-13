# Last modified 07/05/2024

# This R script can be used to simulate individual-based datasets with an intraindividual 
# tradeoff between a mother's fecundity and its growth. The expression of this tradeoff is 
# dependent on the environmental quality at each time step (continuous variation
# of the among-individual correlation depends on the environment).

# The code simulate 4 datasets (only the first one is analysed in the manuscript):
# 1 - fast life history, 20 years of data, 25 new individuals per year
# 2 - fast life history, 20 years of data, 50 new individuals per year
# 3 - slow life history, 20 years of data, 25 new individuals per year
# 4 - slow life history, 20 years of data, 50 new individuals per year

# These datasets are generated and stored in the list "inds_hist_storage"

# variables you will need to run the hybrid CRN are :
# "year" = year, from 1 to 20
# "id" = unique ID for each individual offspring
# "productivity" = number of offspring produced in the year (fecundity)
# "growth" = growth of the given individual during the year
# "climate" = environmental variable, predictor of the among-individual correlation

# The list "expected_cor" contain information on the expected correlation for
# each the dataset. Each dataframe contain the expected correlation and the
# environmental quality for each year

#############################################
#### Clean environment and load packages ####
#############################################

rm(list = ls())
library(mgcv)
library(car)
library(MASS)
library(psych)
library(data.table)
library(plyr)
set.seed(105)

###############################
#### Create Life Histories ####
###############################

s2_init <- c(0.6, 0.9)    # adult survival, 2 life histories (fast, slow)

lhs <- expand.grid(s2_init)

f <- c(1,0.5)             # fecundity, 2 life histories (fast, slow)

lhs <- cbind(lhs, f)


##############################
#### Run IBM  for each LH ####
##############################

########################################
#### Vital rates correlation matrix ####
########################################

context_dependent_cor <- function(inds, t, climate){
  
  env <- climate[t]
  
  # parameters for the context dependent correlation value LM  
  intercept_cor_value <- -0.3      
  slope_env_cor_value <- -0.4
  
  correlation_value <- exp(intercept_cor_value + slope_env_cor_value*env) / 
    (1+exp(intercept_cor_value + slope_env_cor_value*env))
  
  X <- cbind(inds$quality_mass, inds$quality_f)
  
  initial_sd_1 <- sqrt(var(inds$quality_mass))  # initial sd in adult survival
  initial_sd_2 <- sqrt(var(inds$quality_f))   # initial sd in fecundity
  
  # get the eigen values/vectors
  eigen_system <- eigen(cor_matrix_tradeoff)
  
  # helper function - build rescaling matrix from eigenvalues and environment (`a`)
  mk_rescale_matrix <- function(eigen_vals, correlation_value) {
    x <- 2 * sqrt(eigen_vals[1]) - 1
    y <- (1 - sqrt(eigen_vals[2])) / (0.5 * sqrt(eigen_vals[2]))
    diag(c(1 / (1 + x * correlation_value), 1 + y * correlation_value))
  }
  
  # 1. matrix formed of eigenvectors
  eigen_vecs <- eigen_system$vectors
  # 2. diagonal scale matrix formed from eigenvalues and environment
  rescale_mat <- mk_rescale_matrix(eigen_system$values, correlation_value) 
  # 3. net linear transformation - depends on three matrices...
  Q <- eigen_vecs %*% rescale_mat %*% t(eigen_vecs)
  # 4. compute the context dependent values
  X2 <- X %*% Q
  # 5. rescale to original variance in each trait
  if (sqrt(var(X2)[1,1]) !=0 ){
    X2[,1] <- X2[,1]/sqrt(var(X2)[1,1]) * initial_sd_1
  }
  if (sqrt(var(X2)[2,2]) !=0 ){
    X2[,2] <- X2[,2]/sqrt(var(X2)[2,2]) * initial_sd_2
  }
  
  return(X2)
}


##############################################
#### Simulate mean vr for each individual ####
##############################################


simulate_mean_vr <- function(inds, t, climate, X2, lh){
  
  env <- climate[t]
  
  # choose amount of context-dependent among-id variance in each trait
  gamma_survival_ad <- inds$realized_sa * 0.7    
  gamma_fecundity_ad <- inds$realized_f * 0.4
  gamma_mass_change <- inds$realized_mass * 0.4
  
  # choose amount of fixed heterogeneity in each trait
  alpha_fecundity <- inds$fixed_h_f * 0.5
  alpha_mass <- inds$fixed_h_m * 0.5
  
  # parameters for the survival adult glm   
  intercept_survival_ad <- as.numeric(logit(lh[1]))        
  slope_env_survival_ad <- 0
  
  # parameters for the fecundity adult glm
  intercept_fecundity_ad <- as.numeric(log(lh[2]))       
  slope_env_fecundity_ad <- 0.3
  
  # parameters for the mass change glm
  intercept_mass_change <- 0     
  slope_env_mass_change <- 0.3
  
  
  ## Survival adult GLM ##
  survival_ad <- function(env) {
    
    exp(intercept_survival_ad + slope_env_survival_ad*env + gamma_survival_ad)/
      (1+exp(intercept_survival_ad + slope_env_survival_ad*env + gamma_survival_ad))
    
  }
  
  mu_sa <- survival_ad(env = env)
  
  
  ## Fecundity adult GLM ##
  fecundity_ad <- function(env) {
    
    exp(intercept_fecundity_ad + slope_env_fecundity_ad*env + alpha_fecundity + gamma_fecundity_ad)
    
  }
  
  mu_fa <- fecundity_ad(env = env)
  
  
  ## Mass change GLM ##
  fun_mass_change <- function(env) {
    
    intercept_mass_change + slope_env_mass_change*env + alpha_mass + gamma_mass_change
    
  }
  
  mu_mass <- fun_mass_change(env = env)
  
  
  vr.mu=cbind(mu_sa, mu_fa, mu_mass)              # vital rates per individual
  return(vr.mu)
}




#####################################
#### Survival process and ageing ####
#####################################

death <- function(inds, phi_ad){
  
  inds$survival <- rbinom(dim(inds)[1], 1, phi_ad[])
  
  inds$age <- ifelse(inds$survival==1, inds$age+1, inds$age)
  
  return(inds)
}




##############################################
#### Reproduction process and mass change ####
##############################################

birth <- function(inds, lambda_ad, growth){
  
  inds$productivity <- rpois(dim(inds)[1], lambda_ad[]) 
  inds$growth <- rnorm(dim(inds)[1], growth[], 0.01)
  
  return(inds)
}


#############################
#### Immigration process ####
#############################


immigration <- function(inds, t, correlation_matrix, n_recruit){
  
  corr.values.new <- mvrnorm(n_recruit, mu=rep(0,nrow(correlation_matrix)), Sigma=correlation_matrix) 
  
  # ---- We now have the total number of new individuals; now add to inds
  new_inds <- data.frame(id=seq(length.out=n_recruit)+max(inds$id), 
                         year=rep(t, n_recruit), 
                         age=rep(1, n_recruit), 
                         survival=rep(1, n_recruit), 
                         productivity=rep(0, n_recruit),
                         growth=rep(0, n_recruit),
                         quality_sa= corr.values.new[,1], 
                         quality_f= corr.values.new[,2],
                         quality_mass= corr.values.new[,3],
                         fixed_h_f = rnorm(n_recruit,0,1),
                         fixed_h_m = rnorm(n_recruit,0,1),
                         realized_sa= 0,
                         realized_f= 0,
                         realized_mass= 0)
  
  # ---- Our new offspring can now be attached in the inds array
  inds <- rbind(inds, new_inds)
  
  return(inds)
}




#####################################
#### Individual based simulation ####
#####################################

#### Set initial conditions

replicat <- 1          # number of replicates
time_steps <- c(30)    # number of time step (e.g. years)
n_recruit <- c(25,50)  # number of new individuals every year


inds_hist_storage <- list()
climate_storage <- list()


for (lh in 1:dim(lhs)[1]) {
  
  for (study_length in 1:length(time_steps)) {
    
    for (sample_size in 1:length(n_recruit)) {
      
      for (i in 1:replicat) {
        
        set.seed(i+20)
        climate <- rnorm(time_steps[study_length], 0, 1)
        
        ts <- 1              # starting step
        
        
        ### create initial population
        inds <- as.data.frame(array(data = 0, dim = c(n_recruit[sample_size], 14)))   # initial population at time step 0
        colnames(inds) <- c("id", "year", "age", "survival", "productivity", "growth",
                            "quality_sa", "quality_f", "quality_mass", "fixed_h_f", "fixed_h_m",
                            "realized_sa", "realized_f", "realized_mass")
        inds$id <- 1:n_recruit[sample_size]
        
        rho <- -0.8    # most negative correlation is -0.6
        
        correlation_matrix <- matrix(c(1, 0, 0,
                                       0, 1, rho,
                                       0, rho, 1),
                                     ncol = 3, nrow = 3)
        
        #corr.values <- matrix(0, nrow = n_recruit[sample_size], ncol = 3)
        corr.values <- mvrnorm(n_recruit[sample_size], mu=rep(0,nrow(correlation_matrix)), Sigma=correlation_matrix) 
        
        inds$quality_sa <- corr.values[,1]  # survival adult individual heterogeneity
        inds$quality_f <- corr.values[,2]  # fecundity individual heterogeneity
        inds$quality_mass <- corr.values[,3]
        
        inds$fixed_h_f <- rnorm(n_recruit[sample_size],0,1)
        inds$fixed_h_m <- rnorm(n_recruit[sample_size],0,1)
        
        cor_matrix_tradeoff <- correlation_matrix[c(2,3), c(2,3)]
        
        inds$survival <- rep(1, n_recruit[sample_size])
        inds$age <- rep(1, n_recruit[sample_size])         
        
        
        inds_hist <- NULL
        
        
        while(ts <= time_steps[study_length]){
          
          inds <- subset(inds, survival==1)
          inds$year          <- ts                  # year
          X2                 <- context_dependent_cor(inds = inds, t = ts, climate = climate)
          inds$realized_sa   <- inds$quality_sa
          inds$realized_f    <- X2[,2]
          inds$realized_mass <- X2[,1]
          vr.mu              <- simulate_mean_vr(inds = inds, t = ts, climate = climate, X2 = X2, lh = lhs[lh,])
          inds               <- birth(inds, lambda_ad = vr.mu[,2], growth = vr.mu[,3])   # reproduction process for all 
          inds               <- death(inds, phi_ad = vr.mu[,1])     # survival process with aging
          
          
          ts                 <- ts + 1                # next time step
          
          
          inds_hist[[ts]]  <- inds # Add to individual history
          
          
          inds <- immigration(inds, t = ts, correlation_matrix = correlation_matrix, 
                              n_recruit = n_recruit[sample_size])    # add new adults to population
        } # population loop
        
        inds_hist <- rbindlist(inds_hist)
        
        # store datasets of individual histories
        inds_hist_storage[[replicat*length(n_recruit)*length(time_steps) * (lh-1) +
                             replicat*length(n_recruit) * (study_length-1) +
                             replicat * (sample_size-1) +
                             i]]  <- inds_hist
        
        # store climate datasets
        climate_storage[[replicat*length(n_recruit)*length(time_steps) * (lh-1) +
                           replicat*length(n_recruit) * (study_length-1) +
                           replicat * (sample_size-1) +
                           i]] <- data.frame(climate=climate, ye=c(1:time_steps[study_length]))
        
      } #replicate
      
    } #sample size
    
  } #study length
  
  print(lh)
} #lh

# list of final datasets containing individual histories and clinate data
inds_hist_storage <- lapply(seq_along(inds_hist_storage), function(x) merge(inds_hist_storage[[x]], climate_storage[[x]],
                                                                            by.x = "year", by.y = "ye"))
# get expected effect size
expected_cor <- lapply(inds_hist_storage,
                       function(x) ddply(x, "year",
                                         summarize,
                                         correlation_year = cor(realized_mass, realized_f))[,2])

expected_cor <- lapply(seq_along(expected_cor), function(x) cbind(expected_cor[[x]], climate_storage[[x]][1]))
expected_cor <- lapply(expected_cor, setNames, c("expected_correlation", "climate"))


# clean environment
rm(list=setdiff(ls(), c("inds_hist_storage", "expected_cor")))




