# Last modified 03/05/2024

# This R script can be used to simulate individual-based datasets with a 
# tradeoff between fecundity and offpsring quality (mass). The expression of this tradeoff is 
# dependent on the environmental quality at each time step (continuous variation
# of the among-individual correlation depends on the environment).

# The code simulate 4 datasets (only the first one is analysed in the manuscript):
# 1 - fast life history, 20 years of data, 25 new individuals per year
# 2 - fast life history, 20 years of data, 50 new individuals per year
# 3 - slow life history, 20 years of data, 25 new individuals per year
# 4 - slow life history, 20 years of data, 50 new individuals per year

# These datasets are generated and stored in the list "offspring_hist_storage"

# variables you will need to run the hybrid CRN are :
# "year" = year, from 1 to 20
# "id" = unique ID for each individual offspring
# "momid" = unique ID of the mother of the offspring
# "clutch_size" = number of offspring in the clutch/litter
# "mass" = mass of the individual offspring
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
set.seed(126)

###############################
#### Create Life Histories ####
###############################

s2_init <- c(0.6, 0.8)    # adult survival, 2 life histories (fast, slow)

lhs <- expand.grid(s2_init)

f <- c(2.5,0.7)             # fecundity, 2 life histories (fast, slow)

lhs <- cbind(lhs, f)

##############################
##############################
#### Run IBM  for each LH ####
##############################
##############################

################################################################################
#### Function to create vector with given correlation to an existing vector ####
################################################################################

#### function needed to set correlation with an existing variable 
set_cor_intergen <- function(new_inds, climate, t) {
  
  env <- climate[t] 
  
  
  correlation_value <- -((2*plogis(env) - 1) * 0.6)
  
  
  # make correlation matrix
  C <- matrix(correlation_value, nrow = 2, ncol = 2)
  diag(C) <- 1
  
  C <- chol(C)
  
  # create new vector (offspring trait)
  X2 <- rnorm(dim(new_inds)[1], 0, 1)
  # store parental and offspring trait 
  X <- cbind(new_inds$quality_f, X2)
  
  # induce correlation (does not change X[,1])
  df <- X %*% C
  
  return(df[,2])  # offspring trait, with specific correlation to parental trait
}


##########################################
#### Simulate mean vr for each parent ####
##########################################


simulate_mean_vr <- function(inds, t, climate, X2, lh){
  
  env <- climate[t]
  
  # choose amount of among-id variance in each trait
  gamma_survival_ad <- inds$quality_sa * 0.1    
  gamma_fecundity_ad <- inds$quality_f * 0.5
  
  # parameters for the survival adult glm   
  intercept_survival_ad <- as.numeric(logit(lh[1]))        
  slope_env_survival_ad <- 0
  
  # parameters for the fecundity adult glm
  intercept_fecundity_ad <- as.numeric(log(lh[2]))       
  slope_env_fecundity_ad <- 0.3
  
  
  ## Survival adult GLM ##
  survival_ad <- function(env) {
    
    exp(intercept_survival_ad + slope_env_survival_ad*env + gamma_survival_ad)/
      (1+exp(intercept_survival_ad + slope_env_survival_ad*env + gamma_survival_ad))
    
  }
  
  mu_sa <- survival_ad(env = env)
  
  
  ## Fecundity adult GLM ##
  fecundity_ad <- function(env) {
    
    exp(intercept_fecundity_ad + slope_env_fecundity_ad*env + gamma_fecundity_ad)
    
  }
  
  mu_fa <- fecundity_ad(env = env)
  
  
  vr.mu=cbind(mu_sa, mu_fa)              # vital rates per individual
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

birth <- function(inds, lambda_ad){
  
  inds$productivity <- rpois(dim(inds)[1], lambda_ad[]) 
  #inds$productivity <- rnegbin(dim(inds)[1], lambda_ad[], (lambda_ad[]/2)) 
  
  return(inds)
}


##########################
#### Create offspring ####
##########################

create_offspring <- function(inds, climate, t){
  
  total_off <- sum(inds$productivity) # total number of offspring during the time step
  
  # We have the total number of new offspring; now add to inds
  offspring_data <- data.frame(id=seq(length.out=total_off),      #ID
                               year=rep(t, total_off),                         #year
                               momid=rep(inds$id,inds$productivity),           #ID of parent
                               clutch_size=rep(inds$productivity,inds$productivity),
                               quality_f_parent = rep(inds$quality_f, inds$productivity),
                               mass=rep(0,total_off), 
                               quality_mass = NA)
  
  average_mass_cluch <- set_cor_intergen(inds, climate, t) 
  
  offspring_data$quality_mass <-  rep(average_mass_cluch, inds$productivity)
  
  #### GLM offspring mass ####
  env <- climate[t]
  intercept_mass <- 0  
  slope_env_mass <- 0.3
  clutch_quality <- offspring_data$quality_mass
  mass_offspring <- function(env) {
    intercept_mass + slope_env_mass*env + clutch_quality
  }
  offspring_data$mass <- rnorm(dim(offspring_data)[1], mass_offspring(env = env), 0.4)
  
  return(offspring_data)
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
                         quality_sa= corr.values.new[,1], 
                         quality_f= corr.values.new[,2])
  
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
offspring_hist_storage <- list()


for (lh in 1:dim(lhs)[1]) {
  
  for (study_length in 1:length(time_steps)) {
    
    for (sample_size in 1:length(n_recruit)) {
      
      for (i in 1:replicat) {
        
        set.seed(i+125)
        climate <- rnorm(time_steps[study_length], 0, 1)
        
        ts <- 1              # starting step
        
        
        ### create initial population
        inds <- as.data.frame(array(data = 0, dim = c(n_recruit[sample_size], 7)))   # initial population at time step 0
        colnames(inds) <- c("id", "year", "age", "survival", "productivity",
                            "quality_sa", "quality_f")
        inds$id <- 1:n_recruit[sample_size]
        
        correlation_matrix <- matrix(c(1, 0,
                                       0, 1),
                                     ncol = 2, nrow = 2)
        
        corr.values <- matrix(0, nrow = n_recruit[sample_size], ncol = 2)
        corr.values <- mvrnorm(n_recruit[sample_size], mu=rep(0,nrow(correlation_matrix)), Sigma=correlation_matrix) 
        
        inds$quality_sa <- corr.values[,1]  # survival adult individual heterogeneity
        inds$quality_f <- corr.values[,2]  # fecundity individual heterogeneity
        
        inds$survival <- rep(1, n_recruit[sample_size])
        inds$age <- rep(1, n_recruit[sample_size])         
        
        
        inds_hist <- NULL
        offspring_hist <- NULL
        
        while(ts <= time_steps[study_length]){
          
          inds <- subset(inds, survival==1)
          inds$year          <- ts                  # year
          vr.mu              <- simulate_mean_vr(inds = inds, t = ts, climate = climate, X2 = X2, lh = lhs[lh,])
          inds               <- birth(inds, lambda_ad = vr.mu[,2])   # reproduction process for all 
          inds               <- death(inds, phi_ad = vr.mu[,1])     # survival process with aging
          offspring_data    <- create_offspring(inds = inds, climate = climate, t = ts)
          
          
          ts                 <- ts + 1                # next time step
          
          
          inds_hist[[ts]]  <- inds # Add to individual history
          offspring_hist[[ts]] <- offspring_data
          
          inds <- immigration(inds, t = ts, correlation_matrix = correlation_matrix, 
                              n_recruit = n_recruit[sample_size])    # add new adults to population
        } # population loop
        
        inds_hist <- rbindlist(inds_hist)
        offspring_hist <- rbindlist(offspring_hist)
        
        # store datasets of individual histories
        inds_hist_storage[[replicat*length(n_recruit)*length(time_steps) * (lh-1) +
                             replicat*length(n_recruit) * (study_length-1) +
                             replicat * (sample_size-1) +
                             i]]  <- inds_hist
        offspring_hist_storage[[replicat*length(n_recruit)*length(time_steps) * (lh-1) +
                                  replicat*length(n_recruit) * (study_length-1) +
                                  replicat * (sample_size-1) +
                                  i]]  <- offspring_hist
        
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
offspring_hist_storage <- lapply(seq_along(offspring_hist_storage), function(x) merge(offspring_hist_storage[[x]], climate_storage[[x]],
                                                                                      by.x = "year", by.y = "ye"))
# get expected effect size
expected_cor <- lapply(
                       lapply(offspring_hist_storage, 
                       function(x) x[!duplicated(paste(x$year, x$momid, sep = ".")),]),
                function(x) ddply(x, "year",
                    summarize,
                    correlation_year = cor(quality_f_parent, quality_mass))[,2])

expected_cor <- lapply(seq_along(expected_cor), function(x) cbind(expected_cor[[x]], climate_storage[[x]][1]))
expected_cor <- lapply(expected_cor, setNames, c("expected_correlation", "climate"))


# clean environment
rm(list=setdiff(ls(), c("inds_hist_storage", "offspring_hist_storage", "expected_cor")))

