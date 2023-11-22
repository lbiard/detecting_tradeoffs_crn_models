# Last modified 20/11/23

#####################################################
#### clean working environment and load packages ####
#####################################################

rm(list = ls())

library(tidyr)
library(dplyr)
library(plyr)
library(tidyverse)
library(shinystan)
library(cmdstanr)
library(ggdist)
library(ggplot2)
library(patchwork)


####################################
#### load formatted marmot data ####
####################################

# set your own path
final_data <- read.delim("~/My Drive/phd/phd_simulation/Chapter 2/data_code_for_submission/marmot/marmot_data_for_analysis.txt")

#########################################
#### get data for fecundity submodel ####
#########################################

data_mom <- final_data[match(unique(final_data$litter_id), final_data$litter_id),]

######################################
#### scale variables for analysis ####
######################################

final_data$total_snow <- scale(final_data$total_snow)[,1]
final_data$summer_max_t <- scale(final_data$summer_max_tj)[,1]
final_data$age <- scale(final_data$age)[,1]
final_data$massjun <- scale(final_data$massjun)[,1]
final_data$massjun_mom <- scale(final_data$massjun_mom)[,1]

data_mom$total_snow <- scale(data_mom$total_snow)[,1]
data_mom$summer_max_t <- scale(data_mom$summer_max_tj)[,1]
data_mom$age <- scale(data_mom$age)[,1]
data_mom$massjun_mom <- scale(data_mom$massjun_mom)[,1]


#############################
#### data for stan model ####
#############################

#data of growth submodel
df <- final_data
#data of fecundity submodel
df_fecundity <- data_mom


stan.df =
  list(N = nrow(df),                      # number of observation offpring mass
       M = nrow(df_fecundity),            # number of observarion fecundity
       I = length(unique(df$dam)),        # number of mothers
       Y = length(unique(df$yrborn)),     # number of years
       D = 2,                             # numver of traits
       id_g = as.numeric(as.factor(df$dam)),                 # mom id index for offpsring mass model
       id_f = as.numeric(as.factor(df_fecundity$dam)),       # mom id index for fecundity model (matching id with id_g)
       years_g = as.numeric(as.factor(df$yrborn)),           # year index for offpsring mass model
       years_f = as.numeric(as.factor(df_fecundity$yrborn)), # year index for fecundity model
       climate_g = df$total_snow,                   # total overwinter snow offpsring model
       climate_f = df_fecundity$total_snow,         # total overwinter snow fecundity model
       climate_y = scale(unique(df[,c("yrborn","total_snow")])$total_snow)[,1],      # total overwinter snow correlation model
       climate_g2 = df$summer_max_t,                                                 # june average daily max temperature offpsring model
       climate_y2 = scale(unique(df[,c("yrborn","summer_max_t")])$summer_max_t)[,1], # june average daily max tempeeature correlation model
       growth = as.numeric(df$massjun),             # offspring mass (response variable)
       productivity = df_fecundity$nboffspring,     # fecundity (response variable)
       age_g = df$age,                              # age of mom offspring mass model
       age_f = df_fecundity$age,                    # age of mom fecundity model
       mass_mom_g = df$massjun_mom,                 # june mass of mom offspring mass model
       mass_mom_f = df_fecundity$massjun_mom        # june mass of mom fecundity model
  )

########################
#### Run stan model ####
########################

# set your own path, where the stan model is
setwd("~/My Drive/phd/phd_simulation/Chapter 2/data_code_for_submission/marmot")

# Compile model
mod <- cmdstan_model("marmot_model.stan"
                     , stanc_options = list("O1")
)

# Fit model (takes a few hours)
fit <- mod$sample(
  data = stan.df, 
  seed = 1567, 
  chains = 3, 
  parallel_chains = 3,
  iter_warmup = 1000,
  iter_sampling = 3000,
  adapt_delta = 0.975,
  refresh = 40 # print update every 20 iters
)

# transform format of the output file foir convenience (takes a few minutes at least)
stanfit <- rstan::read_stan_csv(fit$output_files())

# Save and load output
#saveRDS(stanfit, "model_marmot_output.RDS")
#stanfit <- read_rds("~/model_marmot_output.RDS")

# Open shiny app with diagnostic plot (visualize the chains to check for convergence)
launch_shinystan(stanfit)



######################
#### make figures ####
######################

fit <- stanfit


chain_1 <- fit@sim[["samples"]][[1]][,c(1:13)]
chain_2 <- fit@sim[["samples"]][[2]][,c(1:13)]
chain_3 <- fit@sim[["samples"]][[3]][,c(1:13)]

dat_plot <- rbind(chain_1, chain_2, chain_3)


############################
#### Figure winter snow ####
############################

# predictions across range of winter snow values
x2.sim <- seq(min(stan.df$climate_y), max(stan.df$climate_y), by = 0.02)

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- (exp(dat_plot$mu_r + dat_plot$beta_r * (x2.sim[i])) /
                     (1+exp(dat_plot$mu_r + dat_plot$beta_r * (x2.sim[i])))) * 2 - 1
}

# calculate quantiles of predictions
bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.045)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.945)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean, bayes.c.eff.lower, bayes.c.eff.upper, bayes.c.eff.lower.bis, bayes.c.eff.upper.bis)

# Plot 
p <- ggplot(plot.dat, aes(x = x2.sim, y = bayes.c.eff.mean)) +
  geom_line(color = "black", alpha = 0.8, size = 1.8)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "black", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "black", alpha = 0.1)+
  ylim(c(-1,1))+
  ggtitle("")+
  xlab("Amount of winter snow")+
  ylab("Observation-level correlation")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size=16),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13))
p

###################################
#### Figure summer temperature ####
###################################

# predictions across range of temperature values
x2.sim <- seq(min(stan.df$climate_y2), max(stan.df$climate_y2), by = 0.02)

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- (exp(dat_plot$mu_r + dat_plot$beta_rt * (x2.sim[i])) /
                     (1+exp(dat_plot$mu_r + dat_plot$beta_rt * (x2.sim[i])))) * 2 - 1
}

# calculate quantiles of predictions
bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.045)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.945)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean, bayes.c.eff.lower, bayes.c.eff.upper, bayes.c.eff.lower.bis, bayes.c.eff.upper.bis)

# plot
q <- ggplot(plot.dat, aes(x = x2.sim, y = bayes.c.eff.mean)) +
  geom_line(color = "black", alpha = 0.8, size = 1.8)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "black", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "black", alpha = 0.1)+
  ylim(c(-1,1))+
  ggtitle("")+
  xlab("Maximum daily summer temperature")+
  ylab("Observation-level correlation")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size=16),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13))
q


###########################
#### Figure posteriors ####
###########################

# change the "9000" to another value if you run chains for less or more iterations than I did (3000 iterations post burn-in * 3 = 9000)

df.posteriors <- data_frame(Submodel = c(rep("Offspring mass", 9000*4), rep("Fecundity", 9000*3), rep("Correlation", 9000*2))
                            , parameter = c(rep("Max Temperature", 9000), rep("Total Snow", 9000), rep("Age", 9000), rep("Mass", 9000)
                                            , rep("Total Snow", 9000), rep("Age", 9000), rep("Mass", 9000)
                                            , rep("Max Temperature", 9000), rep("Total Snow", 9000))
                            , Posterior = c(dat_plot$beta_gt, dat_plot$beta_g, dat_plot$beta_g2, dat_plot$beta_g3
                                            , dat_plot$beta_f, dat_plot$beta_f2, dat_plot$beta_f3
                                            , dat_plot$beta_rt, dat_plot$beta_r))



df.posteriors$Submodel <- factor(df.posteriors$Submodel, levels=c("Offspring mass", "Fecundity", "Correlation"))
df.posteriors$parameter <- as.ordered(df.posteriors$parameter)

submodels <- unique(df.posteriors$Submodel)
parameters <- unique(df.posteriors$parameter)

# get summaries of posteriors
df.summary <- expand.grid(Submodel = factor(x = submodels, levels = submodels, ordered = TRUE)
                          , parameter = factor(x = parameters, levels = parameters, ordered = TRUE)
                          , mean = as.numeric(NA)
                          , BCI50_upper = as.numeric(NA)
                          , BCI50_lower = as.numeric(NA)
                          , BCI89_upper = as.numeric(NA)
                          , BCI89_lower = as.numeric(NA)
                          , significant = as.logical(NA))

# remove empty columns (some variables are not predictors in all the submodels)
df.summary <- df.summary[-c(2,9,12),]

# get quantiles
for(i in 1:nrow(df.summary)){
  row <- which(df.posteriors$Submodel==as.character(df.summary$Submodel[i]) & df.posteriors$parameter==as.character(df.summary$parameter[i]))
  df.summary$mean[i] <- mean(df.posteriors$Posterior[row])
  df.summary$BCI50_upper[i] <- quantile(df.posteriors$Posterior[row], 0.75, na.rm = T)
  df.summary$BCI50_lower[i] <- quantile(df.posteriors$Posterior[row], 0.25, na.rm = T)
  df.summary$BCI89_upper[i] <- quantile(df.posteriors$Posterior[row], 0.945, na.rm = T)
  df.summary$BCI89_lower[i] <- quantile(df.posteriors$Posterior[row], 0.055, na.rm = T)
  df.summary$significant[i] <- ifelse(test = (df.summary$BCI89_lower[i]>0 & df.summary$BCI89_upper[i]>0) || (df.summary$BCI89_lower[i]<0 & df.summary$BCI89_upper[i]<0), yes = TRUE, no = FALSE)
}



# some plot settings
dodge.width <- 0.7
colT1 <- "cornflowerblue"
colT2 <- "orange"
colCov <- "seagreen4"
colPlot <- "black"

# plot
p2 <-ggplot()+
  stat_halfeye(data = df.posteriors, aes(x = Posterior, y = parameter, fill = Submodel), color = NA, alpha = 0.2, position = position_dodge(width = dodge.width), normalize="xy", scale=0.55)+
  geom_point(data = df.summary, aes(x = mean, y = parameter, color = Submodel), size = 2, position = position_dodge(width = dodge.width))+
  geom_linerange(data = df.summary, aes(xmin = BCI89_lower, xmax = BCI89_upper, y = parameter, color = Submodel), size=0.6, linetype="solid", position = position_dodge(width = dodge.width))+
  geom_linerange(data = df.summary, aes(xmin = BCI50_lower, xmax = BCI50_upper, y = parameter, color = Submodel), size=1.5, linetype="solid", position = position_dodge(width = dodge.width))+
  scale_color_manual(values = c(colT1, colT2, colCov))+
  scale_fill_manual(values = c(colT1, colT2, colCov))+
  scale_alpha_manual(values = c(0.3,1))+
  geom_vline(xintercept = 0, linetype = "dashed", color = colPlot, size = 0.4) +
  scale_x_continuous(breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5))+
  coord_cartesian(xlim=c(-1.8,1.8), clip = "off")+
  ylab("")+
  theme_minimal()+
  theme(axis.line.y = element_blank()
        , axis.ticks.y = element_blank()
        , panel.grid.major = element_blank() 
        , panel.grid.minor = element_blank()
  )

# pdf("marmot_plot.pdf", width = 9, height = 10)
wrap_elements(full= (p | q)) / p2
# dev.off()

