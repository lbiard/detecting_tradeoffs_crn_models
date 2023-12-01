# Last modified 22/11/23

#####################################################
#### clean working environment and load packages ####
#####################################################

rm(list = ls())

library(tidyverse)
library(shinystan)
library(cmdstanr)
library(ggplot2)
library(ggdist)
library(patchwork)

###################################
#### load formatted sheep data ####
###################################

# set your own path
df <- read.delim("~/My Drive/phd/phd_simulation/Chapter 2/data_code_for_submission/sheep/sheep_data_for_analysis.txt")

################################
#### dataset for Stan model ####
################################

stan.df =
  list(N = nrow(df), #number of obs
       I = length(unique(df$id)), #number of individuals
       Y = length(unique(df$obsY)), #number of years
       D = 2, #number of traits
       id = as.numeric(as.factor(df$id)), #vector of IDs
       years = as.numeric(as.factor(df$obsY)), #vector of years
       age = scale(df$ageY)[,1], #age vector
       density = scale(df$abundance)[,1], #population density vector
       climate = scale(df$nao_winter)[,1], #nao vector
       weight = scale(log(df$capWgt))[,1], #weight at time t vector
       climate_y = scale(unique(df[,c("obsY", "nao_winter")])$nao_winter)[,1], #nao vector of year length
       density_y = scale(unique(df[,c("obsY", "abundance")])$abundance)[,1], #density vector of year length
       growth = scale(log(df$growth))[,1], #weight at time t+1 
       productivity = df$lambNum1 #fecundity
  )


########################
#### run Stan model ####
########################

# set your own path, where the stan model is
setwd("~/My Drive/phd/phd_simulation/Chapter 2/data_code_for_submission/sheep")

# Compile model
mod <- cmdstan_model("sheep_model.stan"
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
  adapt_delta = 0.985,
  refresh = 20 # print update every 20 iters
)

# summary for variables of interest
fit$summary(variables = c("mu_growth", "mu_productivity", "mu_r", 
                          "beta_g_nao", "beta_g_age", "beta_g_den", "beta_g_wei", 
                          "beta_p_den", "beta_p_age", "beta_p_wei", 
                          "beta_r_nao", "beta_r_den"))

stanfit <- rstan::read_stan_csv(fit$output_files())

# Save and load output
#saveRDS(stanfit, "model_sheep_output.RDS")
#stanfit <- read_rds("~/model_marmot_output.RDS.RDS")

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


####################
#### Figure NAO ####
####################

# predictions across range of NAO values
x2.sim <- seq(min(stan.df$climate_y), max(stan.df$climate_y), by = 0.02)

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- (exp(dat_plot$mu_r + dat_plot$beta_r_nao * (x2.sim[i])) /
                     (1+exp(dat_plot$mu_r + dat_plot$beta_r_nao * (x2.sim[i])))) * 2 - 1
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
  xlab("Winter NAO")+
  ylab("Observation-level correlation")+
  theme_bw() +
  coord_cartesian(clip = "off") + 
  annotation_custom(grid::textGrob(expression(~ italic("(-) Environmental harshness (+)")),
                                   gp=grid::gpar(fontsize=10)),
                    xmin = mean(x2.sim), xmax = mean(x2.sim), ymin = -1.48, ymax = -1.48) +
  theme(plot.margin = margin(5.5,5.5,14,5.5))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))
p

########################
#### Figure density ####
########################

# predictions across range of density values
x2.sim <- seq(min(stan.df$density_y), max(stan.df$density_y), by = 0.02)

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- (exp(dat_plot$mu_r + dat_plot$beta_r_den * (x2.sim[i])) /
                     (1+exp(dat_plot$mu_r + dat_plot$beta_r_den * (x2.sim[i])))) * 2 - 1
}

# calculate quantiles of predictions
bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.045)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.945)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean, bayes.c.eff.lower, bayes.c.eff.upper, bayes.c.eff.lower.bis, bayes.c.eff.upper.bis)

# Plot
q <- ggplot(plot.dat, aes(x = x2.sim, y = bayes.c.eff.mean)) +
  geom_line(color = "black", alpha = 0.8, size = 1.8)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "black", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "black", alpha = 0.1)+
  ylim(c(-1,1))+
  ggtitle("")+
  xlab("Population density")+
  ylab("Observation-level correlation")+
  theme_bw() +
  coord_cartesian(clip = "off") + 
  annotation_custom(grid::textGrob(expression(~ italic("(-) Environmental harshness (+)")),
                                   gp=grid::gpar(fontsize=10)),
                    xmin = mean(x2.sim), xmax = mean(x2.sim), ymin = -1.48, ymax = -1.48) +
  theme(plot.margin = margin(5.5,5.5,14,5.5))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))
q


###########################
#### Figure posteriors ####
###########################

# change the "9000" to another value if you run chains for less or more iterations than I did (3000 iterations post burn-in * 3 = 9000)

df.posteriors <- data_frame(Submodel = c(rep("Mass t+1", 9000*4), rep("Fecundity", 9000*3), rep("Correlation", 9000*2))
                            , parameter = c(rep("Winter NAO", 9000), rep("Density", 9000), rep("Age", 9000), rep("Mass", 9000)
                                            , rep("Density", 9000), rep("Age", 9000), rep("Mass", 9000)
                                            , rep("Winter NAO", 9000), rep("Density", 9000))
                            , Posterior = c(dat_plot$beta_g_nao, dat_plot$beta_g_den, dat_plot$beta_g_age, dat_plot$beta_g_wei
                                            , dat_plot$beta_p_den, dat_plot$beta_p_age, dat_plot$beta_p_wei
                                            , dat_plot$beta_r_nao, dat_plot$beta_r_den))


df.posteriors$Submodel <- factor(df.posteriors$Submodel, levels=c("Mass t+1", "Fecundity", "Correlation"))
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


# some plot settings:
dodge.width <- 0.7
colT1 <- "cornflowerblue"
colT2 <- "orange"
colCov <- "seagreen4"
colPlot <- "black"


# Plot
p2 <- ggplot()+
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
  theme(plot.margin = margin(20,5.5,5.5,5.5))+
  theme(axis.line.y = element_blank()
        , axis.ticks.y = element_blank()
        , panel.grid.major = element_blank() 
        , panel.grid.minor = element_blank()
  )

# pdf("sheep_plot.pdf", width = 9, height = 10)
wrap_elements(full= (p | q)) / p2
# dev.off()
