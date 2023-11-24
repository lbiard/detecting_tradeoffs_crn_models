# Last modified 24/11/2023

# This R script can be used to run the Hybrid CRN model on simulated data
# First you need to run the R script "simulate_data_hybrid_crn.R" to create the simulated datasets.
# The current script will run the model on the first dataset only.

###########################
#### Prepare workspace ####
###########################

library(shinystan)
library(cmdstanr)
library(patchwork)
library(ggdist)
library(tibble)
library(ggplot2)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

df1 = offspring_hist_storage[[1]] # select first dataset

########################################
#### Prepare dataset for Stan model ####
########################################

df <- df1

# make the fecundity dataset (i.e. nrow = number of breeding attempts)
df_fecundity <- df[!duplicated(df[,c(1,3)]),]

stan.df =
  list(N = nrow(df),
       M = nrow(df_fecundity),
       I = length(unique(df$momid)),
       Y = length(unique(df$year)),
       D = 2,
       id_g = as.numeric(as.factor(df$momid)),
       id_f = as.numeric(as.factor(df_fecundity$momid)),
       years_g = df$year, 
       years_f = df_fecundity$year, 
       climate_g = df$climate,
       climate_f = df_fecundity$climate,
       climate_y = unique(df$climate),
       growth = df$mass,
       productivity = df_fecundity$clutch_size
  )

####################################
#### Compile and run Stan model ####
####################################

mod <- cmdstan_model("hybrid_crn_model.stan"
                     , stanc_options = list("O1")
)


fit <- mod$sample(
  data = stan.df, 
  seed = 123, 
  chains = 3, 
  parallel_chains = 3,
  adapt_delta = 0.975,
  refresh = 20 # print update every 20 iters
)

stanfit <- rstan::read_stan_csv(fit$output_files())
#saveRDS(stanfit, "model_hybrid_crn_output.RDS")
#stanfit <- read_rds("model_hybrid_crn_output.RDS")

# Open shiny app with diagnostic plot (visualize the chains to check for convergence)
launch_shinystan(stanfit)


######################
#### make figures ####
######################

fit <- stanfit

chain_1 <- fit@sim[["samples"]][[1]][,c(1:7)]
chain_2 <- fit@sim[["samples"]][[2]][,c(1:7)]
chain_3 <- fit@sim[["samples"]][[3]][,c(1:7)]

dat_plot <- rbind(chain_1, chain_2, chain_3)

########################
#### Figure climate ####
########################

# predictions across range of climatic values
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
  ylim(-1,1)+
  ggtitle("")+
  xlab("Climate")+
  ylab("Observation-level correlation")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))
p <- p + geom_point(data=expected_cor[[1]], aes(x=climate, y=expected_correlation), size=2) 

p

###########################
#### Figure posteriors ####
###########################

# change the "3000" to another value if you run chains for less or more iterations than I did (1000 iterations post burn-in * 3 = 3000)

df.posteriors <- data_frame(Submodel = c(rep("Offspring mass", 3000*3), rep("Litter size", 3000*2))
                            , parameter = c(rep("Intercept", 3000), rep("Climate", 3000), rep("Within-litter variance", 3000)
                                            , rep("Intercept", 3000), rep("Climate", 3000))
                            , Posterior = c(dat_plot$mu_growth, dat_plot$beta_g, dat_plot$sigma
                                            , dat_plot$mu_productivity, dat_plot$beta_p))



df.posteriors$Submodel <- factor(df.posteriors$Submodel, levels=c("Offspring mass", "Litter size"))
df.posteriors$parameter <- ordered(df.posteriors$parameter, levels=c("Within-litter variance","Climate","Intercept"))

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
df.summary <- df.summary[-6,]

df.expected <- data_frame(Submodel = df.summary$Submodel
                          , parameter = df.summary$parameter
                          , mean = c(0, log(2), 0.3, 0.3, 0.4))

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
dodge.width <- 0.8
colT1 <- "cornflowerblue"
colT2 <- "orange"
colCov <- "seagreen4"
colPlot <- "black"


# plot
q <- ggplot()+
  stat_halfeye(data = df.posteriors, aes(x = Posterior, y = parameter, fill = Submodel), color = NA, alpha = 0.2, position = position_dodge(width = dodge.width), normalize="xy", scale=0.7)+
  geom_point(data = df.summary, aes(x = mean, y = parameter, color = Submodel), size = 2, position = position_dodge(width = dodge.width))+
  geom_linerange(data = df.summary, aes(xmin = BCI89_lower, xmax = BCI89_upper, y = parameter, color = Submodel), size=0.6, linetype="solid", position = position_dodge(width = dodge.width))+
  geom_linerange(data = df.summary, aes(xmin = BCI50_lower, xmax = BCI50_upper, y = parameter, color = Submodel), size=1.5, linetype="solid", position = position_dodge(width = dodge.width))+
  scale_color_manual(values = c(colT1, colT2, colCov))+
  scale_fill_manual(values = c(colT1, colT2, colCov))+
  scale_alpha_manual(values = c(1,1))+
  geom_segment(aes(x = log(2.723) , y = 3.1, xend = log(2.723), yend = 3.63), linetype=2) + #2.723 is the mean of a poisson with lambda 2.5 once zero truncation is applied
  geom_segment(aes(x = 0 , y = 2.68, xend = 0, yend = 3.25), linetype=2) +
  geom_segment(aes(x = 0.3 , y = 1.6, xend = 0.3, yend = 2.65), linetype=2) +
  geom_segment(aes(x = 0.4 , y = 0.8, xend = 0.4, yend = 1.5), linetype=2) +
  ylab("")+
  theme_minimal()+
  theme(axis.line.y = element_blank()
        , axis.ticks.y = element_blank()
        , panel.grid.major = element_blank() 
        , panel.grid.minor = element_blank()
  )
q

# pdf("hybrid_crn_plot.pdf", width = 9, height = 6)
wrap_elements(full= (p | q))
# dev.off()
