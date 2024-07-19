# Last modified 19/07/2024

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
library(bayesplot)


####################################
#### load formatted marmot data ####
####################################

# set your own path
final_data <- read.delim("~/My Drive/phd/phd_simulation/Chapter 2/data_code_for_submission/marmot/marmot_data_for_analysis.txt")

#########################################
#### get data for fecundity submodel ####
#########################################

data_mom <- final_data[match(unique(final_data$litter_id), final_data$litter_id),]

final_data$summer_max_t <- final_data$summer_max_tj # use June max daily temperature
data_mom$summer_max_t <- data_mom$summer_max_tj

#############################
#### data for stan model ####
#############################

#data of growth submodel
df <- final_data
#data of fecundity submodel
df_fecundity <- data_mom

# matrix of predictors for correlation submodel
X <- as.matrix(cbind(rep(1, length(unique(df$yrborn))), 
                scale(unique(df[,c("yrborn","total_snow")])$total_snow)[,1],
                scale(unique(df[,c("yrborn","summer_max_t")])$summer_max_t)[,1],
                scale(unique(df[,c("yrborn","total_snow")])$total_snow)[,1] * scale(unique(df[,c("yrborn","summer_max_t")])$summer_max_t)[,1]))

# matrix of predictors for growth submodel
X_g = as.matrix(cbind(rep(1, dim(df)[1]), 
                      scale(df$total_snow)[,1],
                      scale(df$summer_max_t)[,1],
                      scale(df$age)[,1],
                      (scale(df$age)[,1])^2,
                      scale(df$massjun_mom)[,1],
                      scale(df$total_snow)[,1] * scale(df$summer_max_t)[,1]))

# matrix of predictors for fecundity submodel
X_f = as.matrix(cbind(rep(1, dim(df_fecundity)[1]), 
                      scale(df_fecundity$total_snow)[,1],
                      scale(df_fecundity$age)[,1],
                      (scale(df_fecundity$age)[,1])^2,
                      scale(df_fecundity$massjun_mom)[,1]))

cn <- as.numeric(table(as.numeric(as.factor(df_fecundity$yrborn))))

#create empty cmat - matrix of context-specific IDs
cmat <- matrix(NA, 
               nrow = length(unique(df$yrborn)), 
               ncol =  max(as.numeric(table(as.numeric(as.factor(df_fecundity$yrborn))))))

#fill cmat
temporary <- as.data.frame(cbind(as.numeric(as.factor(df_fecundity$yrborn)),
                                 as.numeric(as.factor(df_fecundity$dam))))
for (i in 1:length(unique(df$yrborn))) {
  cmat[i, 1:cn[i]] <- temporary$V2[temporary$V1 == i]
}
cmat_n = apply(cmat, 1, FUN = function(x) sum(!is.na(x)) )
cmat[is.na(cmat)] = 0 #remove NAs

temp = t(cmat)
corder = data.frame(id = temp[temp>0], c = rep(seq(1:nrow(cmat)), times = cmat_n))
idc_f = match(paste0(as.numeric(as.factor(df_fecundity$dam)),
                   as.numeric(as.factor(df_fecundity$yrborn)),
                   sep="."), 
            paste0(corder$id,corder$c,sep="."))

idc_g = match(paste0(as.numeric(as.factor(df$dam)),
                   as.numeric(as.factor(df$yrborn)),
                   sep="."), 
            paste0(corder$id,corder$c,sep="."))

rownames(df_fecundity) <- NULL


stan.df =
  list(N = nrow(df),                      # number of observation offspring mass
       M = nrow(df_fecundity),            # number of observation fecundity
       C = length(unique(df$yrborn)),     # number of years
       I = length(unique(df$dam)),        # number of mothers
       D = 2,                             # number of traits
       P_y = 4,                           # number of predictor on correlation
       P_g = 7,                           # number of predictor on growth
       P_f = 5,                           # number of predictor on fecundity
       
       id_g = as.numeric(as.factor(df$dam)),           # index of mother id
       c_id_g = as.numeric(as.factor(df$yrborn)),      # index of year id
       idc_g = idc_g,                                  # index of reproduction event id
       id_g_lm = as.numeric(rownames(df)),             # index of row names
       
       id_f = as.numeric(as.factor(df_fecundity$dam)),       # index of mother id
       c_id_f = as.numeric(as.factor(df_fecundity$yrborn)),  # index of year id
       idc_f = idc_f,                                        # index of reproduction event id
       id_f_lm = as.numeric(rownames(df_fecundity)),         # index of row names
       
       X = X,                               # matrix of predictors for correlation submodel
       X_g = X_g,                           # matrix of predictors for growth submodel
       X_f = X_f,                           # matrix of predictors for fecundity submodel
       A = diag(length(unique(df$dam))),    # identity matrix 
       
       cm = max(as.numeric(table(as.numeric(as.factor(df_fecundity$yrborn))))),
       cmat = cmat,
       cn = cn,
       cnt = length(unique(paste(df$dam, df$yrborn))), # number of reproduction events
       
       growth = scale(as.numeric(df$massjun))[,1], # offspring mass (response variable)
       productivity = df_fecundity$nboffspring     # fecundity (response variable)
  )
       
       

########################
#### Run stan model ####
########################

# set your own path, where the stan model file is
setwd("~/My Drive/phd/phd_simulation/Chapter 2/models_revision/marmot")

# Compile model
mod <- cmdstan_model("updated_model_marmot.stan"
                     , stanc_options = list("O1")
)

# Fit model 
fit <- mod$sample(
  data = stan.df, 
  seed = 1567, 
  chains = 3, 
  parallel_chains = 3,
  iter_warmup = 1000,
  iter_sampling = 3000,
  adapt_delta = 0.975,
  max_treedepth = 10,
  refresh = 40 # print update every XX iterations
)

#####################################
#### Posterior Predictive Checks ####
#####################################

# Extract y_rep for posterior predictive check
# This is to check if fecundity (litter size) predicted by the model concur with the real data
draws_f <- fit$draws(
  variables = "y_rep_f",
  inc_warmup = FALSE,
  format = "matrix"
)

#Plot for posterior predictive checks
plot_ppc_f <- ppc_dens_overlay(
  stan.df$productivity,
  draws_f[1:1000,],
  size = 0.25,
  alpha = 0.7,
  trim = FALSE
)
plot_ppc_f <- plot_ppc_f + 
  xlab("Litter size") + 
  ylab("Density") + 
  ggtitle("Posterior predictive check \u2014 Marmot litter size")
plot_ppc_f

# This is to check if offspring mass predicted by the model concur with the real data
draws_g <- fit$draws(
  variables = "y_rep_g",
  inc_warmup = FALSE,
  format = "matrix"
)

#Plot for posterior predictive checks
plot_ppc_g <- ppc_dens_overlay(
  (stan.df$growth * sd(df$massjun)) + mean(df$massjun),
  (draws_g[1:1000,]  * sd(df$massjun)) + mean(df$massjun),
  size = 0.25,
  alpha = 0.7,
  trim = FALSE
)
plot_ppc_g <- plot_ppc_g + 
  xlab("Offspring mass") + 
  ylab("Density") + 
  ggtitle("Posterior predictive check \u2014 Marmot offpsring mass")
plot_ppc_g

wrap_elements(full= (plot_ppc_f | plot_ppc_g))

# remove large elements from memory
rm(draws_f, draws_g)

#################################################################
#### Inspect chains and extract posterior draws for plotting ####
#################################################################

# transform format of the output file for convenience
CRN_fit <- rstan::read_stan_csv(fit$output_files())

# Save and load output
#saveRDS(CRN_fit, "model_marmot_output.RDS")
#CRN_fit <- read_rds("model_marmot_output.RDS")

# Open shiny app with diagnostic plot (visualize the chains to check for convergence)
launch_shinystan(CRN_fit)

post <- rstan::extract(CRN_fit)

traits = c("growth", "productivity")

# For the following plots, need to first run R script "functions_postprocessing.R", that contains some functions needed to make predictions using the model output

##########################
#### Plot winter snow ####
##########################

# sequence from scaled minimum to maximum value in x 
seq = seq(min(scale(unique(df[,c("yrborn","total_snow")])$total_snow)[,1]),
          max(scale(unique(df[,c("yrborn","total_snow")])$total_snow)[,1]),
          by =  (max(scale(unique(df[,c("yrborn","total_snow")])$total_snow)[,1]) - min(scale(unique(df[,c("yrborn","total_snow")])$total_snow)[,1])) / 25 ) 

# high temperature (set to +1sd)
X_pred = data.frame(int = 1, X1 = seq, X2 = 1, X3=1*seq)
mv_cpc = cpc_f(x = X_pred, cpc_b = post$B_cpc)
mv_cor = cor_f(x = X_pred, traits = traits, cpc = mv_cpc)

# low temperature (set to -1sd)
X_pred2 = data.frame(int = 1, X1 = seq, X2 = -1, X3=-1*seq)
mv_cpc2 = cpc_f(x = X_pred2, cpc_b = post$B_cpc)
mv_cor2 = cor_f(x = X_pred2, traits = traits, cpc = mv_cpc2)

# function to unscale predictions (because predictors were standardized)
inv_fun_p1 <- function(x){(x*sd(unique(df[,c("yrborn","total_snow")])$total_snow))+mean(unique(df[,c("yrborn","total_snow")])$total_snow)}

plot_1 <- ggplot(mv_cor,
       aes(x = inv_fun_p1(X1), y = value))+
  stat_lineribbon(linewidth=1.8, .width = 0.89, alpha=0.15, fill="#009E73")+
  stat_lineribbon(linewidth=1.8, .width = 0.5, alpha=0.30, fill="#009E73")+
  stat_lineribbon(linewidth=1.8, .width = 0.001, color="#009E73")+
  ylim(c(-1,1))+
  ggtitle("")+
  xlab("Amount of winter snow")+
  ylab("Observation-level correlation")+
  theme_bw() +
  coord_cartesian(clip = "off") + 
  annotation_custom(grid::textGrob(expression(~ italic("(+) Environmental harshness (-)")),
                                   gp=grid::gpar(fontsize=10)),
                    xmin = inv_fun_p1(mean(range(scale(unique(df[,c("yrborn","total_snow")])$total_snow)[,1]))),
                    xmax = inv_fun_p1(mean(range(scale(unique(df[,c("yrborn","total_snow")])$total_snow)[,1]))),
                    ymin = -1.48, 
                    ymax = -1.48) +
  theme(plot.margin = margin(5.5,5.5,14,5.5))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size=16),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13))+
  theme(legend.position="none")

plot_1 <- plot_1 + stat_lineribbon(data=mv_cor2, aes(x = inv_fun_p1(X1), y = value), linewidth=1.8, .width = 0.89, alpha=0.15, fill="#E69F00")+
  stat_lineribbon(data=mv_cor2, aes(x = inv_fun_p1(X1), y = value), linewidth=1.8, .width = 0.5, alpha=0.30, fill="#E69F00")+
  stat_lineribbon(data=mv_cor2, aes(x = inv_fun_p1(X1), y = value), linewidth=1.8, .width = 0.001, color="#E69F00")

plot_1


# annotation
mesure_pol <- data.frame(x1 = c(-2.7,0.3),
                         x2 = c(-2.3,0.7),
                         cat = c(1:2),
                         catnames = c("High temperature","Low temperature"))

ann <- ggplot(data = mesure_pol) +
  geom_rect(aes(xmin = x1,
                xmax = x2,
                ymin = 0.095,
                ymax = 0.105,
                fill = as.factor(cat)),
            fill = c("#009E73", "#E69F00"),
            color = "black",
            size = 0.3) +
  geom_text(aes(x = x2+1.15, y = 0.1, label = catnames),
            #vjust = .8, 
            fontface = "bold", color = "black") +
  coord_cartesian(xlim = c(-2.7085, 2.861)
                  , ylim = c(0.065, 0.14)) +
  theme_void()
ann

snow_plot <- ann / plot_1  + 
  plot_layout(heights = c(1, 8))
snow_plot

###############################
#### Plot June temperature ####
###############################

# sequence from scaled minimum to maximum value in x , high snow cover (set to +1sd)
seq = seq(min(scale(unique(df[,c("yrborn","summer_max_t")])$summer_max_t)[,1]),
          max(scale(unique(df[,c("yrborn","summer_max_t")])$summer_max_t)[,1]),
          by =  (max(scale(unique(df[,c("yrborn","summer_max_t")])$summer_max_t)[,1]) - min(scale(unique(df[,c("yrborn","summer_max_t")])$summer_max_t)[,1])) / 25 ) 

# high snow cover (set to +1sd)
X_pred = data.frame(int = 1, X1 = 1, X2 = seq, X3=1*seq)
mv_cpc = cpc_f(x = X_pred, cpc_b = post$B_cpc)
mv_cor = cor_f(x = X_pred, traits = traits, cpc = mv_cpc)

# low snow cover (set to -1sd)
X_pred2 = data.frame(int = 1, X1 = -1, X2 = seq, X3=-1*seq)
mv_cpc2 = cpc_f(x = X_pred2, cpc_b = post$B_cpc)
mv_cor2 = cor_f(x = X_pred2, traits = traits, cpc = mv_cpc2)

# function to unscale predictions (because predictors were standardized)
inv_fun_p2 <- function(x){(x*sd(unique(df[,c("yrborn","summer_max_t")])$summer_max_t))+mean(unique(df[,c("yrborn","summer_max_t")])$summer_max_t)}

plot_2 <- ggplot(mv_cor,
                 aes(x = inv_fun_p2(X2), y = value))+
  stat_lineribbon(linewidth=1.8, .width = 0.89, alpha=0.15, fill="#009E73")+
  stat_lineribbon(linewidth=1.8, .width = 0.5, alpha=0.30, fill="#009E73")+
  stat_lineribbon(linewidth=1.8, .width = 0.001, color="#009E73")+
  ylim(c(-1,1))+
  ggtitle("")+
  xlab("Maximum daily summer temperature")+
  ylab("Observation-level correlation")+
  theme_bw() +
  coord_cartesian(clip = "off") + 
  annotation_custom(grid::textGrob(expression(~ italic("(-) Environmental harshness (+)")),
                                   gp=grid::gpar(fontsize=10)),
                    xmin = inv_fun_p2(mean(range(scale(unique(df[,c("yrborn","summer_max_t")])$summer_max_t)[,1]))),
                    xmax = inv_fun_p2(mean(range(scale(unique(df[,c("yrborn","summer_max_t")])$summer_max_t)[,1]))),
                    ymin = -1.48,
                    ymax = -1.48) +
  theme(plot.margin = margin(5.5,5.5,14,5.5))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size=16),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13))+
  theme(legend.position="none")

plot_2 <- plot_2 + stat_lineribbon(data=mv_cor2, aes(x = inv_fun_p2(X2), y = value), linewidth=1.8, .width = 0.89, alpha=0.15, fill="#E69F00")+
  stat_lineribbon(data=mv_cor2, aes(x = inv_fun_p2(X2), y = value), linewidth=1.8, .width = 0.5, alpha=0.30, fill="#E69F00")+
  stat_lineribbon(data=mv_cor2, aes(x = inv_fun_p2(X2), y = value), linewidth=1.8, .width = 0.001, color="#E69F00")

plot_2


#annotation
mesure_pol2 <- data.frame(x1 = c(-2.7,0.3),
                          x2 = c(-2.3,0.7),
                          cat = c(1:2),
                          catnames = c("High snow cover","Low snow cover"))

ann2 <- ggplot(data = mesure_pol2) +
  geom_rect(aes(xmin = x1,
                xmax = x2,
                ymin = 0.095,
                ymax = 0.105,
                fill = as.factor(cat)),
            fill = c("#009E73", "#E69F00"),
            color = "black",
            size = 0.3) +
  geom_text(aes(x = x2+1.15, y = 0.1, label = catnames),
            #vjust = .8, 
            fontface = "bold", color = "black") +
  coord_cartesian(xlim = c(-2.7085, 2.861)
                  , ylim = c(0.065, 0.14)) +
  theme_void()
ann2

temperature_plot <- ann2 / plot_2  + 
  plot_layout(heights = c(1, 8))
temperature_plot

wrap_elements(full= (snow_plot | temperature_plot))

###########################
#### Figure posteriors ####
###########################

#make sure you know the order of the variables in the predictor matrices
g_effects <- as.data.frame(post$B_m_g)
f_effects <- as.data.frame(post$B_m_f)
cor_effects <- as.data.frame(post$B_cpc)

# change the "9000" to another value if you run chains for less or more iterations than I did (3000 iterations post burn-in * 3 = 9000)

df.posteriors <- data_frame(Submodel = c(rep("Offspring mass", 9000*6), rep("Fecundity", 9000*4), rep("Correlation", 9000*3))
                            , parameter = c(rep("Max Temperature", 9000), rep("Total Snow", 9000), rep("Age", 9000), rep("Age^2", 9000), rep("Mass", 9000), rep("Max Temperature * Total Snow", 9000)
                                            , rep("Total Snow", 9000), rep("Age", 9000), rep("Age^2", 9000), rep("Mass", 9000)
                                            , rep("Max Temperature", 9000), rep("Total Snow", 9000), rep("Max Temperature * Total Snow", 9000))
                            , Posterior = c(g_effects$V3, g_effects$V2, g_effects$V4, g_effects$V5, g_effects$V6, g_effects$V7
                                            , f_effects$V2, f_effects$V3, f_effects$V4, f_effects$V5
                                            , cor_effects[,3], cor_effects[,2], cor_effects[,4]))



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
df.summary <- df.summary[-c(2,9,12,15,17),]

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



# some plot settings, choose some pretty colors
dodge.width <- 0.7
colT1 <- "cornflowerblue"
colT2 <- "tomato3"
colCov <- "olivedrab3"
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
  theme(plot.margin = margin(20,5.5,5.5,5.5))+
  theme(axis.line.y = element_blank()
        , axis.ticks.y = element_blank()
        , panel.grid.major = element_blank() 
        , panel.grid.minor = element_blank()
  )


#ragg::agg_tiff("Figure 4.tiff", width = 9, height = 10, units = "in", res = 300)
wrap_elements(full= (snow_plot | temperature_plot)) / p2
#dev.off()


#################################################
#### Plot individual relationships on traits ####
#################################################

#### Offspring mass ####

## age ##

x2.sim <- seq(min(scale(df$age)[,1]), max(scale(df$age)[,1]), by = 0.01)

int.sim <- matrix(rep(NA, nrow(g_effects)*length(x2.sim)), nrow = nrow(g_effects))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- ((g_effects$V1 + g_effects$V4 * (x2.sim[i]) + g_effects$V5 * ((x2.sim[i])^2)) * sd(df$massjun)) + mean(df$massjun)
}

bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.055)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.945)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean, bayes.c.eff.lower, bayes.c.eff.upper, bayes.c.eff.lower.bis, bayes.c.eff.upper.bis)

inv_fun1 <- function(x){(x*sd(df$age))+mean(df$age)}

plot_age_g <- ggplot(plot.dat, aes(x = inv_fun1(x2.sim), y = bayes.c.eff.mean)) +
  geom_line(color = "black", alpha = 0.8, size = 1.8)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "black", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "black", alpha = 0.1)+
  ggtitle("")+
  xlab("Age")+
  ylab("Offspring mass")+
  ylim(c(150,800))+
  theme_bw() +
  coord_cartesian(clip = "off") + 
  theme(plot.margin = margin(5.5,5.5,14,5.5))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size=16),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13))
plot_age_g <- plot_age_g + geom_point(data=df, aes(x = age, y=massjun), alpha=0.1, color = "black",
                                      position = position_jitter(w=0.15, h = 0), inherit.aes = F)
plot_age_g


## mass mom ##

x2.sim <- seq(min(scale(df$massjun_mom)[,1]), max(scale(df$massjun_mom)[,1]), by = 0.01)

int.sim <- matrix(rep(NA, nrow(g_effects)*length(x2.sim)), nrow = nrow(g_effects))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- ((g_effects$V1 + g_effects$V6 * (x2.sim[i])) * sd(df$massjun)) + mean(df$massjun)
}

bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.055)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.945)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean, bayes.c.eff.lower, bayes.c.eff.upper, bayes.c.eff.lower.bis, bayes.c.eff.upper.bis)

inv_fun2 <- function(x){(x*sd(df$massjun_mom))+mean(df$massjun_mom)}

plot_mass_g <- ggplot(plot.dat, aes(x = inv_fun2(x2.sim), y = bayes.c.eff.mean)) +
  geom_line(color = "black", alpha = 0.8, size = 1.8)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "black", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "black", alpha = 0.1)+
  ggtitle("")+
  xlab("Mass")+
  ylab("Offspring mass")+
  ylim(c(150,800))+
  theme_bw() +
  coord_cartesian(clip = "off") + 
  theme(plot.margin = margin(5.5,5.5,14,5.5))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size=16),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13))
plot_mass_g <- plot_mass_g + geom_point(data=df, aes(x = massjun_mom, y=massjun), alpha=0.1, color = "black",
                                      position = position_jitter(w=0, h = 0), inherit.aes = F)
plot_mass_g


## snow ##

x2.sim <- seq(min(scale(df$total_snow)[,1]), max(scale(df$total_snow)[,1]), by = 0.01)

int.sim <- matrix(rep(NA, nrow(g_effects)*length(x2.sim)), nrow = nrow(g_effects))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- ((g_effects$V1 + g_effects$V2 * (x2.sim[i])) * sd(df$massjun)) + mean(df$massjun)
}

bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.055)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.945)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean, bayes.c.eff.lower, bayes.c.eff.upper, bayes.c.eff.lower.bis, bayes.c.eff.upper.bis)

inv_fun3 <- function(x){(x*sd(df$total_snow))+mean(df$total_snow)}

plot_snow_g <- ggplot(plot.dat, aes(x = inv_fun3(x2.sim), y = bayes.c.eff.mean)) +
  geom_line(color = "black", alpha = 0.8, size = 1.8)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "black", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "black", alpha = 0.1)+
  ggtitle("")+
  xlab("Amount of winter snow")+
  ylab("Offspring mass")+
  ylim(c(150,800))+
  theme_bw() +
  coord_cartesian(clip = "off") + 
  theme(plot.margin = margin(5.5,5.5,14,5.5))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size=16),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13))
plot_snow_g <- plot_snow_g + geom_point(data=df, aes(x = total_snow, y=massjun), alpha=0.1, color = "black",
                                      position = position_jitter(w=5, h = 0), inherit.aes = F)
plot_snow_g


## temperature ##

x2.sim <- seq(min(scale(df$summer_max_t)[,1]), max(scale(df$summer_max_t)[,1]), by = 0.01)

int.sim <- matrix(rep(NA, nrow(g_effects)*length(x2.sim)), nrow = nrow(g_effects))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- ((g_effects$V1 + g_effects$V3 * (x2.sim[i])) * sd(df$massjun)) + mean(df$massjun)
}

bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.055)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.945)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean, bayes.c.eff.lower, bayes.c.eff.upper, bayes.c.eff.lower.bis, bayes.c.eff.upper.bis)

inv_fun4 <- function(x){(x*sd(df$summer_max_t))+mean(df$summer_max_t)}

plot_temp_g <- ggplot(plot.dat, aes(x = inv_fun4(x2.sim), y = bayes.c.eff.mean)) +
  geom_line(color = "black", alpha = 0.8, size = 1.8)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "black", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "black", alpha = 0.1)+
  ggtitle("")+
  xlab("Summer temperature")+
  ylab("Offspring mass")+
  ylim(c(150,800))+
  theme_bw() +
  coord_cartesian(clip = "off") + 
  theme(plot.margin = margin(5.5,5.5,14,5.5))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size=16),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13))
plot_temp_g <- plot_temp_g + geom_point(data=df, aes(x = summer_max_t, y=massjun), alpha=0.1, color = "black",
                                        position = position_jitter(w=0.07, h = 0), inherit.aes = F)
plot_temp_g


#### fecundity ####

## age ##

x2.sim <- seq(min(scale(df_fecundity$age)[,1]), max(scale(df_fecundity$age)[,1]), by = 0.01)

int.sim <- matrix(rep(NA, nrow(f_effects)*length(x2.sim)), nrow = nrow(f_effects))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- exp(f_effects$V1 + f_effects$V3 * (x2.sim[i]) + f_effects$V4 * ((x2.sim[i])^2))
}

bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.055)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.945)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean, bayes.c.eff.lower, bayes.c.eff.upper, bayes.c.eff.lower.bis, bayes.c.eff.upper.bis)

inv_fun5 <- function(x){(x*sd(df_fecundity$age))+mean(df_fecundity$age)}

plot_age_f <- ggplot(plot.dat, aes(x = inv_fun5(x2.sim), y = bayes.c.eff.mean)) +
  geom_line(color = "black", alpha = 0.8, size = 1.8)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "black", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "black", alpha = 0.1)+
  ggtitle("")+
  xlab("Age")+
  ylab("Litter size")+
  theme_bw() +
  coord_cartesian(clip = "off") + 
  theme(plot.margin = margin(5.5,5.5,14,5.5))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size=16),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13))
plot_age_f <- plot_age_f + geom_point(data=df_fecundity, aes(x = age, y=nboffspring), alpha=0.1, color = "black",
                                      position = position_jitter(w=0.15, h = 0.15), inherit.aes = F)
plot_age_f


## mass ##

x2.sim <- seq(min(scale(df_fecundity$massjun_mom)[,1]), max(scale(df_fecundity$massjun_mom)[,1]), by = 0.01)

int.sim <- matrix(rep(NA, nrow(f_effects)*length(x2.sim)), nrow = nrow(f_effects))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- exp(f_effects$V1 + f_effects$V5 * (x2.sim[i]))
}

bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.055)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.945)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean, bayes.c.eff.lower, bayes.c.eff.upper, bayes.c.eff.lower.bis, bayes.c.eff.upper.bis)

inv_fun6 <- function(x){(x*sd(df_fecundity$massjun_mom))+mean(df_fecundity$massjun_mom)}

plot_mass_f <- ggplot(plot.dat, aes(x = inv_fun6(x2.sim), y = bayes.c.eff.mean)) +
  geom_line(color = "black", alpha = 0.8, size = 1.8)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "black", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "black", alpha = 0.1)+
  ggtitle("")+
  xlab("Mass")+
  ylab("Litter size")+
  theme_bw() +
  coord_cartesian(clip = "off") + 
  theme(plot.margin = margin(5.5,5.5,14,5.5))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size=16),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13))
plot_mass_f <- plot_mass_f + geom_point(data=df_fecundity, aes(x = massjun_mom, y=nboffspring), alpha=0.1, color = "black",
                                      position = position_jitter(w=0, h = 0.15), inherit.aes = F)
plot_mass_f


## snow ##

x2.sim <- seq(min(scale(df_fecundity$total_snow)[,1]), max(scale(df_fecundity$total_snow)[,1]), by = 0.01)

int.sim <- matrix(rep(NA, nrow(f_effects)*length(x2.sim)), nrow = nrow(f_effects))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- exp(f_effects$V1 + f_effects$V2 * (x2.sim[i]))
}

bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.055)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.945)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean, bayes.c.eff.lower, bayes.c.eff.upper, bayes.c.eff.lower.bis, bayes.c.eff.upper.bis)

inv_fun7 <- function(x){(x*sd(df_fecundity$total_snow))+mean(df_fecundity$total_snow)}

plot_snow_f <- ggplot(plot.dat, aes(x = inv_fun7(x2.sim), y = bayes.c.eff.mean)) +
  geom_line(color = "black", alpha = 0.8, size = 1.8)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "black", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "black", alpha = 0.1)+
  ggtitle("")+
  xlab("Amount of winter snow")+
  ylab("Litter size")+
  theme_bw() +
  coord_cartesian(clip = "off") + 
  theme(plot.margin = margin(5.5,5.5,14,5.5))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size=16),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13))
plot_snow_f <- plot_snow_f + geom_point(data=df_fecundity, aes(x = total_snow, y=nboffspring), alpha=0.1, color = "black",
                                        position = position_jitter(w=5, h = 0.15), inherit.aes = F)
plot_snow_f


## combine plots 

#ragg::agg_tiff("marmot_supplementary_plot.tiff", width = 12, height = 10, units = "in", res = 300)
wrap_elements(full= (plot_snow_f | plot_age_f | plot_mass_f) / ((plot_temp_g | plot_snow_g | plot_age_g | plot_mass_g)))
#dev.off()
