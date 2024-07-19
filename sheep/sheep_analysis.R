# Last modified 19/07/2024

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
library(bayesplot)

###################################
#### load formatted sheep data ####
###################################

# set your own path
df <- read.delim("~/My Drive/phd/phd_simulation/Chapter 2/data_code_for_submission/sheep/sheep_data_for_analysis.txt")

################################
#### dataset for Stan model ####
################################

# matrix of predictors for correlation submodel
X <- as.matrix(cbind(rep(1, length(unique(df$obsY))), 
                     scale(unique(df[,c("obsY","nao_winter")])$nao_winter)[,1],
                     scale(unique(df[,c("obsY","abundance")])$abundance)[,1], 
                     scale(unique(df[,c("obsY","nao_winter")])$nao_winter)[,1] * scale(unique(df[,c("obsY","abundance")])$abundance)[,1]))

# matrix of predictors for growth submodel
X_g = as.matrix(cbind(rep(1, dim(df)[1]), 
                      scale(df$nao_winter)[,1],
                      scale(df$abundance)[,1],
                      scale(df$ageY)[,1],
                      (scale(df$ageY)[,1])^2,
                      scale(df$capWgt)[,1],
                      scale(df$abundance)[,1] * scale(df$nao_winter)[,1]))

# matrix of predictors for fecundity submodel
X_f = as.matrix(cbind(scale(df$abundance)[,1],
                      scale(df$ageY)[,1],
                      (scale(df$ageY)[,1])^2,
                      scale(df$capWgt)[,1]))


cn <- as.numeric(table(as.numeric(as.factor(df$obsY))))

#create empty cmat - matrix of context-specific IDs
cmat <- matrix(NA, 
               nrow = length(unique(df$obsY)), 
               ncol =  max(as.numeric(table(as.numeric(as.factor(df$obsY))))))

#fill cmat
temporary <- as.data.frame(cbind(as.numeric(as.factor(df$obsY)),
                                 as.numeric(as.factor(df$id))))
for (i in 1:length(unique(df$obsY))) {
  cmat[i, 1:cn[i]] <- temporary$V2[temporary$V1 == i]
}
cmat_n = apply(cmat, 1, FUN = function(x) sum(!is.na(x)) )
cmat[is.na(cmat)] = 0 #remove NAs

temp = t(cmat)
corder = data.frame(id = temp[temp>0], c = rep(seq(1:nrow(cmat)), times = cmat_n))

idc = match(paste0(as.numeric(as.factor(df$id)),
                     as.numeric(as.factor(df$obsY)),
                     sep="."), 
              paste0(corder$id,corder$c,sep="."))

rownames(df) <- NULL


stan.df =
  list(N = nrow(df),                       # number of observation
       C = length(unique(df$obsY)),        # number of years
       I = length(unique(df$id)),          # number of mothers
       D = 2,                              # number of traits
       P_y = 4,                            # number of predictor on correlation
       P_g = 7,                            # number of predictor on growth
       P_f = 4,                            # number of predictor on fecundity
       
       id = as.numeric(as.factor(df$id)),       # index of mother id  
       c_id = as.numeric(as.factor(df$obsY)),   # index of year id
       idc = idc,                               # index of reproduction event id
       id_lm = as.numeric(rownames(df)),        # index of row names
       
       X = X,                                   # matrix of predictors for correlation submodel
       X_g = X_g,                               # matrix of predictors for growth submodel
       X_f = X_f,                               # matrix of predictors for fecundity submodel
       A = diag(length(unique(df$id))),
       
       cm = max(as.numeric(table(as.numeric(as.factor(df$obsY))))),
       cmat = cmat,
       cn = cn,
       cnt = length(unique(paste(df$id, df$obsY))), # number of reproduction events
       
       growth = scale((as.numeric(df$growth)))[,1],    # offspring mass (response variable)
       productivity = df$lambNum1+1    # fecundity (response variable), [0,2] but coded as [1,3] for the ordinal logistic regression
  )



########################
#### run Stan model ####
########################

# set your own path, where the stan model is
setwd("~/sheep")

# Compile model
mod <- cmdstan_model("sheep_model.stan"
                     , stanc_options = list("O1")
)

# Fit model 
fit <- mod$sample(
  data = stan.df, 
  seed = 156789, 
  chains = 3, 
  parallel_chains = 3,
  iter_warmup = 1000,
  iter_sampling = 3000,
  adapt_delta = 0.975,
  refresh = 40 # print update every XX iters
)

#####################################
#### Posterior Predictive Checks ####
#####################################

# Extract y_rep_f for posterior predictive check
# This is to check if fecundity (litter size) predicted by the model concur with the real data
draws_f <- fit$draws(
  variables = "y_rep_f",
  inc_warmup = FALSE,
  format = "matrix"
)

#Plot for posterior predictive checks
plot_ppc_f <- ppc_dens_overlay(
  stan.df$productivity-1,
  draws_f[1:1000,]-1,
  size = 0.25,
  alpha = 0.7,
  trim = FALSE
)
plot_ppc_f <- plot_ppc_f + 
  xlab("Litter size") + 
  ylab("Density") + 
  ggtitle("Posterior predictive check \u2014 Sheep fecundity")
plot_ppc_f

# Extract y_rep_g for posterior predictive check
# This is to check if growth predicted by the model concur with the real data
draws_g <- fit$draws(
  variables = "y_rep_g",
  inc_warmup = FALSE,
  format = "matrix"
)

#Plot for posterior predictive checks
plot_ppc_g <- ppc_dens_overlay(
  stan.df$growth ,
  draws_g[1:1000,],
  size = 0.25,
  alpha = 0.7,
  trim = FALSE
)
plot_ppc_g <- plot_ppc_g + 
  xlab("Mass t+1") + 
  ylab("Density") + 
  ggtitle("Posterior predictive check \u2014 Sheep mass t+1")
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
#saveRDS(CRN_fit, "model_sheep_output.RDS")
#CRN_fit <- read_rds("model_sheep_output.RDS")

# Open shiny app with diagnostic plot (visualize the chains to check for convergence)
launch_shinystan(CRN_fit)

post <- rstan::extract(CRN_fit)

traits = c("growth", "productivity")

# For the following plots, need to first run R script "functions_postprocessing.R", that contains some functions needed to make predictions using the model output

#########################
#### Plot winter NAO ####
#########################

# sequence from scaled minimum to maximum value in x
seq = seq(min(scale(unique(df[,c("obsY","nao_winter")])$nao_winter)[,1]),
          max(scale(unique(df[,c("obsY","nao_winter")])$nao_winter)[,1]),
          by =  (max(scale(unique(df[,c("obsY","nao_winter")])$nao_winter)[,1]) - min(scale(unique(df[,c("obsY","nao_winter")])$nao_winter)[,1])) / 25 ) 

# high pop density (set to +1sd)
X_pred = data.frame(int = 1, X1 = seq, X2 = 1, X3=1*seq)
mv_cpc = cpc_f(x = X_pred, cpc_b = post$B_cpc)
mv_cor = cor_f(x = X_pred, traits = traits, cpc = mv_cpc)

# low pop density (set to -1sd)
X_pred2 = data.frame(int = 1, X1 = seq, X2 = -1, X3=-1*seq)
mv_cpc2 = cpc_f(x = X_pred2, cpc_b = post$B_cpc)
mv_cor2 = cor_f(x = X_pred2, traits = traits, cpc = mv_cpc2)

# function to unscale predictions (because predictors were standardized)
inv_fun_p1 <- function(x){(x*sd(unique(df[,c("obsY","nao_winter")])$nao_winter))+mean(unique(df[,c("obsY","nao_winter")])$nao_winter)}
  
plot_1 <- ggplot(mv_cor,
                 aes(x = inv_fun_p1(X1), y = value))+
  stat_lineribbon(linewidth=1.8, .width = 0.89, alpha=0.15, fill="#009E73")+
  stat_lineribbon(linewidth=1.8, .width = 0.5, alpha=0.30, fill="#009E73")+
  stat_lineribbon(linewidth=1.8, .width = 0.001, color="#009E73")+
  ylim(c(-1,1))+
  ggtitle("")+
  xlab("Winter NAO")+
  ylab("Observation-level correlation")+
  theme_bw() +
  coord_cartesian(clip = "off") + 
  annotation_custom(grid::textGrob(expression(~ italic("(-) Environmental harshness (+)")),
                                   gp=grid::gpar(fontsize=10)),
                    xmin = inv_fun_p1(mean(range(scale(unique(df[,c("obsY","nao_winter")])$nao_winter)[,1]))),
                    xmax = inv_fun_p1(mean(range(scale(unique(df[,c("obsY","nao_winter")])$nao_winter)[,1]))),
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
mesure_pol <- data.frame(x1 = c(-2.2,0.2),
                         x2 = c(-1.8,0.6),
                         cat = c(1:2),
                         catnames = c("High density","Low density"))

ann <- ggplot(data = mesure_pol) +
  geom_rect(aes(xmin = x1,
                xmax = x2,
                ymin = 0.095,
                ymax = 0.105,
                fill = as.factor(cat)),
            fill = c("#009E73", "#E69F00"),
            color = "black",
            size = 0.3) +
  geom_text(aes(x = x2+0.85, y = 0.1, label = catnames),
            #vjust = .8, 
            fontface = "bold", color = "black") +
  coord_cartesian(xlim = c(-2.7085, 2.861)
                  , ylim = c(0.065, 0.14)) +
  theme_void()
ann

nao_plot <- ann / plot_1  + 
            plot_layout(heights = c(1, 8))
nao_plot

#################################
#### Plot population density ####
#################################

# sequence from scaled minimum to maximum value in x
seq = seq(min(scale(unique(df[,c("obsY","abundance")])$abundance)[,1]),
          max(scale(unique(df[,c("obsY","abundance")])$abundance)[,1]),
          by = (max(scale(unique(df[,c("obsY","abundance")])$abundance)[,1]) - min(scale(unique(df[,c("obsY","abundance")])$abundance)[,1])) / 25 ) 

#high nao (+1sd)
X_pred = data.frame(int = 1, X1 = 1, X2 = seq, X3=1*seq)
mv_cpc = cpc_f(x = X_pred, cpc_b = post$B_cpc)
mv_cor = cor_f(x = X_pred, traits = traits, cpc = mv_cpc)

#low nao (-1sd)
X_pred2 = data.frame(int = 1, X1 = -1, X2 = seq, X3=-1*seq)
mv_cpc2 = cpc_f(x = X_pred2, cpc_b = post$B_cpc)
mv_cor2 = cor_f(x = X_pred2, traits = traits, cpc = mv_cpc2)

# function to unscale predictions (because predictors were standardized)
inv_fun_p2 <- function(x){(x*sd(unique(df[,c("obsY","abundance")])$abundance))+mean(unique(df[,c("obsY","abundance")])$abundance)}

plot_2 <- ggplot(mv_cor,
                 aes(x = inv_fun_p2(X2), y = value))+
  stat_lineribbon(linewidth=1.8, .width = 0.89, alpha=0.15, fill="#009E73")+
  stat_lineribbon(linewidth=1.8, .width = 0.5, alpha=0.15, fill="#009E73")+
  stat_lineribbon(linewidth=1.8, .width = 0.001, color="#009E73")+
  ylim(c(-1,1))+
  ggtitle("")+
  xlab("Population density")+
  ylab("Observation-level correlation")+
  theme_bw() +
  coord_cartesian(clip = "off") + 
  annotation_custom(grid::textGrob(expression(~ italic("(-) Environmental harshness (+)")),
                                   gp=grid::gpar(fontsize=10)),
                    xmin = inv_fun_p2(mean(range(scale(unique(df[,c("obsY","abundance")])$abundance)[,1]))), 
                    xmax = inv_fun_p2(mean(range(scale(unique(df[,c("obsY","abundance")])$abundance)[,1]))), 
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
mesure_pol2 <- data.frame(x1 = c(-2.2,0.2),
                         x2 = c(-1.8,0.6),
                         cat = c(1:2),
                         catnames = c("High NAO","Low NAO"))

ann2 <- ggplot(data = mesure_pol2) +
  geom_rect(aes(xmin = x1,
                xmax = x2,
                ymin = 0.095,
                ymax = 0.105,
                fill = as.factor(cat)),
            fill = c("#009E73", "#E69F00"),
            color = "black",
            size = 0.3) +
  geom_text(aes(x = x2+0.8, y = 0.1, label = catnames),
            #vjust = .8, 
            fontface = "bold", color = "black") +
  coord_cartesian(xlim = c(-2.7085, 2.861)
                  , ylim = c(0.065, 0.14)) +
  theme_void()
ann2

density_plot <- ann2 / plot_2  + 
  plot_layout(heights = c(1, 8))
density_plot

wrap_elements(full= (nao_plot | density_plot))


###########################
#### Figure posteriors ####
###########################

#make sure you know the order of the variables in the predictor matrices
g_effects <- as.data.frame(post$B_m_g)
f_effects <- as.data.frame(post$B_m_f)
cor_effects <- as.data.frame(post$B_cpc)
cutpoints <- as.data.frame(post$cutpoint)

# change the "9000" to another value if you run chains for less or more iterations than I did (3000 iterations post burn-in * 3 = 9000)

df.posteriors <- data_frame(Submodel = c(rep("Mass t+1", 9000*6), rep("Fecundity", 9000*4), rep("Correlation", 9000*3))
                            , parameter = c(rep("Winter NAO", 9000), rep("Density", 9000), rep("Age", 9000), rep("Age^2", 9000), rep("Mass", 9000), rep("Winter NAO * Density", 9000)
                                            , rep("Density", 9000), rep("Age", 9000), rep("Age^2", 9000), rep("Mass", 9000)
                                            , rep("Winter NAO", 9000), rep("Density", 9000), rep("Winter NAO * Density", 9000))
                            , Posterior = c(g_effects$V2, g_effects$V3, g_effects$V4, g_effects$V5
                                            ,g_effects$V6, g_effects$V7
                                            , f_effects$V1, f_effects$V2, f_effects$V3, f_effects$V4
                                            , cor_effects[,2], cor_effects[,3], cor_effects[,4]))


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


#ragg::agg_tiff("Figure 5.tiff", width = 9, height = 10, units = "in", res = 300)
wrap_elements(full= (nao_plot | density_plot)) / p2
#dev.off()



#################################################
#### Plot individual relationships on traits ####
#################################################

#### Fecundity ####

## Age ##

x2.sim <- seq(min(scale(df$ageY)[,1]), max(scale(df$ageY)[,1]), by = 0.01)

int.sim <- matrix(rep(NA, nrow(f_effects)*length(x2.sim)), nrow = nrow(f_effects))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- 
    # 2*prob of having 2 offspring
    (2 * (1 - (exp(cutpoints$V2 - f_effects$V2 * (x2.sim[i]) - f_effects$V3 * ((x2.sim[i])^2))/
                    (1+exp(cutpoints$V2 - f_effects$V2 * (x2.sim[i]) - f_effects$V3 * ((x2.sim[i])^2)))))) +
    # 1*prob of having 1 offspring
    ((exp(cutpoints$V2 - f_effects$V2 * (x2.sim[i]) - f_effects$V3 * ((x2.sim[i])^2))/
        (1+exp(cutpoints$V2 - f_effects$V2 * (x2.sim[i]) - f_effects$V3 * ((x2.sim[i])^2)))) -
       (exp(cutpoints$V1 - f_effects$V2 * (x2.sim[i]) - f_effects$V3 * ((x2.sim[i])^2))/
          (1+exp(cutpoints$V1 - f_effects$V2 * (x2.sim[i]) - f_effects$V3 * ((x2.sim[i])^2)))))
}

bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.055)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.945)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean, bayes.c.eff.lower, bayes.c.eff.upper, bayes.c.eff.lower.bis, bayes.c.eff.upper.bis)

inv_fun <- function(x){(x*sd(df$ageY))+mean(df$ageY)}

plot_age_f <- ggplot(plot.dat, aes(x = inv_fun(x2.sim), y = bayes.c.eff.mean)) +
  geom_line(color = "black", alpha = 0.8, size = 1.8)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "black", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "black", alpha = 0.1)+
  ggtitle("")+
  xlab("Age")+
  ylab("Fecundity")+
  theme_bw() +
  coord_cartesian(clip = "off") + 
  theme(plot.margin = margin(5.5,5.5,14,5.5))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size=16),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13))
plot_age_f <- plot_age_f + geom_point(data=df, aes(x = ageY, y=lambNum1), alpha=0.1, color = "black",
                                      position = position_jitter(w=0.15, h = 0.25), inherit.aes = F)
plot_age_f


## mass ##

x2.sim <- seq(min(scale(df$capWgt)[,1]), max(scale(df$capWgt)[,1]), by = 0.01)

int.sim <- matrix(rep(NA, nrow(f_effects)*length(x2.sim)), nrow = nrow(f_effects))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- 
    # 2*prob of having 2 offspring
    (2 * (1 - (exp(cutpoints$V2 - f_effects$V4 * (x2.sim[i]) )/
                 (1+exp(cutpoints$V2 - f_effects$V4 * (x2.sim[i]) ))))) +
    # 1*prob of having 1 offspring
    ((exp(cutpoints$V2 - f_effects$V4 * (x2.sim[i]) )/
        (1+exp(cutpoints$V2 - f_effects$V4 * (x2.sim[i]) ))) -
       (exp(cutpoints$V1 - f_effects$V4 * (x2.sim[i]) )/
          (1+exp(cutpoints$V1 - f_effects$V4 * (x2.sim[i]) ))))
}

bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.055)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.945)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean, bayes.c.eff.lower, bayes.c.eff.upper, bayes.c.eff.lower.bis, bayes.c.eff.upper.bis)

inv_fun2 <- function(x){(x*sd(df$capWgt))+mean(df$capWgt)}

plot_mass_f <- ggplot(plot.dat, aes(x = inv_fun2(x2.sim), y = bayes.c.eff.mean)) +
  geom_line(color = "black", alpha = 0.8, size = 1.8)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "black", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "black", alpha = 0.1)+
  ggtitle("")+
  xlab("Mass")+
  ylab("Fecundity")+
  theme_bw() +
  coord_cartesian(clip = "off") + 
  theme(plot.margin = margin(5.5,5.5,14,5.5))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size=16),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13))
plot_mass_f <- plot_mass_f + geom_point(data=df, aes(x = capWgt, y=lambNum1), alpha=0.1, color = "black",
                                      position = position_jitter(w=0, h = 0.25), inherit.aes = F)
plot_mass_f


## abundance ##

x2.sim <- seq(min(scale(df$abundance)[,1]), max(scale(df$abundance)[,1]), by = 0.01)

int.sim <- matrix(rep(NA, nrow(f_effects)*length(x2.sim)), nrow = nrow(f_effects))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- 
    # 2*prob of having 2 offspring
    (2 * (1 - (exp(cutpoints$V2 - f_effects$V1 * (x2.sim[i]) )/
                 (1+exp(cutpoints$V2 - f_effects$V1 * (x2.sim[i]) ))))) +
    # 1*prob of having 1 offspring
    ((exp(cutpoints$V2 - f_effects$V1 * (x2.sim[i]) )/
        (1+exp(cutpoints$V2 - f_effects$V1 * (x2.sim[i]) ))) -
       (exp(cutpoints$V1 - f_effects$V1 * (x2.sim[i]) )/
          (1+exp(cutpoints$V1 - f_effects$V1 * (x2.sim[i]) ))))
}

bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.055)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.945)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean, bayes.c.eff.lower, bayes.c.eff.upper, bayes.c.eff.lower.bis, bayes.c.eff.upper.bis)

inv_fun3 <- function(x){(x*sd(df$abundance))+mean(df$abundance)}

plot_abundance_f <- ggplot(plot.dat, aes(x = inv_fun3(x2.sim), y = bayes.c.eff.mean)) +
  geom_line(color = "black", alpha = 0.8, size = 1.8)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "black", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "black", alpha = 0.1)+
  ggtitle("")+
  xlab("Density")+
  ylab("Fecundity")+
  theme_bw() +
  coord_cartesian(clip = "off") + 
  theme(plot.margin = margin(5.5,5.5,14,5.5))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size=16),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13))
plot_abundance_f <- plot_abundance_f + geom_point(data=df, aes(x = abundance, y=lambNum1), alpha=0.1, color = "black",
                                        position = position_jitter(w=0.1, h = 0.25), inherit.aes = F)
plot_abundance_f




#### growth ####

## Age ##

x2.sim <- seq(min(scale(df$ageY)[,1]), max(scale(df$ageY)[,1]), by = 0.01)

int.sim <- matrix(rep(NA, nrow(g_effects)*length(x2.sim)), nrow = nrow(g_effects))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- ((g_effects$V1 + g_effects$V4 * (x2.sim[i]) + g_effects$V5 * ((x2.sim[i])^2)) * sd(as.numeric(df$growth))) + mean(as.numeric(df$growth))
}

bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.055)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.945)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean, bayes.c.eff.lower, bayes.c.eff.upper, bayes.c.eff.lower.bis, bayes.c.eff.upper.bis)

inv_fun4 <- function(x){(x*sd(df$ageY))+mean(df$ageY)}

plot_age_g <- ggplot(plot.dat, aes(x = inv_fun4(x2.sim), y = bayes.c.eff.mean)) +
  geom_line(color = "black", alpha = 0.8, size = 1.8)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "black", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "black", alpha = 0.1)+
  ggtitle("")+
  xlab("Age")+
  ylab("Mass t+1")+
  theme_bw() +
  coord_cartesian(clip = "off") + 
  theme(plot.margin = margin(5.5,5.5,14,5.5))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size=16),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13))
plot_age_g <- plot_age_g + geom_point(data=df, aes(x = ageY, y=growth), alpha=0.1, color = "black",
                                      position = position_jitter(w=0.15, h = 0), inherit.aes = F)
plot_age_g


## abundance ##

x2.sim <- seq(min(scale(df$abundance)[,1]), max(scale(df$abundance)[,1]), by = 0.01)

int.sim <- matrix(rep(NA, nrow(g_effects)*length(x2.sim)), nrow = nrow(g_effects))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- ((g_effects$V1 + g_effects$V3 * (x2.sim[i]) ) * sd(as.numeric(df$growth))) + mean(as.numeric(df$growth))
}

bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.055)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.945)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean, bayes.c.eff.lower, bayes.c.eff.upper, bayes.c.eff.lower.bis, bayes.c.eff.upper.bis)

inv_fun5 <- function(x){(x*sd(df$abundance))+mean(df$abundance)}

plot_abundance_g <- ggplot(plot.dat, aes(x = inv_fun5(x2.sim), y = bayes.c.eff.mean)) +
  geom_line(color = "black", alpha = 0.8, size = 1.8)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "black", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "black", alpha = 0.1)+
  ggtitle("")+
  xlab("Density")+
  ylab("Mass t+1")+
  theme_bw() +
  coord_cartesian(clip = "off") + 
  theme(plot.margin = margin(5.5,5.5,14,5.5))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size=16),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13))
plot_abundance_g <- plot_abundance_g + geom_point(data=df, aes(x = abundance, y=growth), alpha=0.1, color = "black",
                                                  position = position_jitter(w=0, h = 0), inherit.aes = F)
plot_abundance_g


## winter nao ##

x2.sim <- seq(min(scale(df$nao_winter)[,1]), max(scale(df$nao_winter)[,1]), by = 0.01)

int.sim <- matrix(rep(NA, nrow(g_effects)*length(x2.sim)), nrow = nrow(g_effects))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- ((g_effects$V1 + g_effects$V2 * (x2.sim[i]) ) * sd(as.numeric(df$growth))) + mean(as.numeric(df$growth))
}

bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.055)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.945)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean, bayes.c.eff.lower, bayes.c.eff.upper, bayes.c.eff.lower.bis, bayes.c.eff.upper.bis)

inv_fun6 <- function(x){(x*sd(df$nao_winter))+mean(df$nao_winter)}

plot_nao_g <- ggplot(plot.dat, aes(x = inv_fun6(x2.sim), y = bayes.c.eff.mean)) +
  geom_line(color = "black", alpha = 0.8, size = 1.8)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "black", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "black", alpha = 0.1)+
  ggtitle("")+
  xlab("Winter NAO")+
  ylab("Mass t+1")+
  theme_bw() +
  coord_cartesian(clip = "off") + 
  theme(plot.margin = margin(5.5,5.5,14,5.5))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size=16),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13))
plot_nao_g <- plot_nao_g + geom_point(data=df, aes(x = nao_winter, y=growth), alpha=0.1, color = "black",
                                                  position = position_jitter(w=0, h = 0), inherit.aes = F)
plot_nao_g


## mass ##

x2.sim <- seq(min(scale(df$capWgt)[,1]), max(scale(df$capWgt)[,1]), by = 0.01)

int.sim <- matrix(rep(NA, nrow(g_effects)*length(x2.sim)), nrow = nrow(g_effects))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- ((g_effects$V1 + g_effects$V6 * (x2.sim[i]) ) * sd(as.numeric(df$growth))) + mean(as.numeric(df$growth))
}

bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.055)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.945)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean, bayes.c.eff.lower, bayes.c.eff.upper, bayes.c.eff.lower.bis, bayes.c.eff.upper.bis)

inv_fun7 <- function(x){(x*sd(df$capWgt))+mean(df$capWgt)}

plot_mass_g <- ggplot(plot.dat, aes(x = inv_fun7(x2.sim), y = bayes.c.eff.mean)) +
  geom_line(color = "black", alpha = 0.8, size = 1.8)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "black", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "black", alpha = 0.1)+
  ggtitle("")+
  xlab("Mass t")+
  ylab("Mass t+1")+
  theme_bw() +
  coord_cartesian(clip = "off") + 
  theme(plot.margin = margin(5.5,5.5,14,5.5))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size=16),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13))
plot_mass_g <- plot_mass_g + geom_point(data=df, aes(x = capWgt, y=growth), alpha=0.1, color = "black",
                                        position = position_jitter(w=0, h = 0), inherit.aes = F)
plot_mass_g


## combine plots 

# ragg::agg_tiff("sheep_association_plot.tiff", width = 12, height = 10, units = "in", res = 300)
wrap_elements(full= (plot_abundance_f | plot_age_f | plot_mass_f) / ((plot_nao_g | plot_abundance_g | plot_age_g | plot_mass_g)))
# dev.off()
