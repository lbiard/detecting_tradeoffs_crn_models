# Last modified 19/07/2024

# This R script can be used to run the non-repeated measures CRN model on simulated data
# First you need to run the R script "simulate_data_non_repeated_crn.R" to create the simulated datasets.
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

df1 = inds_hist_storage[[1]] # select first dataset

########################################
#### Prepare dataset for Stan model ####
########################################

df <- df1


################################
#### dataset for Stan model ####
################################

X <- as.matrix(cbind(rep(1, length(unique(df$year))), 
                     scale(unique(df[,c("year","climate")])$climate)[,1]))

X_g = as.matrix(cbind(rep(1, dim(df)[1]), 
                      scale(df$climate)[,1]))

X_f = as.matrix(cbind(rep(1, dim(df)[1]),
                      scale(df$climate)[,1]))


cn <- as.numeric(table(as.numeric(as.factor(df$year))))

#create empty cmat
cmat <- matrix(NA, 
               nrow = length(unique(df$year)), 
               ncol =  max(as.numeric(table(as.numeric(as.factor(df$year))))))

#fill cmat
temporary <- as.data.frame(cbind(as.numeric(as.factor(df$year)),
                                 as.numeric(as.factor(df$id))))
for (i in 1:length(unique(df$year))) {
  cmat[i, 1:cn[i]] <- temporary$V2[temporary$V1 == i]
}
cmat_n = apply(cmat, 1, FUN = function(x) sum(!is.na(x)) )
cmat[is.na(cmat)] = 0 #remove NAs

temp = t(cmat)
corder = data.frame(id = temp[temp>0], c = rep(seq(1:nrow(cmat)), times = cmat_n))

idc = match(paste0(as.numeric(as.factor(df$id)),
                   as.numeric(as.factor(df$year)),
                   sep="."), 
            paste0(corder$id,corder$c,sep="."))

rownames(df) <- NULL


stan.df =
  list(N = nrow(df),                      # number of observation offpring mass
       C = length(unique(df$year)),     # number of years
       I = length(unique(df$id)),        # number of mothers
       D = 2,                             # numver of traits
       P_y = 2,
       P_g = 2,
       P_f = 2,
       
       id = as.numeric(as.factor(df$id)),           
       c_id = as.numeric(as.factor(df$year)),
       idc = idc,
       id_lm = as.numeric(rownames(df)),
       
       X = X,
       X_g = X_g,
       X_f = X_f,
       A = diag(length(unique(df$id))),
       
       cm = max(as.numeric(table(as.numeric(as.factor(df$year))))),
       cmat = cmat,
       cn = cn,
       cnt = length(unique(paste(df$id, df$year))), # number of repro events
       
       growth = scale(as.numeric(df$growth))[,1],             # offspring mass (response variable)
       productivity = df$productivity    # fecundity (response variable)
  )



####################################
#### Compile and run Stan model ####
####################################

mod <- cmdstan_model("model_non_repeated_crn.stan"
                     , stanc_options = list("O1")
)


fit <- mod$sample(
  data = stan.df, 
  seed = 1234, 
  chains = 3, 
  parallel_chains = 3,
  adapt_delta = 0.99,
  refresh = 20 # print update every 20 iters
)

stanfit <- rstan::read_stan_csv(fit$output_files())
#saveRDS(stanfit, "model_non_repeated_crn_output.RDS")
#stanfit <- read_rds("model_non_repeated_crn_output.RDS")

# Open shiny app with diagnostic plot (visualize the chains to check for convergence)
#launch_shinystan(stanfit)


######################
#### make figures ####
######################

fit <- stanfit

chain_1 <- fit@sim[["samples"]][[1]][,c("B_m_g.1", "B_m_g.2", "B_m_f.1", "B_m_f.2", "B_cpc.1.1", "B_cpc.2.1")]
chain_2 <- fit@sim[["samples"]][[2]][,c("B_m_g.1", "B_m_g.2", "B_m_f.1", "B_m_f.2", "B_cpc.1.1", "B_cpc.2.1")]
chain_3 <- fit@sim[["samples"]][[3]][,c("B_m_g.1", "B_m_g.2", "B_m_f.1", "B_m_f.2", "B_cpc.1.1", "B_cpc.2.1")]

dat_plot <- rbind(chain_1, chain_2, chain_3)

#################
#### climate ####
#################

# predictions across range of climatic values
x2.sim <- seq(min(scale(unique(df[,c("year","climate")])$climate)[,1]),
              max(scale(unique(df[,c("year","climate")])$climate)[,1]), by = 0.02)

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- tanh(dat_plot$B_cpc.1.1 + dat_plot$B_cpc.2.1 * (x2.sim[i])) 
}

# calculate quantiles of predictions
bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.045)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.945)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean, bayes.c.eff.lower, bayes.c.eff.upper, bayes.c.eff.lower.bis, bayes.c.eff.upper.bis)

inv_fun_p <- function(x){(x*sd(unique(df[,c("year","climate")])$climate))+mean(unique(df[,c("year","climate")])$climate)}


p <- ggplot(plot.dat, aes(x = inv_fun_p(x2.sim), y = bayes.c.eff.mean)) +
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

df.posteriors <- data_frame(Submodel = c(rep("Growth", 3000*2), rep("Fecundity", 3000*2))
                            , parameter = c(rep("Intercept", 3000), rep("Climate", 3000)
                                            , rep("Intercept", 3000), rep("Climate", 3000))
                            , Posterior = c(dat_plot$B_m_g.1, dat_plot$B_m_g.2
                                            , dat_plot$B_m_f.1, dat_plot$B_m_f.2))



df.posteriors$Submodel <- factor(df.posteriors$Submodel, levels=c("Growth", "Fecundity"))
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

df.expected <- data_frame(Submodel = df.summary$Submodel
                          , parameter = df.summary$parameter
                          , mean = c(0, log(1), 0.3, 0.3))

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
#colCov <- "seagreen4"
colPlot <- "black"


# plot
q <- ggplot()+
  stat_halfeye(data = df.posteriors, aes(x = Posterior, y = parameter, fill = Submodel), color = NA, alpha = 0.2, position = position_dodge(width = dodge.width), normalize="xy", scale=0.7)+
  geom_point(data = df.summary, aes(x = mean, y = parameter, color = Submodel), size = 2, position = position_dodge(width = dodge.width))+
  geom_linerange(data = df.summary, aes(xmin = BCI89_lower, xmax = BCI89_upper, y = parameter, color = Submodel), size=0.6, linetype="solid", position = position_dodge(width = dodge.width))+
  geom_linerange(data = df.summary, aes(xmin = BCI50_lower, xmax = BCI50_upper, y = parameter, color = Submodel), size=1.5, linetype="solid", position = position_dodge(width = dodge.width))+
  scale_color_manual(values = c(colT1, colT2))+
  scale_fill_manual(values = c(colT1, colT2))+
  scale_alpha_manual(values = c(1,1))+
  geom_segment(aes(x = 0 , y = 1.6, xend = 0, yend = 2.7), linetype=2) +
  geom_segment(aes(x = 0.3 , y = 0.6, xend = 0.3, yend = 1.7), linetype=2) +
  ylab("")+
  theme_minimal()+
  theme(axis.line.y = element_blank()
        , axis.ticks.y = element_blank()
        , panel.grid.major = element_blank() 
        , panel.grid.minor = element_blank()
  )
q

# ragg::agg_tiff("Figure 3.tiff", width = 9, height = 6, units = "in", res = 300)
wrap_elements(full= (p | q))
# dev.off()
