# Last modified 07/05/2024

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
launch_shinystan(stanfit)

post <- rstan::extract(stanfit)

traits = c("growth", "productivity")


######################
#### make figures ####
######################

#################
#### climate ####
#################

seq = seq(-2.2,2.2, by = 0.2) #standardized values from -2 to +2
X_pred = data.frame(int = 1, X1 = seq)
mv_cpc = cpc_f(x = X_pred, cpc_b = post$B_cpc)
mv_cor = cor_f(x = X_pred, traits = traits, cpc = mv_cpc)

plot_1 <- ggplot(mv_cor,
                 aes(x = X1, y = value))+
  stat_lineribbon(linewidth=1.8, .width = 0.89, alpha=0.15, fill="black")+
  stat_lineribbon(linewidth=1.8, .width = 0.5, alpha=0.15, fill="black")+
  stat_lineribbon(linewidth=1.8, .width = 0.001, fill="black")+
  ylim(c(-1,1))+
  ggtitle("")+
  xlab("Climate")+
  ylab("Observation-level correlation")+
  theme_bw() +
  coord_cartesian(clip = "off") + 
  theme(plot.margin = margin(5.5,5.5,14,5.5))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size=16),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13))
plot_1 <- plot_1 + geom_point(data=expected_cor[[1]], aes(x=climate, y=expected_correlation), size=2) 
plot_1



###########################
#### Figure posteriors ####
###########################

g_effects <- as.data.frame(post$B_mq_g)
f_effects <- as.data.frame(post$B_mq_f)
cor_effects <- as.data.frame(post$B_cpcq)

# change the "3000" to another value if you run chains for less or more iterations than I did (1000 iterations post burn-in * 3 = 3000)

df.posteriors <- data_frame(Submodel = c(rep("Growth", 3000*2), rep("Fecundity", 3000*2))
                            , parameter = c(rep("Intercept", 3000), rep("Climate", 3000)
                                            , rep("Intercept", 3000), rep("Climate", 3000))
                            , Posterior = c(g_effects$V1, g_effects$V2
                                            , f_effects$V1, f_effects$V2))



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

# pdf("non_repeated_crn_plot.pdf", width = 9, height = 6)
wrap_elements(full= (p | q))
# dev.off()


