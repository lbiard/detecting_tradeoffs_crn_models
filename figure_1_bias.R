###################################
#### clean working environment ####
###################################
rm(list = ls())

######################
#### load package ####
######################
library(ggplot2)
library(patchwork)

###########################
#### low repeatability ####
###########################

ra <- seq(-1,1,by=0.001) # among-individual correlation
re <- seq(-1,1,by=0.001) # within-individual correlation
r1 <- 0.2 # repeatability trait 1
r2 <- 0.2 # repeatability trait 2
dat.plot <- expand.grid(ra,re,r1,r2)
colnames(dat.plot) <- c("ra","re","r1","r2")

# observation-level correlation for each combination of trait repeatability, and among and within individual correlation
observation_cor <- function(ra,re,r1,r2){
            ra*sqrt(r1*r2) + re*sqrt((1-r1)*(1-r2)) 
}


dat.plot$outcome <- observation_cor(ra=dat.plot$ra,
                re=dat.plot$re,
                r1=dat.plot$r1,
                r2=dat.plot$r2)
dat.plot$bias <- dat.plot$outcome-dat.plot$ra
dat.plot$too_much_bias <- as.factor(as.numeric(sign(dat.plot$outcome) != sign(dat.plot$ra)))

# plot
p.low <- ggplot(dat.plot,
       aes(x = ra,
           y = re, 
           fill = bias)) +
  ggtitle("Repeatability = 0.2")+
  xlab("Among-individual correlation")+
  ylab("Within-individual correlation")+
  labs(fill = "Bias")+
  geom_raster()+
  scale_fill_gradient2(low = "cornflowerblue", mid = "white", high = "orange",
                       midpoint = 0,
                       limits = c(-1.6, 1.6),
                       breaks = c(-1,-0.5,-0.25,0,0.25,0.5,1))+
  scale_y_continuous(expand = c(0, 0) )+
  scale_x_continuous(expand = c(0, 0) )+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.key.height= unit(1.3, 'cm'),
        legend.key.width= unit(0.8, 'cm'))+
  annotate(geom="raster", 
           x=dat.plot$ra[dat.plot$too_much_bias==1], 
           y=dat.plot$re[dat.plot$too_much_bias==1], alpha=.1)



###########################
#### med repeatability ####
###########################

ra <- seq(-1,1,by=0.001)
re <- seq(-1,1,by=0.001)
r1 <- 0.5
r2 <- 0.5
dat.plot <- expand.grid(ra,re,r1,r2)
colnames(dat.plot) <- c("ra","re","r1","r2")

observation_cor <- function(ra,re,r1,r2){
  ra*sqrt(r1*r2) + re*sqrt((1-r1)*(1-r2)) 
}


dat.plot$outcome <- observation_cor(ra=dat.plot$ra,
                                    re=dat.plot$re,
                                    r1=dat.plot$r1,
                                    r2=dat.plot$r2)
dat.plot$bias <- dat.plot$outcome-dat.plot$ra
dat.plot$too_much_bias <- as.factor(as.numeric(sign(dat.plot$outcome) != sign(dat.plot$ra)))



p.med <- ggplot(dat.plot,
                aes(x = ra,
                    y = re, 
                    fill = bias)) +
  ggtitle("Repeatability = 0.5")+
  xlab("Among-individual correlation")+
  ylab("Within-individual correlation")+
  labs(fill = "Bias")+
  geom_raster()+
  scale_fill_gradient2(low = "cornflowerblue", mid = "white", high = "orange",
                       midpoint = 0,
                       limits = c(-1.6, 1.6),
                       breaks = c(-1,-0.5,-0.25,0,0.25,0.5,1))+
  scale_y_continuous(expand = c(0, 0) )+
  scale_x_continuous(expand = c(0, 0) )+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.position = "none")+
  annotate(geom="raster", 
           x=dat.plot$ra[dat.plot$too_much_bias==1], 
           y=dat.plot$re[dat.plot$too_much_bias==1], alpha=.1)



############################
#### high repeatability ####
############################

ra <- seq(-1,1,by=0.001)
re <- seq(-1,1,by=0.001)
r1 <- 0.8
r2 <- 0.8
dat.plot <- expand.grid(ra,re,r1,r2)
colnames(dat.plot) <- c("ra","re","r1","r2")

observation_cor <- function(ra,re,r1,r2){
  ra*sqrt(r1*r2) + re*sqrt((1-r1)*(1-r2)) 
}


dat.plot$outcome <- observation_cor(ra=dat.plot$ra,
                                    re=dat.plot$re,
                                    r1=dat.plot$r1,
                                    r2=dat.plot$r2)
dat.plot$bias <- dat.plot$outcome-dat.plot$ra
dat.plot$too_much_bias <- as.factor(as.numeric(sign(dat.plot$outcome) != sign(dat.plot$ra)))


p.high <- ggplot(dat.plot,
                aes(x = ra,
                    y = re, 
                    fill = bias)) +
  ggtitle("Repeatability = 0.8")+
  xlab("Among-individual correlation")+
  ylab("Within-individual correlation")+
  labs(fill = "Bias")+
  geom_raster()+
  scale_fill_gradient2(low = "cornflowerblue", mid = "white", high = "orange",
                       midpoint = 0,
                       limits = c(-1.6, 1.6),
                       breaks = c(-1,-0.5,-0.25,0,0.25,0.5,1))+
  scale_y_continuous(expand = c(0, 0) )+
  scale_x_continuous(expand = c(0, 0) )+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.position = "none")+
  annotate(geom="raster", 
           x=dat.plot$ra[dat.plot$too_much_bias==1], 
           y=dat.plot$re[dat.plot$too_much_bias==1], alpha=.1)

  
##############
#### plot ####
##############

p.high + p.med + p.low
