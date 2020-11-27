#==========================================
#
# Koala IPM model for BBNP
#
# Dave Ramsey 27/11/2020 
#==========================================

library(tidyverse)
library(jagsUI)
library(mcmcOutput)
library(grid)
library(gridExtra)

load("Data/BBNP_IPM_data.Rdata")

#--------------------------------------------------------------
# Integrated population model JAGS code
#
sink("koala_IPM.jags")
cat("
model {
# Priors and constraints
 
  for (t in 1:(n.occasions-1)){
    logit(phi.juv[t]) <- mu.phijuv + eps[t]
    logit(phi.ads[t]) <- mu.phiads + eps[t]
		logit(phi.adf[t]) <- mu.phiadf + eps[t]
    logit(p[t]) <- alpha[t]
    eps[t] ~ dnorm(0, tau)
  }
    
     mu.phijuv<- log(mean.phijuv/(1-mean.phijuv))
     mu.phiads<- log(mean.phiads/(1-mean.phiads))
		 mu.phiadf<- log(mean.phiadf/(1-mean.phiadf))
     mean.phijuv ~ dunif(0, 1)          # Prior for mean juv. survival
     mean.phiads ~ dunif(0, 1)           # Prior for mean ad. survival
		 mean.phiadf ~ dunif(0, 1)           # Prior for mean ad. survival
     tau<- pow(sigma, -2)
     sigma ~ dunif(0,20)
    
     phi.py ~ dbeta(71, 29) # prior for pouch young survival to 1yo.
     sr ~ dbeta(100,100) # prior for sex ratio at birth
     
  for(j in 1:(n.occasions-1)) {
    alpha[j] ~ dnorm(mu.p, tau.p)
    }
    mu.p ~ dnorm(0, 0.1)
    tau.p<- pow(sigma.p, -2)
    sigma.p ~ dunif(0,20)
    
  for(j in 1:n.occasions) {
    Br[j] ~ dunif(0,1)      # prior for breeding rate
    rhoA[j] ~ dunif(0,1)  # prior for sterilisation rate
    eta[j] ~ dunif(0,1) # prior for reversion rate (sterile to fertile)
  }
    
    IF2[1]<- 0
    IM2[1]<- 0
    muIF[1]<- 0
    muIM[1]<- 0
    
  for(j in 2:n.occasions) {
    muIF[j] ~ dgamma(1, 0.05) # no. of female immigrants
    muIM[j] ~ dgamma(1, 0.05) # no. of male immigrants
#		muIF[j] ~ dunif(0, 100) # no. of female immigrants
#   muIM[j] ~ dunif(0, 100) # no. of male immigrants
  }
	
		f3 ~ dpois(Sinit)
    f2  ~ dpois(Ninit * 0.42)
		f1 ~ dpois(Ninit * 0.08)
		
		m2  ~ dpois(Ninit * 0.42)
		m1  ~ dpois(Ninit * 0.08)
		
		F2[1] <- f2
		M2[1] <- m2
		F3[1]<- f3
		M1[1]<- m1
		F1[1]<- f1
    #-----------------------------------
    # Likelihoods
    # Abundance data model

  for (i in 1:n) {
       for(j in 1:3) {
          mu[i,j] <- Ntran[year[i]]*psight[year[i],j]*3
          y[i,j] ~ dpois(mu[i,j])
        }
    }
    
    for(j in 1:n.occasions) {
       psight[j,1] <- pA[j]*pB[j]
       psight[j,2] <- pA[j]*(1-pB[j])
       psight[j,3] <- (1-pA[j])*pB[j]
       pA[j] ~ dbeta(1,1)
       pB[j] ~ dbeta(1,1)
    }
   
   for(i in 1:n.occasions) {
        Ntran[i]<- Ntot[i] / 100  # Convert Ntot to koalas/km2
        Ntot[i]<- F1[i] + F2[i] + F3[i] + M1[i] + M2[i]  # total pop
    }
    
    #----------------------------------
    # Breeding & sterilistaion rates
    
    for(i in 1:n.occasions) {
      ny2[i] ~ dbin(Br[i], R2[i]) # Overall breeding rate (of non-sterile females)
      nre[i] ~ dbin(eta[i], NRE[i])  # Probability of sterile reversion to fertile
    }

    for(i in 2:n.occasions) {
      nsA[i] ~ dbin(phi.adf[i-1]*rhoA[i],F2[i-1]) # adult sterilisation rate for the park
    }
      nsA[1] ~ dbin(rhoA[1], F2[1]+Sinit)
   
    #-------------------------------------
      # populatioon growth rate
        for(j in 1:(n.occasions-1)) {
          grate[j]<- Ntot[j+1]/Ntot[j]
          lgrate[j]<- log(grate[j])
          }
        mgrate<- exp((1/(n.occasions-1))*sum(lgrate))
    #-------------------------------------
    # Koala age structured matrix model
    
  for(i in 2:n.occasions) {
    # Females
    F1[i] ~ dbin(Br[i-1]*sr*phi.adf[i-1]*phi.py, F2[i-1])
    F21[i-1] ~ dbin(phi.juv[i-1], F1[i-1])
    F22[i-1] ~ dbin(phi.adf[i-1]*(1-rhoA[i]), F2[i-1])
    F23[i-1] ~ dbin(phi.ads[i-1]*eta[i], F3[i-1])
    F2[i] <- F21[i-1] + F22[i-1] + F23[i-1] + IF2[i]
    
    F32[i-1] ~ dbin(phi.adf[i-1]*rhoA[i], F2[i-1])  
    F33[i-1] ~ dbin(phi.ads[i-1]*(1-eta[i]), F3[i-1])
    F3[i]<- F32[i-1] + F33[i-1]
    # Males
    M1[i]  ~ dbin(Br[i-1]*(1-sr)*phi.adf[i-1]*phi.py, F2[i-1])
    M21[i-1] ~ dbin(phi.juv[i-1], M1[i-1])
    M22[i-1] ~ dbin(phi.adf[i-1], M2[i-1])
    M2[i]<- M21[i-1] + M22[i-1] + IM2[i]
    # immigrants
    IF2[i] ~ dpois(muIF[i])
    IM2[i] ~ dpois(muIM[i])
  }
    #-------------------------
    # Koala mark recapture model (CJS) 
    # Define the multinomial likelihood
  for (t in 1:(n.occasions-1)){
    marr.js[t,1:n.occasions] ~ dmulti(pr.js[t,], rel.js[t])
		marr.jf[t,1:n.occasions] ~ dmulti(pr.jf[t,], rel.jf[t])
    marr.as[t,1:n.occasions] ~ dmulti(pr.as[t,], rel.as[t])
		marr.af[t,1:n.occasions] ~ dmulti(pr.af[t,], rel.af[t])
  }
    # Define the cell probabilities of the m-arrays
    # Main diagonal
  for (t in 1:(n.occasions-1)){
    q[t] <- 1-p[t]            # Probability of non-recapture
    pr.js[t,t] <- phi.juv[t]*p[t]
		pr.jf[t,t] <- phi.juv[t]*p[t]
    pr.as[t,t] <- phi.ads[t]*p[t]
		pr.af[t,t] <- phi.adf[t]*p[t]
    # Above main diagonal
  for (j in (t+1):(n.occasions-1)){
    pr.js[t,j] <- phi.juv[t]*prod(phi.ads[(t+1):j])*prod(q[t:(j-1)])*p[j]
	  pr.jf[t,j] <- phi.juv[t]*prod(phi.adf[(t+1):j])*prod(q[t:(j-1)])*p[j]
    pr.as[t,j] <- prod(phi.ads[t:j])*prod(q[t:(j-1)])*p[j]
		pr.af[t,j] <- prod(phi.adf[t:j])*prod(q[t:(j-1)])*p[j]
  } #j
    # Below main diagonal
  for (j in 1:(t-1)){
    pr.js[t,j] <- 0
	  pr.jf[t,j] <- 0
    pr.as[t,j] <- 0
		pr.af[t,j] <- 0
    } #j
  } #t
    # Last column: probability of non-recapture
  for (t in 1:(n.occasions-1)){
    pr.js[t,n.occasions] <- 1-sum(pr.js[t,1:(n.occasions-1)])
		pr.jf[t,n.occasions] <- 1-sum(pr.jf[t,1:(n.occasions-1)])
    pr.as[t,n.occasions] <- 1-sum(pr.as[t,1:(n.occasions-1)])
	  pr.af[t,n.occasions] <- 1-sum(pr.af[t,1:(n.occasions-1)])
  } #t
    
}
    ",fill = TRUE)
sink()


n.occasions<- ipm.data$n.occasions
iN1<- 20
iN2<- 100

inits <- function(){list(mean.phijuv=runif(1), mean.phiads=runif(1),mean.phiadf=runif(1), sigma=1, mu.p= -2, sigma.p=1, Br=rep(0.3,n.occasions),eta=rep(0.03,n.occasions),rhoA=rep(0.5,n.occasions), f1=iN1,f2=iN2,m1=10,m2=iN2,f3=iN1,muIF=c(NA,runif(9,0,100)),muIM=c(NA,runif(9,0,100)),pA=runif(n.occasions), 
                             pB=runif(n.occasions))}  

# Parameters monitored
parameters<- c("mean.phijuv","mean.phiads","mean.phiadf","phi.ads","phi.adf","phi.juv","p","Ntot","F1","F2","F3","M1","M2","IF2","IM2","rhoA","sigma","eta","Br","F32","mgrate")

# MCMC settings - for testing.  Increase number of samples for inference
ni <- 2000
nt <- 1
nb <- 1000
nc <- 5

BBNP_ipm <- jags(ipm.data, inits, parameters, "koala_IPM.jags", n.chains = nc, 
             n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE)

print(BBNP_ipm, digits = 3)

#-------------------------------------------------
pars<- mcmcOutput(BBNP_ipm)

Year<- 2004:2013

N<- apply(pars$Ntot,2,mean)
lcl<- apply(pars$Ntot,2,quantile, 0.025)
ucl<- apply(pars$Ntot,2,quantile, 0.975)
df<- tibble(Year=Year, Stage="Total population", N=N, lcl=lcl, ucl=ucl)

N<- apply(pars$IF2,2,mean)
lcl<- apply(pars$IF2,2,quantile, 0.025)
ucl<- apply(pars$IF2,2,quantile, 0.975)
df<- bind_rows(df, tibble(Year=Year, Stage="Female immigrants", N=N, lcl=lcl, ucl=ucl))

N<- apply(pars$IM2,2,mean)
lcl<- apply(pars$IM2,2,quantile, 0.025)
ucl<- apply(pars$IM2,2,quantile, 0.975)
df<- bind_rows(df, tibble(Year=Year, Stage="Male immigrants", N=N, lcl=lcl, ucl=ucl))

tmp<- pars$M1 + pars$F1
N<- apply(tmp,2,mean)
lcl<- apply(tmp,2,quantile, 0.025)
ucl<- apply(tmp,2,quantile, 0.975)
df<- bind_rows(df, tibble(Year=Year, Stage="Immature", N=N, lcl=lcl, ucl=ucl))

N<- apply(pars$F2,2,mean)
lcl<- apply(pars$F2,2,quantile, 0.025)
ucl<- apply(pars$F2,2,quantile, 0.975)
df<- bind_rows(df, tibble(Year=Year, Stage="Non-Sterile females", N=N, lcl=lcl, ucl=ucl))

N<- apply(pars$F3,2,mean)
lcl<- apply(pars$F3,2,quantile, 0.025)
ucl<- apply(pars$F3,2,quantile, 0.975)
df<- bind_rows(df, tibble(Year=Year, Stage="Sterile females", N=N, lcl=lcl, ucl=ucl))

N<- apply(pars$M2,2,mean)
lcl<- apply(pars$M2,2,quantile, 0.025)
ucl<- apply(pars$M2,2,quantile, 0.975)
df<- bind_rows(df, tibble(Year=Year, Stage="Mature males", N=N, lcl=lcl, ucl=ucl))

#-------------------------
pd <- position_dodge(0.3)
# The palette with grey:
cbPalette <- c("#999999", "#CC6666", "#9999CC", "#66CC99")

p1<- df %>% filter(Stage %in% c("Total population","Female immigrants","Male immigrants")) %>% 
  ggplot(aes(Year, N)) +
  geom_line(aes(color=Stage), size=1.2, position = pd) +
  geom_linerange(aes(Year, ymin = lcl, ymax=ucl, color=Stage), size=1, position=pd) +
  geom_point(aes(fill=Stage), shape=21, color="black", size=3, position=pd) +
  scale_x_continuous(limits = c(2003.8, 2013.2),breaks = seq(2004,2013, 1)) +
  scale_y_continuous(limits = c(0, 200),breaks = seq(0,200, 50)) +
  labs(y=expression(paste("Density (koalas/",km^2,")")), x="Year", subtitle="A") +
  theme_bw() +
  theme(legend.position=c(0.8,0.85),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        plot.subtitle = element_text(face="bold", size=15))


pd <- position_dodge(0.3)

p2<- df %>% filter(!(Stage %in% c("Total population","Female immigrants","Male immigrants"))) %>% 
  ggplot(aes(Year, N)) +
  geom_line(aes(color=Stage), size=1.2, position = pd) +
  geom_linerange(aes(Year, ymin = lcl, ymax=ucl, color=Stage), size=1, position=pd) +
  geom_point(aes(fill=Stage), shape=21, color="black", size=3, position=pd) +
  scale_color_manual(values = cbPalette, aesthetics = c("colour", "fill")) +
  scale_x_continuous(limits = c(2003.8, 2013.2),breaks = seq(2004,2013, 1)) +
  scale_y_continuous(limits = c(0, 100),breaks = seq(0,100, 20)) +
  labs(y=expression(paste("Density (koalas/",km^2,")")), x="Year", subtitle = "B") +
  theme_bw() +
  theme(legend.position=c(0.8,0.85),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        plot.subtitle = element_text(face="bold", size=15))



win.graph(15,8)
grid.arrange(p1, p2, ncol=2,top=textGrob("Budj Bim NP",gp=gpar(fontsize=20,face="bold")))
