#==========================================
#
# Tree condition model for BBNP
#
# Dave Ramsey 27/11/2020 
#==========================================

library(tidyverse)
library(jagsUI)
library(mcmcOutput)
library(grid)
library(gridExtra)

load("Data/BBNP_tree_data.Rdata")

#----------------------------------------
sink("tree-multinom.jags")
cat("
		model {
		
		# -------------------------------------------------
		# Parameters:
		# s: survival probability
		# psi1:psi5 transitions from defclass 1:5
		# -------------------------------------------------
		# States (S):
		# 1 high
		# 2 mod
		# 3 low
		# 4 vlow
		# 5 dead
		
		# -------------------------------------------------
		
		# Priors and constraints
		# Survival: uniform
		for(i in 1:nind) {
		for (t in first[i]:(last[i]-1)){ 
			logit(s[i,t]) <- alpha + beta*Rain[t] + gam*pfc[i,t] + zeta*log(tm[t])
			}
		}
		
		
		# Transitions:
		for (i in 1:4){
		for(t in 1:(n.occasions-1)) {
		lpsi1[i,t] <- eta1[i] + delta1[i]*Ko[t] + kappa1[i]*Rain[t]
		lpsi2[i,t] <- eta2[i] + delta2[i]*Ko[t] + kappa2[i]*Rain[t]
		lpsi3[i,t] <- eta3[i] + delta3[i]*Ko[t] + kappa3[i]*Rain[t]
		lpsi4[i,t] <- eta4[i] + delta4[i]*Ko[t] + kappa4[i]*Rain[t]
		
		psi1[i,t]<- epsi1[i,t]/sum(epsi1[,t])
		epsi1[i,t]<- exp(lpsi1[i,t])
		psi2[i,t]<- epsi2[i,t]/sum(epsi2[,t])
		epsi2[i,t]<- exp(lpsi2[i,t])
		psi3[i,t]<- epsi3[i,t]/sum(epsi3[,t])
		epsi3[i,t]<- exp(lpsi3[i,t])
		psi4[i,t]<- epsi4[i,t]/sum(epsi4[,t])
		epsi4[i,t]<- exp(lpsi4[i,t])
		}
		}
		
		alpha ~ dnorm(0, 0.1)
		beta ~ dnorm(0, 0.1)
		gam ~ dnorm(0, 0.1)
		zeta ~ dnorm(0, 0.1)

		# Constrain baseline transitions to 0: 
		# baseline is self transition 1->1; 2->2, 3->3,4->4
		
		eta1[1] <- 0
		for(i in 2:4) {eta1[i] ~ dnorm(0, 0.1)}
		eta2[2]<- 0
		for(i in c(1,3,4)) {eta2[i] ~ dnorm(0, 0.1)}
		eta3[3]<- 0
		for(i in c(1,2,4)) {eta3[i] ~ dnorm(0, 0.1)}
		eta4[4]<- 0
		for(i in 1:3) {eta4[i] ~ dnorm(0, 0.1)}
		
		
		delta1[1] <- 0
		for(i in 2:4) {delta1[i] ~ dnorm(0, 0.1)}
		delta2[2]<- 0
		for(i in c(1,3,4)) {delta2[i] ~ dnorm(0, 0.1)}
		delta3[3]<- 0
		for(i in c(1,2,4)) {delta3[i] ~ dnorm(0, 0.1)}
		delta4[4]<- 0
		for(i in 1:3) {delta4[i] ~ dnorm(0, 0.1)}
		
		kappa1[1] <- 0
		for(i in 2:4) {kappa1[i] ~ dnorm(0, 0.1)}
		kappa2[2]<- 0
		for(i in c(1,3,4)) {kappa2[i] ~ dnorm(0, 0.1)}
		kappa3[3]<- 0
		for(i in c(1,2,4)) {kappa3[i] ~ dnorm(0, 0.1)}
		kappa4[4]<- 0
		for(i in 1:3) {kappa4[i] ~ dnorm(0, 0.1)}
		
		# Define transition matrix 	
		for (i in 1:nind){
		for (t in first[i]:(last[i]-1)){
		ps[1,i,t,1] <- s[i,t] * psi1[1,t]
		ps[1,i,t,2] <- s[i,t] * psi1[2,t]
		ps[1,i,t,3] <- s[i,t] * psi1[3,t]
		ps[1,i,t,4] <- s[i,t] * psi1[4,t]
		ps[1,i,t,5] <- 1-s[i,t]
		ps[2,i,t,1] <- s[i,t] * psi2[1,t] 
		ps[2,i,t,2] <- s[i,t] * psi2[2,t]
		ps[2,i,t,3] <- s[i,t] * psi2[3,t]
		ps[2,i,t,4] <- s[i,t] * psi2[4,t]
		ps[2,i,t,5] <- 1-s[i,t] 
		ps[3,i,t,1] <- s[i,t] * psi3[1,t] 
		ps[3,i,t,2] <- s[i,t] * psi3[2,t]
		ps[3,i,t,3] <- s[i,t] * psi3[3,t]
		ps[3,i,t,4] <- s[i,t] * psi3[4,t]
		ps[3,i,t,5] <- 1-s[i,t] 
		ps[4,i,t,1] <- s[i,t] * psi4[1,t] 
		ps[4,i,t,2] <- s[i,t] * psi4[2,t]
		ps[4,i,t,3] <- s[i,t] * psi4[3,t]
		ps[4,i,t,4] <- s[i,t] * psi4[4,t]
		ps[4,i,t,5] <- 1-s[i,t] 
		ps[5,i,t,1] <- 0 
		ps[5,i,t,2] <- 0
		ps[5,i,t,3] <- 0
		ps[5,i,t,4] <- 0
		ps[5,i,t,5] <- 1
		} #t
		} #i
		
		# Likelihood 
		for (i in 1:nind){
		for (t in (first[i]+1):last[i]){
		# transition process
		y[i,t] ~ dcat(ps[y[i,t-1], i, t-1,])
		} #t
		} #i
		}
		",fill=TRUE)
sink()

# Initial values 
inits <- function(){list(alpha=2,beta=0,gam=0,eta1=c(NA,rep(0,3)),zeta=0,
												 eta2=c(0,NA,0,0),eta3=c(0,0,NA,0),eta4=c(0,0,0,NA),
												 delta1=c(NA,rep(0,3)),delta2=c(0,NA,0,0),delta3=c(0,0,NA,0),
												 delta4=c(0,0,0,NA),kappa1=c(NA,rep(0,3)),
												 kappa2=c(0,NA,0,0),kappa3=c(0,0,NA,0),kappa4=c(0,0,0,NA))}
												   
# Parameters monitored
parameters <- c("psi1", "psi2", "psi3", "psi4","alpha","beta","gam","zeta",
								"eta1","eta2","eta3","eta4","delta1","delta2","delta3","delta4",
								"kappa1","kappa2","kappa3","kappa4")


# MCMC settings - for testing.  Increase number of samples for inference
ni <- 2000
nt <- 1
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 56 min)
BBNP_trees<- jags(tree.data, inits, parameters, "tree-multinom.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE)

print(BBNP_trees, digits = 3)


#===========================================
# Tree condition
#===========================================
# Tree survival
#===========================================

pars<- mcmcOutput(BBNP_trees)

fc<- seq(0,1,0.01)
n<- length(fc)
nsims<- 1000
mat<- matrix(NA,nrow=nsims,ncol=n)
for(i in 1:nsims){
for(j in 1:n) {
	mat[i,j]<- plogis(pars$alpha[i] + pars$gam[i]*fc[j] + pars$zeta[i]*log(1))
	}
}

xbar<- apply(mat,2,mean)
lcl<- apply(mat,2, quantile, 0.05)
ucl<- apply(mat,2, quantile, 0.95)


df<- as.data.frame(xbar) %>% mutate(PFC = fc, lcl=lcl, ucl=ucl)

# Survival rate over time
fc<- 0.1
tt<- 1:11
n<- length(tt)
nsims<- 1000
mat<- matrix(NA,nrow=nsims,ncol=n)

for(i in 1:nsims) {
	for(t in tt) {
	mat[i,t]<- plogis(pars$alpha[i] + pars$gam[i]*fc + pars$zeta[i]*log(t))
	}
}

xbar<- apply(mat,2,mean)
lcl<- apply(mat,2, quantile, 0.025)
ucl<- apply(mat,2, quantile, 0.975)

df2<- as.data.frame(xbar) %>% mutate(Year = tt, lcl=lcl, ucl=ucl)


p1<- df %>% ggplot(aes(PFC,1-xbar)) +
	geom_line() +
	geom_ribbon(aes(ymin=1-lcl, ymax=1-ucl), fill="grey80", alpha=0.5) +
	scale_y_continuous(limits = c(0, 1),breaks = seq(0,1,0.2)) +
	labs(x="Projected Foliage Cover",y="Tree mortality rate", subtitle="A") +
	theme_bw() +
	theme(axis.title.x = element_text(face="bold", size=20),
				axis.title.y = element_text(face="bold", size=20),
				axis.text.x = element_text(size=15),
				axis.text.y = element_text(size=15),
				plot.subtitle = element_text(face="bold", size=15))

p2<- df2 %>% ggplot(aes(Year,1-xbar)) +
	geom_line() +
	geom_ribbon(aes(ymin=1-lcl, ymax=1-ucl), fill="grey80", alpha=0.5) +
	scale_y_continuous(limits = c(0, 1),breaks = seq(0,1,0.2)) +
	scale_x_continuous(limits = c(1, 11),breaks = seq(1,11,2), labels=seq(2004, 2014, 2)) +
	labs(x="Year",y="Tree mortality rate",subtitle="B") +
	theme_bw() +
	theme(axis.title.x = element_text(face="bold", size=20),
				axis.title.y = element_text(face="bold", size=20),
				axis.text.x = element_text(size=15),
				axis.text.y = element_text(size=15),
				plot.subtitle = element_text(face="bold", size=15))


win.graph(15,8)
grid.arrange(p1, p2, ncol=2,top=textGrob("Budj Bim NP",gp=gpar(fontsize=20,face="bold")))

#===========================================
# Tree condition
#===========================================
# Covariates - koala density
#===========================================

pars<- mcmcOutput(BBNP_trees)
nsims<- 1000
kdens<- seq(0,2,0.1)

n<- length(kdens)
lmat<- array(NA,c(nsims,n,4))
mat<- lmat

for(k in 1:nsims) {
  for(j in 1:4) {
    for(i in 1:n) {
      lmat[k,i,j]<- pars$eta1[k,j] + pars$delta1[k,j]*kdens[i] 
    }
  }
}
for(k in 1:nsims) {
  for(i in 1:n) {
    mat[k,i,]<- exp(lmat[k,i,])/sum(exp(lmat[k,i,]))
  }
}


xbar<- apply(mat, c(2,3), mean)
lcl<- apply(mat, c(2,3), quantile, 0.025)
ucl<- apply(mat, c(2,3), quantile, 0.975)

xbar<- as.data.frame(xbar) %>% mutate(Density = kdens)
lcl<- as.data.frame(lcl) %>% mutate(Density = kdens)
ucl<- as.data.frame(ucl) %>% mutate(Density = kdens)

xbar<- xbar %>% pivot_longer(-Density, names_to = "strata", values_to ="mean")
lcl<- lcl %>% pivot_longer(-Density, names_to = "strata", values_to ="lcl")
ucl<- ucl %>% pivot_longer(-Density, names_to = "strata", values_to ="ucl")

df<- full_join(xbar, lcl, by = c("Density", "strata"))
df<- full_join(df, ucl, by = c("Density", "strata"))
df<- df %>% mutate(strata = factor(strata, labels=c("to vlow","to low","to mod","to high")))

#--------
lmat<- array(NA,c(nsims,n,4))
mat<- lmat

for(k in 1:nsims) {
  for(j in 1:4) {
    for(i in 1:n) {
      lmat[k,i,j]<- pars$eta2[k,j] + pars$delta2[k,j]*kdens[i]
    }
  }
}
for(k in 1:nsims) {
  for(i in 1:n) {
    mat[k,i,]<- exp(lmat[k,i,])/sum(exp(lmat[k,i,]))
  }
}

xbar<- apply(mat, c(2,3), mean)
lcl<- apply(mat, c(2,3), quantile, 0.025)
ucl<- apply(mat, c(2,3), quantile, 0.975)

xbar<- as.data.frame(xbar) %>% mutate(Density = kdens)
lcl<- as.data.frame(lcl) %>% mutate(Density = kdens)
ucl<- as.data.frame(ucl) %>% mutate(Density = kdens)

xbar<- xbar %>% pivot_longer(-Density, names_to = "strata", values_to ="mean")
lcl<- lcl %>% pivot_longer(-Density, names_to = "strata", values_to ="lcl")
ucl<- ucl %>% pivot_longer(-Density, names_to = "strata", values_to ="ucl")

df2<- full_join(xbar, lcl, by = c("Density", "strata"))
df2<- full_join(df2, ucl, by = c("Density", "strata"))
df2<- df2 %>% mutate(strata = factor(strata, labels=c("to vlow","to low","to mod","to high")))

p1<- df %>% ggplot(aes(Density, mean, color=strata, fill=strata)) +
  geom_ribbon(aes(ymin=lcl, ymax=ucl), fill="grey70", alpha=0.5, show.legend = F, linetype=0) +
  geom_line(size=1) +
  scale_color_manual(values=c("green","blue","orange","red")) +
  scale_y_continuous(limits = c(0, 1),breaks = seq(0,1,0.2)) +
  labs(x = "Koala density (koalas/ha)", y = "Transition rate", color = "from vlow defoliation", subtitle="A") +
  theme_bw() +
  theme(legend.position=c(0.5,0.8),
        legend.box.just = "center",
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        plot.subtitle = element_text(face="bold", size=15))

p2<- df2 %>% ggplot(aes(Density, mean, color=strata, fill=strata)) +
  geom_ribbon(aes(ymin=lcl, ymax=ucl), fill="grey70", alpha=0.5, show.legend = F, linetype=0) +
  geom_line(size=1) +
  scale_color_manual(values=c("green","blue","orange","red")) +
  scale_y_continuous(limits = c(0, 1),breaks = seq(0,1,0.2)) +
  labs(x = "Koala density (koalas/ha)", y = "Transition rate", color = "from low defoliation",subtitle="B") +
  theme_bw() +
  theme(legend.position=c(0.5,0.8),
        legend.box.just = "center",
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        plot.subtitle = element_text(face="bold", size=15))

win.graph(15,8)
grid.arrange(p1, p2, ncol=2,top=textGrob("Budj Bim NP",gp=gpar(fontsize=20,face="bold")))
