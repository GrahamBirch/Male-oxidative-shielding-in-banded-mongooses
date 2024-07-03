library(rjags)
library(jagsUI)
library(tidyverse)
library(ggmcmc)
library(stringr)
library(MCMCvis)

#Key
#betaRS= Reproductive state of male
#betaP = Mating period
#Bcells= cells present (1), no sample(2), cells absent(3)
#intaps= interaction between reproductive state and mating period
#epsG = Group ID
#epsM = Male ID
#epsO = oestrus event ID

#Model oupt and random effect variance
Oms.OS.justoes.basecsum.olap
variancejustOS

Oms.OS.baseintcsum.olap
varianceOS

#See bottom of script for diagnostic plots

spermdamage2.ageweightc<-read_csv("ej.mda.data.csv")
jags.data.OS<- list( y = (
  sqrt(spermdamage2.ageweightc$meanMDA)),n = length(spermdamage2.ageweightc$meanMDA),
  cells=as.numeric(as.factor(spermdamage2.ageweightc$cells.est)),
  rs=(spermdamage2.ageweightc$rsn),
  p=(spermdamage2.ageweightc$pn),
  code = as.numeric(as.factor(spermdamage2.ageweightc$tempcode)),
  group = as.numeric(as.factor(spermdamage2.ageweightc$pack)),
  male = as.numeric(as.factor(spermdamage2.ageweightc$male)),
  malel= length(as.numeric(as.factor((spermdamage2.ageweightc$male)))),
  codel= length(levels(as.factor((spermdamage2.ageweightc$tempcode)))),
  groupl= length(as.numeric(as.factor((spermdamage2.ageweightc$pack)))))

jags.data.OS.justoes<- list( y = (
  sqrt(spermdamage2.ageweightc.justoes$meanMDA)),n = length(spermdamage2.ageweightc.justoes$meanMDA),
  cells=as.numeric(as.factor(spermdamage2.ageweightc.justoes$cells.est)),
  rs=(spermdamage2.ageweightc.justoes$rsn),
  p=(spermdamage2.ageweightc.justoes$pn),
  code = as.numeric(as.factor(spermdamage2.ageweightc.justoes$tempcode)),
  group = as.numeric(as.factor(spermdamage2.ageweightc.justoes$pack)),
  male = as.numeric(as.factor(spermdamage2.ageweightc.justoes$male)),
  malel= length(as.numeric(as.factor((spermdamage2.ageweightc.justoes$male)))),
  codel= length(levels(as.factor((spermdamage2.ageweightc.justoes$tempcode)))),
  groupl= length(as.numeric(as.factor((spermdamage2.ageweightc.justoes$pack)))))



nc<-3
nt<-1000
ni<-100000
nb<-10000


sink("weightOS.base.justoes.modc")
cat(" 
  model{ 
    # likelihood
    for (i in 1:n){
            y[i] ~ dnorm(mu[i], tau)
            mu[i] <- alpha 
            +betaRS*rs[i]
            +Bcells[cells[i]]
           +epsG[group[i]]
            +epsO[code[i]]
            +epsM[male[i]]
                        #+bState*staten[i]
    }
    
Bcells[3] <- 0       # cells absent as the reference category
for(i in 1:2) {
  Bcells[i] ~ dunif(-5,5) # Difference between cells present/no cells and reference category
}
    # priors
alpha ~ dnorm(0, .001)
bState~ dnorm(0, .001)
betaRS~ dnorm(0, .001)

	sigma ~ dunif(0, 100) # standard deviation
	tau <- 1 / (sigma * sigma) # sigma^2 doesn't work in JAGS

  #bState[1] <- 0       # First level is reference category
  #for(i in 2:3) {
   # bState[i] ~ dunif(-5,5) # Difference between each player and the reference player
  #}

tau.oc <- 1 / (sd.oc*sd.oc)
sd.oc ~ dunif(0, 1)
tau.g <- 1 / (sd.g*sd.g)
sd.g ~ dunif(0, 1)
tau.m <- 1 / (sd.m*sd.m)
sd.m ~ dunif(0, 1)

for(i in 1:codel){
epsO[i]~ dnorm(0,tau.oc)}
for(i in 1:groupl){
epsG[i]~ dnorm(0,tau.g)}
for(i in 1:malel){
epsM[i]~ dnorm(0,tau.m)}}")
sink()


sink("weightOS.baseint.modc")
cat(" 
  model{ 
    # likelihood
    for (i in 1:n){
            y[i] ~ dnorm(mu[i], tau)
            mu[i] <- alpha 
            +betaRS*rs[i]
            +betaP*p[i]
            +Bcells[cells[i]]
            +intaps*p[i]*rs[i]
           +epsG[group[i]]
            +epsO[code[i]]
            +epsM[male[i]]
    }
    # priors
alpha ~ dnorm(0, .001)
bState~ dnorm(0, .001)
betaRS~ dnorm(0, .001)
betaP~ dnorm(0, .001)
intaps~ dnorm(0, .001)

	sigma ~ dunif(0, 100) # standard deviation
	tau <- 1 / (sigma * sigma) # sigma^2 doesn't work in JAGS

    
Bcells[3] <- 0       # cells absent as the reference category
for(i in 1:2) {
  Bcells[i] ~ dunif(-5,5) # Difference between cells present/no cells and reference category
}

tau.oc <- 1 / (sd.oc*sd.oc)
sd.oc ~ dunif(0, 1)
tau.g <- 1 / (sd.g*sd.g)
sd.g ~ dunif(0, 1)
tau.m <- 1 / (sd.m*sd.m)
sd.m ~ dunif(0, 1)

for(i in 1:codel){
epsO[i]~ dnorm(0,tau.oc)}
for(i in 1:groupl){
epsG[i]~ dnorm(0,tau.g)}
for(i in 1:malel){
epsM[i]~ dnorm(0,tau.m)}}")
sink()



ms.OS.justoes.basec<-jagsUI::jags(jags.data.OS.justoes, inits.OS, parameters.OS, "weightOS.base.justoes.modc", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=TRUE,DIC=T)#Beh?



ms.OS.baseintc<-jagsUI::jags(jags.data.OS, inits.OS, parameters.OS, "weightOS.baseint.modc", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=TRUE,DIC=T)#Beh?





ms.OS.justoes.basecsum<-MCMCsummary(ms.OS.justoes.basec, round = 2)%>%rename(low = "2.5%", high = "97.5%")
ms.OS.justoes.basecsum.olap<-ms.OS.justoes.basecsum%>%mutate(Effect = case_when(low > 0 & high> 0~ "+",
                                                                      low < 0 & high< 0~"-",
                                                                      TRUE ~"NotSignif"))



write.csv(ms.OS.justoes.basecsum.olap, "Oms.OS.justoes.basecsum.olap.csv")#to make ghost name column a real column
Oms.OS.justoes.basecsum.olap<-read.csv("Oms.OS.justoes.basecsum.olap.csv")%>%filter( !grepl('eps', X))%>% mutate(X = gsub("\\[|\\]","",X))%>%#removing square brackets so can join to f output
  left_join(cbind(names(unlist(ms.OS.justoes.basec$f)),as.data.frame(matrix(unlist(ms.OS.justoes.basec$f ))))%>%rename("X"="names(unlist(ms.OS.justoes.basec$f))", "f"="V1")%>%filter( !grepl('eps', X)),
            by=c("X"))%>%
  rename(Parameter=X)#mcmc doesn't give f output so need to rejoin it



ms.OS.baseintcsum<-MCMCsummary(ms.OS.baseintc, round = 2)%>%rename(low = "2.5%", high = "97.5%")
ms.OS.baseintcsum.olap<-ms.OS.baseintcsum%>%mutate(Effect = case_when(low > 0 & high> 0~ "+",
                                                                      low < 0 & high< 0~"-",
                                                                      TRUE ~"NotSignif"))



write.csv(ms.OS.baseintcsum.olap, "Oms.OS.baseintcsum.olap.csv")#to make ghost name column a real column
Oms.OS.baseintcsum.olap<-read.csv("Oms.OS.baseintcsum.olap.csv")%>%filter( !grepl('eps', X))%>% mutate(X = gsub("\\[|\\]","",X))%>%#removing square brackets so can join to f output
  left_join(cbind(names(unlist(ms.OS.baseintc$f)),as.data.frame(matrix(unlist(ms.OS.baseintc$f ))))%>%rename("X"="names(unlist(ms.OS.baseintc$f))", "f"="V1")%>%filter( !grepl('eps', X)),
            by=c("X"))%>%
  rename(Parameter=X)#mcmc doesn't give f output so need to rejoin it




variancejustOS<-ggs(ms.OS.justoes.basec$samples)%>%filter(grepl("eps", Parameter))%>%group_by(Parameter)%>%mutate(meanrand=mean(value))%>%ungroup()%>%select(-c("value", "Iteration", "Chain"))%>%distinct()%>%
  mutate(Paramgroup = case_when(grepl("epsM", Parameter)~"Male.id",grepl("epsG", Parameter)~"Group",grepl("epsO", Parameter)~"Event"))%>%group_by(Paramgroup)%>%mutate(vargroup=sd(meanrand))%>%ungroup()%>%
  select(Paramgroup,vargroup)%>%distinct()

varianceOS<-ggs(ms.OS.baseintc$samples)%>%filter(grepl("eps", Parameter))%>%group_by(Parameter)%>%mutate(meanrand=mean(value))%>%ungroup()%>%select(-c("value", "Iteration", "Chain"))%>%distinct()%>%
  mutate(Paramgroup = case_when(grepl("epsM", Parameter)~"Male.id",grepl("epsG", Parameter)~"Group",grepl("epsO", Parameter)~"Event"))%>%group_by(Paramgroup)%>%mutate(vargroup=sd(meanrand))%>%ungroup()%>%
  select(Paramgroup,vargroup)%>%distinct()



#betaRS= Reproductive state of male
#betaP = Mating period
#Bcells= cells present (1), no sample(2), cells absent(3)
#intaps= interaction between reproductive state and mating period
#epsG = Group ID
#epsM = Male ID
#epsO = oestrus event ID

ggms.OS.baseintc<-ggs(ms.OS.baseintc$samples)%>%filter(!grepl("deviance", Parameter))%>%filter(!grepl("eps", Parameter),!grepl("alpha", Parameter))%>%mutate(Parameter=case_when(
  #Parameter=="betac" ~ "GroupSexRatio",
  Parameter=="betaRS"~"Reproductive state of male",
  Parameter=="betaP"~"Mating period",
  Parameter=="Bcells[1]"~"No cells",
  Parameter=="Bcells[2]"~"Cells present",
  Parameter=="intaps"~"State:Period"))%>%arrange(desc(Parameter))%>%na.omit()
unique(ggms.OS.baseintc$Parameter)
orderassort <- c( "Reproductive state of male",
                  "Mating period",
                  "State:Period",
                  "Cells present",
                  "No cells")
ggms.OS.baseintc$Parameter<-factor(ggms.OS.baseintc$Parameter, levels=orderassort)

ggms.OS.justoes.basec<-ggs(ms.OS.justoes.basec$samples)%>%filter(!grepl("deviance", Parameter))%>%filter(!grepl("eps", Parameter),!grepl("alpha", Parameter))%>%mutate(Parameter=case_when(
  Parameter=="betaRS"~"Reproductive state of male",
  Parameter=="Bcells[1]"~"No cells",
  Parameter=="Bcells[2]"~"Cells present"))%>%arrange(desc(Parameter))%>%na.omit()



#significance and diagnostic plots

ggs_caterpillar(ggms.OS.baseintc%>%drop_na(Parameter)%>%filter(!grepl('deviance', Parameter)), line=0, sort=FALSE)+ theme_classic()+ xlim(-1,1)+
  labs(x="Posterior probability (HPD)", y="")+scale_y_discrete(limits=rev)

ggs_density(ggms.OS.baseintc, hpd=TRUE)+xlim(-1,1)+ geom_vline(xintercept = 0)+theme_classic()+facet_wrap(~ Parameter, nrow = 4)+labs(y="Posteior samples")


ggs_compare_partial(ggms.OS.baseintc)#compares whole chain with the last value
ggs_autocorrelation(ggms.OS.baseintc)
ggs_traceplot(ggms.OS.baseintc)#traceplot of convergence
ggs_crosscorrelation(ggms.OS.baseintc)#check for highly correlated dependent variables (deep blue or red)



ggs_caterpillar(ggms.OS.justoes.basec%>%drop_na(Parameter)%>%filter(!grepl('deviance', Parameter)), line=0, sort=FALSE)+ theme_classic()+ xlim(-1,1)+
  labs(x="Posterior probability (HPD)", y="")+scale_y_discrete(limits=rev)

ggs_density(ggms.OS.justoes.basec, hpd=TRUE)+xlim(-1,1)+ geom_vline(xintercept = 0)+theme_classic()+facet_wrap(~ Parameter, nrow = 4)+labs(y="Posteior samples")


ggs_compare_partial(ggms.OS.justoes.basec)#compares whole chain with the last value
ggs_autocorrelation(ggms.OS.justoes.basec)
ggs_traceplot(ggms.OS.justoes.basec)#traceplot of convergence
ggs_crosscorrelation(ggms.OS.justoes.basec)#check for highly correlated dependent variables (deep blue or red)

