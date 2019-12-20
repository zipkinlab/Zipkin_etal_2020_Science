#####################################################################
# El Cope snake community model to assess snake species richness 
# before and after the decimation of amphibians from Bd
#####################################################################

# Name model file and set working directory 
# choose "JAGS_MODEL.R"
# rm(list=ls()) #clear workspace if needed
model<-file.choose()
setwd(dirname(model))
#load(file.choose()) # or load a previous workspace
###############################
# Load and format the data
###############################

# Detection data
cap<-read.csv("snake 1997 2012 dataset for model_June2018.csv", stringsAsFactors = FALSE)
colnames(cap)[1]<-c("Survey_ID")

# Effort data
eff<-read.csv("effort 1997 to 2012 transect surveys for model_June2018.csv", stringsAsFactors = FALSE)
colnames(eff)[1]<-c("Survey_ID")

# Format the data for JAGS and create the necessary variables 
nPeriod<-2 # time periods (before and after Bd epizootic)
nSamples<-nrow(eff) # total detections
cap$Gs<-paste(cap$Genus,cap$Species) # combine genus and species names
nSpecies<-length(unique(cap$Gs)) # number of observed species
nZero<-130 # Number of species to augment the data with all zero encounter histories 
nSum<-nZero + nSpecies
species<-data.frame(sort(unique(cap$Gs)), 1:nSpecies,stringsAsFactors = FALSE) # species list
nTran<-length(unique(cap$Transect)) # number of transects
tran<-data.frame(levels(as.factor(cap$Transect)), 1:nTran,stringsAsFactors = FALSE) # transect names

# Replace species and transect names with numbers
cap$sp_num<-NA
for (i in 1:nSpecies){
  cap[cap$Gs==species[i,1],"sp_num"] <-species[i,2] 
}
cap$tran_num<-NA
eff$tran_num<-NA
for (i in 1:nTran){
  cap[cap$Transect==tran[i,1],"tran_num"] <-tran[i,2]
  eff[eff$Transect==tran[i,1],"tran_num"] <-tran[i,2] 
}


# Make time variable for pre and post Bd
eff$bd_time<-NA
eff[eff$Year<2005,"bd_time"]<-1
eff[eff$Year>2005,"bd_time"]<-2
cap$bd_time<-NA
cap[cap$Year<2005,"bd_time"]<-1
cap[cap$Year>2005,"bd_time"]<-2

#Create species array and fill in zeros for non-detections
sp.array<-array(0, dim=c(nSamples, nSum))
for (i in 1:nSamples){
  sp.array[cap[i,"Survey_ID"],cap[i,"sp_num"]]<-1
}

#Create the necessary arguments to run the jags() command
transects<-eff$tran_num
period<-eff$bd_time

#############################################
#Create the elements to run the JAGS model
############################################

# Package data for JAGS
sp.data <- list(nSum=nSum,
                X=sp.array,
                nTran=nTran,
                nPeriod=nPeriod,
                nSamples=nSamples,
                transects=transects,
                period=period)

# create inits function (only needed for Z, w, and psi)
Z.init<-array(NA,dim=c(nSum, nTran, nPeriod))
for (i in 1:nSum){
  for (j in 1:nTran){
    for (t in 1:nPeriod){
      Z.init[i,j,t]<-max(sp.array[transects==j&period==t,i])
    }}}

w.init<-apply(Z.init,c(1,3),max)
inits <- function(){
  list(Z=Z.init,w=w.init,omega=runif(2,0.4,0.9), rho.init=runif(1,0.4,0.9))
} 


# Parameters to save
params<-c("u","u.sigma","omega","w","u.mu","v.mu","u.slope","v.slope",
          "rho", "p.diff")

# Test JAGS run
ni=10
nt=1
nb=0
nc=3
na=0
library(jagsUI)
jags.out<- jags(data=sp.data, inits=inits, parameters.to.save=params, model.file=model, store.data = TRUE,
                n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, n.adapt=na, parallel = TRUE)


# Run the JAGS model to convergence
# In previous simulation studies the autojags function in jagsUI has been unreliable
# This may take ca. 3 days
counter=0
while(max(unlist(jags.out$Rhat),na.rm = TRUE)>1.1){
  counter<-counter+1
  if(counter==1){
    na=10000
    nt=25
    ni=(5000*nt)/nc
  } else{
    na=0
    nt=25
    ni=(5000*nt)/nc
  }
  jags.out<- update(jags.out, parameters.to.save=params,
                    n.iter=ni, n.thin=nt, n.adapt=na, parallel=TRUE)
  
  save.image(paste0("MODEL_OUT",s,".RData"))
  if(length(grep("STOP.txt",list.files()))>0){break}
}

# Check for model convergence
traceplot(jags.out)
head(sort(unlist(jags.out$Rhat),decreasing = TRUE))


################################################################
# Species Richness Estiamtes and Posterior Distribution
################################################################

# Load output if needed
# load("OUTPUT_EL_COPE_SNAKE.RData")

# Check for model fit using the Bayesian p-value by taking the mean of the p.diff vector of 0 and 1 values
mean(jags.out$sims.list$p.diff)

# Create needed objects
nIter<-length(jags.out$sims.list$u.sigma) # number if MCMC iterations in posterior
infec<-grep(2, period)[1] #locate first sample from time peroid 2

# Estimate richness
w<-jags.out$sims.list$w # pull out w matrix (species estimated to be exist and be present)
N<-apply(w,c(1,3),sum)  # Calcualte posterior for species abundance

# Estiamte difference in richness
N.diff<-N[,2] - N[,1]
quantile(N.diff, probs=c(0.025, 0.5, 0.975))
quantile(N[,1], probs=c(0.025, 0.5, 0.975))
quantile(N[,2], probs=c(0.025, 0.5, 0.975))

# Calulate probabilty richness declined
sum(N.diff<0)/length(N.diff)

# Density plots for species richness
seen.pre<-sum(apply(sp.array[period==1,],2,sum)!=0)
seen.post<-sum(apply(sp.array[period==2,],2,sum)!=0)
N.data<-as.data.frame(rbind(cbind(rep("pre",nrow(N)), N[,1]),cbind(rep("post",nrow(N)), N[,2])))
N.data[,2]<-as.numeric(paste(N.data[,2]))
colnames(N.data)<-c("Bd", "N")
pre.color<-"#3399ff"
post.color<-"#cc6633"

# Calcuate mean and mode
d1<-density(N[,1])
rich.mode<-d1$x[grep(max(d1$y), d1$y)]
rich.mean<-mean(N[,1])
d2<-density(N[,2])
rich.mode[2]<-d2$x[grep(max(d2$y), d2$y)]
rich.mean[2]<-mean(N[,2])

# Plot posterior densities for species richness
library(ggplot2)
sp.richness<-ggplot() + geom_density(data=N.data[N.data$Bd=="pre",], aes(x=N), fill=pre.color, alpha=0.6,linetype=0,adjust=2.0) +
  geom_density(data=N.data[N.data$Bd=="post",], aes(x=N), fill=post.color, alpha=0.6,linetype=0,adjust=2.0) +
  geom_segment(aes(x = seen.pre, y = 0, xend = seen.pre, yend = 0.01), size=1.5, color="darkblue", linetype="solid",lineend="butt")  +
  geom_segment(aes(x = seen.post, y = 0, xend = seen.post, yend = 0.01), size=1.5, color="darkorange", linetype="solid",lineend="butt")  +
  geom_point(aes(x=rich.mean,y=c(0,0), shape="Mean"), color=c("darkblue", "darkorange"), size=3)+
  geom_point(aes(x=rich.mode,y=c(0,0), shape="Mode"), color=c("darkblue", "darkorange"), size=3)+
  theme(legend.position = c(0.75, 0.7),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +  
  labs(fill = expression(paste(italic("B. dendrobatidis "))), y="Posterior Density", shape=element_blank()) + 
  scale_x_continuous(name = "Species Richness",
                     breaks = seq(20, 120, 20), limits=c(17, 120))
sp.richness


# Table for 95% CI on occupancy probabilities
psi.all<-plogis(jags.out$sims.list$u)[,1:36,]
psi.all.q<-as.data.frame(species[,1])
psi.all.q[,2]<-apply(psi.all,c(2,3),mean)[,1]
psi.all.q[,3:4]<-t(apply(psi.all,c(2,3),quantile,probs=c(0.025,0.975))[,,1])
psi.all.q[,5]<-apply(psi.all,c(2,3),mean)[,2]
psi.all.q[,6:7]<-t(apply(psi.all,c(2,3),quantile,probs=c(0.025,0.975))[,,2])
colnames(psi.all.q)<-c("Species","pre_mean", "pre_2.5", "pre_97.5","post_mean", "post_2.5", "post_97.5")
print(psi.all.q, digits=2)
write.table(psi.all.q, file="psi_table.csv")


######################################################################
# Occupany statistics, proability psi.post< psi.pre, and graphs
######################################################################
# Graph only species with > 5 detections
colnames(species)<-c("sp", "number")
use.species<-c("Bothriechis schlegelii", "Bothrops asper", "Chironius grandisquamis", 
               "Clelia clelia", "Dipsas sp.", "Imantodes cenchoa", "Imantodes inornatus",
               "Leptodeira septentrionalis", "Leptophis depressirostris", "Oxybelis brevirostris", "Pliocercus euryzonus", 
               "Rhadinaea decorata", "Rhadinaea vermiculaticeps", "Sibon annulatus", "Sibon argus", 
               "Sibon longifrenis", "Sibon nebulatus")
use.species.num<-c(species[species$sp%in%use.species,2])

#psi indexed by iteration, species, period
num<-c(use.species.num)

# Pull out occupancy probabilities
psi.pre<-plogis(jags.out$sims.list$u[,num,1])[,1:length(num)]
psi.post<-plogis(jags.out$sims.list$u[,num,2])[,1:length(num)]

# Calculate probability occupancy decline pre to post Bd
psi.delta<-psi.post-psi.pre
psi.p<-apply(psi.delta<0, 2, sum)/nIter

# Values for graphs
values.psi.pre<-data.frame(use.species, colMeans(psi.pre), t(apply(psi.pre,2,quantile,probs=c(0.025,0.25,0.5,0.75,0.975))))
values.psi.post<-data.frame(use.species, colMeans(psi.post), t(apply(psi.post,2,quantile,probs=c(0.025,0.25,0.5,0.75,0.975))))
colnames(values.psi.post)<-colnames(values.psi.pre)<-c("species","mean","q2.5","q25","q50","q75","q97.5")

#Graph
offset<-0.125
occ.g <- ggplot() + geom_point(data=values.psi.pre, aes(x = mean, y = species), size = 2, col="dodgerblue", position=position_nudge(y=offset,x=0)) + 
  geom_errorbarh(data=values.psi.pre, aes( y = species, xmin = q2.5, xmax = q97.5), col="dodgerblue",position=position_nudge(y=offset,x=0), 
                 size = 0.3, height = 0, alpha = 0.8) +
  geom_errorbarh(data = values.psi.pre, aes( y = species, xmin = q25, xmax = q75), col="dodgerblue", position=position_nudge(y=offset,x=0),
                 size = 1.5, height = 0, alpha = 0.8) +
  geom_point(data=values.psi.post, aes(x = mean, y = species), size = 2, col="chocolate", position=position_nudge(y=-offset,x=0)) + 
  geom_errorbarh(data=values.psi.post, aes( y = species, xmin = q2.5, xmax = q97.5), col="chocolate", position=position_nudge(y=-offset,x=0), 
                 size = 0.3, height = 0, alpha = 0.8) + 
  geom_errorbarh(data = values.psi.post, aes( y = species, xmin = q25, xmax = q75), col="chocolate", position=position_nudge(y=-offset,x=0),
                 size = 1.5, height = 0, alpha = 0.8)   +    
  # geom_point( aes(x = psi.p, y = use.species), size = 3, col="darkblue", position=position_nudge(y=0,x=0)) + 
  coord_cartesian(xlim = c(0, 1)) + scale_y_discrete(limits=rev(use.species)) +  
  xlab("Occupacy Probability") + theme_bw() +
  theme(axis.text.y =element_text(face="italic"), axis.title.y = element_blank(), axis.title  = element_text(colour="black", size = 14),
        axis.text  = element_text(colour="black", size = 14)) #+ scale_y_discrete(limits=rev(use.species))
occ.g

# Probability of decline
species.order
psi.p.g <- ggplot() +    geom_point( aes(x = psi.p, y = use.species), size = 3, col="darkblue", position=position_nudge(y=0,x=0)) +
  #geom_point( aes(x = psi.p[,2], y = species$sp), size = 3, col="darkgreen", shape=17, position=position_nudge(y=0,x=0)) +
  coord_cartesian(xlim = c(0, 1)) + scale_y_discrete(limits=use.species[order(psi.p)]) + 
  theme(axis.text.y =element_text(face="italic"), axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text = element_text(colour="black", size = 14)) 
psi.p.g


###############################
# Power Analysis
###############################
# GOAL: determine the probability of detecting the species which we know occured in study area across all sites/reps between time periods
# approach:  1 - ((1-P[1,1,1])*(1-P[1,1,2]....) or 1 - The probability of not seeing species x at any sites/reps for pre-Bd or post-Bd time peroid
# generate posterior for detection probability at each sampling event
# generate eta, the covariate on detection to account for more detections of species with higher mean occupancy

num<-c(use.species.num)
#num<-1:nSpecies
eta<-array(NA, dim=c(nIter, length(num), nPeriod))
for (i in 1:length(num)){
  for (t in 1:nPeriod){
    eta[, i, t]<- jags.out$sims.list$v.mu + jags.out$sims.list$rho*(jags.out$sims.list$u[,num[i],t] - jags.out$sims.list$u.mu[,t])
  }}

#calculate per transect detection
p<- plogis(eta)

# Calculate transects in each time peroid
pre.surv<-sum(period==1)
post.surv<-sum(period==2)
site.reps<-matrix(NA, ncol=nPeriod, nrow=nTran)
for (j in 1:nTran){
  for (t in 1:nPeriod){
    site.reps[j,t]<-sum(transects==j&period==t)
  }
} 
colnames(site.reps)<-c("pre","post")
rownames(site.reps)<-levels(transect.key[,1])
mean.surv<-colMeans(site.reps)
power.pre<-1-(1-p[,,1])^mean.surv[1]
power.post<-1-(1-p[,,2])^mean.surv[2]
power.mean.pre<-colMeans(power.pre)
power.mean.post<-colMeans(power.post)
power.pre.q<-apply(power.pre,2,quantile,probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
power.post.q<-apply(power.post,2,quantile,probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
#use.species<-species[,1]
plot.names.g<-as.factor(use.species)
values.power.pre<-data.frame(use.species, power.mean.pre, t(power.pre.q))
values.power.post<-data.frame(use.species, power.mean.post, t(power.post.q))
colnames(values.power.pre)<-colnames(values.power.post)<-c("species","mean", "q2.5","q25", "q50", "q75", "q97.5")
mean(values.power.post[,2])
mean(values.power.pre[,2])

#plot
offset<-0.25
plot.names<-use.species  
power.g <- ggplot() + geom_point(data=values.power.pre, aes(x = mean, y = species), size = 2, col="dodgerblue", position=position_nudge(y=offset,x=0)) + 
  geom_errorbarh(data=values.power.pre, aes( y = species, xmin = q2.5, xmax = q97.5), col="dodgerblue",position=position_nudge(y=offset,x=0), 
                 size = 0.3, height = 0, alpha = 0.8) +
  geom_errorbarh(data = values.power.pre, aes( y = species, xmin = q25, xmax = q75), col="dodgerblue", position=position_nudge(y=offset,x=0),
                 size = 1.5, height = 0, alpha = 0.8) +
  geom_point(data=values.power.post, aes(x = mean, y = species), size = 2, col="chocolate", position=position_nudge(y=-offset,x=0)) + 
  geom_errorbarh(data=values.power.post, aes( y = species, xmin = q2.5, xmax = q97.5), col="chocolate", position=position_nudge(y=-offset,x=0), 
                 size = 0.3, height = 0, alpha = 0.8) +
  geom_errorbarh(data = values.power.post, aes( y = species, xmin = q25, xmax = q75), col="chocolate", position=position_nudge(y=-offset,x=0),
                 size = 1.5, height = 0, alpha = 0.8) +
  coord_cartesian(xlim = c(0, 1)) + 
  xlab("Probability to Detect") + theme_bw() +
  theme(axis.text.y =element_text(face="italic"), axis.title.y = element_blank(), axis.title  = element_text(colour="black", size = 14),
        axis.text  = element_text(colour="black", size = 14)) + scale_y_discrete(limits=rev(use.species))

#generate plot
power.g


#########################################################################
# Detection Rate (per survey)
##########################################################################

# Pull dectection rates and sumarize
p.pre<-p.post<-as.data.frame(use.species)
p.pre[,2]<-apply(p,c(2,3),mean)[,1]
p.post[,2]<-apply(p,c(2,3),mean)[,2]
p.pre[,3:4]<-t(apply(p,c(2,3),quantile,probs=c(0.025,0.975))[,,1])
p.post[,3:4]<-t(apply(p,c(2,3),quantile,probs=c(0.025,0.975))[,,2])
colnames(p.pre)<-colnames(p.post)<-c("species", "means", "q2.5", "q97.5")
# plot
offset=0.25
detection<- ggplot() + geom_point(data=p.pre, aes(x=means, y=species), color="dodgerblue",position=position_nudge(y=offset,x=0), size=3) + 
                       geom_point(data=p.post, aes(x=means, y=species),color="chocolate",position=position_nudge(y=-offset,x=0),size=3)  +
  geom_errorbarh(data=p.pre, aes( y = species, xmin = q2.5, xmax = q97.5),color="dodgerblue",
                 position=position_nudge(y=offset,x=0), size = 0.5, height = 0, alpha = 0.8)  +
  geom_errorbarh(data=p.post, aes( y = species, xmin = q2.5, xmax = q97.5),color="chocolate",
                 position=position_nudge(y=-offset,x=0), size = 0.5, height = 0, alpha = 0.8)  +
  xlab("Detection Rate") + theme_bw() +
  theme(axis.text.y =element_text(face="italic"), axis.title.y = element_blank(), axis.title  = element_text(colour="black", size = 14),
        axis.text  = element_text(colour="black", size = 14)) + scale_y_discrete(limits=rev(use.species))
  
detection

# Post-processing of figures in Adobe Photoshop by Dr. Sam Rossman

