# SNAKE BODY CONDITION ANALYSIS 
rm(list=ls())

#Read in data
bc<-read.csv("~/Zipkin_etal_InPress/data/snake_body_condition_data.csv", 
             stringsAsFactors = FALSE)

#Generate list of species 
#Will only be analyzing the 6 species with the most data
species.list<-unique(paste(bc$Genus, bc$Species))
bc$Gs<-paste(bc$Genus, bc$Species)

#Generate body condition metric
#Mass divided by SVL^2 
#We standardized SVL by dividing by 10 so that body conditions values were not so close to zero
bc$condition<-bc$Mass/((bc$SVL/10)^2)

#Omit NAs and make a new data frame bc1
bc1<-bc[!is.na(bc$condition),]

#Add time frame
bc1$period<-NA
bc1$period[bc1$Year<2005]<-1
bc1$period[bc1$Year>2005]<-2

#Create a list of species included in the analysis
#(e.g., all species with at least 5 samples from both pre- and post-Bd time periods)
species = c("Dipsas sp.", "Imantodes cenchoa", "Leptodeira septentrionalis", "Oxybelis brevirostris",
            "Sibon annulatus", "Sibon argus")

#Make a table with with number of observations in each period by species
#Include mean and sd of the body condition pre and post
# Summary table
data.summary<-matrix(nrow=length(species), ncol=6)
rownames(data.summary)=species
colnames(data.summary)= c("#samples.pre", "mean.bd.pre", "sd.bd.pre", 
                          "#samples.post", "mean.bd.post", "sd.bd.post")

for (i in 1:length(species)){
  data.summary[i,1]<-sum(bc1$Gs==species[i]&bc1$period==1)
  data.summary[i,2]<-mean(bc1[bc1$Gs==species[i]&bc1$period==1,"condition"])
  data.summary[i,3]<-sd(bc1[bc1$Gs==species[i]&bc1$period==1,"condition"])
  
  data.summary[i,4]<-sum(bc1$Gs==species[i]&bc1$period==2)
  data.summary[i,5]<-mean(bc1[bc1$Gs==species[i]&bc1$period==2,"condition"])
  data.summary[i,6]<-sd(bc1[bc1$Gs==species[i]&bc1$period==2,"condition"])
}

#Sample sizes for each species pre (n1) and post (n2)
n1<-as.vector(data.summary[,1])
n2<-as.vector(data.summary[,4])

#Create data vectors containing all species pre (y1) and post (y2)
y1 <- matrix(NA, nrow=length(species), ncol=max(n1))
y2 <- matrix(NA, nrow=length(species), ncol=max(n2))

for (i in 1:length(species)){
  y1[i,1:n1[i]] = bc1$condition[bc1$Gs==species[i]&bc1$period==1]
  y2[i,1:n2[i]] = bc1$condition[bc1$Gs==species[i]&bc1$period==2]
}
nSpecies<-length(species)

#Write JAGS model for species' t-tests and save into working directory
sink("bcttest.txt")
cat("
    model {
    #Priors
    #species = i, sample = j
    
    for (i in 1:nSpecies) {

    mu1[i] ~ dnorm(0,0.01)
    mu2[i] ~ dnorm(0,0.01)
    sigma[i] ~ dunif(0, 10)  
    tau[i] <- 1 / ( sigma[i] * sigma[i])
    
    #Likelihood
    for (j in 1:n1[i]) {
    y1[i,j] ~ dnorm(mu1[i], tau[i]) 
    }
    
    for (j in 1:n2[i]) {
    y2[i,j] ~ dnorm(mu2[i], tau[i]) 
    }
    
    #Derived quantities
    delta[i] <- mu1[i] - mu2[i]
    P.delta[i] <- step(delta[i])
  
    }
}
    ",fill=TRUE)
sink()

#Package the components for JAGS
#Bundle data
bc.data <- list("y1", "y2", "n1", "n2", "nSpecies")

#Function that draws random inits
inits <- function(){ list(mu1=rnorm(length(species)), mu2=rnorm(length(species)), 
                          sigma = rlnorm(length(species)) ) }

#Parameters to estimate
params <- c("mu1","mu2", "delta", "sigma", "tau", "P.delta") 

#JAGS settings
ni=10000
nt=5
nb=5000
nc=3
na=10000

#Load the JAGS library and run the analysis
library(jagsUI)
jags.out<- jags(data=bc.data, inits=inits, parameters.to.save=params, 
                model.file="bcttest.txt", n.chains=nc, n.iter=ni, n.burnin=nb, 
                n.thin=nt)

#View the results
print(jags.out, digits=3)

#Calculate the probabilities that body_condition_pre > body_condition_post
#Note that probabilities will differ slightly with each model run because of MCMC error
probabilities = print(apply(jags.out$sims.list$delta>0,2,sum)/3000, digits=3)
names(probabilities)=species

probabilities

##########################################################
##########################################################

##Code to create body condition figure
##Note that the boxplots are created below then combined using illustrator

#Load ggplot library
library(ggplot2)

offset<-0.25
body.pre<-data.frame(species,jags.out$mean$mu1, jags.out$sd$mu1, jags.out$q2.5$mu1, jags.out$q25$mu1, jags.out$q75$mu1, jags.out$q97.5$mu1, 
                     stringsAsFactors = FALSE)  
colnames(body.pre)<-c("species", "mean","sd", "q2.5", "q25", "q75", "q97.5")
body.post<-data.frame(species,jags.out$mean$mu2, jags.out$sd$mu2, jags.out$q2.5$mu2, jags.out$q25$mu2, jags.out$q75$mu2, jags.out$q97.5$mu2,
                      stringsAsFactors = FALSE)  
colnames(body.post)<-c("species", "mean","sd", "q2.5", "q25", "q75", "q97.5")

######################################
# Individual Graphs
######################################
point.size=1
line.size=1.5
alpha.set=1
offset=-0.25


break.seq<-matrix(c(0.4, 0.7, 0.1,
                    0.32, 0.42, 0.02,
                    0.8, 1.1, 0.1,
                    0.4, 1.3, 0.2,
                    0.3, 0.8, 0.1,
                    0.55, 0.75, 0.05), ncol=3, byrow = TRUE)
ylims<-matrix(c(0.4, 0.7,
         0.31, 0.43, 
         0.75, 1.1,
         0.4, 1.3,
         0.35, 0.75,
         0.55, 0.75), ncol=2, byrow = TRUE)

  for (i in 1:6){
  break.seq[i,]<-c(ylims[i,], (ylims[i,2]-ylims[i,1])/4)
  body.condition.g <- ggplot() + geom_point(data=body.pre[i,], aes(y = mean, x = species), size = point.size, col="black", position=position_nudge(x=offset,y=0)) +
    geom_pointrange(data=body.pre[i,],na.rm=TRUE, aes(x=species,y=mean, ymin = q2.5, ymax =q97.5),size=line.size, col="dodgerblue",position=position_nudge(x=offset, y=0)) +
    geom_point(data=body.post[i,], aes(y = mean, x = species), size = point.size, col="black", position=position_nudge(x=-offset,y=0)) + 
    geom_pointrange(data=body.post[i,],na.rm=TRUE, aes(x=species, y=mean, ymin = q2.5, ymax =q97.5),size=line.size, col="chocolate",position=position_nudge(x=-offset,y=0)) +
    theme(axis.text.x = element_text(face = "italic"), 
          panel.border= element_rect(fill = NA, color="black"),
          panel.background = element_blank(),  
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +   
    labs(y="Mass/SVL", x="")  +
    scale_y_continuous(limits = ylims[i,], breaks=seq(break.seq[i,1],break.seq[i,2],break.seq[i,3]))
   
    plot(body.condition.g)
  }


# Post-processing of figure in Adobe Photoshop by Dr. Sam Rossman and Mollie Newman.

