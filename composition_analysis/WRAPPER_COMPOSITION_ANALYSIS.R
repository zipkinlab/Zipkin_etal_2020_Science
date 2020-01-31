###############################################################
# El Cope NMDS snake community composition before and after Bd #
###############################################################
# Remove any objects in the environment
rm(list=ls())

# Read in the data
setwd(dirname(file.choose())) # choose any file in the folder where all code and data are stored

# Detection data
cap<-read.csv("~/Zipkin_etal_InPress/data/snake_occurrence_data_from_transects.csv", stringsAsFactors = FALSE)

# Effort data
eff<-read.csv("~/Zipkin_etal_InPress/data/transect_survey_effort.csv", stringsAsFactors = FALSE)


# First column header reads in strange, so we rename it to "Survey_ID"
colnames(cap)[1]<-c("Survey_ID")
colnames(eff)[1]<-c("Survey_ID")

# Formatting the data
# Make reference objects
# Calculate the number of surveys by calculating the number of rows in the effort data file
nSamples<-nrow(eff)
# Create a new column where the Genus and species names are put together
cap$Gs<-paste(cap$Genus,cap$Species)
# Calculate the number of unique species
nSpecies<-length(unique(cap$Gs))
# Create species object with unique species names ever encountered on the surveys
species<-data.frame(sort(unique(cap$Gs)), 1:nSpecies,stringsAsFactors = FALSE)
# Calcualte the number of transects (should be 7)
nTran<-length(unique(cap$Transect))
# Create transect object with unique names for each transect
tran<-data.frame(levels(as.factor(cap$Transect)), 1:nTran,stringsAsFactors = FALSE)

# Replace species names with numbers
# Create a new column called "sp_num"
cap$sp_num<-NA
# Use a for loop to assign species number for each species
for (i in 1:nSpecies){ # For each species
  cap[cap$Gs==species[i,1],"sp_num"] <-species[i,2] 
}
# Create a new column called "tran_num" to include a unique number for transect ID
cap$tran_num<-NA
eff$tran_num<-NA
# Use a for loop to assign transect number for each transect
for (i in 1:nTran){ # For each transect
  cap[cap$Transect==tran[i,1],"tran_num"] <-tran[i,2]
  eff[eff$Transect==tran[i,1],"tran_num"] <-tran[i,2] 
}

# Make time variable for pre and post Bd
# Create an empty column called "bd_time"
eff$bd_time<-NA
# When year is < 2005 assign the # 1 = pre
eff[eff$Year<2005,"bd_time"]<-1
# When year is > 2005 assign the # 2 = post
eff[eff$Year>2005,"bd_time"]<-2
# Create an empty column called "bd_time"
cap$bd_time<-NA
# When year is < 2005 assign the # 1 = pre
cap[cap$Year<2005,"bd_time"]<-1
# When year is > 2005 assign the # 2 = post
cap[cap$Year>2005,"bd_time"]<-2

# Create species array where we have each survey as a cul
sp.array<-array(0, dim=c(nSamples, nSpecies))
for (i in 1:nSamples){
  sp.array[cap[i,"Survey_ID"],cap[i,"sp_num"]]<-1
}

# Assign columns to specific object names
transects<-eff$tran_num
period<-eff$bd_time

# Create a dataframe with the species array (1:36), the transect ID, and time period
dat2<-data.frame(sp.array[,1:nSpecies], transects, period)

# Add column names
colnames(dat2)<-c(species[,1], "transect", "period") 

# Create an empty data frame where there are 14 rows (7 transects * 2 time periods) and 36 columns (1 for each observed species)
nMDS.data<-as.data.frame(matrix(NA, ncol=nSpecies, nrow=14))

# Fill in the nMDS.data dataframe
for (i in 1:nTran){
  nMDS.data[i,]<-apply(sp.array[period==1&transects==i,1:36],2,sum)
  nMDS.data[(7+i),]<-apply(sp.array[period==2&transects==i,1:36],2,sum)
}

# Convert numbers to 0 or 1- the dissimilarity index takes in binary data
nMDS.data[nMDS.data > 0] <- 1

# Load the vegan package
library(vegan)
# Run the Ordination analysis
nMDS <- metaMDS(comm= nMDS.data, distance = "jaccard", binary = T, k = 2, trymax = 1000000, autotransform=TRUE, expand = FALSE, trace = 0, plot = FALSE)

# Examine the stress (a value <)
nMDS$stress

## Determine if there is a difference between pre- and post- Bd
  # In terms of Centroids
  # And area of the ellipticals

# Write JAGS model description into working directory
sink("ordin_mod.txt")
cat("
model {
# transect = i

# Flat priors on all parameters
    # Pre
    sigma.pre[1] ~ dgamma(0.01, 0.01)
    sigma.pre[2] ~ dgamma(0.01, 0.01)
    rho.pre ~ dunif(-1, 1)
    mu.pre[1] ~ dnorm(0, 0.01)
    mu.pre[2] ~ dnorm(0, 0.01)

    # Post
    sigma.post[1] ~ dgamma(0.01, 0.01)
    sigma.post[2] ~ dgamma(0.01, 0.01)
    rho.post ~ dunif(-1, 1)
    mu.post[1] ~ dnorm(0, 0.01)
    mu.post[2] ~ dnorm(0, 0.01)
    

# Likelihood
for(i in 1:7){
  # Pre-Bd
    y.pre[i,1:2] ~ dmnorm(mu.pre[1:2], prec.pre[1:2, 1:2]) 

  # Post-Bd
    y.post[i,1:2] ~ dmnorm(mu.post[1:2], prec.post[1:2, 1:2]) 
}

# Constructing the covariance matrix and the corresponding precision matrix.
    # Pre
    prec.pre[1:2,1:2] <- inverse(cov.pre[,])
cov.pre[1,1] <- sigma.pre[1] * sigma.pre[1]
cov.pre[1,2] <- sigma.pre[1] * sigma.pre[2] * rho.pre
cov.pre[2,1] <- sigma.pre[1] * sigma.pre[2] * rho.pre
cov.pre[2,2] <- sigma.pre[2] * sigma.pre[2]

    # Post
prec.post[1:2,1:2] <- inverse(cov.post[,])
cov.post[1,1] <- sigma.post[1] * sigma.post[1]
cov.post[1,2] <- sigma.post[1] * sigma.post[2] * rho.post
cov.post[2,1] <- sigma.post[1] * sigma.post[2] * rho.post
cov.post[2,2] <- sigma.post[2] * sigma.post[2]

}

",fill=TRUE)
sink()

# Package the data data
comm.data <- list(y.pre  = nMDS$points[1:7,], 
                  y.post = nMDS$points[8:14,])

# Parameters to monitor
params <- c("mu.post", "sigma.post", "cov.pre", "rho.pre",
            "mu.pre", "sigma.pre", "cov.post", "rho.post") 

# MCMC settings
ni=15000
nt=5
nb=5000
nc=3
na=10000

# Load the package
library(jagsUI)

# Run the model
jags.out<- jags(data=comm.data, parameters.to.save=params, model.file="ordin_mod.txt", n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, parallel = TRUE)

# Look at the model output
print(jags.out, digits=3)

# Save the model output
# save(jags.out, file = "~/Dropbox/Snake El Cope/Snakes_3/ordination.RData")

# Function to calcuate the standard ellipse area
SEA.function<-function(sigma) {
  axes<-sqrt(eigen(sigma)$values)
  return(prod(axes,pi))
}

# Calculate the elliptical surface area pre
SEA.pre<-apply(jags.out$sims.list$cov.pre, 1, SEA.function)
# Calculate the elliptical surface area post
SEA.post<-apply(jags.out$sims.list$cov.post, 1, SEA.function)

# Calculate the probability that the pre-Bd elliptical is larger than the post-Bd elliptical
# Value may differ slighly with each run because of MCMC error
sum(SEA.pre>SEA.post)/length(SEA.post)
  # 0.987

# Calculate the probability that the centroid pre-Bd is different than the post-Bd centroid
# pull posterior centroid estimates
mu.pre<-jags.out$sims.list$mu.pre[1:3000,]
mu.post<-jags.out$sims.list$mu.post
# use half of posterior centroid estimates for null distrobution
null.pre<-jags.out$sims.list$mu.pre[3001:6000,]
# create distance function
distance<-function(x,y){
  d <- sqrt((x[1]-y[1])^2 + (x[2]-y[2])^2)
  return(d)
}
# Calculate distance between centroids
dist<-c()
dist.null<-c()
for (i in 1:nrow(mu.pre)){
  dist[i]<-distance(mu.pre[i,], mu.post[i,])
  dist.null[i]<-distance(mu.pre[i,], null.pre[i,])
}
# sum instances where distance between pre and post centroid is 
# larger than that between pre and the null pre centroids
# May differ with each run because of MCMC error
mean(dist>dist.null)
#0.93

# Make the figure
color1<-"dodgerblue"
color2<-"chocolate"
colors<-c(rep(color1,7),rep(color2,7))

# Load the package
library(ellipse)
# Create the pre-Bd elliptical
e.pre<-ellipse(x=jags.out$mean$cov.pre, centre =jags.out$mean$mu.pre, level=0.6827)
# Create the post-Bd elliptical
e.post<-ellipse(x=jags.out$mean$cov.post, centre =jags.out$mean$mu.post, level=0.6827)

# Plot the nMDS points
plot(NA, ylab="Axis 2", xlab="Axis 1",  
     xlim = c(-1.2, 1.1), ylim = c(-1.1, 1.1))
# Add the ellipticals
polygon(x=e.pre[,1], y=e.pre[,2], col=NA, border=rgb(6,86,112,maxColorValue = 255))
polygon(x=e.post[,1], y=e.post[,2], col=NA, border=rgb(181,102,18,maxColorValue = 255))
points(nMDS$points[1:7,1],nMDS$points[1:7,2],pch = 19, col = rgb(6,86,112,maxColorValue = 255))
points(nMDS$points[8:14,1],nMDS$points[8:14,2],pch = 19, col = rgb(181,102,18,maxColorValue = 255))

# Post-processing of figure in Adobe Photoshop by Dr. Sam Rossman
