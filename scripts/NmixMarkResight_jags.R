##---------------------------------------------------------------------##
## 
## Simulate and analyze repeated count and mark-resight data
##
## Data: Spatially and temporally repeated surveys collecting
##  (i)   counts of known unmarked individuals (n[s,t], e.g., bird detected, legs visible, and not bands) 
##  (ii)  counts of individuals with unknown marked status (h[s,t], e.g., bird detected but legs not visible) 
##  (iii) counts of marked but unidentified individuals (mPart[s,t], e.g. partial band reads)
##  (iv)  counts of marked and identified individuals (mFull[s,t], e.g. copmplete band reads)
##  (v)  capture histories of marked and identified individuals (y[i,t], e.g., complete band reads)
##
## Three modeling approaches: 
##  1. Integrated model 
##  2. N-mix model using only count data (counts[s,t])
##
##  updated: 30-May-2018
##
##---------------------------------------------------------------------##

rm(list=ls())
#setwd("C:\\Nathan\\pwrc\\projects\\American Oystercatchers\\R files\\simulation\\workspaces\\unknown banding status")
#setwd("~/AMOY/simulation/workspace")

library(jagsUI)
nAdapt=10000; nc=3; nb=25000; ni=200000+nb; nt=10  #saved iterations will be nc*((ni-nb)/nt)
#nAdapt=50; nc=3; nb=50; ni=50+nb; nt=1
JAGSparallel <-TRUE  #run JAGS in parallel
nsim<-100

Scenario <- 1
S<-50           #sites
lambda<-15      #abundance rate
K<-4            #occasions
phi<-0.50       #probability individual is previously marked (i.e. proportion of population that is marked)
p<-0.50         #probability of detecting an individual
delta<-0.75     #probability band marked was observed given the individual was detected and previously marked
theta<-0.75     #probability tag was read given the individual was detected, previously marked, and band status was identified
pStar<-1-(1-(p*delta*theta))^K    #probability mark was uniquely identified at least once, conditional on individual being marked 

#placeholders for stuff to save
LambdaSim<-lamBCI<-NtotSim<-pSim<-pBCI<-matrix(NA, nsim,2)                 
thetaSim<-phiSim<-deltaSim<-thetaBCI<-phiBCI<-deltaBCI<-matrix(NA, nsim,1)
mean_counts<-mean_n<-mean_m<-mean_mFull<-mean_mPart<-mean_h<-mean_y<-sum_mFull<-totMarked<-totN<-sum_m<-NA
RhatFULL<-matrix(NA, nr=nsim, nc=10)
RhatNMIX<-matrix(NA, nr=nsim, nc=4)


# start simulation
for(sim in 1:nsim){
  set.seed(2018*sim)  #set.seed inside of simulation loop so datasets are identicial across scenarios (i.e. scenarios sample identical N) 
  
  #Generate abundance and count data
  N<-rpois(S,lambda)
  Ntot<-sum(N)
  site<-rep(1:S,N)                  #site associated with each individual.
  marked<-NA
  for(i in 1:Ntot){
    marked[i]<-rbinom(1,1,phi)       #is individual marked (Y/N)
  }
  NmarkedTot<-sum(marked)
  markedSite<-site[which(marked==1)]
  Nmarked <- rep(0,S)  #placehoder for marked individuals per site
  Nmarked[as.numeric(names(table(markedSite)))]<-table(markedSite) #marked individuals per site
  
  #Individuals have some true (latent) encounter history (yTrue) (this could be one loop, but for reproducability it is two loops)
  DetTrue<-yTrue<-matrix(NA, nr=Ntot, nc=K)      #unique row for each individual
  for(i in 1:Ntot){
    for(k in 1:K){
      DetTrue[i,k]<-rbinom(1,1,p)
    }
  }
  for(i in 1:Ntot){  
    for(k in 1:K){
      if(DetTrue[i,k]==0) yTrue[i,k]<-5   #not detected = 5 in multinomial (placeholder)
      if(DetTrue[i,k]==1) {               #if detected, how is it classified
        yTrue[i,k]<-which(rmultinom(1,1,c((1-marked[i])*delta,   #1 = observed, known unmarked
                                          (1-delta),                         #2 = observed, unknown marking status
                                          (marked[i]*delta*(1-theta)),       #3 = observed, known marking status, no band ID
                                          (marked[i]*delta*theta)))==1)      #4 = observed, known marking status, band ID
      }#if
    }#K
  }#i
  
  #for individuals without complete band reads, encounter histories are collapsed to count data
  #site and survey specific count data (all individuals ('counts'), known unmarked ('n'), unknown marked status ('h'), partial band reads 'r' )
  c_obs<-n_obs<-mPart_obs<-h_obs<-mFull_obs<-list()
  counts<-n<-mPart<-h<-mFull<-matrix(0, nr=S, nc=K)
  for(k in 1:K){
    #counts of unmarked+marked individuals
    c_obs[[k]] <-which(yTrue[,k]<5)     
    allObs<-table(site[c_obs[[k]]])
    counts[as.numeric(names(allObs)),k]<-allObs   
    
    #counts of known unmarked individuals
    n_obs[[k]] <-which(yTrue[,k]==1)    
    nObs<-table(site[n_obs[[k]]])
    n[as.numeric(names(nObs)),k]<-nObs         
    
    #counts of unknown marking status individuals
    h_obs[[k]] <-which(yTrue[,k]==2)    
    hObs<-table(site[h_obs[[k]]])
    h[as.numeric(names(hObs)),k]<-hObs     
    
    #counts of full band read individuals
    mFull_obs[[k]] <-which(yTrue[,k]==4)    
    mFullObs<-table(site[mFull_obs[[k]]])
    mFull[as.numeric(names(mFullObs)),k]<-mFullObs   
    
    #counts of partial band read individuals
    mPart_obs[[k]] <-which(yTrue[,k]==3)    
    mPartObs<-table(site[mPart_obs[[k]]])
    mPart[as.numeric(names(mPartObs)),k]<-mPartObs   
  }
  
  m<-mFull+mPart  #detection of known marked individuals 
  mc<-array(NA, c(2,S,K))
  mc[1,,]<-mFull
  mc[2,,]<-mPart
  
  
  #Mark-resight data
  MR<- which(apply(yTrue,1,function(x) max(x==4))==1)     #marked individuals identified at least once during the study
  mTot<-length(MR)                                        #total number of detected marked individuals
  y<-yTrue[MR,]                                           #MR encounter histories for observed individuals
  y[y!=4]<-0                                              #keep only observations where the band was identified
  y[y==4]<-1
  mSite<-site[MR]                                          #Site of detection 
  
  
  ## END DATA SIMULATION
  
  ## ANALYSES in JAGS
  
  
  ####################################
  #Integrated Nmix & mark-resight with unknown banding status 
  
  # BUGS model, Conditional 4-part model
  cat("
      model {
      # Priors
      lambda ~ dunif(0,100)         #abundance rate
      phi~dunif(0,1)                #probability an individual is marked
      p ~ dunif(0, 1)               #probability of detecting an individual  
      delta ~ dunif(0,1)            #probability marked status observed
      theta~dunif(0,1)              #probability tag was read given the individual was detected
      
      #conditional on detecting a marked individual, probability of classification type
      pmc[1]<-theta            #known marked and bands read
      pmc[2]<-(1-theta)        #known marked but bands not read
      
      #conditional mark-resight parameters
      pMR<-p*theta*delta            #probability of reading bands of a marked individual
      pStar<-(1-(1-pMR)^K)          #probability band is read at least once given individual is marked
      pMRc<-(pMR/pStar)             #conditional probability of recording bands of a marked individual
      
      #Likelihoods
      for (s in 1:S) {
      #Total abundance
      N[s] ~ dpois(lambda)
      Nmarked[s]~dbin(phi, N[s])    #Nmarked[s] could be provided as data if known
      Nunmarked[s]<-N[s]-Nmarked[s] 
      
      #count data
      for (k in 1:K) {
      n[s,k] ~ dbin(p*delta, Nunmarked[s])   #probability of detecting and observing no-bands of an unmarked individual
      m[s,k] ~ dbin(p*delta, Nmarked[s])     #probability of detecting and observing bands of a marked individual
      h[s,k]~dbin((1-delta), counts[s,k])    #probability marked status is unobserved given detection
      
      #classification of count types, condtional on detecting marked individual
      mc[1:2,s,k]~dmulti(pmc[1:2], m[s,k])
      }#j
      }#s
      
      #mark-resight observation data
      for(i in 1:mTot){
      for (k in 1:K) {
      y[i,k]~dbern(pMRc)
      }#j
      }#i
      
      #derived total abundances
      Ntot<-sum(N[1:S])
      NmarkedTot<-sum(Nmarked[1:S])
      NunmarkedTot<-sum(Nunmarked[1:S])
      
      }#model
      ", fill=TRUE, file="NmixIntegrated_UnknownMarkingStatus.txt")
  
  
  # JAGS data
  #JAGS throws an error for matrices with 1 row.
  if(mTot==1) y<-rbind(y,NA)  #add row of NA. Does not alter estimation, but prevents error.
  
  jags.data <- list(S = S, K = K, n=n, m=m, h=h, counts=counts, 
                    y=y, mTot=mTot)
  
  # Initial values
  Nm<-apply(m,1,max)+5                          #intital values for the number of marked individuals
  Ns<- pmax(Nm, apply(counts, 1, max))+10       #intital values for the number of unmarked individuals  
  inits <- function() list(N = Ns, Nmarked=Nm, lambda=runif(1,10,18), phi=runif(1,.25,.50),p=runif(1,.25,.50), theta=runif(1,0.25,0.75), delta=runif(1,.25,0.75))
  
  # Parameters monitored
  params <- c("p", "theta","pStar","phi","delta", "lambda", "Ntot", "NmarkedTot", "NunmarkedTot")
  
  # Call JAGS 
  out <- jags(jags.data, inits, params, "NmixIntegrated_UnknownMarkingStatus.txt", n.adapt=nAdapt, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=JAGSparallel)  
  #print(out,2)
  
  
  ####################################
  ##Nmix only
  
  # BUGS model
  cat("
      model {
      # Priors
      #total abundance
      lambda ~ dunif(0,100) 
      p ~ dunif(0, 1)               #probability of detecting an individual 
      
      #Likelihoods
      for (s in 1:S) {
      #Total abundance
      N[s] ~ dpois(lambda)
      #count data
      for (k in 1:K) {
      #N-mix detection
      counts[s,k] ~ dbin(p, N[s])
      }#j
      }#s
      
      Ntot<-sum(N[1:S])
      }#model", fill=TRUE, file="Nmix.txt")
  
  # JAGS data
  jags.data <- list(S = S, K = K, counts=counts)
  
  Ns<- apply(counts, 1, max)+10
  inits <- function() list(N = Ns, lambda=runif(1,10,18), p=runif(1,.25,.50))
  
  # Parameters monitored
  params <- c("p", "lambda", "Ntot")
  
  # Call JAGS 
  Nmix <- jags(jags.data, inits, params, "Nmix.txt", n.adapt=nAdapt, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=JAGSparallel)  
  #print(Nmix,2)
  
  
  
  ## Harvest results
  LambdaSim[sim,]<-c(out$q50$lambda,Nmix$q50$lambda)
  NtotSim[sim,]<-c(out$q50$Ntot,Nmix$q50$Ntot)
  pSim[sim,]<-c(out$q50$p,Nmix$q50$p)
  thetaSim[sim,]<-c(out$q50$theta)
  phiSim[sim,]<-c(out$q50$phi)
  deltaSim[sim,]<-c(out$q50$delta)
  lamBCI[sim,]<-c(out$q2.5$lambda<lambda&lambda<out$q97.5$lambda,  Nmix$q2.5$lambda<lambda&lambda<Nmix$q97.5$lambda)
  pBCI[sim,]<-c(out$q2.5$p<p&p<out$q97.5$p,  Nmix$q2.5$p<p&p<Nmix$q97.5$p)
  thetaBCI[sim,]<-c(out$q2.5$theta<theta&theta<out$q97.5$theta)
  phiBCI[sim,]<-c(out$q2.5$phi<phi&phi<out$q97.5$phi)
  deltaBCI[sim,]<-c(out$q2.5$delta<delta&delta<out$q97.5$delta)
  
  RhatFULL[sim,]<-out$summary[,"Rhat"]
  RhatNMIX[sim,]<-Nmix$summary[,"Rhat"]
  
  
  #what data look like
  mean_counts[sim]<- mean(counts)
  mean_n[sim]<-mean(n)
  mean_m[sim]<-mean(m)
  mean_mFull[sim]<-mean(mFull)
  mean_mPart[sim]<-mean(mPart)
  mean_h[sim]<-mean(h)
  sum_mFull<-mTot
  mean_y[sim]<-ifelse(is.integer(dim(y)), mean(apply(y,1,sum), mean(y)))  #mean number of occasions an individual was uniquely identified
  totMarked[sim]<-NmarkedTot
  totN[sim]<-Ntot
  
  print(sim)  #to help keep track
  save.image(paste("Scenario",Scenario,"_S",S,"_K",K,"_lambda",lambda,"_phi",phi,"_p",p,"_delta",delta,"_theta",theta,"_",sim,".Rdata",sep=""))
}#sim


#Plot results
dev.new(width=5, height=10)
par(mar=c(2,4,0,0)+.1, mfrow=c(5,1))
boxplot(LambdaSim, names=NA, col="grey", ylab=expression(lambda), las=1, ylim=c(5,30));abline(h=lambda, lwd=2, col="red")
boxplot(cbind(phiSim[,1],NA), names=NA, col="grey", ylab=expression(phi), las=1, ylim=c(0,.6));abline(h=phi, lwd=2, col="red"); text(2,phi,"NA")
boxplot(pSim, names=NA, col="grey", ylab="p", las=1, ylim=c(.05,.75));abline(h=p, lwd=2, col="red")
boxplot(cbind(deltaSim[,1],NA), names=NA, col="grey", ylab=expression(delta), las=1, ylim=c(.5,1.0));abline(h=delta, lwd=2, col="red"); text(2,delta,"NA")
boxplot(cbind(thetaSim[,1],NA), names=c("Integrated", "N-mix"), col="grey", ylab=expression(theta), las=1, ylim=c(.1,1.0));abline(h=theta, lwd=2, col="red"); text(2,theta,"NA")

#credible interval coverage, bias, and RMSE
BCI<- rbind(apply(lamBCI,2,mean), c(mean(phiBCI[,1]),NA), apply(pBCI,2,mean), c(mean(thetaBCI[,1]),NA), c(mean(deltaBCI[,1]),NA))
rownames(BCI)<-c("lambda","phi","p","theta","delta")
colnames(BCI)<-c("Integrated","Nmix")

BiasRMSE<-matrix(NA,6,nr=5)
BiasRMSE[1,]<-c(apply(LambdaSim,2,mean), apply((LambdaSim-lambda)/lambda,2,mean), sqrt(apply((LambdaSim-lambda)^2,2,mean)))
BiasRMSE[2,c(1,3,5)]<-c(apply(phiSim,2,mean), apply((phiSim-phi)/phi,2,mean), sqrt(apply((phiSim-phi)^2,2,mean)))
BiasRMSE[3,]<-c(apply(pSim,2,mean),apply((pSim-p)/p,2,mean), sqrt(apply((pSim-p)^2,2,mean)))
BiasRMSE[4,c(1,3,5)]<-c(apply(deltaSim,2,mean), apply((deltaSim-delta)/delta,2,mean), sqrt(apply((deltaSim-theta)^2,2,mean)))
BiasRMSE[5,c(1,3,5)]<-c(apply(thetaSim,2,mean), apply((thetaSim-theta)/theta,2,mean), sqrt(apply((thetaSim-theta)^2,2,mean)))

simTable<-cbind(c(lambda,phi,p,theta,delta),BiasRMSE, BCI)
colnames(simTable)<-c("true",rep(colnames(BCI),4))
#round(simTable,3)

#what data look like
datTable<-matrix(NA, nr=10, nc=4)
colnames(datTable)<-c("mean","sd","min","max")
rownames(datTable)<-c("Ntot","Nmarked","counts","n","m","mFull","mPart","h","mTot","y")

resultsFun<-function(x){c(mean(x), sd(x), min(x), max(x))}
datTable[1,]<-resultsFun(totN)
datTable[2,]<-resultsFun(totMarked)
datTable[3,]<-resultsFun(mean_counts)
datTable[4,]<-resultsFun(mean_n)
datTable[5,]<-resultsFun(mean_m)
datTable[6,]<-resultsFun(mean_mFull)
datTable[7,]<-resultsFun(mean_mPart)
datTable[8,]<-resultsFun(mean_h)
datTable[9,]<-resultsFun(sum_mFull)
datTable[10,]<-resultsFun(mean_y)
#round(datTable,2)


#convergence
#non-converged simulations
# which(apply(RhatFULL,1,max)>1.1)  
# which(apply(RhatNMIX,1,max)>1.1)  
# which(apply(RhatMISS,1,max, na.rm=T)>1.1) 
#average Rhat
#apply(RhatFULL,2,mean, na.rm=T)
#apply(RhatNMIX,2,mean, na.rm=T)
#apply(RhatMISS,2,mean, na.rm=T)

save.image(paste("FINAL_Scenario",Scenario,"_S",S,"_K",K,"_lambda",lambda,"_phi",phi,"_p",p,"_delta",delta,"_theta",theta,"_",sim,".Rdata",sep=""))


