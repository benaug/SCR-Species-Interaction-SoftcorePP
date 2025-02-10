#this version includes a baseline interaction parameter, h0
#h0 and omega weakly identifiable

library(nimble)
library(coda)
nimbleOptions(determinePredictiveNodesInModel = FALSE)
source("sim.SCR.2species.softcore.R")
source("init.data.2species.softcore.R")
source("NimbleModel 2species Softcore Baseline.R")
source("Nimble Functions 2species Softcore Baseline.R") #nimble functions used in data simulator
source("sSampler Dcov 2species Softcore Discrete.R") #need to modify this if switching to continuous state space
#make state space grid
#get some colors
library(RColorBrewer)
cols1 <- brewer.pal(9,"Greens")

#trapping array, based on sigmas
sigma1 <- 4 #sps1 spatial scale parameter
sigma2 <- 2 #sps2 spatial scale parameter

mean.sigma <- mean(c(sigma1,sigma2))
max.sigma <- max(c(sigma1,sigma2))
spacing <- mean.sigma*2
res <- spacing/3 #resolution, cell width/height
tmp <- seq(max.sigma*3+res/2,by=spacing,length.out=10)
X <- as.matrix(expand.grid(tmp,tmp))

#state space. Must start at (0,0)
xlim <- c(0,max(X[,1]+3*max.sigma+res/2))
ylim <- c(0,max(X[,2]+3*max.sigma+res/2))
if(xlim[1]!=0|ylim[1]!=0)stop("xlim and ylim must start at 0.")
if((diff(range(xlim))/res)%%1!=0)stop("The range of xlim must be divisible by 'res'")
if((diff(range(ylim))/res)%%1!=0)stop("The range of ylim must be divisible by 'res'")

#make discrete state space objects
x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res)
y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res)
dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
n.cells <- nrow(dSS)
n.cells.x <- length(x.vals)
n.cells.y <- length(y.vals)
InSS <- rep(1,n.cells) #can exclude state space cells with InSS. Not excluding any here


#simulate a D.cov, higher cov.pars for large scale cov. takes a while if n.cells is large, say over 10,000
set.seed(1320567)
library(geoR)
D.cov <- grf(n.cells,grid=dSS,cov.pars=c(100,100),messages=FALSE)[[2]]
D.cov <- as.numeric(scale(D.cov)) #scale
par(mfrow=c(1,1),ask=FALSE)

image(x.vals,y.vals,matrix(D.cov*InSS,n.cells.x,n.cells.y),main="D.cov",xlab="X",ylab="Y",col=cols1)
# grid(n.cells.x,n.cells.y,lwd=1,lty=1,col="grey80") #helps see how traps are placed in each
points(X,pch=4,lwd=2)

D.cov1 <- D.cov2 <- D.cov #using same D.cov for both species

K <- 10

#species 1
D.beta01 <- -5.5 #baseline D
D.beta11 <- 1 #density coefficient 
p01 <- 0.125 #baseline detection prob

#species 2
D.beta02 <- -4 #baseline D
D.beta12 <- 1 #density coefficient
h0 <- 0.5 #baseline repulsion parameter
omega <- 4 #repulsion scale parameter
p02 <- 0.25 #baseline detection prob

set.seed(323413) #change this for new data set

data <- sim.SCR.2species.softcore(D.beta01=D.beta01,D.beta11=D.beta11,
               p01=p01,sigma1=sigma1,D.cov1=D.cov1,h0=h0,omega=omega,
               D.beta02=D.beta02,D.beta12=D.beta12,D.beta22=D.beta22,
               p02=p02,sigma2=sigma2,D.cov2=D.cov2,
               X=X,xlim=xlim,ylim=ylim,res=res,InSS=InSS,K=K)

#sps 1
data$truth$lambda1 #expected abundance from D cov inputs
data$truth$N1 #simulated realized abundance
data$truth$n1 #number of inds captured
table(rowSums(data$capture$y1)) #number of inds captures X times

#sps 2
data$truth$lambda2 #expected abundance from D cov inputs
data$truth$N2 #simulated realized abundance
data$truth$n2 #number of inds captured
table(rowSums(data$capture$y2)) #number of inds captures X times

##Fit Model
#data augmentation: Must be larger than simulated N. If N posterior hits M, need to raise M and try again.
M1 <- 75 #data augmentation limit sps 1. 
M2 <- 175 #data augmentation limit sps 2.
if(M1<=data$truth$N1)stop("Raise M1 to be larger than simulated N1.")
if(M2<=data$truth$N2)stop("Raise M2 to be larger than simulated N2.")

inits <- list(sigma1=5,sigma2=5) #needs to be set somewhere in the ballpark of truth
nimbuild <- init.data.2species.softcore(data=data,inits=inits,M1=M1,M2=M2)

Niminits <- list(z1=nimbuild$z1,N1=nimbuild$N1, #must init N to be sum(z.init)
                 s1.cell=nimbuild$s1.cell,
                 D01=sum(nimbuild$z1)/(sum(data$constants$InSS)*data$constants$res^2),D.beta11=0,
                 sigma1=inits$sigma1,
                 p01=0.25,
                 z2=nimbuild$z2,N2=nimbuild$N2, #must init N to be sum(z.init)
                 s2.cell=nimbuild$s2.cell,
                 D02=sum(nimbuild$z2)/(sum(data$constants$InSS)*data$constants$res^2),D.beta12=0,omega=2,
                 sigma2=inits$sigma2,
                 p02=0.25,h0=0.5)

#constants for Nimble
J <- nrow(data$constants$X)
K1D <- rep(data$constants$K,J)
constants <- list(M1=M1,M2=M2,J=J,K1D=K1D,n.cells=data$constants$n.cells,
                  cellArea=data$constants$res^2,D.cov1=data$constants$D.cov1,D.cov2=data$constants$D.cov2)

#supply data to nimble
Nimdata <- list(y1=nimbuild$y2D1,y2=nimbuild$y2D2,X=data$constants$X,
                dSS=data$constants$dSS,InSS=data$constants$InSS)  #fix these to 0 to turn of D.covs D.beta11=0,D.beta12=0

# set parameters to monitor
parameters <- c('p01','D.beta11','sigma1','N1','D01','lambda1',
              'p02','D.beta12','omega','h0','sigma2','N2','D02','lambda2')

#can also monitor a different set of parameters with a different thinning rate
nt <- 2 #thinning rate

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
#tell nimble which nodes to configure so we don't waste time for samplers we will replace below
config.nodes <- c("p01","sigma1","p02","sigma2")
# config.nodes <- c()
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,useConjugacy = FALSE,
                      nodes=config.nodes) 

##Sps 1 z/N sampler
z.ups <- round(M1*0.25) # how many N/z proposals per iteration? Not sure what is optimal, setting to 25% of M here.
y.nodes <- Rmodel$expandNodeNames(paste("y1[1:",M1,",1:",J,"]"))
pd.nodes <- Rmodel$expandNodeNames(paste("pd1[1:",M2,",1:",J,"]"))
N.node <- Rmodel$expandNodeNames(paste("N1"))
z.nodes <- Rmodel$expandNodeNames(paste("z1[1:",M1,"]"))
sps2.nodes <- Rmodel$getDependencies("h")[2:5] #lambda2.cell, pi2.denom, pi2.cell, lambda2 (do not include h)
sps2.lik.nodes <- c(Rmodel$expandNodeNames("s2.cell"),Rmodel$expandNodeNames("N2"))
d2.nodes <- Rmodel$expandNodeNames("d2")
IntKern.nodes <- Rmodel$expandNodeNames("IntKern")
calcNodes <- c(N.node,y.nodes,IntKern.nodes,sps2.lik.nodes)
ind.detected <- 1*(rowSums(nimbuild$y1)>0)
conf$addSampler(target = c("N1"),
                type = 'zSampler1',control = list(z.ups=z.ups,M=M1,J=data$constants$J,ind.detected=ind.detected,
                                                  dSS=data$constants$dSS,n.cells=data$constants$n.cells,
                                                  y.nodes=y.nodes,pd.nodes=pd.nodes,N.node=N.node,z.nodes=z.nodes,
                                                  d2.nodes=d2.nodes,IntKern.nodes=IntKern.nodes,
                                                  sps2.nodes=sps2.nodes,sps2.lik.nodes=sps2.lik.nodes,
                                                  calcNodes=calcNodes),
                silent = TRUE)


##Sps 2 z/N sampler
z.ups <- round(M2*0.25) # how many N/z proposals per iteration? Not sure what is optimal, setting to 25% of M here.
y.nodes <- Rmodel$expandNodeNames(paste("y2[1:",M2,",1:",J,"]"))
pd.nodes <- Rmodel$expandNodeNames(paste("pd2[1:",M2,",1:",J,"]"))
N.node <- Rmodel$expandNodeNames(paste("N2"))
z.nodes <- Rmodel$expandNodeNames(paste("z2[1:",M2,"]"))
calcNodes <- c(N.node,y.nodes)
ind.detected <- 1*(rowSums(nimbuild$y2)>0)
conf$addSampler(target = c("N2"),
                type = 'zSampler2',control = list(z.ups=z.ups,M=M2,J=data$constants$J,ind.detected=ind.detected,
                                                  y.nodes=y.nodes,pd.nodes=pd.nodes,N.node=N.node,z.nodes=z.nodes,
                                                  calcNodes=calcNodes),
                silent = TRUE)

# Sps 1: activity center sampler
for(i in 1:M1){
  calcNodes <- Rmodel$getDependencies(paste("s1.cell[",i,"]", sep=""))
  sps2.nodes <- Rmodel$getDependencies("h")[2:5] #lambda2.cell, pi2.denom, pi2.cell, lambda2 (do not include h)
  sps2.lik.nodes <- c(Rmodel$expandNodeNames("s2.cell"),Rmodel$expandNodeNames("N2"))
  d2.nodes <- Rmodel$expandNodeNames(paste("d2[",i,"]"))
  IntKern.nodes <- Rmodel$expandNodeNames(paste("IntKern[",i,"]"))
  s1.nodes <- setdiff(calcNodes,c(Rmodel$getDependencies("h")[1:5],sps2.lik.nodes,IntKern.nodes,d2.nodes))
  conf$addSampler(target = paste("s1.cell[",i,"]", sep=""),
                  type = 'sSamplerDcov1',control = list(i=i,xlim=data$constants$xlim,
                                                          ylim=data$constants$ylim,
                                                          cells=data$constants$cells,
                                                          res=data$constants$res,
                                                          s1.nodes=s1.nodes,sps2.nodes=sps2.nodes,
                                                          sps2.lik.nodes=sps2.lik.nodes,
                                                          d2.nodes=d2.nodes,IntKern.nodes=IntKern.nodes,
                                                          calcNodes=calcNodes), silent = TRUE)
}
# Sps 2: activity center sampler
for(i in 1:M2){
  calcNodes <- Rmodel$getDependencies(paste("s2.cell[",i,"]", sep=""))
  conf$addSampler(target = paste("s2.cell[",i,"]", sep=""),
                  type = 'sSamplerDcov2',control = list(i=i,xlim=data$constants$xlim,
                                                           ylim=data$constants$ylim,
                                                           cells=data$constants$cells,
                                                           res=data$constants$res,
                                                           calcNodes=calcNodes), silent = TRUE)
}

#use block samplers for these, AF_slice works best, particularly when using D0 instead of D.beta0.
conf$addSampler(target = c("D01","D.beta11"),
                type = 'AF_slice',control=list(adaptive=TRUE),silent = TRUE)
conf$addSampler(target = c("D02","D.beta12",'omega','h0'),
                type = 'AF_slice',control=list(adaptive=TRUE),silent = TRUE)
#adding block samplers for p and sigma (retaining independent samplers, too). 
#posterior correlation can be high, but AF_slice slow for obsmod. 
#RW_block not great either bc parameters on different scales
#was worth the cost in an n=1 comparison I did
conf$addSampler(target = c("p01","sigma1"),
                type = 'RW_block',control=list(adaptive=TRUE),silent = TRUE)
conf$addSampler(target = c("p02","sigma2"),
                type = 'RW_block',control=list(adaptive=TRUE),silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Cmodel$calculate()

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(1500,reset=FALSE) #can keep running this line to extend sampler
end.time <- Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

library(coda)
mvSamples <-  as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[25:nrow(mvSamples),])) #discarding some burnin here. Can't plot 1st sample which is all NA

#sps 1
data$truth$lambda1 #target expected abundance
data$truth$N1 #target realized abundance

#sps 2
data$truth$lambda2 #target expected abundance
data$truth$N2 #target realized abundance

#looking at posterior correlation here
burnin <- 25
rem.idx <- which(colnames(mvSamples)%in%c("N1","N2","lambda1","lambda2"))
tmp <- cor(mcmc(mvSamples[burnin:nrow(mvSamples),-rem.idx]))
diag(tmp) <- NA
which(abs(tmp)>0.6,arr.ind=TRUE)

#plot interaction function for final iteration
par(mfrow=c(1,1),ask=FALSE)
image(x.vals,y.vals,matrix(Cmodel$h,n.cells.x,n.cells.y),main="D.cov",xlab="X",ylab="Y",col=cols1)
points(X,pch=4)
points(Cmodel$s1[Cmodel$z1==1,],pch=16,col="darkred")
points(Cmodel$s2[Cmodel$z2==1,],pch=16,col="goldenrod")
