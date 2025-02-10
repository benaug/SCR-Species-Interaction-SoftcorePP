NimModel <- nimbleCode({
  #--------------------------------------------------------------
  # priors
  #--------------------------------------------------------------
  #Density covariates
  # D.beta01 ~ dnorm(0,sd=?)
  # D.beta02 ~ dnorm(0,sd=?)
  # D01 <- exp(D.beta01)
  # D02 <- exp(D.beta02)
  D01 ~ dunif(0,100) #uninformative, diffuse dnorm on log scale can cause neg bias
  D02 ~ dunif(0,100) #uninformative, diffuse dnorm on log scale can cause neg bias
  D.beta11 ~ dnorm(0,sd=10) #sps 1 response to D.cov1
  D.beta12 ~ dnorm(0,sd=10) #sps 2 response to D.cov2
  #interaction function
  h0 ~ dunif(0,1)
  omega ~ dunif(0,20)
  #baseline detection prob
  p01 ~ dunif(0,1)
  p02 ~ dunif(0,1)
  #detection spatial scale
  sigma1 ~ dunif(0,20)
  sigma2 ~ dunif(0,20)

  #Density model
  D1.intercept <- D01*cellArea
  D2.intercept <- D02*cellArea
  lambda1.cell[1:n.cells] <- InSS[1:n.cells]*exp(D.beta11*D.cov1[1:n.cells])
  h[1:n.cells] <- getH(IntKern=IntKern[1:M1,1:n.cells],z=z1[1:M1]) #Interaction function
  lambda2.cell[1:n.cells] <- InSS[1:n.cells]*exp(D.beta12*D.cov2[1:n.cells])*h[1:n.cells]
  pi1.cell[1:n.cells] <- lambda1.cell[1:n.cells]/pi1.denom #expected proportion of total N in cell c
  pi2.cell[1:n.cells] <- lambda2.cell[1:n.cells]/pi2.denom #expected proportion of total N in cell c
  pi1.denom <- sum(lambda1.cell[1:n.cells])
  pi2.denom <- sum(lambda2.cell[1:n.cells])
  lambda1 <- D1.intercept*pi1.denom #Expected N, sps 1
  lambda2 <- D2.intercept*pi2.denom #Expected N, sps 2
  N1 ~ dpois(lambda1) #Realized N, sps 1
  N2 ~ dpois(lambda2) #Realized N, sps2
  for(i in 1:M1){#sps 1
    s1.cell[i] ~ dcat(pi1.cell[1:n.cells])
    s1[i,1:2] <- dSS[s1.cell[i],1:2] #pull out X,Y for this cell
    d2[i,1:n.cells] <- getd2(s=s1[i,1:2],dSS[1:n.cells,1:2],z=z1[i]) #distance between s1 and dSS
    IntKern[i,1:n.cells] <- getIntKern(h0=h0,d2=d2[i,1:n.cells],z=z1[i],omega=omega) #compute for z1=1 inds only
    pd1[i,1:J] <- GetDetectionProb(s = s1[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma1, p0=p01, z=z1[i])
    y1[i,1:J] ~ dBinomialVector(size=K1D[1:J],prob=pd1[i,1:J],z=z1[i])
  } #custom N1/z1 update, maintains sum(z1)=N1
  for(i in 1:M2){#sps 2
    s2.cell[i] ~ dcat(pi2.cell[1:n.cells])
    s2[i,1:2] <- dSS[s2.cell[i],1:2] #pull out X,Y for this cell
    pd2[i,1:J] <- GetDetectionProb(s = s2[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma2, p0=p02, z=z2[i])
    y2[i,1:J] ~ dBinomialVector(size=K1D[1:J],prob=pd2[i,1:J],z=z2[i])
  } #custom N2/z2 update, maintains sum(z2)=N2
})# end model