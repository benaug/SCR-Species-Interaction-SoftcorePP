e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.SCR.2species.softcore <-
  function(D.beta01=NA,D.beta11=NA,p01=NA,sigma1=NA, D.cov1=NA,
           D.beta02=NA,D.beta12=NA,D.beta22=NA,p02=NA,sigma2=NA,h0=1,omega=NA,
           D.cov2=NA,rsf.cov2=NA,K=NA,X=X,xlim=NA,ylim=NA,res=NA,InSS=NA,seed=NA){
    
    if(!is.na(seed)){
      set.seed(seed)
    }
    par(mfrow=c(1,1),ask=FALSE)
    if(xlim[1]!=0|ylim[1]!=0)stop("xlim and ylim must start at 0.")
    if((diff(range(xlim))/res)%%1!=0)stop("The range of xlim must be divisible by 'res'")
    if((diff(range(ylim))/res)%%1!=0)stop("The range of ylim must be divisible by 'res'")
    # make discrete state space
    x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res)
    y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res)
    dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
    cells <- matrix(1:nrow(dSS),nrow=length(x.vals),ncol=length(y.vals))
    n.cells <- nrow(dSS)
    n.cells.x <- length(x.vals)
    n.cells.y <- length(y.vals)
    cellArea <- res^2
    library(RColorBrewer)
    cols1 <- brewer.pal(9,"Greens")
    
    #Density Model, sps 1
    D1.intercept <- exp(D.beta01)*cellArea
    lambda1.cell <- InSS*exp(D.beta11*D.cov1)
    pi1.denom <- sum(lambda1.cell)
    pi1.cell <- lambda1.cell/pi1.denom
    lambda1 <- D1.intercept*pi1.denom
    lambda1.cell.plot <- lambda1.cell
    lambda1.cell.plot[lambda1.cell==0] <- -Inf
    image(x.vals,y.vals,matrix(lambda1.cell.plot,length(x.vals),length(y.vals)),
          main="Sps 1: Spatially Explicit D and Realized ACs",xlab="X",ylab="Y",col=cols1)
    points(X,pch=4)
    N1 <- rpois(1,lambda1)
    
    #Activity centers, sps 1
    s1.cell <- sample(1:n.cells,N1,replace=TRUE,prob=pi1.cell)
    # s1 <- matrix(NA,N1,2)
    # for(i in 1:N1){
    #   s.xlim <- dSS[s1.cell[i],1] + c(-res,res)/2
    #   s.ylim <- dSS[s1.cell[i],2] + c(-res,res)/2
    #   s1[i,1] <- runif(1,s.xlim[1],s.xlim[2])
    #   s1[i,2] <- runif(1,s.ylim[1],s.ylim[2])
    # }
    s1 <- dSS[s1.cell,]
    points(s1[,1],s1[,2],pch=16)
    
    #repulsion
    dists <- e2dist(s1,dSS)
    kern <- h0*exp(-dists^2/(2*omega^2))
    # plot(log(1-kern[1,])~dists[1,]) #log(1-kern) exponentially approaches 0 from below with distance
    log.h <- colSums(log(1-kern)) 
    h <- exp(log.h)
    image(x.vals,y.vals,matrix(h,length(x.vals),length(y.vals)),
          main="Interaction Function",xlab="X",ylab="Y",col=cols1)
    points(s1[,1],s1[,2],pch=16)
    
    #Density Model, sps 2
    D2.intercept <- exp(D.beta02)*cellArea
    # lambda2.cell <- InSS*exp(D.beta12*D.cov2)
    lambda2.cell <- InSS*exp(D.beta12*D.cov2)*h
    pi2.denom <- sum(lambda2.cell)
    pi2.cell <- lambda2.cell/pi2.denom
    lambda2 <- D2.intercept*pi2.denom
    lambda2.cell.plot <- lambda2.cell
    lambda2.cell.plot[lambda2.cell==0] <- -Inf
    image(x.vals,y.vals,matrix(lambda2.cell.plot,length(x.vals),length(y.vals)),
          main="Sps 2: Spatially Explicit D and Realized ACs",xlab="X",ylab="Y",col=cols1)
    points(X,pch=4)
    N2 <- rpois(1,lambda2)
    
    #Activity centers, sps 2
    s2.cell <- sample(1:n.cells,N2,replace=TRUE,prob=pi2.cell)
    # s2 <- matrix(NA,N2,2)
    # for(i in 1:N2){
    #   s.xlim <- dSS[s2.cell[i],1] + c(-res,res)/2
    #   s.ylim <- dSS[s2.cell[i],2] + c(-res,res)/2
    #   s2[i,1] <- runif(1,s.xlim[1],s.xlim[2])
    #   s2[i,2] <- runif(1,s.ylim[1],s.ylim[2])
    # }
    s2 <- dSS[s2.cell,]
    points(s2[,1],s2[,2],pch=16)
    
    #plot sps2 over sps1
    image(x.vals,y.vals,matrix(h,length(x.vals),length(y.vals)),
          main="Interaction Function + Both Species ACs",xlab="X",ylab="Y",col=cols1)
    points(X,pch=4)
    points(s1[,1],s1[,2],col="darkred",pch=16)
    points(s2[,1],s2[,2],col="goldenrod",pch=16)
    
    #observation models
    J <- nrow(X)
    
    #simulate detection, sps 1
    D1 <- e2dist(s1,X)
    pd1 <- p01*exp(-D1*D1/(2*sigma1*sigma1))
    y1 <- array(NA,dim=c(N1,J,K))
    for(i in 1:N1){
      for(j in 1:J){
        for(k in 1:K){
          y1[i,j,k] <- rbinom(1,1,pd1[i,j])
        }
      }
    }
    
    #simulate detection, sps 2
    D2 <- e2dist(s2,X)
    pd2 <- p02*exp(-D2*D2/(2*sigma2*sigma2))
    y2 <- array(NA,dim=c(N2,J,K))
    for(i in 1:N2){
      for(j in 1:J){
        for(k in 1:K){
          y2[i,j,k] <- rbinom(1,1,pd2[i,j])
        }
      }
    }
    
    image(x.vals,y.vals,matrix(lambda1.cell.plot,length(x.vals),length(y.vals)),
          main="Sps 1: Expected D + Spatial Recaps",xlab="X",ylab="Y",col=cols1)
    grid(n.cells.x,n.cells.y,lwd=1,lty=1,col="grey80")
    detected1 <- which(rowSums(y1)>0)
    points(s1[,1],s1[,2],pch=16,col="grey30",cex=1.25)
    points(X,pch=4)

    y2D1 <- apply(y1,c(1,2),sum)
    for(i in detected1){
      for(j in 1:J){
        if(y2D1[i,j]>0){
          lines(x=c(s1[i,1],X[j,1]),y=c(s1[i,2],X[j,2]))
        }
      }
    }
    
    image(x.vals,y.vals,matrix(lambda2.cell.plot,length(x.vals),length(y.vals)),
          main="Sps 2: Expected D + Spatial Recaps",xlab="X",ylab="Y",col=cols1)
    grid(n.cells.x,n.cells.y,lwd=1,lty=1,col="grey80")
    detected2 <- which(rowSums(y2)>0)
    points(s2[,1],s2[,2],pch=16,col="goldenrod",cex=1.25)
    points(X,pch=4)
    
    y2D2 <- apply(y2,c(1,2),sum)
    for(i in detected2){
      for(j in 1:J){
        if(y2D2[i,j]>0){
          lines(x=c(s2[i,1],X[j,1]),y=c(s2[i,2],X[j,2]))
        }
      }
    }

    #discard uncaptured inds
    caught1 <- which(rowSums(y1)>0)
    n1 <- length(caught1)
    y1 <- y1[caught1,,]
    
    caught2 <- which(rowSums(y2)>0)
    n2 <- length(caught2)
    y2 <- y2[caught2,,]
    
    constants <- list(D.cov1=D.cov1,D.cov2=D.cov2,
                      X=X,J=J,K=K,xlim=xlim,ylim=ylim,dSS=dSS,res=res,cells=cells,x.vals=x.vals,y.vals=y.vals,
                      InSS=InSS,n.cells=n.cells,n.cells.x=n.cells.x,n.cells.y=n.cells.y)
    truth <- list(lambda1=lambda1,lambda1.cell=lambda1.cell,N1=N1,s1=s1,s1.cell=s1.cell,n1=n1,
                  lambda2=lambda2,lambda2.cell=lambda2.cell,N2=N2,s2=s2,s2.cell=s2.cell,n2=n2)
    capture <- list(y1=y1,y2=y2)
    out <- list(constants=constants,truth=truth,capture=capture,seed=seed)
    return(out)
  }
