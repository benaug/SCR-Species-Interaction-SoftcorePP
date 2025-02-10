getCellR <- function(u,res,cells,xlim,ylim){
  inout <- 1*(u[1]>xlim[1]&u[1]<xlim[2]&u[2]>ylim[1]&u[2]<ylim[2])
  if(inout==1){
    this.cell.init <- cells[trunc(u[1]/res)+1,trunc(u[2]/res)+1]
  }else{
    this.cell.init <- 0
  }
  return(this.cell.init)
}

init.data.2species.softcore <- function(data=NA,M1=NA,M2=NA,inits=NA,initTrue=FALSE){
  if(initTrue){
    truth <- data$truth
  }
  data <- c(data$constants,data$capture) #restructure data list
  J <- data$J
  K <- data$K
  xlim <- data$xlim
  ylim <- data$ylim
  cells <- data$cells
  dSS <- data$dSS
  InSS <- data$InSS
  res <- data$res
  X <- data$X
  
  #get some inits, actually sigma is all we need
  sigma1 <- inits$sigma1
  sigma2 <- inits$sigma2
  
  #sps 1
  n1 <- nrow(data$y1)
  y1 <- array(0,dim=c(M1,J,K))
  y1[1:n1,,] <- data$y1
  y2D1 <- apply(y1,c(1,2),sum)
  #sps 2
  n2 <- nrow(data$y2)
  y2 <- array(0,dim=c(M2,J,K))
  y2[1:n2,,] <- data$y2
  y2D2 <- apply(y2,c(1,2),sum)
  
  if(initTrue){
    z1.init <- rep(0,M1)
    z1.init[1:truth$N1] <- 1
    s1.init <- cbind(runif(M1,xlim[1],xlim[2]),runif(M1,ylim[1],ylim[2]))
    s1.init[1:truth$N1,] <- truth$s1
    s1.cell.init <- rep(0,M1)
    for(i in 1:M1){
      s1.cell.init[i] <- getCellR(s1.init[i,],res,cells,xlim,ylim)
    }
    
    z2.init <- rep(0,M2)
    z2.init[1:truth$N2] <- 1
    s2.init <- cbind(runif(M2,xlim[1],xlim[2]),runif(M2,ylim[1],ylim[2]))
    s2.init[1:truth$N2,] <- truth$s2
  }else{
    #Initialize z, just using observed z's
    z1.init <- c(rep(1,n1),rep(0,M1-n1))
    s1.init <- cbind(runif(M1,xlim[1],xlim[2]),runif(M1,ylim[1],ylim[2]))
    for(i in 1:M1){
      if(sum(y1[i,,])>0){#if captured
        trapcaps <- which(y2D1[i,]>0)
        if(length(trapcaps)>1){
          s1.init[i,] <- colMeans(X[trapcaps,]) + 0.000001
        }else{
          s1.init[i,] <- X[trapcaps,] + 0.000001
        }
      }
    }
    #move any initialized outside state space
    s1.cell.init <- rep(NA,M1)
    for(i in 1:M1){
      s1.cell.init[i] <- getCellR(s1.init[i,],res,cells,xlim,ylim)
      if(InSS[s1.cell.init[i]]==0){#not in SS, move to nearest cell
        dists1 <- sqrt((dSS[s1.cell.init[i],1]-dSS[,1])^2+(dSS[s1.cell.init[i],2]-dSS[,2])^2)
        dists1[InSS==0] <- Inf
        pick <- which(dists1==min(dists1))[1] #if more than 1, just use first
        s1.init[i,] <- dSS[pick,]
        s1.cell.init[i] <- getCellR(s1.init[i,],res,cells,xlim,ylim)
      }
    }
    #sps 2
    #Initialize z, just using observed z's
    z2.init <- c(rep(1,n2),rep(0,M2-n2))
    s2.init <- cbind(runif(M2,xlim[1],xlim[2]),runif(M2,ylim[1],ylim[2]))
    for(i in 1:M2){
      if(sum(y2[i,,])>0){#if captured
        trapcaps <- which(y2D2[i,]>0)
        if(length(trapcaps)>1){
          s2.init[i,] <- colMeans(X[trapcaps,])
        }else{
          s2.init[i,] <- X[trapcaps,]
        }
      }
    }
  }
  #move any initialized outside state space
  #also cannot be in s1 cells with z1=1!
  s2.cell.init <- rep(NA,M2)
  InSS2 <- InSS
  InSS2[unique(s1.cell.init[z1.init==1])] <- 0
  for(i in 1:M2){
    s2.cell.init[i] <- getCellR(s2.init[i,],res,cells,xlim,ylim)
    if(InSS2[s2.cell.init[i]]==0){#not in SS, move to nearest cell
      dists2 <- sqrt((dSS[s2.cell.init[i],1]-dSS[,1])^2+(dSS[s2.cell.init[i],2]-dSS[,2])^2)
      dists2[InSS2==0] <- Inf
      pick <- which(dists2==min(dists2))[1] #if more than 1, just use first
      s2.init[i,] <- dSS[pick,]
      s2.cell.init[i] <- getCellR(s2.init[i,],res,cells,xlim,ylim)
    }
  }
  return(list(y1=y1,y2D1=y2D1,z1=z1.init,s1=s1.init,s1.cell=s1.cell.init,N1=sum(z1.init),
              y2=y2,y2D2=y2D2,z2=z2.init,s2=s2.init,s2.cell=s2.cell.init,N2=sum(z2.init)))
}