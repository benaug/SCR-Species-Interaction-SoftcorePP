getd2 <- nimbleFunction(
  run = function(s = double(1),dSS = double(2),z = double(0)){ 
    returnType(double(1))
    n.cells <- nimDim(dSS)[1]
    if(z==1){
      d2 <- (s[1]-dSS[1:n.cells,1])^2 + (s[2]-dSS[1:n.cells,2])^2
    }else{
      d2 <- rep(0,n.cells)
    }
    return(d2)
  }
)

getIntKern <- nimbleFunction(
  run = function(d2 = double(1),z = double(0),omega = double(0)){ 
    returnType(double(1))
    n.cells <- nimDim(d2)[1]
    if(z==1){
      IntKern <- exp(-d2/(2*omega^2)) #impossible for species 2 to live in exact location as species 1
    }else{
      IntKern <- rep(0,n.cells)
    }
    return(IntKern)
  }
)

getH <- nimbleFunction(
  run = function(IntKern = double(2), z = double(1)){ 
    returnType(double(1))
    M <- nimDim(IntKern)[1]
    n.cells <- nimDim(IntKern)[2]
    log_h <- rep(0,n.cells)
    for(i in 1:M){
      if(z[i]==1){
        log_h <- log_h + log(1-IntKern[i,1:n.cells])
      }
    }
    return(exp(log_h))
  }
)

GetDetectionProb <- nimbleFunction(
  run = function(s = double(1), p0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
      d2 <- (s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2
      ans <- p0*exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)

dBinomialVector <- nimbleFunction(
  run = function(x = double(1), size = double(1), prob = double(1), z = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual (never occurs with all known IDs)
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      logProb <- sum(dbinom(x, size=size, prob=prob, log = TRUE))
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rBinomialVector <- nimbleFunction(
  run = function(n = integer(0),size = double(1), prob = double(1), z = double(0)) {
    returnType(double(1))
    J <- nimDim(prob)[1]
    out <- numeric(J,value=0)
    return(out)
  }
)

#Required custom update for N/z
zSampler1 <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    M <- control$M
    J <- control$J
    dSS <- control$dSS
    n.cells <- control$n.cells
    ind.detected <- control$ind.detected
    z.ups <- control$z.ups
    sps2.nodes <- control$sps2.nodes
    sps2.lik.nodes <- control$sps2.lik.nodes
    d2.nodes <- control$d2.nodes
    IntKern.nodes <- control$IntKern.nodes
    y.nodes <- control$y.nodes
    pd.nodes <- control$pd.nodes
    N.node <- control$N.node
    z.nodes <- control$z.nodes
    calcNodes <- control$calcNodes
  },
  run = function() {
    #track these "manually" so computations faster than nimble will do them
    h.initial <- model$h
    log.h.initial <- log(h.initial)
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected call
      if(updown==0){#subtract
        #find all z's currently on
        z.on <- which(model$z1==1)
        n.z.on <- length(z.on)
        pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
        pick <- z.on[pick]
        #prereject turning off individuals currently allocated samples
        if(ind.detected[pick]==1){#is this an individual with captured?
          reject <- TRUE
        }
        if(!reject){
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.sps2 <- model$getLogProb(sps2.lik.nodes)
          lp.initial.y <- model$getLogProb(y.nodes[pick])

          #propose new N/z
          model$N1[1] <<-  model$N1[1] - 1
          model$z1[pick] <<- 0
          
          #turn pd off
          model$calculate(pd.nodes[pick])

          #sps2 stuff
          IntKern.initial <- model$IntKern[pick,] #pull out initial
          model$calculate(d2.nodes[pick])
          model$calculate(IntKern.nodes[pick])
          #need to deal with underflow and subtracting -Inf from -Inf that produces NA
          #this happens for s1.cell[pick] because distance between AC and this grid cell is 0
          #this guy (always?) zeros out a cell. need to recompute it if we turn him off
          tmp <- log(1 - IntKern.initial)
          # tmp[model$s1.cell[pick]] <- sum(log(1-model$IntKern[,model$s1.cell[pick]]))
          # log.h.proposed <- log.h.initial - log(1 - tmp) #subtract this individual out
          log.h.proposed <- log.h.initial - tmp #subtract this individual out
          #recompute this entire cell
          log.h.proposed[model$s1.cell[pick]] <- sum(log(1-model$IntKern[,model$s1.cell[pick]]))
          h.proposed <- exp(log.h.proposed)
          model$h <<- h.proposed #put into model to skip resumming over all individuals
          model$calculate(sps2.nodes) #update sps2 nodes

          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.sps2 <- model$calculate(sps2.lik.nodes)
          lp.proposed.y <- model$calculate(y.nodes[pick]) #will always be 0

          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y + lp.proposed.sps2) -
            (lp.initial.N + lp.initial.y + lp.initial.sps2)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N1",1][1] <<- model[["N1"]]
            mvSaved["h",1][1:n.cells] <<- model[["h"]][1:n.cells]
            mvSaved["lambda2",1][1] <<- model[["lambda2"]]
            mvSaved["lambda2.cell",1][1:n.cells] <<- model[["lambda2.cell"]][1:n.cells]
            mvSaved["pi2.cell",1][1:n.cells] <<- model[["pi2.cell"]][1:n.cells]
            mvSaved["pi2.denom",1][1] <<- model[["pi2.denom"]]
            for(j in 1:J){
              mvSaved["pd1",1][pick,j] <<- model[["pd1"]][pick,j]
            }
            mvSaved["z1",1][pick] <<- model[["z1"]][pick]
            mvSaved["d2",1][pick,1:n.cells] <<- model[["d2"]][pick,1:n.cells]
            mvSaved["IntKern",1][pick,1:n.cells] <<- model[["IntKern"]][pick,1:n.cells]
            log.h.initial <- log.h.proposed
            h.initial <- h.proposed
          }else{
            model[["N1"]] <<- mvSaved["N1",1][1]
            model[["h"]][1:n.cells] <<- mvSaved["h",1][1:n.cells]
            model[["lambda2"]] <<- mvSaved["lambda2",1][1]
            model[["lambda2.cell"]][1:n.cells] <<- mvSaved["lambda2.cell",1][1:n.cells]
            model[["pi2.cell"]][1:n.cells] <<- mvSaved["pi2.cell",1][1:n.cells]
            model[["pi2.denom"]] <<- mvSaved["pi2.denom",1][1]
            for(j in 1:J){
              model[["pd1"]][pick,j] <<- mvSaved["pd1",1][pick,j]
            }
            model[["z1"]][pick] <<- mvSaved["z1",1][pick]
            model[["d2"]][pick,1:n.cells] <<- mvSaved["d2",1][pick,1:n.cells]
            model[["IntKern"]][pick,1:n.cells] <<- mvSaved["IntKern",1][pick,1:n.cells]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
            model$calculate(sps2.lik.nodes)
          }
        }
      }else{#add
        if(model$N1[1] < M){ #cannot update if z maxed out. Need to raise M
          z.off <- which(model$z1==0)
          n.z.off <- length(z.off)
          pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
          pick <- z.off[pick]
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.sps2 <- model$getLogProb(sps2.lik.nodes)
          lp.initial.y <- model$getLogProb(y.nodes[pick]) #will always be 0

          #propose new N/z
          model$N1[1] <<-  model$N1[1] + 1
          model$z1[pick] <<- 1
          
          #turn pd on
          model$calculate(pd.nodes[pick])

          #sps2 stuff
          model$calculate(d2.nodes[pick])
          model$calculate(IntKern.nodes[pick])
          log.h.proposed <- log.h.initial + log(1 - model$IntKern[pick,]) #add this individual in
          h.proposed <- exp(log.h.proposed)
          model$h <<- h.proposed #put into model to skip resumming over all individuals
          model$calculate(sps2.nodes) #update sps2 nodes

          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.sps2 <- model$calculate(sps2.lik.nodes)
          lp.proposed.y <- model$calculate(y.nodes[pick])

          #MH step
          #if using discrete state space, 
          # will always reject if we propose z1==1 when s1.cell is in any cell that s2.cell occupies
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y + lp.proposed.sps2) -
                                      (lp.initial.N + lp.initial.y + lp.initial.sps2)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N1",1][1] <<- model[["N1"]]
            mvSaved["h",1][1:n.cells] <<- model[["h"]][1:n.cells]
            mvSaved["lambda2",1][1] <<- model[["lambda2"]]
            mvSaved["lambda2.cell",1][1:n.cells] <<- model[["lambda2.cell"]][1:n.cells]
            mvSaved["pi2.cell",1][1:n.cells] <<- model[["pi2.cell"]][1:n.cells]
            mvSaved["pi2.denom",1][1] <<- model[["pi2.denom"]]
            for(j in 1:J){
              mvSaved["pd1",1][pick,j] <<- model[["pd1"]][pick,j]
            }
            mvSaved["z1",1][pick] <<- model[["z1"]][pick]
            mvSaved["d2",1][pick,1:n.cells] <<- model[["d2"]][pick,1:n.cells]
            mvSaved["IntKern",1][pick,1:n.cells] <<- model[["IntKern"]][pick,1:n.cells]
            log.h.initial <- log.h.proposed
            h.initial <- h.proposed
          }else{
            model[["N1"]] <<- mvSaved["N1",1][1]
            model[["h"]][1:n.cells] <<- mvSaved["h",1][1:n.cells]
            model[["lambda2"]] <<- mvSaved["lambda2",1][1]
            model[["lambda2.cell"]][1:n.cells] <<- mvSaved["lambda2.cell",1][1:n.cells]
            model[["pi2.cell"]][1:n.cells] <<- mvSaved["pi2.cell",1][1:n.cells]
            model[["pi2.denom"]] <<- mvSaved["pi2.denom",1][1]
            for(j in 1:J){
              model[["pd1"]][pick,j] <<- mvSaved["pd1",1][pick,j]
            }
            model[["z1"]][pick] <<- mvSaved["z1",1][pick]
            model[["d2"]][pick,1:n.cells] <<- mvSaved["d2",1][pick,1:n.cells]
            model[["IntKern"]][pick,1:n.cells] <<- mvSaved["IntKern",1][pick,1:n.cells]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
            model$calculate(sps2.lik.nodes)
          }
        }
      }
    }
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)

#Required custom update for N/z
zSampler2 <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    M <- control$M
    J <- control$J
    ind.detected <- control$ind.detected
    z.ups <- control$z.ups
    y.nodes <- control$y.nodes
    pd.nodes <- control$pd.nodes
    N.node <- control$N.node
    z.nodes <- control$z.nodes
    calcNodes <- control$calcNodes
  },
  run = function() {
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected call
      if(updown==0){#subtract
        #find all z's currently on
        z.on <- which(model$z2==1)
        n.z.on <- length(z.on)
        pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
        pick <- z.on[pick]
        #prereject turning off individuals currently allocated samples
        if(ind.detected[pick]==1){#is this an individual with captured?
          reject <- TRUE
        }
        if(!reject){
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick])
          
          #propose new N/z
          model$N2[1] <<-  model$N2[1] - 1
          model$z2[pick] <<- 0
          
          #turn pd off
          model$calculate(pd.nodes[pick])
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick]) #will always be 0
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)
          
          if(accept) {
            mvSaved["N2",1][1] <<- model[["N2"]]
            for(j in 1:J){
              mvSaved["pd2",1][pick,j] <<- model[["pd2"]][pick,j]
            }
            mvSaved["z2",1][pick] <<- model[["z2"]][pick]
          }else{
            model[["N2"]] <<- mvSaved["N2",1][1]
            for(j in 1:J){
               model[["pd2"]][pick,j] <<- mvSaved["pd2",1][pick,j]
            }
            model[["z2"]][pick] <<- mvSaved["z2",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }else{#add
        if(model$N2[1] < M){ #cannot update if z maxed out. Need to raise M
          z.off <- which(model$z2==0)
          n.z.off <- length(z.off)
          pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
          pick <- z.off[pick]
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick]) #will always be 0
          
          #propose new N/z
          model$N2[1] <<-  model$N2[1] + 1
          model$z2[pick] <<- 1
          
          #turn pd on
          model$calculate(pd.nodes[pick])
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick])
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N2",1][1] <<- model[["N2"]]
            for(j in 1:J){
              mvSaved["pd2",1][pick,j] <<- model[["pd2"]][pick,j]
            }
            mvSaved["z2",1][pick] <<- model[["z2"]][pick]
          }else{
            model[["N2"]] <<- mvSaved["N2",1][1]
            for(j in 1:J){
              model[["pd2"]][pick,j] <<- mvSaved["pd2",1][pick,j]
            }
            model[["z2"]][pick] <<- mvSaved["z2",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }
    }
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)