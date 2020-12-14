# calculation of MVP from a stochastic projection of a Leslie matrix
# where x = the projection matrix
#       initN = initial population size (females and males)
#       projYrs = years to project into the future
#       projGens = generations to project into the future
#       maxPrQe = maximum probability of quasi-extinction acceptable to conclude persistence at minimum N
#       Nint = population size interval used to decline initN at each step to determine MVP
#       Nmin = minimum population size to test
#       iter = number of stochastic iterations
#       S_SD = vector of survival standard deviations across ages (%)
#       F_SD = vector of fertility standard deviations across ages (%)
#       Qthresh = quasi-extinction threshold
#       DFparams = density-feedback parameters from a logistic power function y = DFa/(1+(x/DFb)^DFc)
#       CatParams = catastrophic mortality parameters (Cmag = magnitude, CMagSD = magnitude standard deviation)

StochMVP <- function (x, initN, projYrs, projGens, maxPrQe, Nint, Nmin, iter, S_SD, F_SD, Qthresh=50, DFparams = c(DFa, DFb, DFc), CatParams = c(CMag, CMagSD))
{	
    
    ## matrix dimension
    dimL <- dim(x)[1]
  
    pyrs <- ifelse(missing(projGens), projYrs, round(projGens * Gval(x,dimL), 0))
    
    ## set N-reduction vector
    Ninit.vec <- seq(from=initN, to=Nmin, by=-Nint) 
    
    ## set time limit for projection in 1-yr increments
    yrNow <- 1
    #************************
    yrEnd <- yrNow + (pyrs) # set projection end date
    #************************
    t <- (yrEnd - yrNow)
      
    # check optionals
    if (missing(DFparams) & missing(CatParams)) {
      
      Qext <- rep(0,length(Ninit.vec)) # Pr(Qe) storage vector
      
      for (n in 1:length(Ninit.vec)) {
          
        initvec <- StableStageDist(x) * Ninit.vec[n] #initial population vector
        nSumsMat <- matrix(data = 0, nrow = iter, ncol = (t+1)) # N storage matrix

        for (e in 1:iter) {
          
          popmat <- x # resets matrix to original
          yrvec <- seq(yrNow,yrEnd)
          
          ## set population storage matrices
          nmat <- matrix(0, nrow=dimL,ncol=(t+1)) # empty matrix
          nmat[,1] <- initvec #fill first matrix column with initial population vector
    
          survOrig <- c(diag(popmat[2:dimL,]))
          survLast <- c(popmat[dimL,dimL])
          fertOrig <- popmat[1,]
          
          ## set up projection loop
            for (i in 1:t) {
              
              survSD <- survOrig * (as.vector(as.numeric(S_SD))[-dimL])/100
              survSD.last <- survLast * (as.vector(as.numeric(S_SD))[dimL])/100
              Salpha <- estBetaParams(survOrig, survSD^2)$alpha
              Sbeta <- estBetaParams(survOrig, survSD^2)$beta
              
              if (survLast == 0) {
                Sstoch <- rbeta(length(Salpha), Salpha, Sbeta)
              }
              if (survSD.last > 0) {
                Salpha.last <- estBetaParams(survLast, survSD.last^2)$alpha
                Sbeta.last <- estBetaParams(survLast, survSD.last^2)$beta
                Sstoch.last <- rbeta(1, Salpha.last, Sbeta.last)
                Sstch1 <- rbeta(length(Salpha), Salpha, Sbeta)
                Sstoch <- c(Sstch1, Sstoch.last)
              }
              
              # stochastic fertility sampler (Gaussian)
              fertSD <- fertOrig * as.numeric(F_SD)/100
              fertStch <- rnorm(length(popmat[,1]), popmat[1,], fertSD)
              fertStoch <- ifelse(fertStch < 0, 0, fertStch)
              
              popmat[1,] <- fertStoch
              if (dimL > length(Sstoch)){
                diag(popmat[2:dimL,]) <- Sstoch
              }
              if (dimL == length(Sstoch)){
                diag(popmat[2:dimL,]) <- Sstoch[-dimL]
                popmat[dimL,dimL] <- Sstoch[dimL]
              }
              nmat[,i+1] <- popmat %*% nmat[,i]
              nmat[,i+1] <- ifelse(nmat[,i+1] < 1, 0, nmat[,i+1])
              
            } # end i loop
          
          nSumsMat[e,] <- as.vector(colSums(nmat))
          
        } # end e loop
        
        minNvec <- apply(nSumsMat, MARGIN=1, min, na.rm=T)
        Qext[n] <- sum(ifelse(round(minNvec, 0) < Qthresh, 1, 0), na.rm=T) / iter
        
      } # end n loop
      
      MVPn <- Ninit.vec[which(Qext <= maxPrQe)[1]]
      MVPdat <- data.frame(Ninit.vec,Qext)
      MVPoutl.dat <- tso(ts(Qext,frequency=1), maxit.iloop = 100000)
      MVPstep <- Ninit.vec[max(MVPoutl.dat$outliers$time) - 2]
      
      return(list(MVPout=MVPdat, MVP=MVPn, MVPst=MVPstep))
      
    } # end missing if
    
    
    if (missing(CatParams)) {
      
      Qext <- rep(0,length(Ninit.vec)) # Pr(Qe) storage vector
      
      for (n in 1:length(Ninit.vec)) {
        
        initvec <- StableStageDist(x) * Ninit.vec[n] #initial population vector
        nSumsMat <- matrix(data = 0, nrow = iter, ncol = (t+1)) # N storage matrix
  
        for (e in 1:iter) {
          
          popmat <- x # resets matrix to original
          yrvec <- seq(yrNow,yrEnd)
          
          ## set population storage matrices
          nmat <- matrix(0, nrow=dimL,ncol=(t+1)) # empty matrix
          nmat[,1] <- initvec # fill first matrix column with initial population vector
          
          survOrig <- c(diag(popmat[2:dimL,]))
          survLast <- c(popmat[dimL,dimL])
          fertOrig <- popmat[1,]
          
          ## set up projection loop
          for (i in 1:t) {
            
            predRed <- DFparams[1]/(1+((sum(nmat[,i]))/DFparams[2])^DFparams[3]) # density feedback on survival
            
            survSD <- survOrig * (as.vector(as.numeric(S_SD))[-dimL])/100
            survSD.last <- survLast * (as.vector(as.numeric(S_SD))[dimL])/100
            Salpha <- estBetaParams(survOrig, survSD^2)$alpha
            Sbeta <- estBetaParams(survOrig, survSD^2)$beta
            
            if (survLast == 0) {
              Sstoch <- rbeta(length(Salpha), Salpha, Sbeta) * predRed
            }
            if (survSD.last > 0) {
              Salpha.last <- estBetaParams(survLast, survSD.last^2)$alpha
              Sbeta.last <- estBetaParams(survLast, survSD.last^2)$beta
              Sstoch.last <- rbeta(1, Salpha.last, Sbeta.last)
              Sstch1 <- rbeta(length(Salpha), Salpha, Sbeta)
              Sstoch <- c(Sstch1, Sstoch.last) * predRed
            }
            
            # stochastic fertility sampler (Gaussian)
            fertSD <- fertOrig * as.numeric(F_SD)/100
            fertStch <- rnorm(length(popmat[,1]), popmat[1,], fertSD)
            fertStoch <- ifelse(fertStch < 0, 0, fertStch)
            
            popmat[1,] <- fertStoch
            if (dimL > length(Sstoch)){
              diag(popmat[2:dimL,]) <- Sstoch
            }
            if (dimL == length(Sstoch)){
              diag(popmat[2:dimL,]) <- Sstoch[-dimL]
              popmat[dimL,dimL] <- Sstoch[dimL]
            }
            nmat[,i+1] <- popmat %*% nmat[,i]
            nmat[,i+1] <- ifelse(nmat[,i+1] < 1, 0, nmat[,i+1])
            
          } # end i loop
          
          nSumsMat[e,] <- as.vector(colSums(nmat))
          
        } # end e loop
        
        minNvec <- apply(nSumsMat, MARGIN=1, min, na.rm=T)
        Qext[n] <- sum(ifelse(round(minNvec, 0) < Qthresh, 1, 0), na.rm=T) / iter

      } # end n loop
      
      MVPn <- Ninit.vec[which(Qext <= maxPrQe)[1]]
      MVPdat <- data.frame(Ninit.vec,Qext)
      MVPoutl.dat <- tso(ts(Qext,frequency=1), maxit.iloop = 100000)
      MVPstep <- Ninit.vec[max(MVPoutl.dat$outliers$time) - 2]
      
      return(list(MVPout=MVPdat, MVP=MVPn, MVPst=MVPstep))
      

    } else { # end missing if
      
      genL <- Gval(x,dimL) # mean generation length
      prCat <- 0.14/genL # probability of catastrophe (Reed et al. 2003)

      Qext <- rep(0,length(Ninit.vec)) # Pr(Qe) storage vector
      
      for (n in 1:length(Ninit.vec)) {
        
        initvec <- StableStageDist(x) * Ninit.vec[n] #initial population vector
        nSumsMat <- matrix(data = 0, nrow = iter, ncol = (t+1)) # N storage matrix

        for (e in 1:iter) {
          
          popmat <- x # resets matrix to original
          yrvec <- seq(yrNow,yrEnd)
          
          ## set population storage matrices
          nmat <- matrix(0, nrow=dimL, ncol=(t+1)) # empty matrix
          nmat[,1] <- initvec #fill first matrix column with initial population vector
          
          survOrig <- c(diag(popmat[2:dimL,]))
          survLast <- c(popmat[dimL,dimL])
          fertOrig <- popmat[1,]
          
          ## set up projection loop
          for (i in 1:t) {
            
            if (missing(DFparams)) {
              predRed <- 1
            } else {
              predRed <- (DFparams[1]/(1+((sum(nmat[,i]))/DFparams[2])^DFparams[3]))
            }
  
            # catastrophic mortality    
            catMortAlpha <- estBetaParams(CatParams[1]/100, (CatParams[2]/100)^2)$alpha
            catMortBeta <- estBetaParams(CatParams[1]/100, (CatParams[2]/100)^2)$beta
            catMortStoch <- rbeta(1, catMortAlpha, catMortBeta)
            catMortStoch <- ifelse(is.na(catMortStoch)==T, 1, catMortStoch)
            catastrophe <- rbinom(1, 1, prCat)
            
            survSD <- survOrig * (as.vector(as.numeric(S_SD))[-dimL])/100
            survSD.last <- survLast * (as.vector(as.numeric(S_SD))[dimL])/100
            Salpha <- estBetaParams(survOrig, survSD^2)$alpha
            Sbeta <- estBetaParams(survOrig, survSD^2)$beta
            
            if (survLast == 0) {
              Sstoch <- rbeta(length(Salpha), Salpha, Sbeta) * predRed * ifelse(catastrophe==1, 1-catMortStoch, 1)
            }
            if (survSD.last > 0) {
              Salpha.last <- estBetaParams(survLast, survSD.last^2)$alpha
              Sbeta.last <- estBetaParams(survLast, survSD.last^2)$beta
              Sstoch.last <- rbeta(1, Salpha.last, Sbeta.last)
              Sstch1 <- rbeta(length(Salpha), Salpha, Sbeta)
              Sstoch <- c(Sstch1, Sstoch.last) * predRed * ifelse(catastrophe==1, 1-catMortStoch, 1)
            }
  
            # stochastic fertility sampler (Gaussian)
            fertSD <- fertOrig * as.vector(as.numeric(F_SD))/100
            fertStch <- rnorm(length(popmat[,1]), popmat[1,], fertSD)
            fertStoch <- ifelse(fertStch < 0, 0, fertStch)
            
            popmat[1,] <- fertStoch
            if (length(Sstoch) < dimL){
              diag(popmat[2:dimL,]) <- Sstoch
            }
            if (dimL == length(Sstoch)){
              diag(popmat[2:dimL,]) <- Sstoch[-dimL]
              popmat[dimL,dimL] <- Sstoch[dimL]
            }
            nmat[,i+1] <- popmat %*% nmat[,i]
            nmat[,i+1] <- ifelse(nmat[,i+1] < 1, 0, nmat[,i+1])
            
          } # end i loop
          
          nSumsMat[e,] <- as.vector(colSums(nmat))

        } # end e loop
        
        minNvec <- apply(nSumsMat, MARGIN=1, min, na.rm=T)
        Qext[n] <- sum(ifelse(round(minNvec, 0) < Qthresh, 1, 0), na.rm=T) / iter

      } # end n loop
      
      MVPQwhich <- which(Qext <= maxPrQe)
      MVPn <- Ninit.vec[MVPQwhich[length(MVPQwhich)]]
      MVPdat <- data.frame(Ninit.vec,Qext)
      MVPoutl.dat <- tso(ts(Qext,frequency=1), maxit.iloop = 100000)
      MVPstep <- Ninit.vec[max(MVPoutl.dat$outliers$time) - 2]
      
      return(list(MVPout=MVPdat, MVP=MVPn, MVPst=MVPstep))
      
  } # end else
} # end function
