# project a Leslie matrix
# where x = the projection matrix
#       initN = initial population size (females and males)
#       projYrs = years to project into the future
#       projGens = generations to project into the future
#       iter = number of stochastic iterations
#       S_SD = vector of survival standard deviations across ages (%)
#       F_SD = vector of fertility standard deviations across ages (%)
#       Qthres = quasi-extinction threshold
#       DFparams = density-feedback parameters from a logistic power function y = DFa/(1+(x/DFb)^DFc)
#       CatParams = catastrophic mortality parameters (Cmag = magnitude, CMagSD = magnitude standard deviation)
#       PercOff = percentage single offtake (% total N)
#       FixOff = fixed # individuals to offtake
#       IntOff = interval (% of entire projection interval) during which annual offtake happens

projStochPress <- function (x, initN, projYrs, projGens, iter, S_SD, F_SD, Qthresh=50, DFparams = c(DFa, DFb, DFc), CatParams = c(CMag, CMagSD),
                            PercOff=50, FixOff=10, IntOff = c(Int1, Int2))
{	
    ## matrix dimension
    dimL <- dim(x)[1]
    
    pyrs <- ifelse(missing(projGens), projYrs, round(projGens * Gval(x,dimL), 0))
    
    ## press interval offtake parameters
    if(missing(IntOff)) {
      IntOfft <- c(1, pyrs)
    } else {
      IntOfft <- c(round((pyrs*(as.numeric(IntOff[1])/100)), 0), round((pyrs*(as.numeric(IntOff[2])/100)),0))
    }
    
    ## initial population vector
    initvec <- StableStageDist(x) * initN #initial population vector
    
    ## set time limit for projection in 1-yr increments
    yrNow <- 1
    #************************
    yrEnd <- yrNow + (pyrs) # set projection end date
    #************************
    t <- (yrEnd - yrNow)
    
    # check optionals
    if (missing(DFparams) & missing(CatParams)) {
    
      nSumsMat <- matrix(data = 0, nrow = iter, ncol = (t+1)) # N storage matrix
      rMat <- matrix(data = 0, nrow = iter, ncol = t) # r storage matrix
        
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
            
            if (missing(PercOff) & ((i) >= IntOfft[1] & (i) <= IntOfft[2])) {
              nmat[,i+1] <- nmat[,i+1] - (StableStageDist(popmat)*FixOff/2)
              nmat[,i+1] <- ifelse(nmat[,i+1] < 0, 0, nmat[,i+1])
            }
            if (missing(FixOff) & ((i) >= IntOfft[1] & (i) <= IntOfft[2])) {
              nmat[,i+1] <- nmat[,i+1] - (nmat[,i+1]*(PercOff/100))
            }
            
          } # end i loop
        
        nSumsMat[e,] <- as.vector(colSums(nmat))
        rMat[e, ] <- nSumsMat[e,2:yrEnd]/nSumsMat[e,1:t]
        rMat[e, ] <- ifelse(is.infinite(rMat[e, ])==T, NA, rMat[e, ])
        
      } # end e loop
      
      minNvec <- apply(nSumsMat, MARGIN=1, min, na.rm=T)
      minMnN <- mean(minNvec, na.rm=T)
      minLoN <- quantile(minNvec, probs=0.025, na.rm=T)
      minUpN <- quantile(minNvec, probs=0.975, na.rm=T)
      Qext <- sum(ifelse(round(minNvec, 0) < Qthresh, 1, 0)) / iter
  
      mnRVec <- apply(rMat, MARGIN=1, mean, na.rm=T)
      mnMnR <- mean(mnRVec, na.rm=T)
      mnLoR <- quantile(mnRVec, probs=0.025, na.rm=T)
      mnUpR <- quantile(mnRVec, probs=0.975, na.rm=T)
      
      Nmd <- apply(nSumsMat, MARGIN=2, mean, na.rm=T) # median over all iterations
      Nup <- apply(nSumsMat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
      Nlo <- apply(nSumsMat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
      rMn <- apply(rMat, MARGIN=2, mean, na.rm=T)
      rUp <- apply(rMat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
      rLo <- apply(rMat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
      
      Nout <- data.frame(yrvec,Nlo,Nmd,Nup)
      Rout <- data.frame(yrvec[2:yrEnd],rLo,rMn,rUp)
      colnames(Rout)[1] <- "yrvec" 
      
      return(list(Nrange=Nout, rRange=Rout, minMnN=minMnN, minLoN=minLoN, minUpN=minUpN, mnRmn=mnMnR, mnRlo=mnLoR, mnRup=mnUpR, PrQext=Qext))
    } # end missing if
    
    if (missing(CatParams)) {

      nSumsMat <- matrix(data = 0, nrow = iter, ncol = (t+1)) # N storage matrix
      rMat <- matrix(data = 0, nrow = iter, ncol = t) # r storage matrix
      
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
          
          if (missing(PercOff) & ((i) >= IntOfft[1] & (i) <= IntOfft[2])) {
            nmat[,i+1] <- nmat[,i+1] - (StableStageDist(popmat)*FixOff/2)
            nmat[,i+1] <- ifelse(nmat[,i+1] < 0, 0, nmat[,i+1])
          }
          if (missing(FixOff) & ((i) >= IntOfft[1] & (i) <= IntOfft[2])) {
            nmat[,i+1] <- nmat[,i+1] - (nmat[,i+1]*(PercOff/100))
          }
          
        } # end i loop
        
        nSumsMat[e,] <- as.vector(colSums(nmat))
        rMat[e, ] <- nSumsMat[e,2:yrEnd]/nSumsMat[e,1:t]
        rMat[e, ] <- ifelse(is.infinite(rMat[e, ])==T, NA, rMat[e, ])
        
      } # end e loop
      
      minNvec <- apply(nSumsMat, MARGIN=1, min, na.rm=T)
      minMnN <- mean(minNvec, na.rm=T)
      minLoN <- quantile(minNvec, probs=0.025, na.rm=T)
      minUpN <- quantile(minNvec, probs=0.975, na.rm=T)
      Qext <- sum(ifelse(round(minNvec, 0) < Qthresh, 1, 0)) / iter

      mnRVec <- apply(rMat, MARGIN=1, mean, na.rm=T)
      mnMnR <- mean(mnRVec, na.rm=T)
      mnLoR <- quantile(mnRVec, probs=0.025, na.rm=T)
      mnUpR <- quantile(mnRVec, probs=0.975, na.rm=T)
      
      Nmd <- apply(nSumsMat, MARGIN=2, mean, na.rm=T) # median over all iterations
      Nup <- apply(nSumsMat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
      Nlo <- apply(nSumsMat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

      rMn <- apply(rMat, MARGIN=2, mean, na.rm=T)
      rUp <- apply(rMat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
      rLo <- apply(rMat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

      Nout <- data.frame(yrvec,Nlo,Nmd,Nup)
      Rout <- data.frame(yrvec[2:yrEnd],rLo,rMn,rUp)
      colnames(Rout)[1] <- "yrvec" 
      
      return(list(Nrange=Nout, rRange=Rout, minMnN=minMnN, minLoN=minLoN, minUpN=minUpN, mnRmn=mnMnR, mnRlo=mnLoR, mnRup=mnUpR, PrQext=Qext))
      
    } else { # end missing if
      
      genL <- Gval(x,dimL) # mean generation length
      prCat <- 0.14/genL # probability of catastrophe (Reed et al. 2003)

      nSumsMat <- matrix(data = 0, nrow = iter, ncol = (t+1)) # N storage matrix
      rMat <- matrix(data = 0, nrow = iter, ncol = t) # r storage matrix
      
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
          catMortAlpha <- estBetaParams((CatParams[1]/100), ((CatParams[2]/100)^2))$alpha
          catMortBeta <- estBetaParams((CatParams[1]/100), ((CatParams[2]/100)^2))$beta
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
          
          if (missing(PercOff) & ((i) >= IntOfft[1] & (i) <= IntOfft[2])) {
            nmat[,i+1] <- nmat[,i+1] - (StableStageDist(popmat)*FixOff/2)
            nmat[,i+1] <- ifelse(nmat[,i+1] < 0, 0, nmat[,i+1])
          }
          if (missing(FixOff) & ((i) >= IntOfft[1] & (i) <= IntOfft[2])) {
            nmat[,i+1] <- nmat[,i+1] - (nmat[,i+1]*(PercOff/100))
          }
          
        } # end i loop
        
        nSumsMat[e,] <- as.vector(colSums(nmat))
        rMat[e, ] <- log(nSumsMat[e,2:yrEnd]/nSumsMat[e,1:t])
        rMat[e, ] <- ifelse(is.infinite(rMat[e, ])==T, NA, rMat[e, ])
        
      } # end e loop
      
      minNvec <- apply(nSumsMat, MARGIN=1, min, na.rm=T)
      minMnN <- mean(minNvec, na.rm=T)
      minLoN <- quantile(minNvec, probs=0.025, na.rm=T)
      minUpN <- quantile(minNvec, probs=0.975, na.rm=T)
      Qext <- sum(ifelse(round(minNvec, 0) < Qthresh, 1, 0)) / iter
      
      mnRVec <- apply(rMat, MARGIN=1, mean, na.rm=T)
      mnMnR <- mean(mnRVec, na.rm=T)
      mnLoR <- quantile(mnRVec, probs=0.025, na.rm=T)
      mnUpR <- quantile(mnRVec, probs=0.975, na.rm=T)
      
      Nmd <- apply(nSumsMat, MARGIN=2, mean, na.rm=T) # median over all iterations
      Nup <- apply(nSumsMat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
      Nlo <- apply(nSumsMat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

      rMn <- apply(rMat, MARGIN=2, mean, na.rm=T)
      rUp <- apply(rMat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
      rLo <- apply(rMat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
      
      Nout <- data.frame(yrvec,Nlo,Nmd,Nup)
      Rout <- data.frame(yrvec[2:yrEnd],rLo,rMn,rUp)
      colnames(Rout)[1] <- "yrvec" 
      
      return(list(Nrange=Nout, rRange=Rout, minMnN=minMnN, minLoN=minLoN, minUpN=minUpN, mnRmn=mnMnR, mnRlo=mnLoR, mnRup=mnUpR, PrQext=Qext))
      
  } # end else
} # end function
