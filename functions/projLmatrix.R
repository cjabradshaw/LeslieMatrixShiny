# project a Leslie matrix
# where x = the projection matrix
#       initN = initial population size (females and males)
#       projYrs = years to project into the future
#       projGens = generations to project into the future
#       DFparams = density-feedback parameters from a logistic power function y = DFa/(1+(x/DFb)^DFc)
projLmatrix <- function (x, initN, projYrs, projGens, DFparams = c(DFa, DFb, DFc))
{	
  
    ## matrix dimension
    dimL <- dim(x)[1]
  
    prjyrs <- ifelse(missing(projGens), projYrs, round(projGens * Gval(x,dimL), 0))
    
    ## initial population vector
    initvec <- StableStageDist(x) * initN #initial population vector
    
    ## set time limit for projection in 1-yr increments
    yrNow <- 0
    #************************
    yrEnd <- yrNow + (prjyrs+1) # set projection end date
    #************************
    t <- (yrEnd - yrNow)
    
    totF <- sum(x[1,])
    popmat <- x # resets matrix to original
    yrvec <- seq(yrNow,yrEnd)
    
    ## set population storage matrices
    nmat <- matrix(0, nrow=dim(x)[1],ncol=(t+1)) # empty matrix
    nmat[,1] <- initvec #fill first matrix column with initial population vector

    # check optionals
    if (missing(DFparams)) {
    
    ## set up projection loop
      for (i in 1:t) {
        nmat[,i+1] <- popmat %*% nmat[,i]
      }
    
  } else {

    ## set up projection loop
    survOrig <- c(diag(popmat[2:(dim(x)[1]),]), popmat[(dim(x)[1]),(dim(x)[1])])
    for (i in 1:t) {
      predRed <- DFparams[1]/(1+((sum(nmat[,i]))/DFparams[2])^DFparams[3])
      survUpdate <- survOrig * predRed
      diag(popmat[2:(dim(x)[1]),]) <- survUpdate[1:(dim(x)[1]-1)]
      popmat[(dim(x)[1]),(dim(x)[1])] <- survUpdate[(dim(x)[1])]
      nmat[,i+1] <- popmat %*% nmat[,i]
      }
  } # end else
  
    npred <- colSums(nmat)
    yrs <- seq(yrNow, yrEnd, 1)
    projOut <- data.frame(yrs,npred)
    
    return(projOut)
}
