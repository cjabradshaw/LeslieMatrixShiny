## number of female offspring produced per female during lifetime
Rval <- function(x,age.max) ## reproductive value (R0) where x = Leslie matrix; age.max = maximum age of females
{		
  ## define the transition matrix
  T <- x[1:age.max,1:age.max]
  T[1,1:(age.max)] <- 0
  
  ## define the fertility matrix
  F <- x[1:age.max,1:age.max]
  diag(F[2:age.max,1:(age.max-1)]) <- 0
  
  ## define the identity matrix
  I <- matrix(data<-0,nrow<-age.max,ncol<-age.max)
  diag(I) <- 1
  
  ## define the fundamental matrix
  library(MASS)
  N.fund <- ginv(I-T)
  
  ## define the reproductive matrix
  R <- F %*% N.fund
  
  ## define R0 (number of female offspring produced per female during lifetime)
  R0 <- Re((eigen(R)$values)[1])
  
  ## output
  return(R0)
}
