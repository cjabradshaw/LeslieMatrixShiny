# create pre-breeding-census Leslie matrix
# where age.max = maximum age (longevity)
#       Svec = vector of age-specific survival probabilities (must have length = age.max)
#       Fvec = vector of age-specific fertilities (must have length = age.max)
#       finalStage = default is final survival probability from Svec, or 'abrupt' if the final stage kills all remaining individuals
createLmatrix <- function (age.max, Svec, Fvec, finalStage="abrupt") ## where x is a Leslie Matrix
{	
  Lmat <- matrix(data = 0, nrow=age.max+1, ncol=age.max+1)
  diag(Lmat[2:(age.max+1),]) <- Svec[-(age.max+1)]
  Lmat[age.max+1,age.max+1] <- ifelse(finalStage == "abrupt", 0, Svec[age.max+1])
  Lmat[1,] <- Fvec

  return(Lmat)
}
