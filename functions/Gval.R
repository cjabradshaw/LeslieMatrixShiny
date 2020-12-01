# mean generation time
Gval <- function (x,age.max) ## where x is a Leslie Matrix
{	
  G <- (log(Rval(x,age.max)))/(log(Re((eigen(x)$values)[1])))
  return(G)
}
