# where x = the projection matrix
#       initN = initial population size (females and males)
#       K = carrying capacity (long-term average)
#       SmodMax = survival modifier at low population densities (generally = 1)
#       SmodMin = survival modifier at high population densities (near K)

DcompFunc <- function (initN, K, SmodMax=1, SmodMin)
{	
  Kvec <- c(1,initN/2, initN, 0.75*K, K)
  redvec <- c(SmodMax, 0.99*SmodMax, 0.95*SmodMax, SmodMin + 0.5*(SmodMax - SmodMin), SmodMin)
  KredDat <- data.frame(Kvec,redvec)
  
  # logistic power function a/(1+(x/b)^c) #fits logistic power function to population relative to carry capacity, K
  paramInit <- c(1, 5*initN, 2.5)
  fitlp <- nls(redvec ~ a/(1+(Kvec/b)^c), 
                data = KredDat,
                algorithm = "port",
                start = c(a = paramInit[1], b = paramInit[2], c = paramInit[3]),
                trace = TRUE,      
                nls.control(maxiter = 10000, tol = 1e-05, minFactor = 1/1024))
  fitlpsumm <- summary(fitlp)
  KvecCont <- seq(1, K, K/100)
  predlpfx <- coef(fitlp)[1]/(1+(KvecCont/coef(fitlp)[2])^coef(fitlp)[3])

  CompOut <- data.frame(KvecCont, predlpfx)
  
  alp <- coef(fitlp)[1]
  blp <- coef(fitlp)[2]
  clp <- coef(fitlp)[3]

  return(list(DcompOut = CompOut, a = alp, b = blp, c = clp))
}
