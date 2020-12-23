# logistic power function
# estimate parameters for a logistic power function
# y = a / (1 + (x/b)^c)
# PrPersist = persistence probability
# N = initial population size
# a, b, c = constants
# alpha = alpha probability to calculate confidence intervals
# number of simulations to estimate uncertainty in predicted relationship

LogPowFunc <- function (PrPersist, N, a=0.95, b=37, c=-17, alpha=0.05, sim=10000)
{	
  dat <- data.frame(N, PrPersist)
  param.init <- c(a, b, c)
  fit <- nls(PrPersist ~ a / (1 + (N/b)^c), 
                  data = dat,
                  algorithm = "default",
                  start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                  trace = TRUE, 
                  nls.control(maxiter = 100000, tol = 1e-05, minFactor = 1/2048))
  Nvec <- round(seq(dat$N[1], dat$N[length(dat$N)], by=(dat$N[length(dat$N)] - dat$N[1])/100), 0)
  Npred <- data.frame(Nvec)
  colnames(Npred) <- "N"
  aLP <- coef(fit)[1]
  aLPse  <- summary(fit)$coefficients[4]
  bLP <- coef(fit)[2]
  bLPse  <- summary(fit)$coefficients[5] 
  cLP <- coef(fit)[3]
  cLPse  <- summary(fit)$coefficients[6] 
  
  pred <- predictNLS(object=fit, newdata=Npred, level=(1-alpha), nsim=sim)
  PerPredMed <- pred[,4]
  PerPredUp <- pred[,7]
  PerPredLo <- pred[,6]
  
  MVPasympUp <- round(Nvec[which(round(PerPredLo, 2) == max(round(PerPredLo,2), na.rm=T))[1]], 0)
  MVPasympMed <- round(Nvec[which(round(PerPredMed, 2) == max(round(PerPredMed,2), na.rm=T))[1]], 0)
  MVPasympLo <- round(Nvec[which(round(PerPredUp, 2) == max(round(PerPredUp,2), na.rm=T))[1]], 0)
  PredDat <- data.frame(Nvec,PerPredMed,PerPredUp,PerPredLo)
  
  return(list(PerPredDat=PredDat, MVPasMed=MVPasympMed, MVPasUp=MVPasympUp, MPVasLo=MVPasympLo, aLP=aLP, bLP=bLP, cLP=cLP))
} # end Func