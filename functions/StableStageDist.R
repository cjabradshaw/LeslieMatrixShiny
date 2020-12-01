## stable stage distribution
StableStageDist <- function(x) ((x %*% (Re((eigen(x)$vectors)[,1])))/(sum((x %*% (Re((eigen(x)$vectors)[,1]))))))[,1]
