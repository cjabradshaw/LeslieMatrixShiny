## maximum r function
maxR <- function(x) log(Re((eigen(x)$values)[1])) ## where 'x' is a Leslie matrix
