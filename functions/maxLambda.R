## maximum lambda function
maxLambda <- function(x) Re((eigen(x)$values)[1]) ## where 'x' is a Leslie matrix
