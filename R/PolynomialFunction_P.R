P <- function(x, deg){
  na <- ncol(x)
  avec <- replicate(deg, 1:na, simplify=FALSE)
  aposna <- as.matrix(expand.grid(avec))
  result <- matrix(1, nrow=nrow(x), ncol=na^deg)
  for (s in 1:(na^deg)) {
    for(m in 1:deg){
      result[,s] <- result[,s] * x[,aposna[s,m]]
    }
  }
  qq <- list(x=deparse(substitute(x)), deg=deg)
  attributes(result) <- c(attributes(result), qq)
  result
}
