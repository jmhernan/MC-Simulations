getICC <- function(y,x) {
  groupmeans <- tapply(y,x,mean)
  nj         <- tapply(y,x,length)
  ngroups    <- length(groupmeans)
  n          <- length(y)
  
  ntilde <- (n-sum(nj^2/n)/(ngroups-1))
  
  s2j       <- tapply(y,x,var)
  s2within  <- sum((nj-1)*s2j)/(n-ngroups)
  s2between <- sum(nj*(groupmeans-mean(y))^2)/(ntilde*(ngroups-1))
  
  sigma2y     <- s2within
  sigma2alpha <- s2between-s2within/ntilde
  
  Fvalue <- ntilde*s2between/s2within
  pvalue <- 1-pf(Fvalue,ngroups-1,n-ngroups)
  
  icc <- ifelse(Fvalue<1,0,(Fvalue-1)/(Fvalue+ntilde-1))
  
  return(list(sigma2y = sigma2y,sigma2alpha = sigma2alpha,icc = icc,Fvalue-Fvalue,
              pvalue = pvalue))
}

r.cor.matrix <- function(a,b,c,d,e,f,g,h,i,j) {
  R <-  matrix(c(1,a,b,c,d,
                 a,1,e,f,g,
                 b,e,1,h,i,
                 c,f,h,1,j,
                 d,g,i,j,1), nrow = 5,byrow = TRUE)
  return(R)
}