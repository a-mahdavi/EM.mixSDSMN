 dSDSSEN <- function(x, mu = rep(0,d), Sigma, nu , la){
  if(missing(Sigma))
    stop("Sigma is missing")
  if(sum(nu < 0)>0)
    stop("nu must be greater than, or equal to, 0")

  if(is.matrix(Sigma)){
    d <- ncol(Sigma)
    q <- nrow(x)
  }

  if(!is.matrix(Sigma)){
    d <- 1
    q <- length(x)
  }

  if(is.vector(x)){
    x <- matrix(x,length(x),1)
    Sigma <- matrix(Sigma,nrow=d,ncol=d)
  }

  PDF <- NULL

  Sigmainv <- solve(Sigma)
  const1   <- 2^d/(sqrt(det(Sigma))) * prod(nu) * prod(exp(nu))

  if(d==1){

    Dt       <- nu

    fpass <- function(z){

      Dy     <- z-mu
      P      <- Dy %*% Sigmainv %*% Dy
      Pt     <- P + 2*Dt
      Ptinv  <- solve(Pt)
      const2 <- sqrt(det(Ptinv))

      M <- as.vector(mnormt::recintab(rep(2,d), a=rep(1,d), b=rep(Inf,d), mu=rep(0,d), S=Ptinv))

      PDFpass <- const1 * const2 * M[length(M)]

      return(2*PDFpass*pnorm(la*Sigma^(-1/2)*(z-mu)))

    }
  }

  if(d>1&&d<7){

    Dt       <- diag(nu)

    fpass <- function(z){

      Dy     <- diag(z-mu)
      P      <- Dy %*% Sigmainv %*% Dy
      Pt     <- P + 2*Dt
      Ptinv  <- solve(Pt)
      const2 <- sqrt(det(Ptinv))

      M <- as.vector(mnormt::recintab(rep(2,d), a=rep(1,d), b=rep(Inf,d), mu=rep(0,d), S=Ptinv))

      PDFpass <- const1 * const2 * M[length(M)]
      if(PDFpass < 1.0e-323) PDFpass <- 1.0e-323
      return(2*PDFpass*pnorm(la%*%diag(diag(Sigma)^(-1/2))%*%(z-mu)))
    }
  }
if (d>6){
  fpass <- function(z){
      Dt       <- diag(nu)
      Dy     <- diag(z-mu)
      P      <- Dy %*% Sigmainv %*% Dy
      Pt     <- P + 2*Dt
      Ptinv  <- solve(Pt)
      const2 <- sqrt(det(Ptinv))
Ptinv[upper.tri(Ptinv)] <- t(Ptinv)[upper.tri(Ptinv)]
      pez     <- prod(diag(relliptical::mvtelliptical(lower=rep(1,d),upper=rep(Inf,d),n=10000,mu=rep(0,d),Sigma=Ptinv,dist="Normal")$EYY))*mvtnorm::pmvnorm(lower=rep(1,d),upper=rep(Inf,d), mean=rep(0,d),sigma=Ptinv)[1]
      PDFpass <- const1 * const2 * pez
      if(PDFpass < 1.0e-323) PDFpass <- 1.0e-323
      return(2*PDFpass*pnorm(la%*%diag(diag(Sigma)^(-1/2))%*%(z-mu)))
    } }
  PDF <- sapply(1:q, function(i) fpass(z=x[i,]) )
  return(PDF)
}



