dSDSTIN <- function(x, mu = rep(0,d), Sigma, nu,la){
  if(missing(Sigma))
    stop("Sigma is missing")
  if(min(nu) <= 0 || max(nu) >= 1)
    stop("each element in nu must be in (0,1)")
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
  const1   <- 2^d/prod(nu)
  if(d==1){
    fpass <- function(z){
      Dy     <- (z-mu)
      P      <- Dy %*% Sigmainv %*% Dy
      Pinv   <- solve(P)
      const2 <- 1/abs(Dy)
      M   <- as.vector(mnormt::recintab(rep(2,d), a=sqrt(1-nu), b=rep(1,d), mu=rep(0,d), S=Pinv))
      pez <- M[length(M)]
      if(pez < 0 || abs(z-mu)<0.0001){
        const2 <- (2*pi)^(-1/2)/sqrt(Sigma)
        pez    <- 1/3*(1-(1-nu)^(3/2))
      }
      PDFpass <- const1 * const2 * pez
      return(2*PDFpass*pnorm(la*Sigma^(-1/2)*(z-mu)))
    }
  }
  if(d>1&&d<6){
    fpass <- function(z){
      Dy     <- diag(z-mu)
      P      <- Dy %*% Sigmainv %*% Dy
      Pinv   <- MASS::ginv(P)
      const2 <- 1/(det(Sigma)^(1/2)*det(P)^(1/2))
      #print(Pinv)
      M      <- as.vector(mnormt::recintab(rep(2,d), a=sqrt(1-nu), b=rep(1,d), mu=rep(0,d), S=Pinv))
      pez    <- M[length(M)]
      #print(pez)
      if(is.nan(M[length(M)]) | M[length(M)]<0 | sum(abs(z-mu)<0.0003)>0 ){
        ind <- which(abs(z-mu)==min(abs(z-mu)))
        #print(ind)
        Dy     <- diag(z-mu)
        Ps     <- Dy[-ind,-ind] %*% Sigmainv[-ind,-ind] %*% Dy[-ind,-ind]
        Psinv  <- MASS::ginv(Ps)
        nus <- nu[-ind]
        const2 <- 1/(det(Sigma)^(1/2)*det(Ps)^(1/2))
        M       <- as.vector(mnormt::recintab(rep(2,(d - length(ind))), a=sqrt(1-nus), b=rep(1,(d - length(ind))), mu=rep(0,(d - length(ind))), S=Psinv))
        pez1    <- M[length(M)]
        pez2    <- (2*pi)^(-length(ind)/2)*prod(1/3*(1-(1-nu[ind])^(3/2)))
        pez     <- pez1*pez2
        #print("Attention M")
        #pez <- 1e-310
      }
      PDFpass <- const1 * const2 * pez

      PDFmax  <- const1 * 1/det(Sigma)^(1/2) * (2*pi)^(-d/2)*prod(1/3*(1-(1-nu)^(3/2)))

      if(PDFpass > PDFmax)   PDFpass <- PDFmax

PDFpass <- 2*PDFpass*pnorm(la%*%diag(diag(Sigma)^(-1/2))%*%(z-mu))

      if(PDFpass < 1.0e-323) PDFpass <- 1.0e-323

      return(PDFpass)
    }
  }
if (d>5){
  fpass <- function(z){
      Dy     <- diag(1/(z-mu))
      Pinv      <- Dy %*% Sigma %*% Dy
Pinv[upper.tri(Pinv)] <- t(Pinv)[upper.tri(Pinv)]
      const2 <- abs(det(Dy))
      pez     <- prod(diag(relliptical::mvtelliptical(lower=sqrt(1-nu),upper=rep(1,d),n=10000,mu=rep(0,d),Sigma=Pinv,dist="Normal")$EYY))*mvtnorm::pmvnorm(lower=sqrt(1-nu),upper=rep(1,d), mean=rep(0,d),sigma=Pinv)[1]
      if(is.nan(pez) | pez <0 | sum(abs(z-mu)<0.0003)>0 ){
        ind <- which(abs(z-mu)==min(abs(z-mu)))
        Dy     <- diag(1/(z-mu))
        Psinv     <- Dy[-ind,-ind] %*% Sigma[-ind,-ind] %*% Dy[-ind,-ind]
 Psinv[upper.tri(Psinv)] <- t(Psinv)[upper.tri(Psinv)]
        nus <- nu[-ind]
        const2 <- abs(det(Dy))
  pez1     <- prod(diag(relliptical::mvtelliptical(lower=sqrt(1-nus),upper=rep(2,(d - length(ind))),n=10000,mu=rep(2,(d - length(ind))),Sigma=Psinv,dist="Normal")$EYY))*mvtnorm::pmvnorm(lower=sqrt(1-nus),upper=rep(2,(d - length(ind))), mean=rep(2,(d - length(ind))),sigma=Psinv)[1]
        pez2    <- (2*pi)^(-length(ind)/2)*prod(1/3*(1-(1-nu[ind])^(3/2)))
        pez     <- pez1*pez2
      }
      PDFpass <- const1 * const2 * pez
      PDFmax  <- const1 * 1/det(Sigma)^(1/2) * (2*pi)^(-d/2)*prod(1/3*(1-(1-nu)^(3/2)))
      if(PDFpass > PDFmax)   PDFpass <- PDFmax
PDFpass <- 2*PDFpass*pnorm(la%*%diag(diag(Sigma)^(-1/2))%*%(z-mu))
      if(PDFpass < 1.0e-323) PDFpass <- 1.0e-323
      return(PDFpass)
    } }
  PDF <- sapply(1:q, function(i) fpass(z=x[i,]) )
  return(PDF)
  }





