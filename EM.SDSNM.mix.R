EM.SDSNM.mix<-function(y, M, family="STIN", CML=T ,w=NULL, xi=NULL, S=NULL, nu=NULL, la=NULL, iter.max=30, tol=10^-6,get.init = TRUE, group=FALSE, eqnu=FALSE){
 begin <- proc.time()[3] ; library(mnormt);library(zipfR)
dSDSTIN.mix=function(y, M, w, xi, S, nu, la){
dens <- 0
for(m in 1:M)
dens <- dens + w[m]*dSDSTIN(y, as.vector(xi[[m]]), as.matrix(S[[m]]), as.vector(nu[[m]]), as.vector(la[[m]]))
return(dens)
}
dSDSSEN.mix=function(y, M, w, xi, S, nu, la){
dens <- 0
for(m in 1:M)
dens <- dens + w[m]*dSDSSEN(y, as.vector(xi[[m]]), as.matrix(S[[m]]), as.vector(nu[[m]]), as.vector(la[[m]]))
return(dens)
}
y <- as.matrix(y) ; n <- dim(y)[1] ; p <- dim(y)[2]
        dif <- 1 ;  count <- 0 ; LL <- 1 
if (get.init == TRUE) {
               init <- kmeans(y, M,  algorithm="Hartigan-Wong")
               w <- init$size/nrow(y)
xi <-  la <- S <- nu <- list() 
for(m in 1:M){
xi [[m]] <- init$centers[m, ]
S[[m]] <- var(y[init$cluster == m, ])
               la[[m]] <-  sign(apply((y[init$cluster == 
                    m, ] - matrix(rep(xi[[m]], nrow(y[init$cluster == 
                    m, ])), nrow = nrow(y[init$cluster == m, 
                    ]), ncol = p, byrow = TRUE))^3, 2, sum))
if(family=="STIN"){
if(eqnu)
nu[[m]] <- rep(runif(1,0,1),p)
else
               nu[[m]] <- runif(p,0,1)}
else {
if(eqnu)
nu[[m]] <- rep(runif(1,0,1),p)
else
               nu[[m]] <- runif(p,0,1)}
          }     }
if(family=="STIN"){
while ((dif > tol) && (count <= iter.max)) {
# E step
gam.hat <- array(0,c(n,p,M)); z.hat <- W.hat <- matrix(0,n,M)
Del.gam <- array(0, c(p,p,n,M)) ; dSDSTIN.total <- dSDSTIN.mix(y,M,w,xi,S,nu,la)
for (m in 1:M){
z.hat[,m] <- w[m]*dSDSTIN(y,as.vector(xi[[m]]), as.matrix(S[[m]]), as.vector(nu[[m]]), as.vector(la[[m]]))/dSDSTIN.total
Sm <- as.matrix(S[[m]]) ; xi.mat <- matrix(xi[[m]],n,p,byrow=T)
y.xi <- y-xi.mat 
alpha.T <- t(la[[m]])%*%diag(diag(S[[m]])^{-1/2}); alpha <- t(alpha.T); nu.m <- nu[[m]]
for (i in 1:n){
W.hat[i,m] <- as.numeric(alpha.T%*%y.xi[i,])
if (pnorm(alpha.T%*%y.xi[i,])!=0)
W.hat[i,m] <- W.hat[i,m]+ dnorm(alpha.T%*%y.xi[i,])/pnorm(alpha.T%*%y.xi[i,])
for(j in 1:p){
etaij <- (y[i,j]-xi.mat[i,j])^2/diag(Sm)[j]
      num  <-  2*(zipfR::Igamma(a=(1/2+2), x=(1-nu.m[j])*etaij/2, lower=FALSE) - zipfR::Igamma(a=(1/2+2), x=etaij/2, lower=FALSE))
      den  <- etaij*(zipfR::Igamma(a=(1/2+1), x=(1-nu.m[j])*etaij/2, lower=FALSE) - zipfR::Igamma(a=(1/2+1), x=etaij/2, lower=FALSE))
gam.h <- num/den
     if(is.nan(gam.h) | gam.h >= 1) # etaij -> 0 menas y_ij -> xi_j
           gam.h <- 0.999
          if(gam.h <= (1-nu.m[j]))
            gam.h <- (1-nu.m[j]) + 0.001
gam.hat[i,j,m] <- gam.h
}
gam.hat[i,is.na(gam.hat[i,,m]),m] <- 0
Del.gam[,,i,m] <- diag(as.vector(1/gam.hat[i,,m]))
}
# MCE steps
w[m] <- sum(z.hat[,m])/n
a <- 0 ; b <- 0
for(i in 1:n){
a <- a+ z.hat[i,m]*(diag(diag(Del.gam[,,i,m])^(-1/2))%*%solve(Sm)%*%diag(diag(Del.gam[,,i,m])^(-1/2))+alpha%*%alpha.T)
b <- b + z.hat[i,m]*(diag(diag(Del.gam[,,i,m])^(-1/2))%*%solve(Sm)%*%diag(diag(Del.gam[,,i,m])^(-1/2))%*%y[i,]+alpha%*%alpha.T%*%y[i,]-alpha*W.hat[i,m])
}
xi[[m]] <- solve(a)%*%b
xi.mat <- matrix(xi[[m]],n,p,byrow=T)
y.xi <- y-xi.mat
c <- 0
for(i in 1:n)
c <- c + z.hat[i,m]*(diag(diag(Del.gam[,,i,m])^(-1/2))%*%y.xi[i,]%*%t(y.xi[i,])%*%diag(diag(Del.gam[,,i,m])^(-1/2)))
S[[m]] <- c/sum(z.hat[,m])
cc <- 0 ; dd <- 0
for(i in 1:n){
cc <- cc + z.hat[i,m]*y.xi[i,]%*%t(y.xi[i,])
dd <- dd + z.hat[i,m]*W.hat[i,m]*y.xi[i,] }
alpha <- solve(cc)%*% dd
alpha.T <- t(alpha)
la[[m]] <- diag(diag(S[[m]])^(1/2))%*%alpha
if(eqnu){
nu0 <- optim(nu.m[1],function(x){
x <- rep(x,p)
-sum(w[[m]]*log(dSDSTIN(y, mu=as.vector(xi[[m]]), Sigma=as.matrix(S[[m]]), nu=as.vector(x), la=as.vector(la[[m]])))) 
},method="L-BFGS-B",lower=.01,upper=.99, control = list(maxit=500,trace = 0,factr=1e100,pgtol=1e1))$par
nu[[m]] <- rep(nu0,p) }else{
nu[[m]] <- optim(nu.m,function(x){
-sum(w[[m]]*log(dSDSTIN(y, mu=as.vector(xi[[m]]), Sigma=as.matrix(S[[m]]), nu=as.vector(x), la=as.vector(la[[m]])))) 
},method="L-BFGS-B",lower=.01,upper=.99, control = list(maxit=500,trace = 0,factr=1e60,pgtol=1e1))$par
 }  }
LL.new <- sum(log(dSDSTIN.mix(y,M, w, xi, S, nu, la))) # log-likelihood function
count <- count +1 
dif <- abs(LL.new/LL-1)
LL <- LL.new
cat('iter =', count, '\tloglike =', LL.new, '\n')
}} else { ####################################### Family=SSEN   #####################
while ((dif > tol) && (count <= iter.max)) {
# E step
gam.hat <- array(0,c(n,p,M)); z.hat <- W.hat <- matrix(0,n,M)
Del.gam <- array(0, c(p,p,n,M)) ; dSDSSEN.total <- dSDSSEN.mix(y,M,w,xi,S,nu,la)
for (m in 1:M){
z.hat[,m] <- w[m]*dSDSSEN(y,as.vector(xi[[m]]), as.matrix(S[[m]]), as.vector(nu[[m]]), as.vector(la[[m]]))/dSDSSEN.total
Sm <- as.matrix(S[[m]]) ; xi.mat <- matrix(xi[[m]],n,p,byrow=T)
y.xi <- y-xi.mat 
alpha.T <- t(la[[m]])%*%diag(diag(S[[m]])^{-1/2}); alpha <- t(alpha.T); nu.m <- nu[[m]]
for (i in 1:n){
W.hat[i,m] <- as.numeric(alpha.T%*%y.xi[i,])
if (pnorm(alpha.T%*%y.xi[i,])!=0)
W.hat[i,m] <- W.hat[i,m]+ dnorm(alpha.T%*%y.xi[i,])/pnorm(alpha.T%*%y.xi[i,])
for(j in 1:p){
etaij <- (y[i,j]-xi.mat[i,j])^2/diag(Sm)[j]
      num    <- expint::gammainc(a = (1/2 + 2), x = (etaij/2 + nu.m[j]))
      den    <- (etaij/2 + nu.m[j])*expint::gammainc(a = (1/2 + 1), x = (etaij/2 + nu.m[j]))
gam.hat[i,j,m] <- num/den
}
gam.hat[i,is.na(gam.hat[i,,m]),m] <- 0
Del.gam[,,i,m] <- diag(as.vector(1/gam.hat[i,,m]))
}
# MCE steps
w[m] <- sum(z.hat[,m])/n
a <- 0 ; b <- 0
for(i in 1:n){
a <- a+ z.hat[i,m]*(diag(diag(Del.gam[,,i,m])^(-1/2))%*%solve(Sm)%*%diag(diag(Del.gam[,,i,m])^(-1/2))+alpha%*%alpha.T)
b <- b + z.hat[i,m]*(diag(diag(Del.gam[,,i,m])^(-1/2))%*%solve(Sm)%*%diag(diag(Del.gam[,,i,m])^(-1/2))%*%y[i,]+alpha%*%alpha.T%*%y[i,]-alpha*W.hat[i,m])
}
xi[[m]] <- solve(a)%*%b
xi.mat <- matrix(xi[[m]],n,p,byrow=T)
y.xi <- y-xi.mat
c <- 0
for(i in 1:n)
c <- c + z.hat[i,m]*(diag(diag(Del.gam[,,i,m])^(-1/2))%*%y.xi[i,]%*%t(y.xi[i,])%*%diag(diag(Del.gam[,,i,m])^(-1/2)))
S[[m]] <- c/sum(z.hat[,m])
cc <- 0 ; dd <- 0
for(i in 1:n){
cc <- cc + z.hat[i,m]*y.xi[i,]%*%t(y.xi[i,])
dd <- dd + z.hat[i,m]*W.hat[i,m]*y.xi[i,] }
alpha <- solve(cc)%*% dd
alpha.T <- t(alpha)
la[[m]] <- diag(diag(S[[m]])^(1/2))%*%alpha
if (CML){
if(eqnu){
nu0 <- optim(nu.m[1],function(x){
x <- rep(x,p)
-sum(w[[m]]*log(dSDSSEN(y, mu=as.vector(xi[[m]]), Sigma=as.matrix(S[[m]]), nu=as.vector(x), la=as.vector(la[[m]])))) 
},method="L-BFGS-B",lower=.01,upper=Inf, control = list(maxit=500,trace = 0,factr=1e100,pgtol=1e1))$par
nu[[m]] <- rep(nu0,p) }else{
nu[[m]] <- optim(nu.m,function(x){
-sum(w[[m]]*log(dSDSSEN(y, mu=as.vector(xi[[m]]), Sigma=as.matrix(S[[m]]), nu=as.vector(x), la=as.vector(la[[m]])))) 
},method="L-BFGS-B",lower=.01,upper=Inf, control = list(maxit=500,trace = 0,factr=1e60,pgtol=1e1))$par
 } } else {
if(eqnu) 
nu[[m]]=rep((sum(z.hat[,m])*p)/(sum(rowSums(z.hat[,m]*gam.hat[,,m]))-sum(z.hat[,m])*p),p)
else {
for(j in 1:p)
nu[[m]][j]=sum(z.hat[,m])/(sum(z.hat[,m]*gam.hat[,j,m])-sum(z.hat[,m]))
} }  } 
LL.new <- sum(log(dSDSSEN.mix(y,M, w, xi, S, nu, la))) # log-likelihood function
count <- count +1 
dif <- abs(LL.new/LL-1)
LL <- LL.new 
cat('iter =', count, '\tloglike =', LL.new, '\n')
} }
if(eqnu){
aic <- -2 * LL.new + 2 * (M*(2*p+p*(p+1)/2+1)+M-1)
bic <- -2 * LL.new + log(n) * (M*(2*p+p*(p+1)/2+1)+M-1)
} else{
aic <- -2 * LL.new + 2 * (M*(3*p+p*(p+1)/2)+M-1)
bic <- -2 * LL.new + log(n) * (M*(3*p+p*(p+1)/2)+M-1) }
end <- proc.time()[3]
time <- end-begin
obj.out <- list(w=w, xi=xi, S=S , nu=nu, la=la, loglik=LL.new, AIC=aic, BIC=bic, iter=count,elapsed=as.numeric(time),group = apply(z.hat, 
                1, which.max))
if (group==FALSE)
obj.out <- obj.out[-length(obj.out)]
obj.out
}
