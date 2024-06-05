
r.SDSTIN.mix <- function(n, w , xi, S, nu, la){
r.SDSTIN <- function(n , xi, S, nu,la){
	p <- length(xi) ;	S <- as.matrix(S) ; Sb=diag(diag(S)^(-1/2))%*%S%*%diag(diag(S)^(-1/2))
	 y<-matrix(0,n,p)
		for(i in 1:n){
		del.tau <- 0
		for(j in 1:p)
		del.tau[j] <- runif(1,1-nu[j],1)^(-1/2)
		Del.tau <- diag(del.tau)
		Z <- t(mvtnorm::rmvnorm(n = 1, mean = rep(0,p), sigma =Sb ))
		Z0 <- rnorm(1)
		if(Z0 < (t(la)%*%Del.tau%*%Z))
		y[i,] <- xi+t(Del.tau%*%diag(diag(S)^(1/2))%*%Z)
		else
		y[i,] <- xi-t(Del.tau%*%diag(diag(S)^(1/2))%*%Z)
		}
		return(y)
		}
	M<- length(w) ; Z <- rmultinom(n,size=1,prob=w); n <- rowSums(Z)
      y<-0
	for( m in 1:M)
	y <- rbind(y,r.SDSTIN(n[m], xi[[m]], S[[m]], nu[[m]],la[[m]]))
	y <- y[-1,]
		return(y)
		}