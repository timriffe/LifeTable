axSchoen <-
function(Mx,n,axsmooth=TRUE){
	N <- length(Mx)
	if (axsmooth==TRUE){
		ages <- cumsum(n)-n
		span <- if (N>30) .15 else .4
		Mx <- log(Mx)
		Mx[2:N] <- predict(loess(Mx~ages,span=span,control=loess.control(surface="interpolate")),newdata=ages[2:N])
		Mx <- exp(Mx)
	}
	ax <- ux <- wx <- lx <- vector(length=N)
	lx[1] <- 1
	for (i in 2:N){
		lx[i] <- lx[i-1]*exp(-n[i-1]*Mx[i-1])
	}
	dx <- -diff(lx)
	dx <- c(dx,1-sum(dx))
	for (i in 2:(N-1)){
		ux[i] <- ((n[i]^2)/240)*(Mx[i+1]+38*Mx[i]+Mx[i-1])
		wx[i] <- ((n[i]^2)/240)*(14*Mx[i+1]+72*Mx[i]-6*Mx[i-1])
		ax[i] <- (ux[i]*lx[i]+wx[i]*lx[i+1])/dx[i]
	}
	ax[1] <- .07+1.7*Mx[1]
	########### linear assumption for last ax ####################
	########### using 3 info points           ####################
	coefs <- lm(ax[(N-3):(N-1)]~c((N-3):(N-1)))$coef
	ax[N] <- coefs[1]+N*coefs[2]
	return(ax)
}

