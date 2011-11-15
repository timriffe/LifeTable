axPreston <-
function(Mx,n,axsmooth=TRUE,sex="female"){
	N <- length(Mx)
	if (axsmooth==TRUE){
		ages <- cumsum(n)-n
		span <- if (N>30) .15 else .4
		Mx <- log(Mx)
		Mx[2:N] <- predict(loess(Mx~ages,span=span,control=loess.control(surface="interpolate")),newdata=ages[2:N])
		Mx <- exp(Mx)
	}
	ax <- n + (1/Mx) - n/(1-exp(-n*Mx))
	
	# Preston rules of thumb:
	if (Mx[1]>=.107){
		if (sex=="female"){
			ax[1] <- .35
		}
		if (sex=="male"){
			ax[1] <- .33
		}
	}
	if (Mx[1]<.107){
		if (sex=="female"){
			ax[1] <- .053+2.8*Mx[1]
		}
		if (sex=="male"){
			ax[1] <- .045+2.684*Mx[1]
		}
	}
	if (n[2]==4){
		if (Mx[1]>=.107){
			if (sex=="female"){
				ax[2] <- 1.352
			}
			if (sex=="male"){
				ax[2] <- 1.361
			}
		}
		if (Mx[1]<.107){
			if (sex=="female"){
				ax[2] <- 1.522-1.581*Mx[1]
			}
			if (sex=="male"){
				ax[2] <- 1.651-2.816*Mx[1]
			}
		}
	}
	# my own adaptation for a1-a8 in the case of single years: otherwise they're *very* close to .5
	if (length(unique(n[1:10]))==1){
		if (unique(n[1:10])==1){
			int <- ax[9]-ax[1]
			jumps <- c(100,30,10,8,6,4,2)
			jumps <- jumps/sum(jumps)
			for (i in 1:7){
				ax[(i+1)] <- ax[i]+jumps[i]*int
			}
		}
	}
	return(ax)
}

