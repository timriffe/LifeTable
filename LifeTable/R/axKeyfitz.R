axKeyfitz <-
function(Mx, n, axsmooth = TRUE){
	# iterative ax-dx process decribed on page 44-45 of 
	# Preston et al, Demography: Measuring and Modelling Population Processes. Blackwell Publishing, 2001
	N <- length(Mx)
	if (axsmooth){
		ages        <- cumsum(n) - n
		span        <- ifelse(N > 30, .15, .4)
		Mx          <- log(Mx)
		Mx[2:N]     <- predict(loess(Mx ~ ages, 
                                     span = span, 
                                     control = loess.control(surface = "interpolate")
                                     ), 
                               newdata = ages[2:N]
                       )
		Mx          <- exp(Mx)
	}
	axit        <- .5 * n
	axit[1] <- .07 + 1.7 * Mx[1]
	for (i in 1:7){
		qx              <- (n * Mx) / (1 + (n - axit) * Mx)
		qx[length(Mx)]  <- 1
		px              <- 1 - qx
		lx              <- 1 # typically radix would go here, but it makes no difference since values don't pass on.
		for (i in 2:length(Mx))	{ 
            lx[i]   <- lx[i - 1] * px[i - 1] 
        }
		dx          <- -diff(lx)
		for (i in 2:(length(Mx) - 1)){
			axit[i] <- (-(n[i - 1] / 24) * dx[i - 1] + (n[i] / 2) * dx[i] + (n[i + 1] / 24) * dx[i + 1]) / dx[i]
        }
		
		# this is just my own way of finishing off the ax's, not sooo creative, 
		# but it doesn't usually make a difference
		axit[N - 1] <- axit[N - 2] - (axit[N - 3] - axit[N - 2]) * 1.5
		axit[N]     <- axit[N - 1] - (axit[N - 2] - axit[N - 1]) * 1.5
		# it assumes continued senescence at the final ages:
		axit[N - 1] <- axit[N - 2] - (axit[N - 3] - axit[N - 2]) * 1.5
		axit[N]     <- axit[N - 1] - (axit[N - 2] - axit[N - 1]) * 1.5
	}
	axit[1]     <- .07 + 1.7 * Mx[1]
	return(axit)
}

