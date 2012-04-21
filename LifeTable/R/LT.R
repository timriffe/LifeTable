LT <-
function(Nx, Dx, Mx, ages = "auto", type = "single-age", axmethod = "keyfitz", sex = "female", 
		mxsmooth = TRUE, axsmooth = TRUE, n = "auto", radix = 1, verbose = TRUE){
	# the verbose function:
	Verb <- function(v, x){
		if (v) {
			cat(paste(x, "\n", sep = ""))
		}
	}
	# first a series of checks, messages and imputations to make sure the given Dx&Nx *or* Mx values are viable
	if (missing(Mx)){
		# two checks that will stop the function
		# 1) in absence of Mx, need both Dx and Nx
		if (missing(Nx) | missing(Dx)){ 
			Verb(verbose,"you're missing the key arguments")
			stop("either specify Mx directly, or else specify both Nx and Dx.")
		}
		# 2) both Nx and Dx must be of equal length
		if (length(Nx) != length(Dx)){
			Nxl <- length(Nx)
			Dxl <- length(Dx)
			Verb(verbose,"Nx and Dx lengths not equal\nplease examine and adjust these vectors")
			stop(paste("Nx length =", Nxl, ", Dx length =", Dxl))
		}
		# try to coerce to numeric vectors
		Nx <- as.numeric(c(unlist(Nx))) 
		Dx <- as.numeric(c(unlist(Dx)))
		
		# safety to avoid zeros in the denominator. will make Mx of 1 in those cases
		if (any(Nx == 0)) {
			Verb(verbose,"there was at least one 0 in your Nx vector.")
			Verb(verbose, (paste("index value(S)=", which(Nx == 0))))
			Verb(verbose, "the function imputed the corresponding Dx value into the Nx 0s\nwhich makes an Mx of 1 in those slots, and it was thus able to continue.\nif this was just the last Nx value(s), then it makes no difference.\nLook over the Nx vector to verify that there are no other 0s that could be more problematic")
			Nx[Nx == 0] <- Dx[Nx == 0]
		}
		Mx <- Dx / Nx
	}
	
	# by this point we should have an Mx vector with no NAs: no more need for Nx,Dx
	# we want to be able to accept 0s...
	
	# N is just used for counting, to save space
	N <- length(Mx) # Mx <- mxmat[1,]
	
	# identify single-age or abridged table
	types <- c("single-age", "abridged", "other")
	typesmenu <- c("single-age (0,1,2,3....110+)", "abridged (0,1-4,5-9...)", "other?")
	if (! type %in% types | missing(type)) {
		type <- types[menu(typesmenu, graphics = TRUE, title = "pick relevant lifetable type")]
	}
	if (type == "other") {
		cat("\nsorry, I only have the first 2 kinds of age-categories implemented so far\n")
		stop("please send suggestions to tim.riffe@gmail.com")
	}
	
	# assign interval widths, including the open interval
	if (is.character(n) & n == "auto" & type == "abridged"){
		Widths <- c(1, 4, rep(5, N - 2))
	}
	if (is.character(n) & n == "auto" & type == "single-age"){
		Widths <- rep(1, N)
	}
	if (is.numeric(n)){	
		Widths <- n
		Verb(verbose, "used user-supplied widths.If there's an error, be sure this vector is of the same length as Nx,Dx")
	}
	
	if (ages == "auto" & type == "abridged") {
		Age <- c(0, paste(c(1, seq(from = 5, by = 5, length = (N - 3))), "-",
						c(4, seq(from = 9, by = 5, length = (N - 3))), sep = ""), paste(((N-2) * 5), "+", sep = ""))
	}
	if (ages == "auto" & type == "single-age") {
		Age <- c(seq(0, N - 2), paste(N - 1, "+", sep = ""))
	}
	# redefine term ages for later use.
	if (ages == "auto"){
		ages <- cumsum(Widths) - Widths
	}
	ages.mids.pre <- ages + Widths/2
	ages.mids.pre[1] <- .1
	
	if(!missing(Nx) & !missing(Dx) & mxsmooth){
		# Giancarlo's package. I recommend either supplying an already-smoothed Mx vector (for complete control)
		# or else supplying Dx and Nx and setting mxsmooth to TRUE. 
		require("MortalitySmooth")
		fitBIC <- Mort1Dsmooth(x = ages.mids.pre, y = Dx, offset = log(Nx))
		Mx[2:N] <- (fitted(fitBIC) / Nx)[2:N]
	}
	
	if(missing(Nx) & missing(Dx) & mxsmooth){
		Verb(verbose,"mxsmooth was specified as TRUE, but since Mx was supplied directly, \nthere are no implicit weights (Nx). Function used a loess smoother \nto smooth out the Mx, but please be wary.")
		span <- ifelse(N > 30, .15, .4)
		logMx <- log(Mx)
		logMx[is.infinite(logMx)] <- NA
		logMx[2:N] <- stats:::predict(stats:::loess(logMx ~ ages.mids.pre, span = span, control = stats:::loess.control(surface = "interpolate")), newdata = ages[2:N])
		Mx <- exp(logMx)
	}

	# the link identity. we assume that the lifetable mx is equal to the central death rate.
	mx <- Mx
	
	# these later 3 imputations should not be needed, but are there just in case.
	mx[is.na(mx)] <- 0
	mx[is.infinite(mx)] <- 0
	
	# we don't want 0s in m(x) for calculating a(x), because a(x) (except midpoint) is iterative
	# and we'd erroneously bring down neighboring age's ax values if with zeros.
	# for later calulations, the zeros are taken 'as-is'
	if (length(axmethod) == 1 & min(mx) == 0){
		Ind <- which(mx == 0)
		
		Verb(verbose, paste("\n\n*there were some ages (", ages[Ind], ") with no mortality.\nValues are imputed for calculating a(x), but the zeros are kept for the rest of the lifetable.\n"))
		span <- ifelse(N > 30,.15,.4)
		logMx <- log(Mx)
		logMx[is.infinite(logMx)] <- NA
		Imp <- exp(stats:::predict(stats:::loess(logMx~ages.mids.pre,span=span,control=stats:::loess.control(surface="interpolate")),newdata=ages))[Ind]
		if (any(is.na(Imp))){
			Imp <- exp(stats:::spline(ages.mids.pre, logMx, xmin = 0, xmax = max(ages))$y[Ind])
		}
		mx[Ind] <- Imp
		
	}
	
	if (length(axmethod) == 1){
		if (axmethod %in% c("keyfitz", "schoen", "midpoint", "preston")){
			# here, the ax iteration is externalized to axEstimate()
			if (mxsmooth){
				axsmooth <- FALSE
			}
			ax <- axEstimate(Mx = mx, n = Widths, axsmooth = axsmooth, method = axmethod, sex = sex)
		}
		if (axmethod == "keyfitz" & type == "abridged"){
			Verb(verbose, "It appears you have an abridged lifetable, but have specified the keyfitz method of ax estimation.\nBe aware that this method presumes equal age intervals, as Preston et. al. (2001)\n warn on page 45. Consider using a different method or else specifying your own ax vector.\n Function continued nonetheless.")
		}
	}
	
	if (is.numeric(axmethod) & length(axmethod) == N) {
		ax <- axmethod
	}
	# last default
	if (!exists("ax")){
		ax <- axEstimate(Mx = mx, n = Widths, axsmooth = axsmooth, method = "midpoint")
		Verb(verbose, "axmethod must be specified either as 'schoen','keyfitz','midpoint'\nor as a numeric vector the same length as Nx.\nThis wasn't the case, so the function defaulted to the midpoint method.")
	}
	
	# if zeros were imputed for ax estimation, then we put them back for the rest of the calculations
	if (exists("Ind")){
		mx[Ind] <- 0
	}
	
	qx <- (Widths * mx) / (1 + (Widths - ax) * mx)
	qx[N] <- 1
	
	# can't have qx > 1, so we impute 1s where necessary: hopefully only at the penultimate, as the case may be
	qx[qx > 1] <- 1
	px <- 1 - qx
	
	lx <- Lx <- vector(length = N)
	lx[1] <- radix
	for (i in 2:N) {
		lx[i] <- lx[i - 1] * px[i - 1] 
	}
	
	dx <- -diff(lx)
	dx[N] <- lx[N]
	
	Lx[1:(N-1)] <- Widths[1:(N-1)] * lx[2:N] + ax[1:(N-1)] * dx[1:(N-1)]
	Lx[N] <- lx[N] / mx[N]
	Lx[is.infinite(Lx)] <- 1
	Lx[is.na(Lx)] <- 0
	
	Tx <- rev(cumsum(rev(Lx)))
	ex <- Tx / lx
	
	# any missing ex values are replaced by the corresponding ax values 
	if (any(is.na(ex))) {
		cat("\n\n*some value(s) of age-specific life expectancy could not be calculated in the conventional way.\nIt has been assumed that this was only the case in the final age groups, and a(x) was imputed here.")
		ex[is.na(ex)] <- ax[is.na(ex)]
	}
	
	# in the case that there was no exposure in the last open group (as was the case in my test population)
	# I decided it made sense to plug in the last ax value for ex, for the hypothetical case tha someone reaches
	# that age. I was only able to get that ax value by extrapolating during the iteration anyway. This wouldn't
	# have much affect, and will have no effect if indeed there are no people in that age category
	
	# Sx is the pertinent output for projections
    Sx <- vector(length = N)
	Sx[1:(N - 1)] <- (Lx[2:N] / Widths[2:N]) / (Lx[1:(N-1)] / Widths[1:(N - 1)])
	Sx[N]       <- Tx[N] / Tx[(N - 1)]
	# two not-very satisfying, and possibly redundant, substitutions:
	Sx[Lx == 0]   <- 0
	Sx[is.na(Sx)] <- 0
	
	
	# another calculation of e0:
	e0estimates <- matrix(nrow = 1, ncol = 3)
	e0estimates[1] <- e0LT <- ex[1]
	e0estimates[2] <- e0dx <- sum((ages + ax) * dx) / radix
	e0estimates[3] <- e0lx <- sum(lx * Widths) / radix - .5
	colnames(e0estimates) <- c("T0/l0", "sum((ages+ax)*dx)", "sum(lx)-.5")
	rownames(e0estimates) <- "e0"
    # several estimates of edagger (the first one is probably the best)

# edagger deprecated for now, to be reexamined
#	edaggerestimate <- function(lx,ex,dx,mx,ax){
#		N <- length(lx)
#		ed1 <- ed2 <- ed3 <- ed4 <- vector(length=N)
#		endterm <- 1/lx[N]*dx[N]*.5*ex[N]
#		for (i in 1:N){
#			ed2[i] <- (1/lx[i])*sum(lx[i:N]*mx[i:N]*ex[i:N])
#			ed3[i] <- (1/(2*lx[i]))*sum(dx[i:(N-1)]*(ex[i:(N-1)]+ex[(i+1):N]))
#			ed4[i] <- (1/lx[i])*sum(dx[i:(N-1)]*((1-ax[i:(N-1)])+ex[(i+1):N]))
#			sumthing <- 0
#			if (i < 111){
#				for (j in i:(N-1)){
#					sumthing <- sumthing+(dx[j]*(ex[j+1]+1-ax[j]))
#				}
#			}
#			ed1[i] <- (1/lx[i])*sumthing+endterm
#		}
#		edagger <- matrix(cbind(ed1,ed2,ed3,ed4),ncol=4,nrow=N)
#		return(edagger)
#	}
#	edagger <- edaggerestimate(lx,ex,dx,mx,ax)
#	colnames(edagger) <- c("formula1","formula2","formula3","formula4")
#	rownames(edagger) <- Age
	# pasting together a lifetable object
	LT <- data.frame(cbind("Age" = Age, "ages" = ages, "mx" = round(mx, 4), "ax" = round(ax, digits = 4),
					"qx" = round(qx, 4), "px" = round(px, 4), "lx" = round(lx, 4),
					"dx" = round(dx, 4), "Lx" = round(Lx, 4), "Tx" = round(Tx, 4), "ex" = round(ex, 4)))
	# both LT as well as the individual pieces (not rounded) can be called
	return(list(LT = LT, Age = Age, ages = ages, mx = mx, ax = ax, qx = qx, lx = lx, dx = dx, Lx = Lx, Tx = Tx, ex = ex, Sx = Sx, Widths = Widths, e0est = e0estimates))
}

