\name{axSchoen}
\alias{axSchoen}

\title{Schoen's MDAT approximation of Chiang's a(x), average number of years lived by persons dying in the interval x,x+n.}
\description{
This estimation technique of Chiang's a(x) was proposed by Schoen (1978) and is called the Mean Duration at Transfer (MDAT) method. Adjacent values of the force of mortality function, M(x), assumed to be the lifetable m(x), are used to estimate the mean age at exit from each interval.
}
\usage{axSchoen(Mx, n, axsmooth = TRUE)}

\arguments{
  \item{Mx}{a numeric vector of the age-specific central death rates, calculated as D(x)/N(x) (deaths/exposure)}
  \item{n}{a numeric vector of age interval widths.}
  \item{axsmooth}{
logical. default = TRUE. Should the a(x) values be calculated from a smoothed M(x) series? In this case, the M(x) series is smoothed within the function for a(x) estimation, but the smoothed M(x) function that was used is not returned. In general, it is better to smooth the M(x) function prior to putting it in this function, because the loess smoother used here has no weights or offset. If this is not possible, loess M(x) smoothing still produces more consistent and less erratic a(x) estimates.
}
}

\details{
The very last a(x) value is imputed using linear extrapolation of the prior 3 points. For most demographic measures, precision at the upper ages has little or no effect on results. In general, this estimation procedure works very well, and produces results very similar to the Keyfitz iterative procedure, but the estimates at the final ages tend to turn upward, where they ought to slope downward instead. 
}

\value{returns \code{ax}, a numeric vector of a(x) values.}
\references{
Chiang C.L.(1968) Introduction to Stochastic Processes in Biostatistics. New York: Wiley.\\
Schoen R. (1978) Calculating lifetables by estimating Chiang's a from observed rates. Demography 15: 625-35.\\
Preston et al (2001) Demography: Measuring and Modeling Population Processes. Malden MA: Blackwell.
}
\author{Tim Riffe}

\seealso{See Also as \code{\link{axEstimate}}, a wrapper function for this and three other a(x) estimation procedures (\code{\link{axMidpoint}}, \code{\link{axKeyfitz}} and \code{\link{axPreston}}).}
\examples{
library(LifeTable)
data(UKRmales1965)
Mx <- UKRmales1965[,4]
n <- rep(1,length(Mx))
axs1 <- axSchoen(Mx,n,axsmooth=F)
axs2 <- axSchoen(Mx,n,axsmooth=T)
plot(0:110,axs1,main="comparing smoothed and unsmoothed a(x)",xlab="age",ylab="a(x)")
lines(0:110,axs2,col="blue")
# note that at old ages the a(x) estimate shoots upward. It should in fact decline.
# I am still not sure whether this is due to the formula or my implementation of it...

########## with 5-year age groups
data(UKR5males1965)
Mx <- UKR5males1965[,4]
n <- c(1,4,rep(5,(length(Mx)-2)))
ages <- cumsum(n)-n
axs1 <- axSchoen(Mx,n,axsmooth=F)
axs2 <- axSchoen(Mx,n,axsmooth=T)
plot(ages,c(axs1[1],axs1[2]/4,axs1[3:24]/5),main="comparing smoothed and unsmoothed a(x)",xlab="age",ylab="a(x)")
lines(ages,c(axs2[1],axs2[2]/4,axs2[3:24]/5),col="blue")



## The function is currently defined as
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
}

