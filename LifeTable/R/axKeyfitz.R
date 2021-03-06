#' @title A Nathan Keyfitz approximation of Chiang's a(x), average number of years lived by persons dying in the interval x,x+n.
#' 
#' @description An iterative approximation of Chiang's a(x), as proposed by Nathan Keyfitz (1966) and described in Preston (2001).
#' 
#' @param Mx numeric vector of the age-specific central death rates, calculated as D(x)/N(x) (deaths/exposure)
#' @param n a numeric vector of age interval widths.
#' @param axsmooth logical. default = \code{TRUE}. Should the a(x) values be calculated from a smoothed M(x) series? In this case, the M(x) series is smoothed within the function for a(x) estimation, but the smoothed M(x) function that was used is not returned. In general, it is better to smooth the M(x) function prior to putting it in this function, because the loess smoother used here has no weights or offset. If this is not possible, loess M(x) smoothing still produces more consistent and less erratic a(x) estimates.
#' 
#' @details This procedure usually converges very quickly. First, dx is estimated from an ax vector of interval midpoints (except age 0, which is ignored in this process. This d(x) is used to estimate another a(x) based on d(x) slopes. Repeat 7 times. The Keyfitz iterative procedure provides no estimates for the final two ages. These I impute (see last part of code, below) based on the a(x) slope for the last estimated interval. In the penultimate interval, this increment is multiplied by 1.5, and again multiplied by 1.5 for the final a(x) value, which has the effect of exaggerating the tendency at the end of the series.
#' 
#' @return returns \code{ax}, a numeric vector of a(x) values.
#' 
#' @references Chiang C.L.(1968) Introduction to Stochastic Processes in Biostatistics. New York: Wiley.
#' Keyfitz, Nathan (1966) A Life Table that Agrees with the Data. Journal of the American Statistical Association, 61 (314):305-12. (As described on page 44-45 of Preston et al (2001). Demography: Measuring and Modelling Population Processes. Blackwell Publishing)
#' 
#' @note Preston warns that this method should only be used when the age intervals are equal, which is not the case with abridged data. The function will still estimate using such data, but be wary.
#' 
#' @seealso \code{\link{axEstimate}}, a wrapper function for this and three other a(x) estimation procedures (\code{\link{axMidpoint}}, \code{\link{axSchoen}} and \code{\link{axPreston}}).
#' 
#' @examples 
#' \dontrun{
#' library(LifeTable)
#' data(UKRmales1965)
#' Mx <- UKRmales1965[, 4]
#' n <- rep(1, length(Mx))
#' axk1 <- axKeyfitz(Mx, n, axsmooth = FALSE)
#' axk2 <- axKeyfitz(Mx, n, axsmooth = TRUE)
#' plot(0:110, axk1, main = "comparing smoothed and unsmoothed a(x)", xlab = "age", ylab = "a(x)")
#' lines(0:110, axk2, col = "blue")
#' 
#' ########## with 5-year age groups (not recommended for adbridged data)
#' data(UKR5males1965)
#' Mx <- UKR5males1965[, 4]
#' n <- c(1, 4, rep(5, (length(Mx) - 2)))
#' ages <- cumsum(n) - n
#' axk1 <- axKeyfitz(Mx, n, axsmooth = FALSE)
#' axk2 <- axKeyfitz(Mx, n, axsmooth = TRUE)
#' plot(ages, c(axk1[1], axk1[2] / 4, axk1[3:24] / 5), main = "not recommended for use with abridged data", xlab = "age", ylab = "a(x)")
#' lines(ages, c(axk2[1], axk2[2] / 4, axk2[3:24] / 5), col = "blue")
#'}
#' 
#' @export

axKeyfitz <-
function(Mx, n, axsmooth = FALSE){
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
	for (j in 1:7){
		qx              <- (n * Mx) / (1 + (n - axit) * Mx)
		qx[N]           <- 1
		px              <- 1 - qx
		lx              <- c(1,cumprod(px))[1:N]
		dx              <- qx * lx
		for (i in 2:(N - 1)){
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




