
#' @title axPreston approximation of a(x) using a few rules of thumb outlined by Preston et al (2001)
#' 
#' @description a(x) estimated assuming death rate constant in each interval. Most a values are estimated using \code{ax = n + (1/Mx) - n/(1-exp(-n*Mx))}, except for a0, which uses a few rules of thumb derived from Coale and Demeny (1983) and displayed in table 3.3 in Preston (2001). In the case of single ages, I found that ages a1-a10 were all estimated very close to .5, whereas the older ages a50+ were all estimated very close to what other methods produce. In order to adjust a(x) values
#' 
#' @param Mx a numeric vector of the age-specific central death rates, calculated as D(x)/N(x) (deaths/exposure).
#' @param n a numeric vector of age interval widths.
#' @param axsmooth logical. default = \code{TRUE}. Should the a(x) values be calculated from a smoothed M(x) series? In this case, the M(x) series is smoothed within the function for a(x) estimation, but the smoothed M(x) function that was used is not returned. In general, it is better to smooth the M(x) function prior to putting it in this function, because the loess smoother used here has no weights or offset. If this is not possible, loess M(x) smoothing still produces more consistent and less erratic a(x) estimates.
#' @param sex \code{"male"} or \code{"female"}. default = \code{"female"}. The Coale Demeny rules of thumb are different for males and females.
#' 
#' @details In the case of single ages, I found that ages a1-a10 were all estimated very close to .5, whereas the older ages a50+ were all estimated very close to what other methods produce. In order to adjust a(x) values to reflect dropping mortality at young ages, I wrote the following rule of thumb, which scales based on the level of mortality at age 0 and age 8. Basically the drop in mortality from age 0 to 8 ought to produce successively larger a(x) values, approaching .5. The increments in a(x) at each age from 1 until 8 are thus applied according to some fixed proportions, contained in the variable 'jumps'. This is the last code chunk in the function, which is displayed below under 'examples'. For more info, look at the code.
#' 
#' @return a numeric vector of a(x) values.
#' 
#' @references 
#' Coale Anseley and Paul Demeny, with B Vaughan (1983). Regional Model Life Tables and Stable Populations. New York Academic Press.
#' 
#' Preston, S. et al (2001) Demography: measuring and modeling population processes. Blackwell Publishing. Malden
#' 
#' @seealso See Also as \code{\link{axEstimate}}, a wrapper function for this and three other a(x) estimation procedures (\code{\link{axMidpoint}}, \code{\link{axKeyfitz}} and \code{\link{axSchoen}}).
#' 
#' @examples 
#' \dontrun{
#' 
#' library(LifeTable)
#' data(UKRmales1965)
#' Mx       <- UKRmales1965[, 4]
#' n        <- rep(1, length(Mx))
#' axp1     <- axPreston(Mx, n, axsmooth = FALSE, sex = "male")
#' axp2     <- axPreston(Mx, n, axsmooth = TRUE, sex = "male")
#' plot(0:110, axp1, main = "comparing smoothed and unsmoothed a(x)", xlab = "age", ylab = "a(x)")
#' lines(0:110, axp2, col = "blue")
#' arrows(24.15, 0.4, 8.6, 0.45, lty = 1, lwd = 1, col = "black")
#' arrows(22.7, 0.38, 4.9, 0.35, lty = 1, lwd = 1, col = "black")
#' text(x = 26, y = 0.4, labels = "This is my hack rule of thumb for a1 to a8", cex = 1, pos = 4, col = 1)
#' text(x = 26, y = 0.38, labels = "see last code chunk to see how it works", cex = 1, pos = 4, col = 1)
#' 
#' ########## with 5-year age groups
#' data(UKR5males1965)
#' Mx       <- UKR5males1965[, 4]
#' n        <- c(1, 4, rep(5, (length(Mx) - 2)))
#' ages     <- cumsum(n) - n
#' axp1     <- axPreston(Mx, n, axsmooth = FALSE, sex = "male")
#' axp2     <- axPreston(Mx, n, axsmooth = TRUE, sex = "male")
#' plot(ages, c(axp1[1], axp1[2] / 4, axp1[3:24] / 5), main = "comparing smoothed and unsmoothed a(x)", xlab = "age", ylab = "a(x)")
#' lines(ages, c(axp2[1], axp2[2] / 4, axp2[3:24] / 5), col = "blue")
#' }
#' 
#' @export

axPreston <-
function(Mx, n, axsmooth = TRUE, sex = "female"){
	N <- length(Mx)
	if (axsmooth){
		ages    <- cumsum(n) - n
		span    <- ifelse(N > 30, .15, .4)
		Mx      <- log(Mx)
		Mx[2:N] <- predict(loess(Mx ~ ages, 
                                 span = span, 
                                 control = loess.control(surface = "interpolate")
                                 ), 
                           newdata = ages[2:N]
                   )
		Mx      <- exp(Mx)
	}
	ax      <- n + (1 / Mx) - n / (1 - exp(-n * Mx))
	
	# Coale-Demeny rules of thumb:
    ax[1] <- ifelse(Mx[1] >= .107, 
                ifelse(sex == "female", .35, .33),
                ifelse(sex == "female", {.053 + 2.8 * Mx[1]}, {.045 + 2.684 * Mx[1]})
             )
	if (n[2] == 4){
        ax[2] <- ifelse(Mx[1] >= .107,
                        ifelse(sex == "female", 1.352, 1.361),
                        ifelse(sex == "female", {1.522 - 1.581 * Mx[1]}, {1.651 - 2.816 * Mx[1]})
                 )
	}
	# my own adaptation for a1-a8 in the case of single years: otherwise they're *very* close to .5
	if (length(unique(n[1:10])) == 1){
		if (unique(n[1:10]) == 1){
			int         <- ax[9] - ax[1]
			jumps       <- c(100, 30, 10, 8, 6, 4, 2)
			jumps       <- jumps / sum(jumps)
			for (i in 1:7){
				ax[(i + 1)]     <- ax[i] + jumps[i] * int
			}
		}
	}
	return(ax)
}

