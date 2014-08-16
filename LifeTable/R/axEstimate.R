#' @title axEstimate A function to estimate Chiang's a(x) values using a variety of methods.
#' 
#' @description This is a wrapper function to estimate the average distance into an age interval lived by those dying within that age interval. It calls 4 different methods: \code{"keyfitz"}, \code{"schoen"}, \code{"midpoint"} or \code{"preston"}. These generally use formulas proposed by their respective namesakes, although all have been modified by the current author in minor ways, usually in order to provide values for the final ages, which are left as NAs using certain methods. The Preston method is called so not because he explicitly invented it, but rather because it follows a series of rules of thumbs drawn from other sources and described so well in Preston et al (2001). See the individual ax estimation functions to see the details of the various methods. a0 is handled using a variant of the Andreev-Kingkade method.
#' 
#' @param Mx a numeric vector of the age-specific central death rates, calculated as D(x)/N(x) (deaths/exposure)
#' @param n a numeric vector of age interval widths.
#' @param axsmooth logical. default = \code{TRUE}. Should the a(x) values be calculated from a smoothed M(x) series? In this case, the M(x) series is smoothed within the function for a(x) estimation, but the smoothed M(x) function that was used is not returned. In general, it is better to smooth the M(x) function prior to putting it in this function, because the loess smoother used here has no weights or offset. If this is not possible, loess M(x) smoothing still produces more consistent and less erratic a(x) estimates.
#' @param method either \code{"keyfitz"}, \code{"schoen"}, \code{"preston"} or \code{"midpoint"}. Default = \code{"keyfitz"}, although this is not recommended for abridged ages. See comparisons in examples below.
#' @param sex \code{"female"} or \code{"male"}. default \code{"female"}. This is only used by the \code{"preston"} method and need no be specified for any other method.
#' 
#' @details This function is a wrapper, and it is called by the lifetable function \code{LT()}, although it can be used independently as well. For fuller explanations, see the descriptions and code of the various methods. Formulas are available in the referenced works.
#' 
#' @return \code{ax}, a numeric vector of a(x) values.
#' 
#' @references 
#' Chiang C.L.(1968) Introduction to Stochastic Processes in Biostatistics. New York: Wiley.
#' 
#' Coale Anseley and Paul Demeny, with B Vaughan (1983). Regional Model Life Tables and Stable Populations. New York Academic Press.
#' 
#' Keyfitz, Nathan (1966) A Life Table that Agrees with the Data. Journal of the American Statistical Association, 61 (314):305-12. (as described on page 44-45 of Preston et al (2001). Demography: Measuring and Modelling Population Processes. Blackwell Publishing)
#' 
#' Schoen R. (1978) Calculating lifetables by estimating Chiang\'s a from observed rates. Demography 15: 625-35.
#' 
#' Andreev, Evgueni M and Kingkade, Ward W (2011) Average age at death in infancy and infant mortality level: reconsidering the Coale-Demeny formulas at current levels of low mortality. MPIDR Working Paper WP-2011-016.
#' 
#' @note Be aware that all of the above methods are in some way a hybrid: In the \code{"schoen"} and \code{"keyfitz"} methods, I added procedures to produce values for the final age(s) in a rudimentary way, and in the \code{"preston"} method I also made a rudimentary estimation procedure for a1 - a8 for single age data. For all methods, a0 is calculated using a modified version of the Andreev-Kingkade a0 rule. It is best not to use \code{"keyfitz"} the default method, with abridged age groups.
#' 
#' @seealso This function dispatches to one of four different a(x) estimation functions \code{\link{axMidpoint}} for the "midpoint" method, \code{\link{axSchoen}} for the \code{"schoen"} method, \code{\link{axPreston}} for the "preston" method and \code{\link{axKeyfitz}} for the \code{"keyfitz"} method. Look to these pages for specifics. Compare using the examples below. This function is called by \code{\link{LT}}, a single decrement lifetable function.
#' 
#' @examples # single age comparisons:
#' \dontrun{
#' library(LifeTable)
#' data(UKRmales1965)
#' Mx       <- UKRmales1965[,4]
#' Widths   <- c(rep(1, length(Mx)))
#' ages     <- cumsum(Widths) - Widths
#' axk      <- axEstimate(Mx, Widths, axsmooth = TRUE, method = "keyfitz")
#' axs      <- axEstimate(Mx, Widths, axsmooth = TRUE, method = "schoen")
#' axm      <- axEstimate(Mx, Widths, axsmooth = TRUE, method = "midpoint")
#' axp      <- axEstimate(Mx, Widths, axsmooth = TRUE, method = "preston")
#' 
#' plot(ages, axk, type = 'l', col = "blue", main = "comparing a(x) methods, single ages (smooth)")
#' lines(ages, axs, col = "red", lty = 2)
#' lines(ages, axm, col = "orange", lty = 3)
#' lines(ages, axp, col = "green", lty = 4)
#' text(55, .54, "data from HMD, Ukraine males. 1965", xpd = TRUE)
#' legend("bottom", col = c("blue", "red", "orange", "green"), lty = c(1, 2, 3, 4), legend =c("keyfitz", "schoen", "midpoint", "preston"))
#' ## set axsmooth to FALSE to compare unsmoothed versions of these.
#' 
#' ## abridged 5-year age comparison:
#' data(UKR5males1965)
#' Mx       <- UKR5males1965[, 4]
#' Widths   <- c(1, 4, rep(5, (length(Mx) - 2)))
#' ages     <- cumsum(Widths) - Widths
#' axk      <- axEstimate(Mx, Widths, axsmooth = TRUE, method = "keyfitz")
#' axs      <- axEstimate(Mx, Widths, axsmooth = TRUE, method = "schoen")
#' axm      <- axEstimate(Mx, Widths, axsmooth = TRUE, method = "midpoint")
#' axp      <- axEstimate(Mx, Widths, axsmooth = TRUE, method = "preston")
#' 
#' plot(ages, c(axk[1], axk[2] / 4, axk[3:length(axk)] / 5), type = 'l', col = "blue", main = "comparing a(x) methods, single ages (smooth)")
#' lines(ages, c(axs[1], axs[2] / 4, axs[3:length(axs)] / 5), col = "red", lty = 2)
#' lines(ages, c(axm[1], axm[2] / 4, axm[3:length(axm)] / 5), col = "orange", lty = 3)
#' lines(ages, c(axp[1], axp[2] / 4, axp[3:length(axp)] / 5), col = "green", lty = 4)
#' text(55, .54, "data from HMD, Ukrain males. 1965", xpd = TRUE)
#' legend("bottom", col = c("blue", "red", "orange", "green"), lty = c(1, 2, 3, 4), legend = c("keyfitz", "schoen", "midpoint", "preston"))
#' ## set axsmooth to FALSE to compare unsmoothed versions of these.
#' ## here, you can see why Preston advises against the Keyfitz when age intervals are not equal. 
#' }
#' 
#' @export

axEstimate <-
function(Mx, n, axsmooth = TRUE, method = "keyfitz", sex,verbose){
    if (! method %in% c("keyfitz", "schoen", "midpoint", "preston")){
        stop("ax method specified is not valid, must be 'keyfitz','schoen','midpoint' or 'preston'")
    }
	if (method == "keyfitz"){
		ax <- axKeyfitz(Mx, n, axsmooth)
	}
	if (method == "schoen"){
		ax <- axSchoen(Mx, n, axsmooth)
	}
	if (method == "midpoint"){
		ax <- axMidpoint(Mx, n)
	}
	if (method == "preston"){
		if (missing(sex)) {
            Verb(verbose,"\nWarning: You didn't specify 'sex', assumed 'female'!\n")
            sex <- "female"
        }
        ax <- axPreston(Mx, n, axsmooth, sex)
	}
    # use Andreev-Kingkade a0 formula.
    ax[1]  <- AKm02a0(Mx[1],sex)
	return(ax)
}

