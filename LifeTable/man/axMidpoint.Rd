\name{axMidpoint}
\alias{axMidpoint}
\title{approximates a(x) using interval midpoints}
\description{
approximates a(x) using interval midpoints (except a0, which is estimated using a rule-of-thumb formula scaled by m0)
}
\usage{
axMidpoint(Mx, n)
}
\arguments{
  \item{Mx}{a numeric vector of central death rates, defined as deaths/exposure.}
  \item{n}{a numeric vector of the age-interval widths. Must be the same length as Mx.}
}
\details{all values will be in the exact middle of the age interval, except for a0, which uses: a0 = .07 + 1.7 * M0}
\value{returns a numeric vector of ax values}
\references{
formula for a0 as suggested by:
Keyfitz, N. (1970), Finding probabilities from observed rates, or how to make a life table. The American Statistician (1970), pp. 28–33
}
\author{
Tim Riffe
}

\seealso{
See Also \code{\link{axEstimate}}, a wrapper function for this and 3 other a(x) estimation methods (\code{\link{axKeyfitz}},\code{\link{axSchoen}} and \code{\link{axPreston}}).
}
\examples{

## The function is currently defined as
function(Mx,n){
	ax <- .5*n		
	ax[1] <- .07+1.7*Mx[1]
	return(ax)
  }
}