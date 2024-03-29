% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/axMidpoint.R
\name{axMidpoint}
\alias{axMidpoint}
\title{axMidpoint approximates a(x) using interval midpoints}
\usage{
axMidpoint(Mx, n)
}
\arguments{
\item{Mx}{a numeric vector of central death rates, defined as deaths/exposure.}

\item{n}{a numeric vector of the age-interval widths. Must be the same length as Mx.}
}
\value{
returns a numeric vector of ax values
}
\description{
approximates a(x) using interval midpoints (except a0, which is estimated using a rule-of-thumb formula scaled by m0)
}
\details{
All values will be in the exact middle of the age interval, except for a0, which uses: \code{a0 = .07 + 1.7 * M0}.
}
\references{
formula for a0 as suggested by:
Keyfitz, N. (1970), Finding probabilities from observed rates, or how to make a life table. The American Statistician (1970), pp. 28-33
}
\seealso{
See Also \code{\link{axEstimate}}, a wrapper function for this and 3 other a(x) estimation methods (\code{\link{axKeyfitz}},\code{\link{axSchoen}} and \code{\link{axPreston}}).
}
