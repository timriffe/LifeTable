\name{LT}
\alias{LT}
\title{A function to build lifetables and extract certain lifetable statistics from data.}
\usage{
  LT(Nx, Dx, Mx, ages = "auto", type = "single-age",
    axmethod = "keyfitz", sex = "female", mxsmooth = TRUE,
    axsmooth = TRUE, n = "auto", radix = 1, verbose = TRUE)
}
\arguments{
  \item{Nx}{numeric vector of population exposures by age.}

  \item{Dx}{numeric vector of death counts by age.}

  \item{Mx}{numeric vector of central death rates (assumed
  in the function to be the lifetable m(x)), calculated as
  deaths/exposure.}

  \item{ages}{optional, default = \code{"auto"}. If
  \code{"auto"}, the function tried to guess whether you
  have single-age data or 5-year abridged age groups. These
  are used in the lifetable labels, and do not enter into
  calculations. Otherwise, the user can specify a vector
  (character or numeric).}

  \item{type}{either \code{"single-age"} (default) or
  \code{"abridged"}. These are the only two options. Five
  year abridged age groups presumes that you have age
  intervals of 1,4,5,5,5,5.... Send requests to the author
  if this isn't flexible enough.}

  \item{axmethod}{either \code{"keyfitz"}, \code{"schoen"},
  \code{"preston"} or \code{"midpoint"}. Default =
  \code{"keyfitz"}, although this is not recommended for
  abridged ages. See comparisons in
  \code{\link{axEstimate}} examples. The user can also
  supply a numeric vector of a(x) values here (e.g. from a
  different estimation procedure or from a different
  population).}

  \item{sex}{either \code{"male"} or \code{"female"}
  (default). It is only necessary to specify this if
  \code{"preston"} is the \code{axmethod}. It does not
  affect any other lifetable calculations.}

  \item{mxsmooth}{logical, default = \code{TRUE}. Should
  the mx vector be smoothed? If \code{TRUE} and both
  \code{Nx} and \code{Dx} vectors are supplied (the ideal
  case), smoothing is done using the function
  \code{Mort1Dsmooth()} from Giancarlo Camarda's
  \code{MortalitySmooth} package. In this case, \code{Dx}
  values are smoothed using \code{log(Nx)} as an offset,
  and all other items are the function defaults. If
  \code{Mx} is provided instead of \code{Nx} and \code{Dx}
  a loess smoother is used, \code{loess}, with span set to
  .15 for single age data and .4 for 5-year abridged data.
  If these smoothing procedures are not satisfactory, the
  user may wish to pre-process the Mx estimate and specify
  \code{mxsmooth = FALSE}, or else leave it rough.}

  \item{axsmooth}{logical, default = \code{TRUE}. Ignored
  if \code{mxsmooth = TRUE}. Should the a(x) values be
  calculated from a smoothed M(x) series? In this case, the
  M(x) series is smoothed within the \code{axEstimate()}
  function for a(x) estimation, but the smoothed M(x)
  function that was used is not returned. In general, it is
  better to smooth the M(x) function prior to putting it in
  this function, because the loess smoother used here has
  no weights or offset. If this is not possible, loess M(x)
  smoothing still produces more consistent and less erratic
  a(x) estimates. If \code{mxsmooth = FALSE} and
  \code{axsmooth = TRUE}, the Mx series is only smoothed
  for use in a(x) estimation, and does not affect any other
  lifetable calculations that are dependent on Mx.}

  \item{n}{numeric vector of interval widths or
  \code{"auto"} (default).}

  \item{radix}{The lifetable starting population at age 0,
  l(0). default = 1. Other common values are 1000 and
  100000, although any value may be given.}

  \item{verbose}{logical, default = \code{TRUE}. Should
  informative but possibly annoying messages be returned
  when the function does something that you might want to
  know about?}
}
\value{
  a list \itemize{ \item \code{LT} A \code{data.frame} with
  11 columns and as many rows as you have ages. The columns
  are "Age" (character age labels, i.e. with "+"), "ages"
  (numeric), \code{"mx", "ax", "qx", "px", "lx", "dx",
  "Lx", "Tx", "ex"}. All the individual components of the
  lifetable can be called and are unrounded when
  individually called. Some additional values are also
  available: \item \code{Age} character vector of ages.
  Denotes intervals in the case of an abridged table.
  \item \code{ages} numeric vector of ages. Left side of
  interval.  \item \code{mx} the lifetable mx (may differ
  from a given Mx).  \item \code{ax} Chiang's ax, either
  given by the user or estimated in \code{axEstimate}.
  \item \code{qx} typical lifetable qx. Death probability
  for interval x, x + n.  \item \code{lx} typical lifetable
  lx. Number of radix individuals entering age x. l(0) =
  the radix population.  \item \code{dx} typical lifetable
  dx. When \code{radix = 1} (default), this is the
  probability density function of deaths.  \item \code{Lx}
  typical lifetable Lx. Lifetable exposure for the interval
  x,x+n.  \item \code{Tx} typical lifetable Tx. Total
  number of years remaining to be lived by the cohort
  entering age x.  \item \code{ex} typical lifetable ex.
  Life remaining life expectancy at age x. e(0) = life
  expectancy at birth. Two other estimates of e(0) are
  given below.  \item \code{Sx} probability of surviving
  from age x until age x + n (l_{x+n}/l_{x}).  \item
  \code{Widths} vector of age intervals (n).  \item
  \code{e0est} a matrix with three elements, each e(0)
  estimated in a different way (formula in column headers).
  }
}
\description{
  Accepts either single age or 5-year abridged data.
  Accepts many optional arguments, such as differing
  methods for a(x) estimation, optional smoothing for M(x)
  or a(x) values, a changeable radix.
}
\details{
  Either \code{Nx} must be specified along with \code{Dx},
  *or* \code{Mx} must be specified directly. If smoothing
  is used, it is better to specify both \code{Nx} and
  \code{Dx}, since the \code{Nx} vector can be used as an
  offset in the \code{MortalitySmooth} smoother.
}
\examples{
\dontrun{
library(LifeTable)
data(UKRmales1965)
head(UKRmales1965)
Nx <- UKRmales1965[, 3]
Dx <- UKRmales1965[, 2]

# comparing different estimates of life expectancy based on whether mx or ax smoothing
# is used  (no major differences). Could also just specify Mx directly instead of both Nx and Dx.
# if you want to smooth the mx, better to use Nx and Dx together instead of Mx.
LT(Nx, Dx, type = "single-age", axsmooth = TRUE)$e0est
LT(Nx, Dx, type = "single-age", axsmooth = FALSE)$e0est
LT(Nx, Dx, type = "single-age", mxsmooth = FALSE, axsmooth = FALSE)$e0est
LT(Nx, Dx, type = "single-age", mxsmooth = FALSE, axsmooth = TRUE)$e0est

# now comparing life expectancy estimates depending on the ax estimation method used (no major differences)
LT(Nx, Dx, type = "single-age", axmethod = "keyfitz", mxsmooth = FALSE)$e0est
LT(Nx, Dx, type = "single-age", axmethod = "schoen")$e0est
LT(Nx, Dx, type = "single-age", axmethod = "preston")$e0est
# there are tons more combinations of with/without smoothing, ax methods...

# here a graph with the major functions:
LT1 <- LT(Nx, Dx, type = "single-age", axmethod = "keyfitz", mxsmooth = TRUE)
plot(LT1$ages, LT1$lx, type = 'l', col = "blue", main = "major lifetable functions, UKR 1965 (HMD)",
     xlab = "age", ylab = "lx ; dx*10 ; mux")
lines(LT1$ages, LT1$mx, col = "red")
lines(LT1$ages, LT1$dx * 10, col = "orange", type = 'l')
legend("topright", lty = 1, col = c("blue", "orange", "red"), legend = c("l(x)", "d(x) * 10", "mux"), bg = "white")

######### with 5-year age groups:
data(UKR5males1965)
head(UKR5males1965)
Nx <- UKR5males1965[, 3]
Dx <- UKR5males1965[, 2]

LT5 <- LT(Nx, Dx, type = "abridged", axmethod = "schoen", mxsmooth = FALSE, axsmooth = FALSE)
plot(LT5$ages, LT5$lx, type = 'l', col = "blue", main = "major lifetable functions, UKR 1965 (HMD)", xlab = "age", ylab = "lx ; dx * 10 ; mux")
lines(LT5$ages, LT5$mx, col = "red")
lines(LT5$ages, c(LT5$dx[1], LT5$dx[2] / 4, LT5$dx[3:24] / 5) * 10, col = "orange", type = 'l')
legend("topright", lty = 1, col = c("blue", "orange", "red"), legend = c("l(x)", "d(x) * 10", "mux"), bg = "white")

######### just want a nice Sx vector for your Leslie matrix?:
Sx <- LT1$Sx
plot(Sx)
# don't worry, the last Sx is already the ratio of the last Tx divided by the penultimate Tx
}
}
\references{
  The main reference for this function has been:

  Preston et al (2001). Demography: Measuring and Modelling
  Population Processes. Blackwell Publishing

  ax estimation also received input from:

  Chiang C.L.(1968) Introduction to Stochastic Processes in
  Biostatistics. New York: Wiley.

  Coale Anseley and Paul Demeny, with B Vaughan (1983).
  Regional Model Life Tables and Stable Populations. New
  York Academic Press.\\ Keyfitz, Nathan (1966) A Life
  Table that Agrees with the Data. Journal of the American
  Statistical Association, 61 (314):305-12. As described on
  page 44-45 of: Schoen R. (1978) Calculating lifetables by
  estimating Chiang's a from observed rates. Demography 15:
  625-35.

  function calls \code{MortalitySmooth}: Carlo G Camarda
  (2009) MortalitySmooth: Smoothing Poisson counts with
  P-splines. (version 2.3 at the time of this writing)
  \url{http://CRAN.R-project.org/package=MortalitySmooth}.
}
\seealso{
  \code{\link{MortalitySmooth}}, \code{\link{axEstimate}}.
}

