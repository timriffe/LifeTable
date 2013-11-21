\docType{package}
\name{LifeTable}
\alias{LifeTable}
\alias{LifeTable-package}
\title{LifeTable, a package with a small set of useful lifetable functions.}
\description{
  The main function, \code{LT()}, is rather versatile,
  accepting single-age or 5-year abridged data, and with a
  variety of optional arguments. There is also a
  sub-function, \code{axEstimate()}, which serves as a
  wrapper for 4 different a(x) estimation methods. All
  functions take as their basic argument the M(x) vector,
  although the \code{LT()} function lets you start simply
  with a vector of exposures (Nx) and a vector of death
  counts (Dx) (better for smoothing). From that, all the
  other lifetable functions are estimated using a series of
  identities. Either an entire composed and rounded table
  can be returned, or else the individual un-rounded
  components. See further examples under individual
  function pages. Example data is also provided, as
  extracted from the Human Mortality Database in May, 2010.
  These data are for Ukrainian males in 1965, and are given
  in single ages and in 5-year abridged ages as well.
}
\examples{
\dontrun{
# here a graph with the major functions:
library(LifeTable)
data(UKRmales1965)
head(UKRmales1965)
Nx   <- UKRmales1965[, 3]
Dx   <- UKRmales1965[, 2]
LT1  <- LT(Nx, Dx, type = "single-age", axmethod = "keyfitz", mxsmooth = TRUE)
plot(LT1$ages, LT1$lx, type = 'l', col = "blue", main = "major lifetable functions, UKR 1965 (HMD)",
     xlab = "age", ylab = "lx ; dx*10 ; mux")
lines(LT1$ages, LT1$mx, col = "red")
lines(LT1$ages, LT1$dx * 10, col = "orange")
legend("topright", lty = 1, col = c("blue", "orange", "red"), legend = c("l(x)", "d(x) * 10", "mux"), bg = "white")
}
}
\references{
  The main reference for this function and its subfunctions
  has been:\\ Preston et al (2001). Demography: Measuring
  and Modelling Population Processes. Blackwell Publishing
  ax estimation also received input from:

  Chiang C.L.(1968) Introduction to Stochastic Processes in
  Biostatistics. New York: Wiley.

  Coale Anseley and Paul Demeny, with B Vaughan (1983).
  Regional Model Life Tables and Stable Populations. New
  York Academic Press.

  Keyfitz, Nathan (1966) A Life Table that Agrees with the
  Data. Journal of the American Statistical Association, 61
  (314):305-12. As described on page 44-45 of: Schoen R.
  (1978) Calculating lifetables by estimating Chiang\'s a
  from observed rates. Demography 15: 625-35.

  funciton calls \code{MortalitySmooth}: Carlo G Camarda
  (2009) MortalitySmooth: Smoothing Poisson counts with
  P-splines. (version 2.3 at the time of this writing)
  \url{http://CRAN.R-project.org/package=MortalitySmooth}
}
\seealso{
  \code{\link[MortalitySmooth]{Mort1Dsmooth}},\code{\link[stats]{loess}}
}

