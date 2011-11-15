\name{UKR5males1965}
\alias{UKR5males1965}
\docType{data}
\title{
Example data for Ukraine males, 1965, extracted from the Human Mortality Database.
}
\description{
5-year adbridged data from the HMD for Deaths (count) Population (exposure) and central death rates, death/exposure. Data
}
\usage{data(UKR5males1965)}
\format{
  The format is:
 num [1:24, 1:4] 0 1 5 10 15 20 25 30 35 40 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:4] "age" "Dx" "Nx" "Mx"
}
\details{
All four columns are numeric. Mx provided is not exactly equal to D(x)/N(x).
}
\source{
Human Mortality Database, as downloaded in May, 2010.
}
\references{
Human Mortality Database. University of California, Berkeley (USA), and Max Planck Institute for Demographic Research (Germany). Available at www.mortality.org or www.humanmortality.de (data downloaded on 10 May, 2010). 
}
\examples{
## see examples using this data in \code{\link{LT()}} and in \code{\link{axEstimate}}.
}
\keyword{datasets}
