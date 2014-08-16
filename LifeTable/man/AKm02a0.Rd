\name{AKm02a0}
\alias{AKm02a0}
\title{\code{AKm02a0} estimates a0 using the Andreev-Kinkade rule of thumb.}
\usage{
  AKm02a0(m0, sex = "male")
}
\arguments{
  \item{m0}{the event / exposure infant mortality rate (not
  IMR)}

  \item{sex}{either "male" or "female"}
}
\value{
  a0, the estimated average age at death of those dying in
  the first year of life, either a single value or a vector
  of a_0 values.
}
\description{
  \code{AKm02a0} is an auxiliary function used by version 6
  of the four HMD lifetable functions, \code{ltper_AxN()},
  \code{ltcoh_AxN()}, \code{ltperBoth_AxN()},
  \code{ltcohBoth_AxN()}. This function calls
  \code{AKm02q0()} to help get the work done, since the HMD
  needed to adapt the Andreev-Kingkade formulas to work
  with the period lifetable flow.
}
\details{
  This is based on an MPIDR Working Paper: Andreev, Evgueni
  M and Kingkade, Ward W (2011) "Average age at death in
  infancy and infant mortality level: reconsidering the
  Coale-Demeny formulas at current levels of low
  mortality". short link: http://goo.gl/b5m5pg.
}
\author{
  Tim Riffe \email{triffe@demog.berkeley.edu}
}

