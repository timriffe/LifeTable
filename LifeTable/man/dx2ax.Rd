% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lx2dx.R
\name{dx2ax}
\alias{dx2ax}
\title{dx2ax creates ax vector given \code{dx}}
\usage{
dx2ax(dx, sex = "m")
}
\arguments{
\item{dx}{vector of single age \code{dx} from the lifetable.}

\item{sex}{\code{"m"} or \code{"f"}. The a0 method is sex-specific.}
}
\value{
ax vector of Chiang's ax.
}
\description{
\code{ax} is just a vector of .5, using single age mid-interval, except for age 0 which uses the HMD Methods Protocol version 6 beta method based on Andreev Kingkade (2011)
}
