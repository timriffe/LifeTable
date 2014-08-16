# a series of identities.
###############################################################################

#' @title dx2lx coverts single-age dx to single age lx.
#' 
#' @description This assumes 1) that dx sums to the radix exactly and 2) that the last element of dx is also the last element of lx.
#' 
#' @param dx vector of single age \code{dx} from the lifetable. This is the density function, and it has various interpretations.
#' 
#' @return lx vector of the survival function, of same radix as \code{dx}
#' 
#' @export
dx2lx <- function(dx){
    rev(cumsum(rev(dx)))
}

#' @title lx2qx coverts single-age lx to single age qx.
#' 
#' @description This assumes that lx contains all ages, and it also closes the \code{qx} vector out with a value of 1, which is common practice. \code{lx} may be of any radix. No checking is done for monotanicity or for 0s in the denominator. Please be aware that 0s in \code{lx}, except for in the last age, will mess this up, and smooth your data if necessary.
#' 
#' @param lx vector of the survival function from age 0 to omega.
#' 
#' @return qx vector of single age conditional death probabilities
#' 
#' @export
lx2qx <- cmpfun(function(lx){
            c(1-(lx[2:length(lx)] / lx[1:(length(lx)-1)]),1)
        })

#' @title lqx2lx coverts single-age qx to single age lx.
#' 
#' @description This assumes that \code{qx} starts at age 0. \code{lx} is created with a radix of 1, so rescale outside of this function if needed. No \code{NA}s or \code{Inf}s!
#' 
#' @param qx vector of single age conditional death probabilities
#' 
#' @return lx vector of the radix-1 survival function from age 0 to omega.
#' 
#' @export
qx2lx <- cmpfun(function(qx){
            c(1,cumprod(1-qx)[-length(qx)])
        })

#' @title lx2dx coverts single-age lx to single age dx.
#' 
#' @description \code{dx} will sum to the \code{lx} radix. Assumes but does not check for monontonicity in \code{lx}, so negative values can be produced, but are actually impossible. 
#' 
#' @param lx vector of the single-age survival function values
#' 
#' @return dx vector of single age \code{dx} from the lifetable. This is the density function, and it has various interpretations. Sums to radix.
#' 
#' @export
lx2dx <- function(lx){
    -diff(c(lx,0))
}
#' @title dx2qx coverts single-age dx to single age qx.
#' 
#' @description This assumes 1) that dx sums to the radix exactly and 2) that the last element of dx is also the last element of lx.
#' 
#' @param dx vector of single age \code{dx} from the lifetable. This is the density function, and it has various interpretations.
#' 
#' @return lx vector of the survival function, of same radix as \code{dx}
#' 
#' @export
dx2qx <- function(dx){
    lx2qx(dx2lx(dx))
}
#' @title qx2dx coverts single-age qx to single age dx.
#'
#' @description The resulting \code{dx} vector will sum to 1 (correspond to a radix-1 lifetable)
#' 
#' @param qx vector of single age conditional death probabilities
#' 
#' @return dx vector of single age \code{dx} from the lifetable.
#' 
#' @export
qx2dx <- function(qx){
    lx2dx(qx2lx(qx))
}
#' @title dx2ax creates ax vector given \code{dx}
#'
#' @description \code{ax} is just a vector of .5, using single age mid-interval, except for age 0 which uses the HMD Methods Protocol version 6 beta method based on Andreev Kingkade (2011)
#' 
#' @param dx vector of single age \code{dx} from the lifetable.
#' @param sex \code{"m"} or \code{"f"}. The a0 method is sex-specific.
#' 
#' @return ax vector of Chiang's ax.
#' 
#' @export
dx2ax <- cmpfun(function(dx,sex= "m"){
            qx <- dx2qx(dx)
            c(AKq02a0(qx[1],sex),rep(.5,length(dx)-1))
        })
#' @title dx2mx creates mx vector given \code{dx}
#'
#' @description Calculates \code{ax} ad \code{qx} and then uses the identity between \code{ax}, \code{qx} and \code{mx} to calculate \code{mx}.
#' 
#' @param dx vector of single age \code{dx} from the lifetable.
#' @param sex \code{"m"} or \code{"f"}. The a0 method is sex-specific.
#' 
#' @return  mx vector of lifetable death rates.
#' 
#' @export
dx2mx <- function(dx,sex="m"){
    qx <- dx2qx(dx)
    ax <- dx2ax(dx, sex)
    qx / (1 - (1 - ax) * qx) 
}
#' @title dx2e0 calculates e0 given \code{dx}
#'
#' @description Calculates life expectancy at birth given \code{dx}.
#' 
#' @param dx vector of single age \code{dx} from the lifetable.
#' @param sex \code{"m"} or \code{"f"}. The a0 method is sex-specific.
#' 
#' @return e0 single value, not a vector.
#' 
#' @export
dx2e0 <- function(dx, sex = "m"){
    lx              <- dx2lx(dx)
    N               <- length(dx)
    ax              <- dx2ax(dx, sex)
    Lx              <- lx - (1 - ax) * dx 
    Lx[N]           <- lx[N] * ax[N]
    sum(Lx)
}

