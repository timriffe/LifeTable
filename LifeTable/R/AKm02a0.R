
#'
#' @title AKm02q0 derive q0 from m0 using the Andreev-Kingkade rule of thumb.
#' 
#' @description Derive m0 from q0 according to the relevant segment of the Andreev-Kingkade formula. This is elegant because it's an analytic solution, but ugly because, man, look at it. Carl Boe got this from Maple I think. This formula is only necessary because AK start with q0 whereas the HMD starts with m0, so we needed to adapt. This is an auxiliary function, and not likely needed for direct use.
#' 
#' @param m0 the event / exposure infant mortality rate (not IMR)
#' @param constant the intercept of the relevant Andreev-Kingkade segment
#' @param slope the slope of the relevant Andreev-Kingkade segment
#' 
#' @return q0 the estimate of q0 according to the identity between a0, m0, q0
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#' 
#' @export
#' 
#' @details This is based on an MPIDR Working Paper: Andreev, Evgueni M and Kingkade, Ward W (2011) "Average age at death in infancy and infant mortality level: reconsidering the Coale-Demeny formulas at current levels of low mortality". short link: http://goo.gl/b5m5pg.

AKm02q0 <- function(m0, constant, slope){
    -1 / slope / m0 * (-m0 +  (m0 * constant) - 0.1e1 + sqrt(((m0 ^ 2) - 2 * constant * (m0 ^ 2) + 2 * m0 + (constant ^ 2) * (m0 ^ 2) - 2 *  (m0 * constant) + 1 - 4 * slope * (m0 ^ 2)))) / 2
}

#' @title \code{AKm02a0} estimates a0 using the Andreev-Kinkade rule of thumb.
#'
#' @description
#' \code{AKm02a0} is an auxiliary function used by version 6 of the four HMD lifetable functions, \code{ltper_AxN()}, \code{ltcoh_AxN()}, \code{ltperBoth_AxN()}, \code{ltcohBoth_AxN()}. This function calls \code{AKm02q0()} to help get the work done, since the HMD needed to adapt the Andreev-Kingkade formulas to work with the period lifetable flow.
#'
#' @param m0 the event / exposure infant mortality rate (not IMR)
#' @param sex either "male" or "female"
#' 
#' @return a0, the estimated average age at death of those dying in the first year of life, either a single value or a vector of a_0 values.
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#' 
#' @export
#' 
#' @details This is based on an MPIDR Working Paper: Andreev, Evgueni M and Kingkade, Ward W (2011) "Average age at death in infancy and infant mortality level: reconsidering the Coale-Demeny formulas at current levels of low mortality". short link: http://goo.gl/b5m5pg.
AKm02a0 <- function(m0,sex="male"){
    sex <- rep(sex,length(m0))
    ifelse(sex == "male",
            ifelse(m0 < 0.02306737, 0.1493 - 2.0367 * AKm02q0(m0, 0.1493, -2.0367),
                    ifelse(m0 < 0.0830706, 0.0244 + 3.4994 * AKm02q0(m0, 0.0244, 3.4994), .2991)),
            ifelse(m0 < 0.01725977, 0.1490 - 2.0867 * AKm02q0(m0, 0.1490, -2.0867),
                    ifelse(m0 < 0.06919348, 0.0438 + 4.1075 * AKm02q0(m0, 0.0438, 4.1075), 0.3141))
    )
}
