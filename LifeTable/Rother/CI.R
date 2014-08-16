# original code submitted by Adrien Remund,
# modified by Tim Riffe

LE <- compiler::cmpfun(function(mx, ax, Widths){
    if (is.null(dim(mx))){
        mx <- as.matrix(mx)
    }
    if (is.null(dim(ax))){
        ax <- as.matrix(ax)
    }
    if (is.null(dim(Widths))){
        Widths <- as.matrix(Widths)
    }
    N             <- nrow(mx)

    qx            <- (Widths * mx) / (1 + (Widths - ax) * mx)
    qx[N, ]       <- 1
    qx[qx > 1]    <- 1
    px            <- 1 - qx
    lx            <- apply(px, 2,function(.px,.N){
                        c(1,cumprod(.px)[-.N])
                     }, .N = N)
    dx            <- lx * qx
    dx[N, ]       <- 1 - colSums(dx[-N, , drop = FALSE])
    Lx            <- rbind(Widths[1:(N - 1), , drop = FALSE] * lx[2:N, , drop = FALSE] + ax[1:(N - 1), , drop = FALSE] * dx[1:(N - 1), , drop = FALSE], 
                           lx[N, , drop = FALSE] * ax[N, , drop = FALSE]
                   )

    Lx[is.infinite(Lx)] <- 1
    Lx[is.na(Lx)] <- 0
    Tx            <- apply(Lx, 2,function(.Lx){
                        rev(cumsum(rev(.Lx)))
                     })
    ex            <- Tx / lx
    return(ex[1, ])
})

# This function estimates the confidence intervals for e0
# there are two available methods: "chiang" and "MC", which are based respectively on 
# (Chiang 1978) and (Andreev & Shkolnikov 2010)
# see also http://www.demogr.mpg.de/en/projects_publications/publications_1904/mpidr_technical_reports/spreadsheet_for_calculation_of_confidence_limits_for_any_life_table_or_healthy_life_table_quantity_3853.htm
# the other parameters are alpha (the confidence level) and iter (the number of iteration for the "MC method)
# NB: Chiang published two variants of his method. I used the second one (the first one is commented in the code).
e0.CI <- function(Nx, Dx, mx, sex, ages = 0:(length(mx)-1), ax = NULL, iter = 1000, alpha = 0.05){
    
    # box to hold results
    
    N                   <- length(mx)
    Widths              <- diff(ages)
    Widths              <- c(Widths, Widths[N - 1])
    
    # if ax not provided, assume midpoint (HMD period method)
    if (is.null(ax)){
        ax  <- axEstimate(Mx = mx, 
                        n = Widths, 
                        axsmooth = axsmooth, 
                        method = "midpoint", 
                        sex = sex, 
                        verbose = verbose)
    }
    e0          <- LE(mx, ax, Widths)
    # MC
    # 1) generate random mx
    rmx         <- t(mapply(function(.mx, .Nx, iter){
                        rpois(iter, .mx * .Nx)
                    }, mx, ceiling(Nx), iter = iter)) / Nx
    re0         <- LE(rmx, 
                      matrix(ax[row(rmx)],ncol=iter),
                      matrix(Widths[row(rmx)],ncol=iter)) # take random mx as matrix
    MCCI    <- c(lower=quantile(re0, probs = alpha / 2), upper = quantile(re0, probs = 1 - alpha / 2))
    
    # Chiang's method
    if(any(Dx < 5)){
        message("Number of deaths too low to use Wald's method.")
    }

    qx      <- (Widths * mx) / (1 + (Widths - ax) * mx)
    qx[N]   <- 1
    qx[qx > 1] <- 1
    px      <- 1 - qx
    lx      <- c(1,cumprod(px)[-N]) # vectorized, was loop
    dx      <- -diff(lx)
    dx[N]   <- lx[N]
 
    Lx      <- c(Widths[1:(N - 1)] * lx[2:N] + ax[1:(N - 1)] * dx[1:(N - 1)], lx[N] * ax[N])
    Lx[is.infinite(Lx)] <- 1
    Lx[is.na(Lx)] <- 0
    Tx      <- rev(cumsum(rev(Lx)))
    ex      <- Tx/lx
    varqx   <- (mx * (1 - ax * mx)) / (Nx * (1 + (1 - ax) * mx) ^3) # Chiang 2
    yx      <- lx^2 * ((1 - ax) + c(ex[2:N], 0)) ^2 * varqx
    
    cumyx   <- rev(cumsum(rev(yx)))
    varex   <- cumyx / (lx^2)
    seex    <- sqrt(varex)
    
    ChiangCI <- c(ex[1] - qnorm(1 - alpha / 2) * seex[1], ex[1] + qnorm(1 - alpha / 2) * seex[1])
    
return(list(e0 = e0,MCCI = MCCI, ChiangCI = ChiangCI))
}

e0.CI(Nx, Dx,mx,ax,sex="male")
