axMidpoint <-
function(Mx,n){
	ax <- .5*n		
	ax[1] <- .07+1.7*Mx[1]
	return(ax)
}

