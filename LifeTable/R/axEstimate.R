axEstimate <-
function(Mx,n,axsmooth=TRUE,method="keyfitz",sex){
	if (method=="keyfitz"){
		ax <- axKeyfitz(Mx,n,axsmooth)
	}
	if (method=="schoen"){
		ax <- axSchoen(Mx,n,axsmooth)
	}
	if (method=="midpoint"){
		ax <- axMidpoint(Mx,n)
	}
	if (method=="preston"){
		if (missing(sex)) sex <- "female"
		ax <- axPreston(Mx,n,axsmooth,sex)
	}
	return(ax)
}

