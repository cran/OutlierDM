.onAttach <- function(libname, pkgname) {
    ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    packageStartupMessage("")
    packageStartupMessage("Package ", pkgname, " (",ver,") loaded.")
}

###########################################################
## show methods for the OutlierDM class
setGeneric("show")
setMethod("show", "OutlierDM", function(object){
			options(warn = -1)
			cat("\nCall:\n")
			print(object@call)
			cat("\n")
			cat("Outlier Detection for Multiplicative High-Throughput Data\n\n")
			if( object@type %in% c("pair","diff","proj")){
				cat(" Method: ")
				switch(object@method, 
					constant = cat("constant quantile regression"),
					linear = cat("linear quantile regression"),
					nonlin = cat("non-linear quantile regression"),
					nonpar = cat("non-parametric quantile regression")
				)
				cat(" (",object@method,")","\n")
			}
			cat(" Type: ")
			switch(object@type, 
				dixon = cat("Dixon's Q-test"),
				grubbs = cat("Grubbs' Test"),
				pair = cat("pairwise approach"),
				diff = cat("difference approach"),
				proj = cat("projection approach")
			)
			cat(" (",object@type,")","\n")
			if( object@type %in% c("pair","diff","proj")){
				cat(" k: ",object@k,"\n")
				cat(" Upper Quantile: ", object@contrl.para$Upper,"\n")
				cat(" Lower Quantile: ",object@contrl.para$Lower, "\n")
			}
			cat(" Number of Rows: ", nrow(object@res),"\n")
			cat(" Number of Outliers: ", object@n.outliers,"\n\n")
			cat(" centering: ",object@contrl.para$centering,"\n")
			cat(" transformation: ", object@contrl.para$trans ,"\n")
			cat("\n Head of Output Data\n" )
			print(object@res[1:10,])
			cat("To see the full output, use a command, 'your_object_name@res'. \n")
			}
)

####################################################################
## plot methods for the OutlierDM class
setGeneric("plot")
setMethod("plot", signature(x = "OutlierDM", y = "missing"), 
	function(x, y = NA, pch = 20, cex = 0.5, xlab = "A", ylab = "M", legend.use = TRUE, ...) {
		options(warn = -1)
		res <- x@res
		#Set parallel computing 
		if (x@type == "proj") {
			plot(res$A, res$M, pch = pch, cex = cex, xlab = xlab, ylab = ylab, 
					col = ifelse(res$Outlier, "red", "black"), ...)
			i <- sort.list(res$A)
			lines(res$A[i], res$UB[i], col = "red", lty = 1)
			lines(res$A[i], res$Q3[i], col = "blue", lty = 2)
			if (legend.use) legend("topright", c("Upper Bound", "Q3"), lty = 1:2, col= c("red", "blue"))
		}
		else if(x@type == "diff") {
			## Need modify t
			## now just work when n = 3
			plot(res$A, res$M1, pch = pch, cex = cex, xlab = xlab, ylab = ylab, 
					col = ifelse(res$Outlier, "red", "black"), ...)		
			points(res$A, res$M2, pch = pch, cex = cex, xlab = xlab, ylab = ylab, 
					col = ifelse(res$Outlier, "red", "black"))		
			points(res$A, res$M3, pch = pch, cex = cex, xlab = xlab, ylab = ylab, 
					col = ifelse(res$Outlier, "red", "black"))		
			i <- sort.list(res$A)
			lines(res$A[i], res$UB[i], col = "red", lty = 1)
			lines(res$A[i], res$Q3[i],, col = "blue", lty = 2)
			lines(res$A[i], res$LB[i],, col = "red", lty = 1)
			lines(res$A[i], res$Q1[i],, col = "blue", lty = 2)
			if (legend.use) legend("topright", c("Upper & Lower Bound", "Q3 and Q1"), lty = 1:2, col= c("red", "blue"))
		}
		else if(x@type == "pair"){
			#predefined functions
			textPanel <- function(x = 0.5, y = 0.5, txt, cex, font)
				text(x, y, txt, cex = cex, font = font)
				
			localAxis <- function(side, x, y){
				if(side %% 2 == 1) Axis(x, side = side)
				else Axis(y, side = side)
			}

			x.pair <- x@x.pair
			# Set dimensions
			n.dim <- ncol(sapply(x@x.pair,dim))
			labels <- paste("N", 1:n.dim, sep="")

			opar <- par(mfrow = c(n.dim, n.dim), mar = rep.int( 1/2, 4), oma = c(4,4,6,4)) 
			on.exit(par(opar))
			dev.hold()
			on.exit(dev.flush(), add = TRUE)
			
			for( i in 1L:n.dim) {
				for( j in 1L:n.dim) {
					## set up the range both xlim and ylim
					## plot the empty box
					k <- 1
					with(x.pair[[k]], plot(A, M, type ="n", axes = FALSE, xlab = "", ylab = ""))
					if(i == j || i < j) {
						box()
						# define axis on matrix cell
						#if(i == 1 &&  !(j %% 2))
						#	localAxis(3, res[,j], res[,i])
						#if(i == n.dim && j %% 2)
						#	localAxis(1, res[,j], res[,i])
						#if(j == 1 &&  !(j %% 2))
						#	localAxis(2, res[,i])
						#if(j == n.dim && j %% 2)
						#	localAxis(4, res[,i])						

						mfg <- par("mfg")
						if(i == j) {
							par(usr = c(0,1,0,1))
							l.wid <- strwidth(labels, "user")
							cex.labels <- max(0.5, min(2, 0.9 / max(l.wid)))
							textPanel(0.5, 0.5 , labels[i],	cex = cex.labels, font = 1)
						}
						else if(i < j) {
							# plotting the combinations of pairs MA's
							# new notation k which means the # of plots
							points(x.pair[[k]]$A, x.pair[[k]]$M, pch = pch, cex = cex, xlab = xlab, ylab = ylab, 
								col = ifelse(x.pair[[k]]$Outlier, "red", "black"), ...)
							kk <- sort.list(x.pair[[k]]$A)
							lines(x.pair[[k]]$A[kk], x.pair[[k]]$UB[kk], col = "red", lty = 1)
							lines(x.pair[[k]]$A[kk], x.pair[[k]]$Q3[kk], col = "blue", lty = 2)
							lines(x.pair[[k]]$A[kk], x.pair[[k]]$LB[kk], col = "red", lty = 1)
							lines(x.pair[[k]]$A[kk], x.pair[[k]]$Q1[kk], col = "blue", lty = 2)
							if (legend.use) legend("topright", c("Upper & Lower Bound", "Q3 and Q1"), lty = 1:2, col= c("red", "blue"))
							k <- k + 1
						}
					}
				}
			}			
			invisible(NULL)
		}
	}
)

####################################################################
dist.two <- function(type, q1, q2, p){
    # distance between q=(q2-q1) vector and p
    # algorithm : (p'q / q'q)q
    # type
    ## Pred : new point, A_{j}
    ## Resp : A_{j}
	# additional assumption
	# A_j =
	## |y*v|/norm, 	 if y*v/v'v >= 0
	## -|y*v|/norm,  o.w.
    norm.two <- sum((q2 - q1)^2)
    #norm.two <- sum((q2 - c(0,0,0))^2)
    if(norm.two != 0){
        u <- sum((p-q1)*(q2-q1)) / norm.two
        new.pt <- q1 + u*(q2-q1)
#		if( u < 0) new.pt <- - new.pt
		A.tmp <- sqrt(sum((new.pt - q1)^2))
        if(type == "Resp") return(new.pt)
        else if(type == "Pred") return(c(M=sqrt(sum((new.pt - p)^2)),A=ifelse(u >= 0,A.tmp, -A.tmp),new.pt))
    }
    else{
        # R^{2} space
        return(sqrt(sum( q2 - p)^2))
    }
}

####################################################################
quant.const <- function(x, Lower, Upper, type){
	if(type == "pair"){
		M <- (x[,1]-x[,2])
		A <- (x[,1]+x[,2])/2

		quant.val <- (Lower + Upper) / 2
		Q2 <- quantile(M, probs= quant.val)
		M <- M-Q2
	}
	else if(type == "proj"){
	    M <- x[,1]
	    A <- x[,2]
	}

		Q3 <- quantile(M, probs= Upper)
		Q3 <- rep(Q3, length(A))

		Q1 <- quantile(M, probs= Lower)
		Q1 <- rep(Q1, length(A))

	return(list(M=M, A=A, Q1=Q1, Q3=Q3))
}

####################################################################
quant.linear <- function(x, Lower, Upper, type){
	if(type == "pair"){
		M <- (x[,1]-x[,2])
		A <- (x[,1]+x[,2])/2

		quant.val <- (Lower + Upper) / 2
		Q2 <- rq(M~A, tau = quant.val )
		Q2 <- Q2$coef[1]+A*Q2$coef[2]
		M <- M-Q2

		Q3 <- rq(abs(M)~A, tau = quant.val )
		Q3 <- Q3$coef[1]+A*Q3$coef[2]

		Q1 <- -Q3		
	}
	else if(type == "proj"){
		M <- x[,1]
		A <- x[,2]

		Q3 <- rq(abs(M)~A, tau = Upper)
		Q3 <- Q3$coef[1]+A*Q3$coef[2]

		Q1 <- rq(abs(M)~A, tau = Lower)
		Q1 <- Q1$coef[1]+A*Q1$coef[2]
	}
    return(list(M=M, A=A, Q1=Q1, Q3=Q3))
}

####################################################################

quant.nonlin <- function(x, Lower, Upper, type, nonlin.method, nonlin.SS, nonlin.Frank){

	FrankModel <- function(x,delta,mu,sigma,df,tau){
		z <- qt(-log(1-(1-exp(-delta))/(1+exp(-delta*pt(x,df))*((1/tau)-1)))/delta,df)
		mu + sigma*z
	}

	if(type == "pair"){
		M <- (x[,1]-x[,2])
		A <- (x[,1]+x[,2])/2
		quant.val <- (Lower + Upper) / 2		

		tmp <- try(nlrq(M ~ SSlogis(A, Asym, mid, scal), tau = quant.val ), silent = TRUE)

		if(class(tmp) != "try-error") {
		   Q2 <- nlrq(M ~ SSlogis(A, Asym, mid, scal), tau= quant.val)
		   Q2 <- predict(Q2, newdata=list(x=A))
		   M <- M-Q2
		   Q3 <-  switch( nonlin.SS,
				Frank = try( nlrq( abs(M) ~ FrankModel(A, delta, mu, sigma, df = nonlin.Frank[1], tau = quant.val),
					tau = quant.val, start=list(delta=nonlin.Frank[2], mu = nonlin.Frank[3], sigma = nonlin.Frank[4]), 
					method = nonlin.method ), silent = TRUE),
				Self = try(nlrq(abs(M) ~ SSlogis(A, Asym, mid, scal), tau = quant.val, method = nonlin.method ), silent = TRUE),
				Asym = try(nlrq(abs(M) ~ SSasymp( A, Asym, R0, lrc), tau = quant.val, method = nonlin.method  ), silent = TRUE),
				Orig = try(nlrq(abs(M) ~ SSasympOrig(A, Asym, lrc), tau = quant.val, method = nonlin.method  ), silent = TRUE),
				biexp = try(nlrq(abs(M) ~ SSbiexp(A, A1, lrc1, A2, lrc2), tau = quant.val, method = nonlin.method  ), silent = TRUE),
				AsymOff = try(nlrq(abs(M) ~ SSasympOff(A, Asym, lrc, c0), tau = quant.val, method = nonlin.method  ), silent = TRUE)
			)
			if(class(Q3)=="try-error") {
			   fit <- quant.linear(x, Lower, Upper, type)
			   return(list(M=fit$M, A=fit$A, Q1=fit$Q1, Q3=fit$Q3))
			}
			Q3 <- predict(Q3, newdata=list(x=A))
			Q1 <- -Q3
		}
		else {
		   fit <- quant.linear(x, Lower, Upper, type)
		   return(list(M=fit$M, A=fit$A, Q1=fit$Q1, Q3=fit$Q3))
		}
	}
	else if(type == "proj" ){
		M <- x[,1]
		A <- x[,2]
		
		Q3 <-  switch( nonlin.SS,
				Frank = try( nlrq( abs(M) ~ FrankModel(A, delta, mu, sigma, df = nonlin.Frank[1], tau = Upper),
					tau = Upper, start=list(delta=nonlin.Frank[2], mu = nonlin.Frank[3], sigma = nonlin.Frank[4]), 
					method = nonlin.method ), silent = TRUE),
				Self = try(nlrq(abs(M) ~ SSlogis(A, Asym, mid, scal), tau = Upper, method = nonlin.method ), silent = TRUE),
				Asym = try(nlrq(abs(M) ~ SSasymp( A, Asym, R0, lrc), tau = Upper, method = nonlin.method  ), silent = TRUE),
				Orig = try(nlrq(abs(M) ~ SSasympOrig(A, Asym, lrc), tau = Upper, method = nonlin.method  ), silent = TRUE),
				biexp = try(nlrq(abs(M) ~ SSbiexp(A, A1, lrc1, A2, lrc2), tau = Upper, method = nonlin.method  ), silent = TRUE),
				AsymOff = try(nlrq(abs(M) ~ SSasympOff(A, Asym, lrc, c0), tau = Upper, method = nonlin.method  ), silent = TRUE)
				)
		if(class(Q3)=="try-error") {
		   fit <- quant.linear(x, Lower, Upper, type)
		   return(list(M=fit$M, A=fit$A, Q1=fit$Q1, Q3=fit$Q3))
		}

		Q3 <- predict(Q3, newdata=list(x=A))

		Q1 <-  switch( nonlin.SS,
				Frank = try( nlrq( abs(M) ~ FrankModel(A, delta, mu, sigma, df = nonlin.Frank[1], tau = Lower),
					tau = Lower, start=list(delta=nonlin.Frank[2], mu = nonlin.Frank[3], sigma = nonlin.Frank[4]), 
					method = nonlin.method ), silent = TRUE),
				Self = try(nlrq(abs(M) ~ SSlogis(A, Asym, mid, scal), tau = Lower, method = nonlin.method ), silent = TRUE),
				Asym = try(nlrq(abs(M) ~ SSasymp( A, Asym, R0, lrc), tau = Lower, method = nonlin.method  ), silent = TRUE),
				Orig = try(nlrq(abs(M) ~ SSasympOrig(A, Asym, lrc), tau = Lower, method = nonlin.method  ), silent = TRUE),
				biexp = try(nlrq(abs(M) ~ SSbiexp(A, A1, lrc1, A2, lrc2), tau = Lower, method = nonlin.method  ), silent = TRUE),
				AsymOff = try(nlrq(abs(M) ~ SSasympOff(A, Asym, lrc, c0), tau = Lower, method = nonlin.method  ), silent = TRUE)
				)
		if(class(Q1)=="try-error") {
		   fit <- quant.linear(x, Lower, Upper, type)
		   return(list(M=fit$M, A=fit$A, Q1=fit$Q1, Q3=fit$Q3))
		}

		Q1 <- predict(Q1, newdata=list(x=A))
	}	
	
    if(length(which(abs(M) < abs(Q1))) < length(M)/4 & type == "pair") {
       fit <- quant.linear(x, Lower, Upper , type)
	   return(list(M=fit$M, A=fit$A, Q1=fit$Q1, Q3=fit$Q3))
    }
    return(list(M=M, A=A, Q1=Q1, Q3=Q3))
}


####################################################################

quant.nonpar <- function(x, Lower, Upper, type, lbda){

	if(type == "pair"){
	    M <- (x[,2]-x[,1])
		A <- (x[,2]+x[,1])/2
	}
	else if(type == "proj"){
		M <- x[,1]
		A <- x[,2]
	}
	#First Quantile
	Q1 <- rqss(M~qss(A, lambda=lbda), tau = Lower )
	i <- match(A, Q1$qss[[1]]$xyz[-1,1])
	y  <- Q1$coef[-1]
	Q1 <- Q1$coef[1] + y[i]
	#Q1 <- -Q1

	#Third Quantile
	Q3 <- rqss(M~qss(A, lambda=lbda), tau = Upper )
	i <- match(A, Q3$qss[[1]]$xyz[-1,1])
	y  <- Q3$coef[-1]
	Q3 <- Q3$coef[1] + y[i]
    return(list(M=M, A=A, Q1=Q1, Q3=Q3))
}

####################################################################

quant.const.D <- function(x, dist.mthd, Lower, Upper, n.obs, n.para){
	if(is.data.frame(x)) x <- as.data.frame(x)
    A.med <- apply(x,1,dist.mthd,na.rm=TRUE)
	M <- x - A.med
	M.long <- reshape(M, direction = "long", varying = list(1:n.obs))[,2]	

	A <- apply(x, 1, mean, na.rm = TRUE)
	
    Q1 <- quantile(M.long, probs=Lower)
    Q3 <- quantile(M.long, probs=Upper)

    Q1 <- rep(Q1, length(A))
    Q3 <- rep(Q3, length(A))

    return(list(M = M, A=A, Q1=Q1, Q3=Q3))
}

####################################################################
quant.linear.D <- function(x, dist.mthd, Lower, Upper, n.obs, n.para){
	if(is.data.frame(x)) x <- as.data.frame(x)
    A.med <- apply(x,1,dist.mthd,na.rm=TRUE)
	M <- x - A.med
	M.long <- reshape(M, direction = "long", varying = list(1:n.obs))[,2]	

	A <- apply(x, 1, mean, na.rm = TRUE)
	A <- rep(A,  n.obs)
	
	quant.val <- (Lower + Upper) / 2
    Q3 <- rq(abs(M.long) ~ A, tau= quant.val)
    Q3 <- Q3$coef[1]+A*Q3$coef[2]

	Q1 <- -Q3
    return(list(M = M, A=A[1:n.para], Q1=Q1[1:n.para], Q3=Q3[1:n.para]))
}


####################################################################

quant.nonlin.D <- function(x, dist.mthd, Lower, Upper, n.obs, n.para){
	if(is.data.frame(x)) x <- as.data.frame(x)
    A.med <- apply(x,1,dist.mthd,na.rm=TRUE)
	M <- x - A.med
	M.long <- reshape(M, direction = "long", varying = list(1:n.obs))[,2]	

	A <- apply(x, 1, mean, na.rm = TRUE)
	A <- rep(A,  n.obs)

	quant.val <- (Lower + Upper) / 2
	
    tmp <- try(nlrq(abs(M.long) ~ SSlogis(A, Asym, mid, scal), tau= quant.val ), silent = TRUE)
    if(class(tmp) =="try-error") {
           fit <- quant.linear.D(x,dist.mthd, Lower, Upper, n.obs, n.para)
           return(list(M=fit$M, A=fit$A, Q1=fit$Q1, Q3=fit$Q3))
    }

    Q3 <- nlrq(abs(M.long) ~ SSlogis(A, Asym, mid, scal), tau= quant.val)
    Q3 <- predict(Q3, newdata=list(x=A))

	Q1 <- -Q3
    return(list(M = M, A=A[1:n.para], Q1=Q1[1:n.para], Q3=Q3[1:n.para]))
}
####################################################################

quant.nonpar.D <- function(x, dist.mthd, Lower, Upper, n.obs, n.para, lbda){
	if(is.data.frame(x)) x <- as.data.frame(x)
    A.med <- apply(x,1,dist.mthd,na.rm=TRUE)
	M <- x - A.med
	M.long <- reshape(M, direction = "long", varying = list(1:n.obs))[,2]	

	A <- apply(x, 1, mean, na.rm = TRUE)
	A <- rep(A,  n.obs)
	
    #First Quantile
    Q1 <- rqss(M.long~qss(A, lambda=lbda), tau= Lower)
    i <- match(A, Q1$qss[[1]]$xyz[-1,1])
    y  <- Q1$coef[-1]
    Q1 <- Q1$coef[1] + y[i]

    #Third Quantile
    Q3 <- rqss(M.long~qss(A, lambda=lbda), tau= Upper )
    i <- match(A, Q3$qss[[1]]$xyz[-1,1])
    y  <- Q3$coef[-1]
    Q3 <- Q3$coef[1] + y[i]

    return(list(M = M, A=A[1:n.para], Q1=Q1[1:n.para], Q3=Q3[1:n.para]))
}
