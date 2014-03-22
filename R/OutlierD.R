##########################################################################
#
#        
#    Outlier Detection for Multiplicative High-throughput Data
#
#             by
#
#      Soo-Heang Eo and HyungJun Cho
#      Deparment of Statistics 
#      Korea University
#
#      Last updated: March 2014
#
##########################################################################

###################################################################################
#
# Main function  for OutlierDM 
#
###################################################################################

setClass("OutlierDM", representation(call = "language",
									 raw.data = "data.frame",
									 res = "data.frame",
									 x.pair = "list",
									 k = "numeric",
									 n.outliers = "integer",
									 method = "character",
									 type = "character",
									 contrl.para = "list"
									 )
	)

odm <- function(x, k=1.5, method= c("nonlin", "constant", "linear", "nonpar"), 
	type = c("proj", "diff", "pair", "grubbs", "dixon"), ...){
# Outlier Detection for Multiplicative High-throughput Data
# Soo-Heang Eo and HyungJun Cho
## Input parameters
# x : Data vectors or matrix
# k : tuning parameter 
# method : fitting method for quantile regression
# type : choose multiplicative detection algorithm
# ... : minor control parameters including pair.cre, dist.mthd, Lower, Upper, trans, and centering.
# pair.cre : creterion for pairwise approach
# dist.mthd : median or mean
# Lower : Lower values
# Upper : Upper values
# trans : a parameter for logarithm transformation.
# centering : Logical parameter for centering. If TRUE, data are centered by column means.
# projection.type : projection type, naive, PCA, LPC.
## output 

    ##########
    #Preparation
	call <- match.call()
	method <- match.arg(method)
	type <- match.arg(type)

	contrl.para = odm.control(...)
	pair.cre <- contrl.para[[1]]
	dist.mthd <- contrl.para[[2]]
	Lower <- contrl.para[[3]]
	Upper <- contrl.para[[4]]
	trans <- contrl.para[[5]]
	centering <- contrl.para[[6]]
	projection.type <- contrl.para[[7]]
	lbda <- contrl.para[[8]]
	nonlin.method <- contrl.para[[9]]
	nonlin.SS <- contrl.para[[10]]
	nonlin.Frank <- contrl.para[[11]]	
	ncl <- contrl.para[[12]]
	
	### FIXME: what type of parallel computing is optimal for our function? ##	
	# Set parallel computing in order to incease computing power

	# Start constructing data matrix
	if(!is.data.frame(x)) x <- as.data.frame(x) 
    n.obs <- ncol(x) 
    n.para <- nrow(x) 
    rownames(x) <- 1:nrow(x) 
	colnames(x) <- paste("N", 1:n.obs, sep = "")
	raw.data <- x
 
	if(n.para < 30) stop('the number of peptides is too few.')
	
	# transformation
	x <- eval(call(trans,x))

    # End constructing data
	
    #cat("Please wait... \n")
	
	if(type %in% c("dixon","grubbs")) {
	# Dixon's Tange Test and Grubbs Test for multiplicative experiments
	# The algorithm is based on the R package, outliers.
		if(type == "dixon")	fit <- apply( x, 1, function(x) outliers::dixon.test(x, opposite = TRUE)$p.value)
		else if(type =="grubbs") { 
			fit <- apply(x, 1, function(x) {
									outliers::grubbs.test(x,type = 10)$p.value
									#tmp[2] <- try(grubbs.test(x,type = 11)$p.value, TRUE)
									#tmp[3] <- try(grubbs.test(x,type = 20)$p.value,  TRUE)
									#return(min(ifelse(class(tmp) == "try-error", NA, tmp), na.rm = TRUE))
									}
							)
		}
		
       #Outputs
        x <- cbind(x, round(fit,4))
        colnames(x) <- c(paste("Rep",1:n.obs, sep = ""),"pvalue")
        i <-  which(x$pvalue <= 0.05)
        n.outliers <- length(i)
        Outlier <- rep(FALSE, n.para)
        if(length(i) > 0) Outlier[i] <- TRUE
        x <- cbind(Outlier, x)
		new("OutlierDM", call = call, raw.data = raw.data, res=x, n.outliers=n.outliers, type = type, contrl.para = contrl.para)
	}
    else if(type == "pair"){
	# Pairwise outlier detection for multiplicative experiments
	# The algorithm is the same as that of OutlierD package when the number of replicaates are 2.
	# Step 1. Apply OutlierD to all possible pairs of experiments.
	# Step 2. Declare peptide j as an outlier if it is declared as an outlier for at least one pair.
		nonlin.SS <- "Asym" # fixed for nonlinear modelling
		predict.num <- combn(n.obs,2)
        fit.tmp <- list()
        x.tmp <- list()
        for(i in 1:ncol(predict.num)){
            fit.tmp[[i]] <- switch(method,
                   constant = quant.const(x[,predict.num[,i]], Lower, Upper, type),
                   linear      = quant.linear(x[,predict.num[,i]], Lower, Upper, type),
                   nonlin     = quant.nonlin(x[,predict.num[,i]], Lower, Upper, type,nonlin.method, nonlin.SS, nonlin.Frank),
                   nonpar   = quant.nonpar(x[,predict.num[,i]], Lower, Upper, type, lbda)
                   )
        #Outputs
        x.tmp[[i]] <- cbind(x[,predict.num[,i]], fit.tmp[[i]]$A, fit.tmp[[i]]$M, fit.tmp[[i]]$Q1, fit.tmp[[i]]$Q3, 
                fit.tmp[[i]]$Q1-k*(fit.tmp[[i]]$Q3-fit.tmp[[i]]$Q1), fit.tmp[[i]]$Q3+k*(fit.tmp[[i]]$Q3-fit.tmp[[i]]$Q1))
		#x.tmp[[i]] <- as.data.frame(x.tmp[[i]])
        colnames(x.tmp[[i]]) <- c(paste("Rep",predict.num[,i], sep = ""),"A","M","Q1","Q3","LB","UB")
        j <-  which((x.tmp[[i]][,"M"] < x.tmp[[i]][,"LB"])|(x.tmp[[i]][,"M"] > x.tmp[[i]][,"UB"]))
        n.outliers <- length(j)
        Outlier <- rep(FALSE, n.para)
        if(length(j) > 0) Outlier[j] <- TRUE
        x.tmp[[i]] <- cbind(Outlier, x.tmp[[i]])
        }
        outlier.logi <- as.data.frame(sapply( x.tmp, base::subset , select = Outlier))
        if(ncol(outlier.logi) > 2) {
			colnames(outlier.logi) <- paste("Outlier.pair",1:n.obs,sep = "")
			outlier.logi <- cbind( outlier = ifelse( apply(outlier.logi, 1, sum) >= pair.cre, TRUE, FALSE), outlier.logi, x)
        }
		else outlier.logi <- cbind( outlier = ifelse( apply(outlier.logi, 1, sum) >= pair.cre, TRUE, FALSE), x)

        #cat(" Done. \n")

		new("OutlierDM", call = call, raw.data = raw.data, res= outlier.logi, x.pair = x.tmp, k= k, n.outliers= n.outliers, method= method, type= type, contrl.para = contrl.para)
    } # 22 Aug 2011 by Soo-Heang Eo.
    else if(type == "diff"){ 
	# Difference approach for multiplicative highthroughput data 
	# Step 1. Compute the difference M = y - median and the average A .
	# Step 2. Obtain the first and third quantile values, Q1(A) and Q3(A), on a MA plot using quantile regression approach.
	# Step 3. Calculate IQR(A) = Q3(A) - Q1(A)
	# Step 4. Construct the lower and upper fences, LB(A) = Q1(A) - kIQR(A) and UB(A) = Q3(A) + kIQR(A), where k is a tuning parameter.
	# Step 5. Declare the i-th replicate of peptide j as an outlier if it locates over the upper fence or under the lower fence.
        fit <- switch(method,
               constant = quant.const.D(x, dist.mthd, Lower, Upper, n.obs, n.para),
               linear      = quant.linear.D(x, dist.mthd, Lower, Upper, n.obs, n.para),
               nonlin     = quant.nonlin.D(x, dist.mthd, Lower, Upper, n.obs, n.para),
               nonpar   = quant.nonpar.D(x, dist.mthd, Lower, Upper, n.obs, n.para, lbda)
               )
        #Outputs
        new.x <- cbind(x, fit$M, fit$A, fit$Q1, fit$Q3, 
                fit$Q1-k*(fit$Q3-fit$Q1), fit$Q3+k*(fit$Q3-fit$Q1))
        colnames(new.x) <- c(paste("Rep",1:n.obs, sep = ""),paste("M",1:n.obs, sep = ""),"A","Q1","Q3","LB","UB")
        i <- (new.x[,(n.obs+1):(2*n.obs)] < new.x$LB)|(new.x[,(n.obs+1):(2*n.obs)] > new.x$UB)
        i <- apply(i,1, any)
        n.outliers <- sum(i,na.rm = TRUE)
        Outlier <- rep(FALSE, n.para)
        if(sum(i,na.rm=TRUE) > 0) Outlier[i] <- TRUE
        new.x <- cbind(Outlier, new.x)
        cat("Done. \n")
		new("OutlierDM", call = call, raw.data = raw.data, res= new.x, k= k, n.outliers= n.outliers, method= method, type= type, contrl.para = contrl.para)
    }
    else{
        if( type == "proj"){ 
		# Outlier Detection using projections for multiplicative experiments
		# Step 1. Shift the sample means to the origin.(in preparation step)
		# Step 2. Find the first PC vector v using PCA on the space of y.(projection.type option)
		# Step 3. Obtain the projection of a vector of each peptide j on v.
		# Step 4. Compute the length of the projection A and the length of the difference between a vector of peptide j and the projection M.
		# Step 5. Obtain the third quantile value Q3(A), on a MA plot using a quantile regression approach, and assume Q1(A) = -Q3(A).
		# Step 6. Calculate IQR.
		# Step 7. Construct the lower and upper fences.
		# Step 8. Declare peptide j as an outlier if it locates over the upper fence or under the lower fence.

		#Centering
		if(centering){
			xMeans <- colMeans(x)
			x <- x - rep(xMeans, each = n.para)
		}
		
		if(projection.type == "naive"){
			pt.naive <- 1000 * rep(1,n.obs)
			pt.tmp <- t(apply( x, 1, dist.two, type = "Pred", q1 = rep(0,n.obs), q2 = pt.naive))
		}
		else if(projection.type %in% c("rPCA","PCA")){
			if(projection.type == "rPCA"){ # Project-Pursuit PCA
				pca1 <- PCAproj(x,method = "sd",CalcMethod = "lincomb")$loadings
				wt.line <- abs(pca1[,1])
			}
			else if(projection.type == "PCA"){ # Classical PCA
				pca1 <- prcomp(x, scale=TRUE)
				wt.line <- abs(pca1$rotation[,1])
			}
			pt.pca <- 1000 * wt.line
			pt.tmp <- t(apply( x, 1,dist.two, type = "Pred", q1 = rep(0, n.obs), q2 = pt.pca))
		}
		else{
			stop("TYPE input does not match in a function.")
		}        
		new.x <- cbind(M = pt.tmp[,1], A = pt.tmp[,2])
		fit <- switch( method,
			   constant = quant.const(new.x, Lower, Upper, type),
			   linear      = quant.linear(new.x,Lower, Upper, type),
			   nonlin     = quant.nonlin(new.x,Lower, Upper, type, nonlin.method, nonlin.SS, nonlin.Frank),
			   nonpar   = quant.nonpar(new.x,Lower, Upper, type, lbda)
			   )
       }
       #Outputs
        x <- cbind(x, fit$A, fit$M, fit$Q1, fit$Q3, 
                fit$Q1-k*(fit$Q3-fit$Q1), fit$Q3+k*(fit$Q3-fit$Q1))
        colnames(x) <- c(paste("Rep",1:n.obs, sep = ""),"A","M","Q1","Q3","LB","UB")
        
        i <-  which(x$M >= x$UB)
        n.outliers <- length(i)
        Outlier <- rep(FALSE, n.para)
        if(length(i) >0) Outlier[i] <- TRUE
        x <- cbind(Outlier, x)
        cat(" Done. \n")
 		new("OutlierDM", call = call, raw.data = raw.data, res=x, k=k, n.outliers=n.outliers, method=method, type= type, contrl.para = contrl.para)
    }
}
#END################################################################s
