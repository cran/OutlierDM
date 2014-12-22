\name{odm}
\alias{odm}
\title{Outlier Dectection for Multi-replicated data}
\description{
	This function provides some routines for detecting outlying observations (peptides) for multi-replicated high-throughput data, especially in LC/MS experiments.
}
\usage{
odm(x, k = 3, 
  quantreg = c("linear", "nonlin", "constant", "nonpar"), 
  method = c("proj", "diff", "pair", "grubbs", "dixon", 
  "iqr",  "siqr", "Zscore"), ...)
}

\arguments{
  \item{x}{
    data vectors or matrices. These can be given as named arguments. If the number of predictors is 2, x1 describes one n-by-1 vector for data and x2 describes the other n-by-1 vector for data (n= number of peptides, proteins, or genes)
  }
  \item{k}{
    non-negative tuning parameter for the outlier detection algorithm. For IQR-based algorithms such as \code{'iqr'}, \code{'siqr'}, \code{'proj'}, \code{'diff'}, and \code{'pair'}, it works in the formula of Q1-k*IQR and Q3+k*IQR, where IQR=Q3-Q1. For \code{'Zscore'}, it works for the \code{'k'} in |Z| > k. A default value is 3.
  }
  \item{quantreg}{type of quantile regression models used for the outlier detection method. You can use one of the \code{'constant'}, \code{'linear'}, \code{'nonlin'}, and \code{'nonpar'} which mean the constant, linear, non-linear, and non-parametric quantile regression in order. For more details, see the \code{quantreg} package.}
  \item{method}{type of outlier detection methods. You can select one of the \code{'Zscore'}, \code{'iqr'}, \code{'dixon'}, \code{'grubbs'}, \code{'pair'}, \code{'diff'}, and \code{'proj'} algorithms as follows. \cr
    Zscore: Z-score based criterion (Cho and Eo, 2015) \cr
    iqr: Interquartile range (IQR) criterion (Cho and Eo, 2015) \cr
    siqr: Semi-interquartile range (IQR) criterion (Cho and Eo, 2015) \cr
  	dixon: Dixon's test (Dixon, 1950; 1951) \cr
  	grubbs: Grubbs test (Grubbs, 1950; 1969) \cr
  	pair: Pariwise OutlierD algorithm (Cho et al., 2008; Eo et al., 2012) \cr
  	proj: Projection-based OutlierD algortihm (Eo et al., 2012) \cr	
  	diff: Difference-based OutlierD algorithm (Eo and Cho, 2015)\cr
    }
  \item{...}{
    minor tuning parameters used in odm.control(). See \code{\link{odm.control}}.
  }
}
\value{
   \item{\code{call}:}{evaluated function call}
    \item{\code{raw.data}:}{raw dataset used in the model fitting}
    \item{\code{res}:}{result matrix of the model fitting. It consists of used data set with some transformation and outlying statistic.}
    \item{\code{x.pair}:}{Object of class \code{"list"} }
    \item{\code{k}:}{threshold parameter for constructing outlier detection methods}
    \item{\code{outlier}:}{matrix including the status of each outlying peptide and sample}
    \item{\code{n.outliers}:}{the number of outlying parameters (peptides) to be detected by the model fitting.}
    \item{\code{quantreg}:}{type of quantile regression used for the model fitting}
    \item{\code{method}:}{type of outlier detection method used for the model fitting}
    \item{\code{contrl.para}:}{a list of minor parameters}
}
\seealso{
  \code{\link{OutlierDM-package}} to provide the general information about the OutlierDC package\cr
	\code{\link{OutlierDM-class}} to provide the information about the \code{"OutlierDM"} class \cr
  \code{\link{odm.control}} to control tuning parameters \cr
}
\references{
Eo, S-H and Cho, H (2015) OutlierDM: More robust outlier detection algorithms for multi-replicated high-throughput data.

Cho, H and Eo, S-H. (2015) Outlier detection for mass-spectrometry data.

Eo, S-H, Pak D, Choi J, Cho H (2012) Outlier detection using projection quantile regression for mass spectrometry data with low replication. BMC Res Notes. 

Cho H, Lee JW, Kim Y-J, et al. (2008) OutlierD: an R package for outlier detection using quantile regression on mass spectrometry data. Bioinformatics 24:882--884. 

Grubbs FE (1969) Procedures for detecting outlying observations in samples. Technometrics 11:1--21. 

Dixon WJ (1951) Ratios involving extreme values. Ann Math Statistics 22:68--78. 

Dixon WJ (1950) Analysis of extreme values. Ann Math Statistics 21:488--506. 

Grubbs FE (1950) Sample criteria for testing outlying observations. Ann Math Statistics 21:27--58. 

}
\examples{
  \dontrun{
    ##############################################################
    #
    #   Outlier Detection for Mass Spectrometry Data
    #   Section 3. Illustration
    #   by HyungJun Cho and Soo-Heang Eo,
    #   Dept of Statistics, Korea University, Seoul, Korea
    #
    ##############################################################

    #####
    # Load a package OutlierDM

    # If an OutlierDM package is not installed on your system, type
    #install.package('OutlierDM', dependency = TRUE)
    library(OutlierDM)

    #####
    # Sec 3.1 When the number of replicates is large enough
    ## Load toy dataset
    data(toy)
    head(toy)
    pairs(log2(toy), pch = 20, cex = .7)

    #####
    # Fit 1. Z-score based criterion
    fit1 = odm(x = toy, method = "Zscore", k = 3)
    fit1
    summary(fit1)
    head(input(fit1))
    head(output(fit1))
    print(outliers(fit1), digits = 3)
    plot(fit1)
    rect(1, -4, 10, 4, col = heat.colors(20,alpha = 0.3), border = heat.colors(20,alpha = 0.5))
    oneplot(object = fit1, i = 4)
    title("Outlier Detection by the Z-score criterion")

    # Add a peptide name on a dot-plot
    #oneplot(fit1, 191,1)
    #title("Outlier Detection by the Z-score criterion")

    #####
    # Fit 2. Grubbs test criteria
    fit2 = odm(x = toy, method ="grubbs", alpha = 0.01)
    fit2
    summary(fit2)
    head(output(fit2))
    print(outliers(fit2), digits = 3)
    oneplot(object = fit2, i = 1)
    title("Outlier Detection by the Grubbs criterion")

    # Add text
    #oneplot(fit2, 191,1)
    #title("Outlier Detection by the Grubbs criterion")

    #####
    # Fit 3. IQR criteria
    fit3 = odm(x = toy, method = "iqr", k = 3)
    fit3
    summary(fit3)
    print(outliers(fit3), digits = 3)
    plot(fit3)
    rect(1, -4, 10, 40, col = heat.colors(20,alpha = 0.3), border = heat.colors(20,alpha = 0.5))
    oneplot(fit3, 1)
    title("Outlier Detection by the IQR criterion")

    # Add a peptide name on a dot-plot
    #oneplot(fit3, 1, 1)
    #title("Outlier Detection by the IQR criterion")

    #####
    # Fit 4. SIQR criteria
    fit4 = odm(x = toy, method = "siqr", k = 3)
    fit4
    summary(fit4)
    print(outliers(fit4), digits = 3)
    plot(fit4)
    rect(1, -4, 10, 4, col = heat.colors(20,alpha = 0.3), border = heat.colors(20,alpha = 0.5))
    oneplot(fit4, 1)
    title("Outlier Detection by the SIQR criterion")

    #####################
    ## Real data example
    #####################
    data(lcms3)
    head(lcms3)
    pairs(log2(lcms3), pch = 20, cex = .7)

    #####
    # Fit 5. OutlierD
    fit5 = odm(lcms3[,1:2], method = "pair", k = 3)
    fit5
    summary(fit5)
    head(output(fit5))
    print(outliers(fit5), digits = 3)
    plot(fit5)
    title("Outlier Detection by the OutlierD algorithm")

    #####
    # Fit 6. OutlierDM 
    fit6 = odm(lcms3, method = "proj", k = 3, center = TRUE)
    fit6
    summary(fit6)
    print(outliers(fit6), digits = 3)
    plot(fit6)
    title("Outlier Detection by the OutlierDM algorithm")
    oneplot(fit6, 18)
    #oneplot(fit6, 18, 1)
    title("The dotplot for the 18th samples of the lcms3 data")
    ### End of the illustration

    
    #####
    # Other OutlierDM algorithms
    data(lcms3)
    
    ## Load 
    ## Fit projection approaches
    fit.proj.const <- odm(lcms3, quantreg="constant")
    fit.proj.linear <- odm(lcms3, quantreg="linear")
    fit.proj.nonlin <- odm(lcms3, quantreg="nonlin")
    fit.proj.nonpara <- odm(lcms3, quantreg="nonpar", lbda = 1)
    
    par(mfrow = c(2,2))
    plot(fit.proj.const, main = "Constant")
    plot(fit.proj.linear, main = "Linear")
    plot(fit.proj.nonlin, main = "NonLinear")
    plot(fit.proj.nonpara, main = "Nonparametric")

  }
}
\keyword{odm}


