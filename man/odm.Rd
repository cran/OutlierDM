\name{odm}
\alias{odm}
\title{Outlier dectection for (multiple) replicated high-throughput data}
\description{
	This function is desinged to detect outlying observations, or peptides, based on a MA plot with boxplot by quantile regression methods.
}
\usage{
     odm(x, k=1.5, method= c("nonlin", "constant", "linear", "nonpar"), 
  type = c("proj", "diff", "pair", "grubbs", "dixon"), ...)
}

\arguments{
   \item{x}{data vectors or matrices. These can be given as named arguments. If the number of predictors is 2, x1 describes one n-by-1 vector for data and x2 describes the other n-by-1 vector for data (n= number of peptides, proteins, or genes)}
   \item{k}{a tuning parameter used in Q1-k*IQR and Q3+k*IQR, IQR=Q3-Q1. Default value is 1.5 }
   \item{method}{Selet a type of quantile regression methods used in an outlier detection algorithm. Use one of "constant", "linear", "nonlin", and "nonpar"m which mean constant, linear, non-linear, and nonparametric quantile regression in order. For more details, see the \code{quantreg} package.}
   \item{type}{input parameter for multiplicative detection algorithm. Use one of "dixon", "grubbs", "pair", "diff", and "proj". \cr
			dixon: Dixon's Q-test \cr
			grubbs: Grubbs' test \cr
			pair: OutlierD algorithm (Cho et al., 2008) \cr
			diff: Difference approach for multiplicative highthroughput data \cr
			proj: Outlier Detection using projections for multiplicative experiments (Eo et al., 2012) \cr	
			}
   \item{...}{tuning parameters used in odm.control(). See \code{\link{odm.control}}. 
   }
}
\details{
	Outlier detection is an important issue in high-throughput data generated from various biological or chemical experiments such as mass spectrometry and microarrays because elaborate pre-processing is essential for successful analysis and outlier detection is one of pre-processing steps. The heterogeneity of variability and low replication are often obstables for successful analysis including outlier detection. An outlier detection algorithm using quantile regression was developed to solve the problems; however, it can be used only for duplicate experiments. We here proposed an outlier detection algorithm using projection for multiplicative experiments. 
}
\value{
   \item{\code{call}:}{evaluated function call}
    \item{\code{raw.data}:}{data to be used in the fitted model}
    \item{\code{res}:}{a data.frame including the information about the fitted model. It consists of several columns including \code{outlier, M, A, Q3, Q1, UB} and \code{LB}.}
    \item{\code{x.pair}:}{Object of class \code{"list"} ~~ }
    \item{\code{k}:}{a scalar parameter for constructing boxplot used in the fitted models}
    \item{\code{n.outliers}:}{a scalar value that denotes the number of outliers to be detected by the fitted model.}
    \item{\code{method}:}{the type of method used in the fitted model}
    \item{\code{type}:}{the type of algorithm used in the fitted model}
    \item{\code{contrl.para}:}{a list including information about tuning parameters}
}
\seealso{
    \code{\link{OutlierDM-package}} to show general information about the package\cr
	\code{\link{OutlierDM-class}} to provide information about the "OutlierDM" class \cr
    \code{\link{odm.control}} to control tuning parameters \cr
}
\references{
S-H Eo, D Pak, J Choi, H Cho (2012). Outlier Detection for Multiplicative High-throughput Data. \emph{BMC Research Notes}, \bold{5}, 1--6.
}
\examples{
\dontrun{
data(lcms3)

## Fit projection approaches
fit.proj.const <- odm(lcms3, method="constant")
fit.proj.linear <- odm(lcms3, method="linear")
fit.proj.nonlin <- odm(lcms3, method="nonlin")
fit.proj.nonpara <- odm(lcms3, method="nonpar", lbda = 1)

par(mfrow = c(2,2))
plot(fit.proj.const, main = "Constant")
plot(fit.proj.linear, main = "Linear")
plot(fit.proj.nonlin, main = "NonLinear")
plot(fit.proj.nonpara, main = "Nonparametric")

## Fit pairwise OutlierD algorithm
fit0 <- odm(lcms3, type = "pair")
plot(fit0)
}
}
\keyword{odm}


