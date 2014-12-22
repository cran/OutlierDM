\name{OutlierDM-package}
\alias{OutlierDM}
\alias{OutlierDM-package}
\docType{package}
\title{
	Functions for detecting outlying parameters (peptides) or observations (samples) in multi-replicated high-throughtput data such as mass spectrometry experiments
}
\description{
  This package provides several outlier detection algorithms for multi-replicated high-throughput data ranged from classical approaches to boxplot approaches based on a MA plot.
}
\details{
\tabular{ll}{
Package: \tab OutlierDM\cr
Type: \tab Package\cr
Version: \tab 1.1.1\cr
Date: \tab 2014-12-23\cr
License: \tab GPL version 3\cr
LazyLoad: \tab no\cr
}
}
\author{
Soo-Heang Eo <eo.sooheang@gmail.com> and HyungJun Cho <hj4cho@korea.ac.kr>

Maintainer: Soo-Heang Eo <eo.sooheang@gmail.com>
}
\references{
Eo, S-H and Cho, H (2015) OutlierDM: More robust outlier detection algorithms for multi-replicated high-throughput data.

Cho, H and Eo, S-H. (2015) Outlier detection for mass-spectrometry data.

Eo, S-H, Pak D, Choi J, Cho H (2012) Outlier detection using projection quantile regression for mass spectrometry data with low replication. BMC Res Notes. 

Cho H, Lee JW, Kim Y-J, et al. (2008) OutlierD: an R package for outlier detection using quantile regression on mass spectrometry data. Bioinformatics 24:882--884. 

Min H-K, Hyung S-W, Shin J-W, et al. (2007). Ultrahigh-pressure dual online solid phase extraction/capillary reverse-phase liquid chromatography/tandem mass spectrometry (DO-SPE / cRPLC / MS / MS): A versatile separation platform for high-throughput and highly sensitive proteomic analyses. Electrophoresis 28:1012--1021. 

Grubbs FE (1969) Procedures for detecting outlying observations in samples. Technometrics 11:1--21. 

Dixon WJ (1951) Ratios involving extreme values. Ann Math Statistics 22:68--78. 

Dixon WJ (1950) Analysis of extreme values. Ann Math Statistics 21:488--506. 

Grubbs FE (1950) Sample criteria for testing outlying observations. Ann Math Statistics 21:27--58. 
}
\seealso{
    \code{\link{odm}}, \code{\link{odm.control}}, \code{\link[quantreg:rq]{quantreg}}
}
\keyword{package}