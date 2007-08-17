\name{paran}
\alias{paran}
\title{Horn's Parallel Analysis of Principal Components}
\description{
\code{paran} performs Horn's 'parallel analysis' to a principal components analysis, to adjust for sample bias in the retention of components.
}
\usage{
paran(x, iterations=0, centile=0, has=FALSE, hasmean=FALSE, quietly=FALSE, status=TRUE, all=FALSE)
}
\arguments{
 \item{x}{a numeric matrix or data frame for principal components analysis}
 \item{iterations}{a whole number representing the number of random data sets to be produced in the analysis. The default, indicated by zero, is 30*P, where P is the number of variables or columns in \code{x}.}
 \item{centile}{a whole number between 1 and 99 indicating the centile used in estimating bias. The default is to use the mean. By selecting a conservative number, such as 95 or 99, and a large number of iterations, paran can be used to perform the modified version of parallel analysis suggested by Glorfeld (1995).}
 \item{has}{this option implements a random data set generator in line with the recommendations by Hayton, Allen and Scarpello (2004), namely, that the random data are normally distributed, with the same midpoint (and, nominally, the same variance) as the simulated variable, with no simulated observations having values outside the observed maximum and minimum for that variable. Selecting this option will add considerably to the computation time of \code{paran}. (see Remarks below)}
 \item{hasmean}{as per \code{has}, except that the data are simulated using the observed mean, rather than the mid-point. (see Remarks below)}
 \item{quietly}{suppresses tabled output of the analysis, and only returns the vector of estimated biases.}
 \item{status}{indicates progress in the computation. Parallel analysis can take some time to complete given a large data set and/or a large number of iterations.}
 \item{all}{outputs the results of the parallel analysis to the table for all components, not only those with unadjusted eigenvalues greater than 1.}
}
\details{
\code{paran} is an implementation of Horn's (1965) technique for evaluating the components retained in a principle components analysis (PCA). According to Horn, a common interpretation of non-correlated data is that they are perfectly non-colinear, and one would expect therefore to see eigenvalues equal to 1 in a PCA of such data. However, Horn notes that multi-colinearity occurs due to sampling error and least-squares "bias," even in uncorrelated data, and therefore actual PCAs of such data will reveal eigenvalues of components greater than and less than 1. His strategy is to contrast eigenvalues produced through a PCA on a number of random data sets (of uncorrelated variables) with the same number of variables and observations as the experimental or observational dataset to produce eigenvalues for components that are adjusted for the sample error-induced inflation. Values greater than zero are retained in the adjustment given by:

\deqn{Ev_{obs, n} - (Ev_{sim, n} - 1)}{Observed data Eigenvalue_n - (Simulated Data Eigenvalue_n - 1)}

\code{paran} is used in place of a \code{princomp(x)} command. The user may also specify how many times to make the contrast with a random dataset (default is 30 per variable). Values less than 1 will be ignored, and the default value assumed. Random datasets are generated using the \code{rnorm()} function. The program returns a vector of length P of the estimated bias for each eigenvector, where P = the number of variables in the analysis. if centile is specified, \code{paran} may be thus be used to conduct parallel analysis following Glorfeld's suggestions to reduce the likelihood of over-retention. (Glorfeld, 1995)
}
\section{Remarks}{
Caveat regarding the \code{has} and \code{hasmean} options: Hayton, et al. urge a parameterization of the random data to approximate the distribution of the observed data with respect to the middle ("mid-point") and observed min and max.  However, the PCA as I understand it is an analysis involving a standardization of all variables. So any linear transformation of one, some or all analyzed variables produces the same eigenvalues used in component or factor retention decisions. This is born by the notable lack of difference between analyses conducted using the \code{has} or \code{hasmean} option and those conducted without it.  Certainly we would expect that residuals will be normally distributed whether the random data are generated using a uniform or normal distribution), and the central limit theorem would seem to make the selection of a distributional form for the random data moot with any sizeable number of iterations.
}
\value{
a vector of estimated sample biases in the eigenvalues of a data set analyzed by PCA.
}
\references{
Horn JL. 1965. "A rationale and a test for the number of factors in factor analysis." \emph{Psychometrika}. 30: 179-185

Zwick WR, Velicer WF. 1986. "Comparison of Five Rules for Determining the Number of Components to Retain." \emph{Psychological Bulletin}. 99: 432-442

Glorfeld, LW. 1995. "An Improvement on Horn's Parallel Analysis Methodology for Selecting the Correct Number of Factors to Retain. \emph{Educational and Psychological Measurement}. 55(3): 377-393

Hayton JC, Allen DG, and Scarpello V. 2004. "Factor Retention Decisions in Exploratory Factor Analysis: A Tutorial on Parallel Analysis" \emph{Organizational Research Methods}. 7(2): 191-205
}
\author{
Alexis Dinno (adinno at post dot harvard dot edu)
}
\seealso{
\code{\link{princomp}}.
}
\examples{
## perform a standard parallel analysis on the US Arrest data
paran(USArrests, iterations=5000)

## perform a conservative parallel analysis on the US Arrest data note the different result!
paran(USArrests, iterations=5000, centile=95)
}
\keyword{multivariate}
