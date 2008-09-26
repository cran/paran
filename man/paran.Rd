\name{paran}
\alias{paran}
\title{Horn's Parallel Analysis of Principal Components/Factors}
\description{
\code{paran} performs Horn's parallel analysis for a principal component or common factor analysis, so as to adjust for finite sample bias in the retention of components.
}
\usage{
paran(x, iterations=0, centile=0, quietly=FALSE, 
      status=TRUE, all=FALSE, cfa=FALSE, graph=FALSE, 
      color=TRUE, col=c("black","red","blue"), 
      lty=c(1,2,3), lwd=1, legend=TRUE, file="", 
      width=640, height=640, grdevice="png", seed=0)
}
\arguments{
 \item{x}{a numeric matrix or data frame for principal component analysis, or common factor analysis.}
 \item{iterations}{sets the number of iterations with a user specified whole number representing the number of random data sets to be produced in the analysis. The default, indicated by \code{iterations=0}, is 30\enc{*}P, where P is the number of variables or columns in \code{x}.}
 \item{centile}{employs Monte Carlo estimates according to the user specified whole number between 1 and 99 indicating the centile used in estimating bias. The default is to use the mean. By selecting a conservative number, such as 95 or 99, and a large number of iterations, \code{paran} can be used to perform the modified version of parallel analysis suggested by Glorfeld (1995).}
 \item{quietly}{suppresses tabled output of the analysis, and only returns the vector of estimated biases.}
 \item{status}{indicates progress in the computation. Parallel analysis can take some time to complete given a large data set and/or a large number of iterations. The \code{cfa} option may noticably increase the computational requirements of \code{paran}.}
 \item{all}{report all eigenvalues (default reports only those components or factors that are retained).}
 \item{cfa}{performs a common factor analysis instead of a principal component analysis. This provides only the unrotated eigenvalues from the common factor model. As of version 1.4.0 \code{paran} performs parallel analysis for common factor analysis using a modified method. See Remarks for details.}
 \item{graph}{requests that a plot of the unadjusted, adjusted, and random eigenvalues in a format similar to that presented by Horn in his 1965 paper. Retained components or factors are indicated by the solid circular markers on the adjusted eigenvaue plot, and non-retained components or factors are indicated with hollow circular markers.}
 \item{color}{renders the graph in color with unadjusted eigenvalues drawn in red, adjusted eigenvalues drawn in black, and random eigenvalues drawn in blue if \code{color=TRUE}, and all lines drawn solid. If \code{color=FALSE}, the graph is rendered in black and white, and the line connecting the unadjusted eigenvalues is dashed, the line connecting the random eigenvalues is dotted, and the line connecting the adjusted eigenvalues is solid.}
 \item{col}{sets the colors using a character vector with the color names of adjusted eigenvalues, unadjusted eigenvalues, and estimated random eigenvalues on the plot. These settings are used, only if \code{color=TRUE}.}
 \item{lty}{sets the line type using an integer vector of the line type codes for adjusted eigenvalues, unadjusted eigenvalues, and estimated random eigenvalues on the plot. These settings are used only if \code{color=FALSE}.}
 \item{lwd}{sets the line width. The default is \code{lwd=1}.}
 \item{legend}{draws a legend in the upper right corner of the plot. The default is \code{legend=TRUE}.}
 \item{file}{the png file in which to save the graph output if the analysis is graphed and \code{file} is given a character string representing a valid path. The default is not to save the graph.}
 \item{width}{the width in pixels of the png file. The default is \code{width=640}.}
 \item{height}{the height in pixels of the png file. The default is \code{height=640}.}
 \item{grdevice}{specifies which graphic \code{\link{device}} to format the graph as, if the user has used the \code{file} option. The default is \code{grdevice=png}.}
 \item{seed}{specifies that the random number is to be seeded with the supplied integer. Each random data set is seeded with the supplied value times the number of the iteration, so that entire parallel analyses may be exactly replicated, while each simulated data set maintains a pseudorandom distinction from each of the others. The default value of \code{seed=0} tells \code{paran} to use R's default timer based seed (see \link{RNG}).}
}
\details{
\code{paran} is an implementation of Horn's (1965) technique for evaluating the components or factors retained in a principle component analysis (PCA) or common factor analysis (FA). According to Horn, a common interpretation of non-correlated data is that they are perfectly non-colinear, and one would expect therefore to see eigenvalues equal to 1 in a PCA (or 0 in an FA) of such random data. However, Horn notes that multi-colinearity occurs due to sampling error and least-squares "bias," even in uncorrelated data, and therefore actual PCAs of random data will reveal eigenvalues of components greater than and less than 1, and FAs will reveal common factors greater than and less than 0. Horn's strategy is to contrast eigenvalues produced through a parallel PCA or FA on a number of random data sets (i.e. uncorrelated variables) with the same number of variables and observations as the experimental or observational dataset to produce eigenvalues for components or factors that are adjusted for the sample error-induced inflation. For PCA, values greater than 1 are retained in the adjustment given by:

\deqn{Ev_{obs, n} - (Ev_{sim, n} - 1)}{Observed data Eigenvalue_n - (Simulated Data Eigenvalue_n - 1)}

For FA, values greater than 0 are retained in the adjustment given by:

\deqn{Ev_{obs, n} - Ev_{sim, n}}{Observed data Eigenvalue_n - Simulated Data Eigenvalue_n}

\code{paran} performs a PCA or FA with no rotation and performs Horn's adjustment. The user may also specify how many times to make the contrast with a random dataset (default is 30 per variable). Values less than 1 will be ignored, and the default value assumed. Random datasets are generated using the \code{rnorm()} function. The program returns a vector of length P of the estimated bias for each eigenvector, where P = the number of variables in the analysis. By specifying a high \code{centile} users may employ \code{paran} to conduct parallel analysis following Glorfeld's suggestions to reduce the likelihood of over-retention. (Glorfeld, 1995)
}
\section{Remarks}{
Hayton, et al. (2004) urge a parameterization of the random data to approximate the distribution of the observed data with respect to the middle ("mid-point") and the observed min and max.  However, the PCA as I understand it is insensitive to standardizing transformations of each variable, and any linear transformation of all variables, and produces the same eigenvalues used in component or factor retention decisions. This is born by the notable lack of difference between analyses conducted using the a variety of simulated distributional assumptions (Dinno, 2008). The central limit theorem would seem to make the selection of a distributional form for the random data moot with any sizeable number of iterations. Former functionality implementing the recommendation by Hayton et al. (2004) has been removed, since parallel analysis is insensitive to it, and it only adds to the computation time required to conduct parallel analysis.

As of \code{paran} version 1.4.0 application of parallel analysis to common factor analysis has been revised. See the accompanying document \emph{Gently Clarifying the Application of Horn's Parallel Analysis to Principal Components Analysis Versus Factor Analysis} downloadable at \code{http://doyenne.com/Software/files/PA_for_PCA_vs_FA.pdf}.
}
\value{
a list of objects relating to the parallel analysis:
 \item{Retained components/factors}{a scalar integer representing the number of components/factors retained}
 \item{Adjusted eigenvalues}{a vector of the estimated eigenvalues adjusted for a finite sample size}
 \item{Unadjusted eigenvalues}{a vector of the eigenvalues of the observed data from either an unrotated principal component analysis or an unrotated common factor analysis}
 \item{Random eigenvalues}{a vector of the estimated (mean or centile) eigenvalues from \code{iterations} number of \code{N} by \code{P} random data sets}
 \item{Bias}{a vector of the estimated bias of the unadjusted eigenvalues (i.e. the difference between the adjusted and unadjusted eigenvalues)}
 \item{Simulated eigenvalues}{an \code{iterations} by \code{P} matrix with each row containing the eigenvalues from an equivalent principal component or common factor analysis on an \code{N} by \code{P} data set of uncorrelated random data}
}
\section{Acknowledgement}{
A big thank you to Ulrich Keller of the University of Luxembourg for his thoughtful suggestions improving the intrerface for \code{paran}, especially the legend, and the \code{invisibly()} method for returning data, and \code{seed} option.
}
\references{

Dinno A. 2008. "Exploring the Sensitivity of Horn's Parallel Analysis to the Distributional Form of Simulated Data" \emph{Multivariate Behavioral Research}. Submitted. (Available upon request)

Glorfeld, L. W. 1995. "An Improvement on Horn's Parallel Analysis Methodology for Selecting the Correct Number of Factors to Retain". \emph{Educational and Psychological Measurement}. 55(3): 377--393

Hayton J. C., Allen D. G., and Scarpello V. 2004. "Factor Retention Decisions in Exploratory Factor Analysis: A Tutorial on Parallel Analysis" \emph{Organizational Research Methods}. 7(2): 191--205

Horn J. L. 1965. "A rationale and a test for the number of factors in factor analysis." \emph{Psychometrika}. 30: 179--185

Zwick W. R., Velicer WF. 1986. "Comparison of Five Rules for Determining the Number of Components to Retain." \emph{Psychological Bulletin}. 99: 432--442
}
\author{
Alexis Dinno (adinno at post dot harvard dot edu)
}
\examples{
## perform a standard parallel analysis on the US Arrest data
paran(USArrests, iterations=5000)

## a conservative analysis with different result!
paran(USArrests, iterations=5000, centile=95)
}
\keyword{multivariate}
