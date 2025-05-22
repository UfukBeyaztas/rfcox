\name{rfcox}
\alias{rfcox}
\title{Robust Functional Cox Regression Model with Scalar and Functional Covariates}
\description{
Performs robust estimation of the functional Cox regression model (RFLCRM) by combining projection-pursuit-based robust functional principal component analysis (RFPCA) with an M-type robust partial likelihood estimator. This approach is designed to handle outliers in both functional predictors and survival responses.
}
\usage{
rfcox(time, status, X, Z, trunc = 0.9, nb = NULL, gp = NULL,
      f.weight = c("linear", "quadratic", "exponential"))
}
\arguments{
  \item{time}{A numeric vector of observed event or censoring times.}
  \item{status}{A numeric vector of event indicators: 1 for event, 0 for censored.}
  \item{X}{A matrix of functional predictors (dimension \code{n} x \code{j}), observed over a common grid.}
  \item{Z}{A matrix of scalar covariates (dimension \code{n} x \code{p}).}
  \item{trunc}{Truncation parameter for the robust M-estimator used in the Cox regression step.}
  \item{nb}{Number of basis functions to be used in RFPCA.}
  \item{gp}{Grid points over the functional domain. If \code{NULL}, defaults to \code{seq(0, 1, length.out = j)}.}
  \item{f.weight}{A character string specifying the weight function used in the robust Cox model. One of \code{"linear"}, \code{"quadratic"}, or \code{"exponential"}.}
}
\details{
This function implements the robust functional linear Cox regression model (RFLCRM), which integrates two stages:

\enumerate{
  \item \emph{Robust Functional Principal Component Analysis (RFPCA)}: The functional predictors are decomposed using a projection-pursuit approach that identifies principal components by maximizing robust dispersion measures.
  \item \emph{Robust Partial Likelihood Estimation}: The Cox model is then fitted using an M-type estimator that incorporates a robust weighting function to downweight extreme observations in both covariates and survival times.
}
}
\value{
A list with components:
  \item{bhat}{Estimated functional coefficient evaluated over the grid points.}
  \item{gammahat}{Estimated scalar coefficients.}
  \item{concordance}{Concordance index for model predictive performance.}
  \item{model}{Fitted robust Cox model object.}
  \item{rfpca}{List of RFPCA outputs: basis evaluations, scores, and eigenfunctions.}
}
\references{
Bali, J. L., Boente, G., Tyler, D. E., and Wang, J.-L. (2011). Robust functional principal components: A projection-pursuit approach. \emph{Annals of Statistics}, \bold{39}(6), 2852-2882. \cr
Bednarski, T. (1993). Robust estimation in Cox's proportional hazards model. \emph{Scandinavian Journal of Statistics}, \bold{20}(2), 189-204.
}
\author{Gizel Bakicierler Sezer and Ufuk Beyaztas}
\note{Requires the \pkg{coxrobust} package for robust estimation and \code{getRPCA()} for robust FPCA.}
\seealso{
\code{\link{data_generation}}, \code{\link{rfcox_predict}}, \code{\link[coxrobust]{coxr}}
}
\examples{
# Simulate data
sim_data <- data_generation(n = 500, out.p = 0.05)

# Fit robust functional Cox model
fit <- rfcox(time = sim_data$time, status = sim_data$event,
             X = sim_data$Xt, Z = sim_data$Z, gp = sim_data$gp)
}
