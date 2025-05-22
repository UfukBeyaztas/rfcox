\name{data_generation}
\alias{data_generation}
\title{
Data Generation for Robust Functional Cox Regression Simulations
}
\description{
Generates synthetic data for evaluating robust functional Cox regression models, including both functional and scalar covariates, censoring, and contamination mechanisms.
}
\usage{
data_generation(n = 1000, j = 101, p = 4, gamma0 = rep(0.5, p),
                h0t = 1, cenrate = 0.5, tau = 1.55, mev = 0.5,
                out.p = 0)
}
\arguments{
  \item{n}{Number of observations to generate.}
  \item{j}{Number of grid points used to discretize the functional domain \([0,1]\).}
  \item{p}{Number of scalar covariates.}
  \item{gamma0}{True regression coefficients for the scalar covariates.}
  \item{h0t}{Baseline hazard rate (assumed constant).}
  \item{cenrate}{Target censoring rate.}
  \item{tau}{Upper bound for censoring time sampling (via uniform distribution).}
  \item{mev}{Measurement error variance added to functional predictors.}
  \item{out.p}{Proportion of observations to be contaminated as outliers.}
}
\details{
This function simulates data under a robust functional Cox regression framework with one functional covariate and multiple scalar covariates. The functional predictor is generated from a truncated Karhunenâ€“Loeve expansion with correlated components, while scalar covariates follow a multivariate Gaussian distribution with specified correlation.

The true coefficient function is a smooth nonlinear combination of sinusoidal and Gaussian terms, reflecting rich temporal effects. Survival times are drawn from an exponential distribution with a hazard dependent on both the functional covariate (via numerical integration) and scalar predictors. Right-censoring is imposed via independent uniform random times.

Measurement error is introduced into the functional predictor to reflect realistic noise. Additionally, contamination can be injected via \code{out.p} to simulate outlier scenarios, affecting both the functional trajectories and survival times. Outlier functional data are perturbed using a magnitude-shift Ornsteinâ€“Uhlenbeck process, while survival times are replaced with extreme early or late values.
}
\value{
A list with the following components:
  \item{time}{Observed time-to-event or censoring (numeric vector of length \code{n}).}
  \item{event}{Event indicator (1 = event, 0 = censored).}
  \item{Xt}{Noisy functional predictor matrix (dimension \code{n} x \code{j}).}
  \item{Xt_true}{True smooth functional predictor matrix (before noise).}
  \item{Z}{Scalar covariate matrix (dimension \code{n} x \code{p}).}
  \item{beta}{True coefficient function evaluated on the grid points.}
  \item{gamma}{True scalar coefficients.}
  \item{gp}{Grid points over domain [0,1].}
}

\author{
Gizel Bakicierler Sezer and Ufuk Beyaztas
}
\note{
This function is primarily intended for simulation studies examining robustness properties in functional Cox regression modeling. For robust estimation methods, see the corresponding functions using RFPCA and M-type partial likelihood.
}
\examples{
# Generate synthetic data with default settings
sim_data <- data_generation()

}

