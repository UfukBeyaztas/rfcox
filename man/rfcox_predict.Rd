\name{rfcox_predict}
\alias{rfcox_predict}
\title{
Prediction for Robust Functional Cox Regression Models
}
\description{
Generates linear predictors and computes the concordance index for new survival data using a fitted robust functional Cox regression model.
}
\usage{
rfcox_predict(object, time, status, X, Z)
}
\arguments{
  \item{object}{An object returned by \code{\link{rfcox}}, containing the fitted robust model and associated RFPCA components.}
  \item{time}{A numeric vector of observed survival or censoring times for the new data.}
  \item{status}{A numeric vector of event indicators for the new data: 1 for event, 0 for censored.}
  \item{X}{A matrix of functional predictors for new observations (dimension: \code{n} x \code{j}).}
  \item{Z}{A matrix of scalar covariates for new observations (dimension: \code{n} x \code{p}).}
}
\details{
This function performs prediction using a robust functional Cox regression model fitted via \code{\link{rfcox}}. The procedure includes:
\enumerate{
  \item Projecting new functional data onto the robust FPCA space from training,
  \item Combining FPCA scores with scalar covariates,
  \item Predicting risk scores (linear predictors) using the fitted robust Cox model,
  \item Calculating the concordance index based on new survival outcomes.
}
This ensures that prediction is robust to outliers in both functional and scalar covariates.
}
\value{
A list with the following components:
  \item{predictions}{A numeric vector of linear predictors for the test observations.}
  \item{concordance}{Concordance index comparing predictions to observed survival times and status.}
}
\author{
Gizel Bakicierler Sezer and Ufuk Beyaztas
}
\note{
Requires a fitted model object from \code{rfcox}. The function uses \code{getRPCA_test} for projecting new functional data.
}
\examples{
# Simulate training and test data
train_data <- data_generation(n = 500, out.p = 0.05)
test_data <- data_generation(n = 1000)

# Fit the model
fit <- rfcox(time = train_data$time, status = train_data$event,
             X = train_data$Xt, Z = train_data$Z, gp = train_data$gp)

# Predict on test data
pred_out <- rfcox_predict(object = fit,
                          time = test_data$time,
                          status = test_data$event,
                          X = test_data$Xt,
                          Z = test_data$Z)
}
