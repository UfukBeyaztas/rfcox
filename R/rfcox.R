rfcox <- function(time, status, X, Z, trunc = 0.9, nb = NULL, gp = NULL,
                  f.weight = c("linear", "quadratic", "exponential"))
{

  n <- dim(X)[1]
  j <- dim(X)[2]

  if(is.null(gp))
    gp <- seq(0, 1, length.out = j)

  if(is.null(nb))
    nb <- 5

  if(is.null(f.weight)){
    f.weight <- "quadratic"
  }

  rfpca <- getRPCA(data = X, nbasis = nb, gp = gp)
  rfscore <- rfpca$PCAscore
  evalbase <- rfpca$evalbase
  rPCAcoef <- rfpca$PCAcoef

  model_mat <- cbind(rfscore, Z)
  model_matrix <- data.frame(cbind(time, status, model_mat))
  for(i in 3:dim(model_matrix)[2])
    colnames(model_matrix)[i] = paste("V", (i-2), sep = "")

  var_name <- paste(colnames(model_matrix)[-(1:2)], collapse = "+")

  cox_formula <- as.formula(paste("Surv(time, status)~", var_name, sep = ""))

  rob_model <- coxr(cox_formula, data = model_matrix , trunc = trunc,
                    f.weight = f.weight)

  rob_coefs <- rob_model$coefficients

  linear_predictors <- coxrobust:::predict.coxr(rob_model, model_mat)

  conc <- compute_c_index(time, status, linear_predictors)

  bhat_t <- evalbase %*% (rPCAcoef$coefs %*% as.matrix(rob_coefs[1:dim(rfscore)[2]]))
  gamma_hat <- rob_coefs[-(1:dim(rfscore)[2])]

  return(list(bhat = bhat_t, gammahat = gamma_hat, concordance = conc,
              model = rob_model, rfpca = rfpca))
  }
