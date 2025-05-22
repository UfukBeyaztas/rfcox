rfcox_predict <- function(object, time, status, X, Z)
{

  rfpca <- object$rfpca
  rfscore <- getRPCA_test(rfpca, X)

  model_mat <- cbind(rfscore, Z)
  for(i in 1:dim(model_mat)[2])
    colnames(model_mat)[i] = paste("V", i, sep = "")

  rcox_model <- object$model
  lpreds <- coxrobust:::predict.coxr(rcox_model, model_mat)
  conc <- compute_c_index(time, status, lpreds)

  return(list(predictions = lpreds, concordance = conc))
}
