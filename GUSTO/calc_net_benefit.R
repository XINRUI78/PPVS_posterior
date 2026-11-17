calc_net_benefit <- function(yval, p_val, thresholds, method) {
  
  n <- length(yval)
  
  nb <- sapply(thresholds, function(pt) {
    pred_positive <- p_val >= pt
    TP <- sum(pred_positive & yval == 1, na.rm = TRUE)
    FP <- sum(pred_positive & yval == 0, na.rm = TRUE)
    
    TP / n - FP / n * (pt / (1 - pt))
  })
  
  data.frame(
    method = method,
    threshold = thresholds,
    net_benefit = nb
  )
}