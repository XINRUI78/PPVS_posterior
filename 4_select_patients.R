# Randomly select 2 low-, 2 medium-, and 2 high-risk patients
# from three equal-sized risk groups

select_patients <- function(p_mean) {
  
  n <- length(p_mean)
  ord <- order(p_mean)
  
  groups <- cut(
    seq_len(n),
    breaks = 3,
    labels = c("low", "medium", "high")
  )
  
  low_ids <- sample(ord[groups == "low"], 2)
  medium_ids <- sample(ord[groups == "medium"], 2)
  high_ids <- sample(ord[groups == "high"], 2)
  
  c(low_ids, medium_ids, high_ids)
}

# Keep posterior probability draws for selected six patients only
# Output has 6 columns: low1, low2, medium1, medium2, high1, high2

make_patient_draws <- function(p_draws, selected_ids, yval) {
  
  if (ncol(p_draws) != length(yval)) {
    p_draws <- t(p_draws)
  }
  
  p_small <- p_draws[, selected_ids, drop = FALSE]
  
  colnames(p_small) <- c(
    "low1",
    "low2",
    "medium1",
    "medium2",
    "high1",
    "high2"
  )
  
  as.data.frame(p_small)
}