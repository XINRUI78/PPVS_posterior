select_patients <- function(p_true) {
  
  p_true <- as.numeric(p_true)
  
  targets <- as.numeric(
    quantile(
      p_true,
      probs = c(0.1, 0.50, 0.9),
      na.rm = TRUE
    )
  )
  
  choose_two <- function(target, excluded = integer(0)) {
    
    available <- setdiff(
      seq_along(p_true),
      excluded
    )
    
    available[
      order(
        abs(p_true[available] - target)
      )[1:2]
    ]
  }
  
  low_ids <- choose_two(
    target = targets[1]
  )
  
  medium_ids <- choose_two(
    target = targets[2],
    excluded = low_ids
  )
  
  high_ids <- choose_two(
    target = targets[3],
    excluded = c(low_ids, medium_ids)
  )
  
  c(
    low_ids,
    medium_ids,
    high_ids
  )
}

# Keep posterior probability draws for selected six patients only
# Output has 6 columns: low1, low2, medium1, medium2, high1, high2
make_patient_draws <- function(
    p_draws,
    selected_ids,
    n_observations
) {
  
  p_draws <- orient_draw_matrix(
    draw_matrix = p_draws,
    n_observations = n_observations
  )
  
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
