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

set.seed(999)

patient_pool <- generate_ss(
  n = 1000,
  n.para = n.para,
  beta0 = beta0,
  beta = beta,
  n.true = n.true
)

patient_pool_x <- patient_pool[
  ,
  -1,
  drop = FALSE
]

patient_pool_eta <- as.numeric(
  beta0 +
    as.matrix(patient_pool_x) %*% beta
)

patient_pool_ptrue <- plogis(
  patient_pool_eta
)

fixed_patient_ids <- select_patients(
  p_true = patient_pool_ptrue
)

patient_names <- c(
  "low1",
  "low2",
  "medium1",
  "medium2",
  "high1",
  "high2"
)

fixed_patients <- patient_pool_x[
  fixed_patient_ids,
  ,
  drop = FALSE
]

rownames(fixed_patients) <- patient_names

fixed_patient_truth <- patient_pool_ptrue[
  fixed_patient_ids
]

names(fixed_patient_truth) <- patient_names

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
