
# ============================================================
# generate_fixed_patients.R
# Generate six fixed patients for each of five scenarios
# ============================================================

source("1_parameters.R")
source("2_data_generation.R")

select_patients <- function(p_true) {

  p_true <- as.numeric(p_true)

  targets <- as.numeric(
    quantile(
      p_true,
      probs = c(0.10, 0.50, 0.90)
    )
  )

  choose_two <- function(target, excluded = integer(0)) {

    available <- setdiff(
      seq_along(p_true),
      excluded
    )

    head(
      available[
        order(abs(p_true[available] - target))
      ],
      2L
    )
  }

  low_ids <- choose_two(targets[1])

  medium_ids <- choose_two(
    targets[2],
    excluded = low_ids
  )

  high_ids <- choose_two(
    targets[3],
    excluded = c(low_ids, medium_ids)
  )

  c(low_ids, medium_ids, high_ids)
}

patient_names <- c(
  "low1", "low2",
  "medium1", "medium2",
  "high1", "high2"
)

for (scenario_id in 1:5) {

  beta0_current <- beta0_list[[scenario_id]]
  beta_current <- beta_list[[scenario_id]]
  sigma_current <- sigma_list[[scenario_id]]

  # Reproducible and distinct pool for each scenario
  set.seed(999 + scenario_id)

  patient_pool <- generate_ss(
    n = 10000,
    n.para = n.para,
    beta0 = beta0_current,
    beta = beta_current,
    sigma = sigma_current
  )

  # Assumes first column is outcome y
  patient_pool_x <- patient_pool[
    , setdiff(names(patient_pool), "y"),
    drop = FALSE
  ]

  patient_pool_eta <- as.numeric(
    beta0_current +
      as.matrix(patient_pool_x) %*%
      as.numeric(beta_current)
  )

  patient_pool_ptrue <- plogis(patient_pool_eta)

  fixed_patient_ids <- select_patients(
    patient_pool_ptrue
  )

  fixed_patients <- patient_pool_x[
    fixed_patient_ids,
    ,
    drop = FALSE
  ]

  rownames(fixed_patients) <- patient_names

  # Store true probabilities alongside predictors
  fixed_patients$p_true <- as.numeric(
    patient_pool_ptrue[fixed_patient_ids]
  )

  saveRDS(
    fixed_patients,
    sprintf(
      "scenario%d_fixed_patients.rds",
      scenario_id
    )
  )
}
