
# ============================================================
# 5_hs_array.R
#
# Myriad SGE array runner
#
# 5 scenarios
# 3 development sample sizes
# 100 replications per scenario/sample size
#
# Task 1: Validation performance and variable selection
# Task 2: Fixed-patient posterior mean and mode
# Task 3: Posterior probability distributions
#
# Final folder structure:
# results/
#   s1/
#     ndev_1505/
#       _replications/
#       task1_validation.csv
#       task2_patients.csv
#       task3_posterior.rds
#
# The final three files are produced separately by
# 6_combine_hs_results.R after all replications finish.
# ============================================================


# ------------------------------------------------------------
# 1. Load parameters and simulation function
# ------------------------------------------------------------

source("1_parameters.R")
source("3_measures.R")
source("run_hs_combined.R")


# ------------------------------------------------------------
# 2. Read Myriad SGE array variables
# ------------------------------------------------------------

sim_id <- as.integer(
  Sys.getenv("SGE_TASK_ID")
)

scenario_id <- as.integer(
  Sys.getenv("SCENARIO_ID")
)

ndev_id <- as.integer(
  Sys.getenv("NDEV_ID")
)

ndev_current <- ndev_list[[ndev_id]]

# ------------------------------------------------------------
# 5. Scenario-specific true coefficients
# ------------------------------------------------------------

beta0_current <- beta0_list[[scenario_id]]

beta_current <- beta_list[[scenario_id]]


# ------------------------------------------------------------
# 6. Load six fixed patients
# ------------------------------------------------------------

fixed_patient_file <- sprintf(
  "scenario%d_fixed_patients.rds",
  scenario_id
)

fixed_patients <- readRDS(
  fixed_patient_file
)

# ------------------------------------------------------------
# 7. Check required simulation dataset exists
# ------------------------------------------------------------

data_file <- simulation_file(
  scenario_id = scenario_id,
  ndev = ndev_current,
  i = sim_id,
  data_dir = "simulation_data"
)

# ------------------------------------------------------------
# 8. Run the three analyses using the same fitted models
# ------------------------------------------------------------

out <- run_hs(
  i = sim_id,
  scenario_id = scenario_id,
  ndev = ndev_current,
  n.para = n.para,
  beta0 = beta0_current,
  beta = beta_current,
  fixed_patients = fixed_patients,
  nterms_max = 30,
  prediction_draws = 400,
  projection_draws = 400,
  data_dir = "simulation_data"
)

# ------------------------------------------------------------
# 9. Define result directories
# ============================================================

result_dir <- file.path(
  "results",
  paste0("s", scenario_id),
  paste0("ndev_", ndev_current)
)

task1_dir <- file.path(result_dir, "measures")
task2_dir <- file.path(result_dir, "fixed patients")
task3_dir <- file.path(result_dir, "posterior")

dir.create(task1_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(task2_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(task3_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 10. Save Task 1: Validation performance
# ============================================================

measures_rds <- file.path(
  task1_dir,
  sprintf("m_rep_%03d.rds", sim_id)
)

measures_csv <- file.path(
  task1_dir,
  sprintf("m_rep_%03d.csv", sim_id)
)

saveRDS(
  out$summary,
  measures_rds
)

write.csv(
  out$summary,
  file = measures_csv,
  row.names = FALSE
)


# ============================================================
# 11. Save Task 2: Fixed-patient posterior mean and mode
# ============================================================

fixed_rds <- file.path(
  task2_dir,
  sprintf("fix_rep_%03d.rds", sim_id)
)

fixed_csv <- file.path(
  task2_dir,
  sprintf("fix_rep_%03d.csv", sim_id)
)

saveRDS(
  out$patient_predictions,
  fixed_rds
)

write.csv(
  out$patient_predictions,
  file = fixed_csv,
  row.names = FALSE
)


# ============================================================
# 12. Save Task 3: Posterior distributions
# Only replication 1 at ndev = 1505
# ============================================================

if (!is.null(out$posterior_distributions)) {

  pd_file <- file.path(
    task3_dir,
    sprintf("pd_rep_%03d.rds", sim_id)
  )

  saveRDS(
    out$posterior_distributions,
    pd_file
  )
}
