source('1_parameters.R')
source('3_measures.R')
source('prune_lasso_predictor.R')
source('unilogit.R')
source('backward_pvalue.R')
source('mod_penal_ave_foreach.R')
source('run_freq.R')


# ============================================================
# 5_freq_array.R
#
# Myriad SGE array runner
#
# 5 scenarios
# 3 development sample sizes
# 100 replications per scenario/sample size
#
# Nine frequentist methods:
#   1. MLE
#   2. UVS-5%
#   3. UVS-15%
#   4. BE-5%
#   5. BE-15%
#   6. LS-min
#   7. step-min
#   8. LS-mod
#   9. step-mod
#
# Task 1: Validation performance and variable selection
# Task 2: Fixed-patient predicted probabilities
#
# Results:
# results_freq/
#   s1/
#     ndev_1505/
#       measures/
#         m_rep_001.rds
#         m_rep_001.csv
#       fixed patients/
#         fix_rep_001.rds
#         fix_rep_001.csv
# ============================================================


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
# 3. Scenario-specific true coefficients
# ------------------------------------------------------------

beta0_current <- beta0_list[[scenario_id]]

beta_current <- beta_list[[scenario_id]]


# ------------------------------------------------------------
# 4. Load six fixed patients
# ------------------------------------------------------------

fixed_patient_file <- sprintf(
  "scenario%d_fixed_patients.rds",
  scenario_id
)

fixed_patients <- readRDS(
  fixed_patient_file
)


# ------------------------------------------------------------
# 5. Check required simulation dataset exists
# ------------------------------------------------------------
simulation_file <- function(
    scenario_id, ndev, i,
    data_dir = "simulation_data"
) {
  file.path(
    data_dir,
    sprintf(
      "s%02d_ndev_%d_rep_%03d.rds",
      scenario_id, ndev, i
    )
  )
}

data_file <- simulation_file(
  scenario_id = scenario_id,
  ndev = ndev_current,
  i = sim_id,
  data_dir = "simulation_data"
)

if (!file.exists(data_file)) {
  stop("Simulation dataset not found: ", data_file)
}


# ------------------------------------------------------------
# 6. Run the nine frequentist methods
# ------------------------------------------------------------

out <- run_freq(
  i = sim_id,
  scenario_id = scenario_id,
  ndev = ndev_current,
  n.para = n.para,
  beta0 = beta0_current,
  beta = beta_current,
  fixed_patients = fixed_patients,
  data_dir = "simulation_data"
)


# ------------------------------------------------------------
# 7. Define result directories
# ------------------------------------------------------------

result_dir <- file.path(
  "results_freq",
  paste0("s", scenario_id),
  paste0("ndev_", ndev_current)
)

task1_dir <- file.path(
  result_dir,
  "measures"
)

task2_dir <- file.path(
  result_dir,
  "fixed patients"
)

dir.create(
  task1_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  task2_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 8. Save Task 1: Validation performance and variable selection
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
  out$performance,
  measures_rds
)

write.csv(
  out$performance,
  file = measures_csv,
  row.names = FALSE
)


# ============================================================
# 9. Save Task 2: Fixed-patient predicted probabilities
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

