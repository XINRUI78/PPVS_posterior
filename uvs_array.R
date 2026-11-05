# ============================================================
# uvs_array.R
#
# Myriad SGE array runner
#
# 5 scenarios
# 3 development sample sizes
# 100 replications per scenario/sample size
#
# Methods:
#   1. UVS-5%
#   2. UVS-15%
#
# Task 2: Corrected fixed-patient predicted probabilities
#
# Results:
# results_freq/
#   s1/
#     ndev_376/
#       measures/
#       fixed patients/
#       UVS_corrected/
#         uvs_rep_001.rds
#         uvs_rep_001.csv
#
# Existing results are NOT overwritten.
# ============================================================


# ------------------------------------------------------------
# 1. Load required functions and parameters
# ------------------------------------------------------------

source("1_parameters.R")
source("unilogit.R")
source("run_uvs_only.R")


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
# 3. Load six fixed patients for this scenario
# ------------------------------------------------------------

fixed_patient_file <- sprintf(
  "scenario%d_fixed_patients.rds",
  scenario_id
)

fixed_patients <- readRDS(
  fixed_patient_file
)


# ------------------------------------------------------------
# 4. Define simulation dataset path
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


# ------------------------------------------------------------
# 5. Check simulation dataset exists
# ------------------------------------------------------------

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
# 6. Run only UVS-5% and UVS-15%
# ------------------------------------------------------------

out <- run_uvs_only(
  i = sim_id,
  scenario_id = scenario_id,
  ndev = ndev_current,
  n.para = n.para,
  fixed_patients = fixed_patients,
  data_dir = "simulation_data"
)


# ------------------------------------------------------------
# 7. Define original result directory
# ------------------------------------------------------------

result_dir <- file.path(
  "results_freq",
  paste0("s", scenario_id),
  paste0("ndev_", ndev_current)
)

uvs_dir <- file.path(
  result_dir,
  "UVS_corrected"
)

dir.create(
  uvs_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 8. Save corrected UVS fixed-patient predictions
# ------------------------------------------------------------

uvs_rds <- file.path(
  uvs_dir,
  sprintf("uvs_rep_%03d.rds", sim_id)
)

uvs_csv <- file.path(
  uvs_dir,
  sprintf("uvs_rep_%03d.csv", sim_id)
)

saveRDS(
  out,
  uvs_rds
)

write.csv(
  out,
  file = uvs_csv,
  row.names = FALSE
)
