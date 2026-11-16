
# ============================================================
# 5_freq_array.R
#
# Myriad SGE array runner
#
# 5 scenarios
# 3 development sample sizes
# 100 replications per scenario/sample size
#
# Frequentist methods:
#   MLE
#   UVS-5%
#   UVS-15%
#   BE-5%
#   BE-15%
#   LS-min
#   LS-mod
#   ALASSO
#   SD-MLE
#   SD-LS-min
#   SD-LS-mod
#
# Outputs:
#   Validation performance and variable selection
#   Fixed-patient predicted probabilities
#
# Results:
# results_freq_new/
#   s1/
#     ndev_1505/
#       measures/
#         m_rep_001.rds
#         m_rep_001.csv
#       fixed patients/
#         fix_rep_001.rds
#         fix_rep_001.csv
# ============================================================


# ============================================================
# Source required functions
# ============================================================

source("1_parameters.R")
source("3_measures.R")
source("unilogit.R")
source("backward_pvalue.R")
source("mod_penal_ave_foreach.R")
source("run_freq_new.R")


# ============================================================
# Read Myriad SGE array variables
# ============================================================

sim_id <- as.integer(Sys.getenv("SGE_TASK_ID"))
scenario_id <- as.integer(Sys.getenv("SCENARIO_ID"))
ndev_id <- as.integer(Sys.getenv("NDEV_ID"))

if (anyNA(c(sim_id, scenario_id, ndev_id))) {
  stop("Missing or invalid SGE_TASK_ID, SCENARIO_ID, or NDEV_ID")
}

if (sim_id < 1L || sim_id > 100L) stop("SGE_TASK_ID must be between 1 and 100")
if (scenario_id < 1L || scenario_id > length(beta_list))
  stop("SCENARIO_ID is outside the available scenarios")
if (ndev_id < 1L || ndev_id > length(ndev_list))
  stop("NDEV_ID is outside the available sample sizes")

ndev_current <- ndev_list[[ndev_id]]


# ============================================================
# Scenario-specific true coefficients
# ============================================================

beta0_current <- beta0_list[[scenario_id]]
beta_current <- beta_list[[scenario_id]]


# ============================================================
# Load six fixed patients
# ============================================================

fixed_patient_file <- sprintf("scenario%d_fixed_patients.rds", scenario_id)

if (!file.exists(fixed_patient_file)) {
  stop("Fixed-patient dataset not found: ", fixed_patient_file)
}

fixed_patients <- readRDS(fixed_patient_file)


# ============================================================
# Locate existing simulation dataset
# ============================================================

simulation_file <- function(scenario_id, ndev, i, data_dir = "simulation_data") {
  file.path(data_dir, sprintf("s%02d_ndev_%d_rep_%03d.rds", scenario_id, ndev, i))
}

data_file <- simulation_file(
  scenario_id = scenario_id, ndev = ndev_current, i = sim_id,
  data_dir = "simulation_data"
)

if (!file.exists(data_file)) stop("Simulation dataset not found: ", data_file)


# ============================================================
# Run frequentist methods
# ============================================================

out <- run_freq_new(
  i = sim_id,
  scenario_id = scenario_id,
  ndev = ndev_current,
  n.para = n.para,
  beta0 = beta0_current,
  beta = beta_current,
  fixed_patients = fixed_patients,
  data_dir = "simulation_data"
)


# ============================================================
# Define result directories
# ============================================================

result_dir <- file.path("results_freq_new", paste0("s", scenario_id),
                        paste0("ndev_", ndev_current))

task1_dir <- file.path(result_dir, "measures")
task2_dir <- file.path(result_dir, "fixed patients")

dir.create(task1_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(task2_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# Save validation performance and variable selection
# ============================================================

measures_rds <- file.path(task1_dir, sprintf("m_rep_%03d.rds", sim_id))
measures_csv <- file.path(task1_dir, sprintf("m_rep_%03d.csv", sim_id))

saveRDS(out$performance, measures_rds)
write.csv(out$performance, file = measures_csv, row.names = FALSE)


# ============================================================
# Save fixed-patient predicted probabilities
# ============================================================

fixed_rds <- file.path(task2_dir, sprintf("fix_rep_%03d.rds", sim_id))
fixed_csv <- file.path(task2_dir, sprintf("fix_rep_%03d.csv", sim_id))

saveRDS(out$patient_predictions, fixed_rds)
write.csv(out$patient_predictions, file = fixed_csv, row.names = FALSE)
