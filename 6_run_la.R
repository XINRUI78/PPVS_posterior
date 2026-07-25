library(tidyverse)
library(rstanarm)
library(projpred)
library(loo)
library(glmnet)
library(foreach)
library(doParallel)
library(MASS)

source("0_install_packages.R")
source("1_parameter.R")
source("2_data_generation.R")
source("5_run_la_p.R")


# ============================================================
# Read scenario number from command line
# ============================================================

args <- commandArgs(
  trailingOnly = TRUE
)

if (length(args) < 1) {
  stop(
    "Supply a scenario number: 1, 2, 3, or 4."
  )
}

scenario_id <- as.integer(
  args[1]
)




Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

options(mc.cores = 1)


# ============================================================
# Scenario parameters
# ============================================================

scenario_parameters <- list(
  
  `1` = list(
    beta0 = beta0_1,
    beta = beta_1,
    n.true = NULL
  ),
  
  `2` = list(
    beta0 = beta0_2,
    beta = beta_2,
    n.true = NULL
  ),
  
  `3` = list(
    beta0 = beta0_3,
    beta = beta_3,
    n.true = NULL
  ),
  
  `4` = list(
    beta0 = beta0_4,
    beta = beta_4,
    n.true = 15
  )
)

scenario <- scenario_parameters[
  [as.character(scenario_id)]
]

ndev_current <- ndev

# ============================================================
# Load fixed patients for this scenario
# ============================================================

fixed_patient_file <- paste0(
  "scenario",
  scenario_id,
  "_fixed_patients.rds"
)


fixed_patients <- as.data.frame(
  readRDS(fixed_patient_file)
)

n_cores <- as.integer(
  Sys.getenv(
    "NSLOTS",
    unset = "32"
  )
)

n_cores <- max(
  1L,
  n_cores
)

cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

on.exit(
  parallel::stopCluster(cl),
  add = TRUE
)

# ============================================================
# Run 100 simulations
# ============================================================

results <- foreach(
  i = seq_len(100),
  
  .packages = c(
    "tidyverse",
    "rstanarm",
    "projpred",
    "loo",
    "glmnet",
    "MASS", 
    "mvtnorm", 
    "speedglm"
  ),
  .export = c(
    "run_la_p",
    "generate_ss",
    "orient_draw_matrix"
  ),
  
  .combine = "c",
  .multicombine = TRUE,
  .errorhandling = "stop"
) %dopar% {
  
  out <- run_la_p(
    i = i,
    ndev = ndev_current,
    n.para = n.para,
    beta0 = scenario$beta0,
    beta = scenario$beta,
    fixed_patients = fixed_patients,
    n.true = scenario$n.true,
    nterms_max = 30,
    prediction_draws = 400,
    projection_draws = 400
  )
  
  list(out)
}


# ============================================================
# Combine all posterior mean predictions
# ============================================================

all_mean_predictions <- bind_rows(
  lapply(
    results,
    function(result) {
      result$patient_mean_predictions
    }
  )
) |>
  arrange(
    simulation,
    method,
    patient
  )


# ============================================================
# Keep full posterior draws for simulations 1 to 10
# ============================================================

all_patient_draws <- lapply(
  results,
  function(result) {
    result$patient_draws
  }
)

all_patient_draws <- Filter(
  Negate(is.null),
  all_patient_draws
)

# ============================================================
# Save one combined CSV and RDS per scenario
# ============================================================

result_dir <- file.path(
  "results",
  paste0(
    "scenario",
    scenario_id,
    "_ndev_la"
  )
)

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

mean_rds_file <- file.path(
  result_dir,
  paste0(
    "scenario",
    scenario_id,
    "_ndev",
    ndev_current,
    "_la_100_simulations.rds"
  )
)

mean_csv_file <- file.path(
  result_dir,
  paste0(
    "scenario",
    scenario_id,
    "_ndev",
    ndev_current,
    "_la_100_simulations.csv"
  )
)

draw_rds_file <- file.path(
  result_dir,
  paste0(
    "scenario",
    scenario_id,
    "_ndev",
    ndev_current,
    "_la_draws_simulations_1_to_10.rds"
  )
)

saveRDS(
  all_mean_predictions,
  file = mean_rds_file
)

write.csv(
  all_mean_predictions,
  file = mean_csv_file,
  row.names = FALSE
)

saveRDS(
  all_patient_draws,
  file = draw_rds_file
)

