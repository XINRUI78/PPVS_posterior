source("0_install_packages.R")
source("1_parameter.R")
source("2_data_generation.R")
source("5_run_hs_p.R")


sim_id <- as.integer(Sys.getenv("SGE_TASK_ID"))

scenario_id <- as.integer(Sys.getenv("SCENARIO_ID"))
ndev_id <- as.integer(Sys.getenv("NDEV_ID"))

ndev_list <- c(ndev, round(ndev / 2), round(ndev / 4))
ndev_current <- ndev_list[ndev_id]

scenarios <- list(
  list(beta0 = beta0_1, beta = beta_1, n.true = NULL),
  list(beta0 = beta0_2, beta = beta_2, n.true = NULL),
  list(beta0 = beta0_3, beta = beta_3, n.true = NULL),
  list(beta0 = beta0_4, beta = beta_4, n.true = n.true)
)

scenario <- scenarios[[scenario_id]]

fixed_patient_file <- paste0(
  "scenario",
  scenario_id,
  "_fixed_patients.rds"
)

fixed_patients <- readRDS(
  fixed_patient_file
)

fixed_patients <- as.data.frame(
  fixed_patients
)


# ============================================================
# Run horseshoe analysis
# ============================================================
out <- run_hs_p(
  i = sim_id,
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

result_dir <- paste0("results/scenario", scenario_id, "_ndev", ndev_id, "_hs")

mean_prediction_dir <- file.path(
  result_dir,
  "patient_mean_predictions"
)

patient_draw_dir <- file.path(
  result_dir,
  "patient_draws"
)

dir.create(mean_prediction_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(patient_draw_dir, recursive = TRUE, showWarnings = FALSE)

base_name <- paste0(
  "scenario", scenario_id,
  "_ndev", ndev_id,
  "_sim", sim_id,
  "_hs"
)

saveRDS(
  out$patient_mean_predictions,
  file = file.path(
    mean_prediction_dir,
    paste0(
      base_name,
      "_patient_mean_predictions.rds"
    )
  )
)


if (!is.null(out$patient_draws)) {
  
  saveRDS(
    out$patient_draws,
    file = file.path(
      patient_draw_dir,
      paste0(
        base_name,
        "_patient_draws.rds"
      )
    )
  )
}


