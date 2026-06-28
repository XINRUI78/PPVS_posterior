source("0_install_packages.R")
source("1_parameter.R")
source("2_data_generation.R")
source("3_measures.R")
source("4_select_patients.R")
source("5_la_draw.R")

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

out <- run_la(
  i = sim_id,
  ndev = ndev_current,
  nval = nval,
  n.para = n.para,
  beta0 = scenario$beta0,
  beta = scenario$beta,
  n.true = scenario$n.true,
  nterms_max = 30,
  ns = 2000,
  save_draw_details = sim_id <= 10
)

result_dir <- paste0("results/scenario", scenario_id, "_ndev", ndev_id, "_nor")

dir.create(file.path(result_dir, "summary"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(result_dir, "draw_measures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(result_dir, "patient_probs"), recursive = TRUE, showWarnings = FALSE)

base_name <- paste0(
  "scenario", scenario_id,
  "_ndev", ndev_id,
  "_sim", sim_id,
  "_nor"
)

saveRDS(
  out$summary,
  file = file.path(result_dir, "summary", paste0(base_name, "_summary.rds"))
)

if (!is.null(out$draw_measures_all)) {
  saveRDS(
    out$draw_measures_all,
    file = file.path(result_dir, "draw_measures", paste0(base_name, "_draw_measures.rds"))
  )
}

if (!is.null(out$patient_prob_draws)) {
  saveRDS(
    out$patient_prob_draws,
    file = file.path(result_dir, "patient_probs", paste0(base_name, "_patient_probs.rds"))
  )
}