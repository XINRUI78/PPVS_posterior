source("1_parameters.R")
source("2_data_generation.R")

# ============================================================
# Simulation settings
# ============================================================

R <- 100

dir.create(
  "simulation_data",
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Generate datasets for all scenarios and sample sizes
# ============================================================

for (s in seq_len(5)) {

  beta0 <- beta0_list[[s]]
  beta  <- beta_list[[s]]
  sigma <- sigma_list[[s]]

  for (n_current in ndev_list) {

    for (i in seq_len(R)) {

      # Unique reproducible seed for each combination
      seed <- 100000+s*100000+n_current*100+i

      set.seed(seed)

      # ------------------------------------------------------
      # Development dataset
      # ------------------------------------------------------

      data.dev <- generate_ss(
        n = n_current,
        n.para = n.para,
        beta0 = beta0,
        beta = beta,
        sigma = sigma
      )

      # ------------------------------------------------------
      # Validation dataset
      # ------------------------------------------------------

      data.val <- generate_ss(
        n = nval,
        n.para = n.para,
        beta0 = beta0,
        beta = beta,
        sigma = sigma
      )

      # ------------------------------------------------------
      # Save both datasets
      # ------------------------------------------------------

      file_name <- sprintf(
        "simulation_data/s%02d_ndev_%d_rep_%03d.rds",
        s,
        n_current,
        i
      )

      saveRDS(
        list(
          scenario = s,
          simulation_id = i,
          ndev = n_current,
          seed = seed,
          data.dev = data.dev,
          data.val = data.val
        ),
        file = file_name
      )
    }
  }
}
