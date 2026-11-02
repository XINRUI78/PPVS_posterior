source("1_parameters.R")
source("2_data_generation.R")

R <- 100

dir.create("simulation_data", showWarnings = FALSE)

for (s in 1:5) {

  # Select scenario-specific parameters
  beta0 <- beta0_list[s]
  beta  <- beta_list[[s]]
  sigma <- sigma_list[[s]]

  for (i in 1:R) {

    set.seed(100000 + s * 1000 + i)

    data.dev <- generate_ss(
      n = ndev,
      n.para = n.para,
      beta0 = beta0,
      beta = beta,
      sigma = sigma
    )

    data.val <- generate_ss(
      n = nval,
      n.para = n.para,
      beta0 = beta0,
      beta = beta,
      sigma = sigma
    )

    saveRDS(
      list(
        scenario = s,
        simulation_id = i,
        data.dev = data.dev,
        data.val = data.val
      ),
      sprintf(
        "simulation_data/s%02d_rep_%03d.rds",
        s, i
      )
    )
  }
}
