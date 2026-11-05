run_uvs_only <- function(i, scenario_id, ndev, n.para,
                         fixed_patients,
                         data_dir = "simulation_data") {

  data_file <- simulation_file(
    scenario_id, ndev, i, data_dir
  )

  sim_data <- readRDS(data_file)
  data.dev <- sim_data$data.dev

  predictor_names <- paste0("X", seq_len(n.para))

  x <- data.dev[, predictor_names, drop = FALSE]
  y <- as.numeric(data.dev$y)

  fixed_x <- fixed_patients[
    , predictor_names, drop = FALSE
  ]

  patient_result <- list()

  for (threshold in c(0.05, 0.15)) {

    uni <- unilogit(x, y, threshold)

    selected <- uni$varsel_uni == 1

    eta <- as.numeric(
      as.matrix(
        cbind(1, fixed_x[, selected, drop = FALSE])
      ) %*% coef(uni$uni)
    )

    p <- plogis(eta)

    method_name <- if (threshold == 0.05) {
      "UVS-5%"
    } else {
      "UVS-15%"
    }

    patient_result[[method_name]] <- data.frame(
      simulation = i,
      method = method_name,
      patient_id = seq_len(nrow(fixed_patients)),
      p_true = as.numeric(fixed_patients$p_true),
      p_estimated = p
    )
  }

  result <- do.call(rbind, patient_result)
  rownames(result) <- NULL

  result
}
