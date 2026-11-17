# Source required scripts

source("perform_gusto.R")
source("measures.R")
source("calc_net_benefit.R")
source("backward_pvalue.R")
source("../method/berank.R")
source("../method/lasso_exact.R")
source("../method/mod_penal_ave_foreach.R")
source("../method/unilogit.R")
source("../method/unirank.R")

library(doParallel)
library(foreach)
library(RcppNumerical)
library(brglm2)

gusto_cut <- read.csv("gusto_cut.csv")
order_2n <- read.csv("order_2n.csv")
order_n <- read.csv("order_n.csv")
order_n_2 <- read.csv("order_n_2.csv")
order_n_4 <- read.csv("order_n_4.csv")

n.para <- 20
ndev <- 2277
ndev1 <- round(ndev/2)
ndev2 <- round(ndev/4)
ndev3 <- 2*ndev

order_2n <- as.numeric(unlist(order_2n))
order_n <- as.numeric(unlist(order_n))
order_n_2 <- as.numeric(unlist(order_n_2))
order_n_4 <- as.numeric(unlist(order_n_4))

ex_n <- perform_gusto(order_n, gusto_cut, n.para)
ex_n_2 <- perform_gusto(order_n_2, gusto_cut, n.para)
ex_n_4 <- perform_gusto(order_n_4, gusto_cut, n.para)
ex_2n <- perform_gusto(order_2n, gusto_cut, n.para)

# Save standard performance results
write.csv(ex_n$method_result, "result_n.csv", row.names = FALSE)
write.csv(ex_n_2$method_result, "result_n_2.csv", row.names = FALSE)
write.csv(ex_n_4$method_result, "result_n_4.csv", row.names = FALSE)
write.csv(ex_2n$method_result, "result_2n.csv", row.names = FALSE)

# Save net benefit results
write.csv(ex_n$net_benefit, "net_benefit_n.csv", row.names = FALSE)
write.csv(ex_n_2$net_benefit, "net_benefit_n_2.csv", row.names = FALSE)
write.csv(ex_n_4$net_benefit, "net_benefit_n_4.csv", row.names = FALSE)
write.csv(ex_2n$net_benefit, "net_benefit_2n.csv", row.names = FALSE)
