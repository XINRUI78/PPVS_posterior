library(rms)
library(data.table)
# gusto data

load(url('https://hbiostat.org/data/repo/gusto.rda'))

# setDT(gusto)
gusto

save(gusto, file="gusto.Rdata")
#############################
# full model with 25 parameters, cstat: 0.826
gusto[,c(19,21:24)] <- gusto[,c(19,21:24)]-1

levels(gusto$miloc)[levels(gusto$miloc) == "Inferior"] <- "Other"
gusto$miloc <- factor(gusto$miloc, levels = c("Anterior", "Other"))
levels(gusto$smk)[levels(gusto$smk) != "never"] <- "Other"
gusto$smk <- factor(gusto$smk, levels = c("never", "Other"))
# Killip: I, Other(II, III, IV)
levels(gusto$Killip)[levels(gusto$Killip) != "I"] <- "Other"
gusto$Killip <- factor(gusto$Killip, levels = c("I", "Other"))
levels(gusto$tx)[levels(gusto$tx) != "tPA"] <- "Other"
gusto$tx <- factor(gusto$tx, levels = c("tPA", "Other"))
# cut down the dataset
vars <- c("day30", "age", "sysbp", "Killip", "pulse", "miloc", "pmi", 
          "height", "smk", "dia", "hrt", "fam", "ste", "ttr", "weight", 
          "prevcabg", "prevcvd", "tx", "htn", "sex", "pan")
gusto_cut <- gusto[, vars]
colnames(gusto_cut)[1] <- "y"
gusto_cut[] <- lapply(gusto_cut, function(x) {
  if (inherits(x, "labelled") && is.factor(x)) {
    factor(x)
  } else if (inherits(x, "labelled")) {
    as.numeric(x)
  } else {
    x
  }
})
##########sample size
fit_all <- lrm(y ~ ., data = gusto_cut)
c = 0.823
n.para = 20
prev = 0.07

install.packages("devtools")
require("devtools")
devtools::install_github("mpavlou/samplesizedev")
require(samplesizedev)
library(samplesizedev)

rss <- samplesizedev(outcome = "Binary", S = 0.9, phi = prev, c = c, p = n.para)
ndev <- rss$sim
ndev1 <- round(ndev/2)
ndev2 <- round(ndev/4)
ndev3 <- 2*ndev
#########################
set.seed(1)   # optional, for reproducibility
order_n <- sample(1:nrow(gusto_cut), ndev)
order_n_2 <- sample(1:nrow(gusto_cut), ndev1)
order_n_4 <- sample(1:nrow(gusto_cut), ndev2)
order_2n <- sample(1:nrow(gusto_cut), ndev3)
