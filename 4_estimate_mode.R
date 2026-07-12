library(stats)
estimate_mode <- function(x) {
  
  x <- as.numeric(x)
  dens <- density(x, na.rm = TRUE)
  dens$x[which.max(dens$y)]
}
