## network distance test

library(igraph)
library(apcluster)
library(microbenchmark)
library(RSpectra)

# Set the seed for reproducibility
set.seed(828)

# Generate a 1000x1000 matrix
matrix_size <- 1000
mat1 <- matrix(rbinom(matrix_size^2, 1, 0.5), nrow = matrix_size)
mat2 <- matrix(rbinom(matrix_size^2, 1, 0.5), nrow = matrix_size)

# Create two correlation networks from the matrix
get_cor_net <- function(x) {
  net <- rcorr(x, type = "spearman")
  mat_cor <- net$r
  mat_p <- net$P
  
  mat_cor[mat_p > 0.05] <- 0
  
  return(mat_cor)
}

net1 <- get_cor_net(mat1)
net2 <- get_cor_net(mat2)

get_jaccard <- function(adj_m, adj_n){
  contrast <- sum(abs(adj_m- adj_n))
  max <- sum(pmax(abs(adj_m), abs(adj_n)))
  
  dis <- contrast / max
}

get_spectra <- function(x,  k = 100){
  x <- as.matrix(x)
  x[is.na(x)] <- 0
  spectra <- eigs_sym(x, k, opts = list(retvec = FALSE))
  y <- spectra$values
}

get_spectra_dis <- function(adj_m, adj_n){
  spectra_m <- get_spectra(adj_m)
  spectra_n <- get_spectra(adj_n)
  return(sqrt(sum((spectra_m - spectra_n)^2)))
}

# Compute the distances 1000 times and measure the computation time
results <- microbenchmark(
  jaccard_distance = get_jaccard(net1, net2),
  spectra_distance = get_spectra_dis(net1, net2),
  #GED = get_GED(net1, net2)
  times = 1000
)

# Print the results
print(results)
