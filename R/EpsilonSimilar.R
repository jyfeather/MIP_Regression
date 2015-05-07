# generate epsilon-similar problems via changing objective coefficients
similar_obj <- function(obj) {
  set.seed(runif(1, 0, 2000))
  v <- obj$v  
  noise <- round(rnorm(length(v), mean = mean(v), sd = sd(v))/100)  
  v <- v + noise
  newobj <- simple_triplet_matrix(obj$i, obj$j, v, obj$nrow, obj$ncol, obj$dimnames)
  return(newobj)
}

# epsilon-similar via constraints
similar_cons <- function(cons) {
  set.seed(runif(1, 0, 2000))
  noise <- runif(1, -0.1, 0.1) * cons
  cons <- cons + noise
  return(cons)
}