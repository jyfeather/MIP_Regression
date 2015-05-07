################################
#   for temporary use
################################
# construct training v: predictor variable vector
v_train <- function(result) {
  return(c(result$sol, result$C$v))  
}

# construct testng v: predictor variable vector
v_test <- function(A, b, c) {
  A <- as.matrix(A)
  c <- as.vector(c)
  d <- cbind(A, b)
  d <- as.vector(t(d))
  v <- c(d, c)
  # d_B <- P^T*A_B)^-1*P^T*b
  v <- c(sample(v, nrow(A)), rep(0, ncol(A)-nrow(A)), c)    
  v <- as.data.frame(t(v))
  return(v)  
}

################################
# testing
################################
A <- matrix(c(1,2,0,5,0, 2,0,0,5,4, 0,2,1,0,0), nrow = 3)
P <- matrix(0, nrow = 3, ncol = 5)
b <- c(10, 2, 5)
sol_train <- c(1,0,0,1,1)

################################
# function definition
################################
predictor_vector <- function(P, A, b, sol_train) {
  n <- ncol(P)
  P <- rule1(P, A)
  res <- rule2a(P, A, sol_train)    
  res <- rule3(res$P, res$A)
  d <- rule4(res$P, res$A, b, n) 
  return(d)
}

rule1 <- function(P, A) {
  m <- nrow(P)  
  for (i in 1:m) {
    loc <- which(A[i,] != 0)
    P[i, loc] <- 1
  }
  return(P)
}

# only for training
rule2 <- function() {
  
}

# only for testing
rule2a <- function(P, A, sol_train) {
  m <- nrow(P)  
  n <- ncol(P)
  
  if (m < n) { # exclude some candidate variables randomly
    loc <- which(sol_train!=0)
    esc <- sample(loc, m)
    P[,-esc] <- 0
    A[,-esc] <- 0
  }

  # ensure each row and column has one unitary element
#   for (i in 1:n) {
#     loc <- which(P[,i] == 1)  
#     sig <- length(loc)
#     if (sig > 1) {
#       esc <- sample(loc, sig-1)
#       P[esc,i] <- 0
#     }
#   }

  return(list(P=P, A=A))
}

# eliminate empty columns in P and A
rule3 <- function(P, A) {
  loc <- which(colSums(A) != 0)
  P <- P[,loc]
  A <- A[,loc]
  return(list(P=P, A=A))
}

# construct d
rule4 <- function(P, A, b, n) {
  m <- nrow(P)    
  #dk <- solve(t(P)%*%A)%*%t(P)%*%b
  #dk <- ginv(t(P)%*%A)%*%t(P)%*%b
  dk <- solve(A) %*% b
  dnk <- rep(0, n-m)
  d <- c(dk, dnk)
  return(d)
}

# choose optimal d from a set of ds for LP_te 
rule5 <- function() {
  
}
