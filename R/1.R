library("Rglpk")
library("slam")
library("gbm")
####################################
#  input
####################################
#path <- paste(getwd(), "/example/air05.mps", sep = "")
path <- paste(getwd(), "/example/mas74.mps", sep = "")
model <- Rglpk_read_file(path, type = "MPS_fixed")
# too time consuming
#res <- Rglpk_solve_LP(model$objective, model$constraints[[1]], model$constraints[[2]],
#                      model$constraints[[3]], model$bounds, model$types, model$maximum)
objective <- model$objective
constraint <- model$constraints
constraint_1 <- constraint[[1]]
constraint_2 <- constraint[[2]]
constraint_3 <- constraint[[3]]
bounds <- model$bounds
types <- model$types
maximum <- model$maximum
nrow <- constraint_1$nrow 
ncol <- constraint_1$ncol 

####################################
#   training process
####################################
numtree <- 500
numnode <- 0 
#reg <- array(NA, c(numtree, ncol*2, ncol)) too large
vs <- matrix(NA, ncol = ncol)
sols <- matrix(NA, ncol = ncol)
# time consuming, about 2~3 mins for 100 trees
for (i in 1:numtree) {
  newobjs <- similar_obj(objective)
  newcons <- similar_cons(constraint_3)
  res <- sub_prob(numnode, newobjs, constraint_1, constraint_2, newcons, bounds, maximum, nrow, ncol)  
  v <- v_train(res)
  sol <- res$sol
  vs <- rbind(vs, v)
  sols <- rbind(sols, sol)
}
vs <- na.omit(vs)
sols <- na.omit(sols)

# ncol prediction models, time consuming, 0.5 seconds per prediction
model_list <- list()
for (i in 1:ncol) {
  reg <- cbind(sols, vs, sols[,i])
  reg <- as.data.frame(reg)
  names(reg)[ncol(reg)] <- "y"
  model <- gbm(y~., data = reg, n.trees = 500, distribution = "gaussian")  
  model_list[[i]] <- model
}

####################################
#   testing process
####################################
testobj <- similar_obj(objective)
testcon <- similar_cons(constraint_3)
v <- v_test(constraint_1, testcon, testobj)
# solution from LP solver
sol_solver <- sub_prob(numnode, testobj, constraint_1, constraint_2, testcon, bounds, maximum, nrow, ncol)
sol_solver <- sol_solver$sol
# solution from regression
names(v) <- names(reg)[-303]
sol <- c()
for (i in 1:ncol) {
  sol <- c(sol, predict.gbm(model_list[[i]], v, n.trees = 50))
}

####################################
#    function
####################################
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

#objective <- similar_prob(objective)

# solve sub problems
sub_prob <- function(numnode, obj, con1, con2, con3, bound, max, nrow, ncol) {
  i <- con1$i
  j <- con1$j
  v <- con1$v
  if (numnode > 0) {
    # branch procedure
    for (k in 1:numnode) {
      i <- c(i, nrow+k)
      j <- c(j, k)
    }
    v <- c(v, rep(1, numnode))
    con1 <- simple_triplet_matrix(i, j, v, nrow = max(i), ncol = max(j), dimnames = NULL)
    con2 <- c(con2, rep("==", numnode))
    con3 <- c(con3, round(runif(numnode)))
  }
  # change to LP problems from IP/BP 
  type <- rep("C", max(j))
  res <- Rglpk_solve_LP(obj, con1, con2, con3, bound, type, max)
  # return A,b,c,x*  
  return(list(A=con1, b=con3, c=obj, sol=res$sol))
}

#tmp <- sub_prob(4, objective, constraint_1, constraint_2, constraint_3, bounds, maximum, nrow, ncol)

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