####################################
# In this example, purpose is to test 
# speed and accuracy of prediction
# method.
# 
# prediction method: nearest neighbour
####################################
rm(list=ls())

library("Rglpk")
library("lpSolveAPI")

source("./R/PredictionVector.R")
source("./R/BAB.R")
source("./R/EpsilonSimilar.R")
source("./R/MIPSolution.R")
source("./R/PredictionModel.R")
####################################
#  input
####################################
path <- paste(getwd(), "/example/Data/Sample/air05", sep = "")
model <- read.lp(filename = path, type = "mps")
ncol <- length(get.type(model))
set.bounds(model, upper = rep(Inf, ncol), columns = 1:ncol) # remove upper bound, currently limited
set.type(model, columns = 1:ncol, type = "real") # change MIP to LP 
nrow <- length(get.rhs(model))

# get A matrix
model_Rglpk <- Rglpk_read_file(path, type = "MPS_fixed")
Amat <- as.matrix(model_Rglpk$constraints[[1]])
Amat <- cbind(diag(nrow), Amat)

# get cvec, bvec
cvec <- model_Rglpk$objective
bvec <- model_Rglpk$constraints[[3]]
dir <- as.vector(model_Rglpk$constraints[[2]])

####################################
#  solving training LPs
#  time consuming work
####################################
first = TRUE
beg_t <- proc.time()
if (first) {
  numtree <- 300
  vs_train <- matrix(NA, ncol = 2*ncol)
  sols <- matrix(NA, ncol = ncol)
  Basis <- matrix(NA, ncol = nrow)
  beg_t <- proc.time()
  for (i in 1:numtree) {
    newobjs <- as.vector(similar_obj(cvec))
    newcons <- as.vector(similar_cons(bvec))
    ## solve LP and generate A_b
    set.objfn(model, obj = newobjs)
    set.rhs(model, b = newcons)
    solve(model)
    B_loc <- get.basis(model)
    Basis <- rbind(Basis, B_loc)  
    sol <- get.variables(model)
    v <- c(sol, newobjs)
    vs_train <- rbind(vs_train, v)
    sols <- rbind(sols, sol)
    cat(i)
  }
  exc_t <- proc.time() - beg_t
  vs_train <- na.omit(vs_train)
  sols <- na.omit(sols)
  Basis <- na.omit(Basis)
  exc_t <- proc.time() - beg_t
  write.csv(sols, file = "./example/temp/train_sols.csv")
  write.csv(vs_train, file = "./example/temp/train_vector.csv")
  write.csv(Basis, file = "./example/temp/basis.csv")
  rm(numtree, beg_t, newobjs, newcons, v, sol, i)
} else {
  sols <- read.csv(file = "./example/temp/train_sols.csv", header = T)  
  vs_train <- read.csv(file = "./example/temp/train_vector.csv", header = T)  
  sols <- sols[,-1]
  vs_train <- vs_train[,-1]
}

####################################
#   testing process
####################################
# test solver
num = 100
vs_test <- matrix(NA, ncol = 2*ncol)
vs_pred <- matrix(NA, ncol = 2*ncol)
beg_t <- proc.time()
while(num > 0) {
  testobj <- as.vector(similar_obj(cvec))
  testcon <- as.vector(similar_cons(bvec))

  # solution from LP solver
  set.objfn(model, obj = testobj)
  set.rhs(model, b = testcon)
  solve(model)
  
  vs_test <- rbind(vs_test, c(get.variables(model), testobj))
  
  # construct test vs
  B_loc <- abs(Basis[sample(300, 1),])
  B <- Amat[, B_loc]
  d_B <- round(solve(B) %*% testcon, 2)
  vs <- rep(0, nrow + ncol)
  vs[B_loc] <- d_B
  vs_pred <- rbind(vs_pred, c(vs[-c(1:nrow)], testobj))

  num <- num - 1
  cat(num)
}
exc_t_solver <- proc.time() - beg_t
vs_test <- na.omit(vs_test)

# test prediction
num = 100

# solution from regression
beg_t <- proc.time()
# sol <- matrix(0, nrow = num, ncol = ncol)
# for (i in 1:ncol) {
#   sol[,i] <- TestPrediction(model_list[[i]], vs, method = "xgboost")
# }

# 3 nearest neighbour
vs_train <- as.matrix(vs_train)
sol <- matrix(NA, ncol = ncol)
for (i in 1:num) {
  diff <- sweep(vs_train, 2, vs_test[i,])
  distance <- sqrt(rowSums(diff^2))
  loc <- which(distance %in% sort(distance)[1:5])[1:5]
  sol <- rbind(sol, colSums(vs_train[loc,1:ncol])/5)
}
sol <- na.omit(sol)
exc_t_prediction <- proc.time() - beg_t

######################################
#   precision
######################################
# compare the objective value
sol_solver <- as.matrix(vs_test[,1:ncol])
sol_prediction <- sol
C <- as.matrix(vs_test[,-c(1:ncol)])
obj_solver <- rowSums(sol_solver * C)
obj_prediction <- rowSums(sol_prediction * C)
error <- round((obj_prediction-obj_solver)/obj_solver,2)
hist(error, main = "Objective Value Error of 5-NN in LP problem")

# compare MIP
BP_sol <- read.table(file = "./example/Data/Sample/air05.sol")
BP_sol <- BP_sol[,2]
BP_sol_pred <- round(sol_prediction)
BP_obj <- rowSums(BP_sol * C)
BP_obj_pred <- rowSums(BP_sol_pred * C)
error <- round((BP_obj_pred-BP_obj)/BP_obj,2)
hist(error, main = "Objective Value Error of 5-NN in BP problem")

# compare correlation
cor <- c()
for (i in 1:num) {
  cor <- c(cor, cor(sol_prediction[i,], sol_solver[i,]))  
}
plot(cor, main = "Solution Correlation Between Solver and Prediction")
