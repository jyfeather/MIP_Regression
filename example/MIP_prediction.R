rm(list=ls())

library("Rglpk")
library("slam")
library("gbm")
library("MASS")

source("./R/PredictionVector.R")
source("./R/BAB.R")
source("./R/EpsilonSimilar.R")
source("./R/MIPSolution.R")
####################################
#  input
####################################
#path <- paste(getwd(), "/example/Data/Sample/mas74", sep = "")
path <- paste(getwd(), "/example/Data/Sample/pack1.mps", sep = "")
#path <- paste(getwd(), "/example/Data/Sample/air04", sep = "")
model <- Rglpk_read_file(path, type = "MPS_fixed")
if (FALSE) {
  # solve MIP 
  # too time consuming
  ###################
  # for mas74.mps
  # R:               more than 20 mins unsolved
  # coin-or cbc:     more than 40 mins unsolved
  # cplex:          
  # boosting-tree:   less than 20 seconds
  ###################
  # for p0033.mps
  # R:               ~1.31s
  # coin-or cbc:     ~0.16s
  # boosting-tree:   ~2s
  beg_t <- proc.time()
  res_solver <- Rglpk_solve_LP(model$objective, model$constraints[[1]], model$constraints[[2]],
                        model$constraints[[3]], model$bounds, model$types, model$maximum)
  exc_t <- proc.time() - beg_t
}
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
numnode <- 2
#reg <- array(NA, c(numtree, ncol*2, ncol)) too large
vs <- matrix(NA, ncol = ncol)
sols <- matrix(NA, ncol = ncol)
# time consuming, about 2~3 mins for 100 trees
for (i in 1:numtree) {
  newobjs <- similar_obj(objective)
  newcons <- similar_cons(constraint_3)
  for (j in 1:numnode) {
    res <- sub_prob(numnode, newobjs, constraint_1, constraint_2, newcons, bounds, maximum, nrow, ncol)  
    v <- v_train(res)
    sol <- res$sol
    vs <- rbind(vs, v)
    sols <- rbind(sols, sol)
  }
}
vs <- na.omit(vs)
sols <- na.omit(sols)

# ncol prediction models, time consuming, 0.5 seconds per prediction
model_list <- list()
for (i in 1:ncol) {
  reg <- cbind(sols, vs, sols[,i])
  reg <- as.data.frame(reg)
  names(reg)[ncol(reg)] <- "y"
  model <- gbm(y~., data = reg, 
               n.trees = 50, distribution = "gaussian",
               verbose = FALSE)  
  model_list[[i]] <- model
}

####################################
#   testing process
####################################
testobj <- similar_obj(objective) 
testcon <- similar_cons(constraint_3)

d <- predictor_vector(matrix(0, nrow = nrow, ncol = ncol), as.matrix(constraint_1), constraint_3, sol)
#v <- v_test(constraint_1, testcon, testobj)

# solution from LP solver
sol_solver <- sub_prob(numnode, testobj, constraint_1, constraint_2, testcon, bounds, maximum, nrow, ncol)
sol_solver <- sol_solver$sol
# solution from regression
names(v) <- names(reg)[-length(reg)]
sol <- c()
for (i in 1:ncol) {
  sol <- c(sol, predict.gbm(model_list[[i]], v, n.trees = 50))
}
sol_regr <- MIPround(sol)
opt_regr <- as.vector(objective) %*% sol_regr
res_solver$optimum
res_solver$solution
