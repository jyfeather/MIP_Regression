####################################
# In this example, purpose is to test 
# speed and accuracy of prediction
# method.
####################################
rm(list=ls())

library("Rglpk")
library("slam")
library("gbm")
library("MASS")
library("xgboost")
library("ggplot2")

source("./R/PredictionVector.R")
source("./R/BAB.R")
source("./R/EpsilonSimilar.R")
source("./R/MIPSolution.R")
source("./R/PredictionModel.R")
####################################
#  input
####################################
#path <- paste(getwd(), "/example/Data/Sample/air05", sep = "")
path <- paste(getwd(), "/example/Data/Sample/p0033.mps", sep = "")
model <- Rglpk_read_file(path, type = "MPS_fixed")

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
#  solving training LPs
#  time consuming work
####################################
first = TRUE
if (first) {
  numtree <- 100
  #reg <- array(NA, c(numtree, ncol*2, ncol)) too large
  vs <- matrix(NA, ncol = 2*ncol)
  sols <- matrix(NA, ncol = ncol)
  beg_t <- proc.time()
  for (i in 1:numtree) {
    newobjs <- similar_obj(objective)
    newcons <- similar_cons(constraint_3)
    res <- Rglpk_solve_LP(obj = newobjs, mat = constraint_1, dir = constraint_2, rhs = newcons, 
                          bounds = rep("C",nrow), max = maximum)  
    v <- c(v_train(res), as.vector(newobjs))
    sol <- res$solution
    vs <- rbind(vs, v)
    sols <- rbind(sols, sol)
    cat(i)
  }
  exc_t <- proc.time() - beg_t
  vs <- na.omit(vs)
  sols <- na.omit(sols)
  write.csv(sols, file = "./example/temp/train_sols.csv")
  write.csv(vs, file = "./example/temp/train_vector.csv")
  rm(numtree, beg_t, newobjs, newcons, res, v, sol, i)
} else {
  sols <- read.csv(file = "./example/temp/train_sols.csv", header = T)  
  vs <- read.csv(file = "./example/temp/train_vector.csv", header = T)  
  sols <- sols[,-1]
  vs <- vs[,-1]
}

####################################
#   training process
####################################
# ncol prediction models, time consuming, 0.5 seconds per prediction
model_list <- list()
for (i in 1:ncol) {
  reg <- cbind(vs, sols[,i])
  reg <- as.data.frame(reg)
  names(reg)[ncol(reg)] <- "y"
  #model <- TrainModel(reg, "gbm")  #  Error: protect(): protection stack overflow 
  #model <- TrainModel(reg, "linear")
  model <- TrainModel(as.matrix(reg), "xgboost")
  model_list[[i]] <- model
}

####################################
#   testing process
####################################
# test solver
num = 100
vs <- matrix(NA, ncol = 2*ncol)
beg_t <- proc.time()
while(num > 0) {
  testobj <- similar_obj(objective) 
  testcon <- similar_cons(constraint_3)

  # solution from LP solver
  sol_solver <- Rglpk_solve_LP(obj = testobj, mat = constraint_1, dir = constraint_2, rhs = testcon, 
                               bounds = rep("C",nrow), max = maximum)
  vs <- rbind(vs, c(sol_solver$sol, as.vector(testobj)))
  num <- num - 1
}
exc_t_solver <- proc.time() - beg_t
vs <- na.omit(vs)

# test prediction
num = 100
# vs <- matrix(NA, ncol = 2*ncol)
# while (num > 0) {
#   testobj <- similar_obj(objective) 
#   testcon <- similar_cons(constraint_3)
#   #d <- predictor_vector(matrix(0, nrow = nrow, ncol = ncol), as.matrix(constraint_1), constraint_3, sols[1,])
#   v <- v_test(constraint_1, testcon, testobj)
#   vs <- rbind(vs, v)
#   num <- num - 1
# }
# vs <- na.omit(vs)

# solution from regression
vs <- as.data.frame(vs)
names(vs) <- names(reg)[-length(reg)]

beg_t <- proc.time()
# sol <- matrix(0, nrow = num, ncol = ncol)
# for (i in 1:ncol) {
#   sol[,i] <- TestPrediction(model_list[[i]], vs, method = "xgboost")
# }

# 3 nearest neighbour
train_vs <- as.matrix(reg[,-ncol(reg)])
test_vs <- as.matrix(vs)
sol <- matrix(NA, ncol = ncol)
for (i in 1:num) {
  diff <- sweep(train_vs, 2, test_vs[i,])
  distance <- sqrt(rowSums(diff^2))
  sol <- rbind(sol, reg[which.min(distance),1:ncol])
}
sol <- na.omit(sol)
exc_t_prediction <- proc.time() - beg_t

######################################
#   precision
######################################
# compare the objective value
sol_solver <- as.matrix(vs[,1:ncol])
sol_prediction <- sol
C <- as.matrix(vs[,-c(1:ncol)])
obj_solver <- rowSums(sol_solver * C)
obj_prediction <- rowSums(sol_prediction * C)
error <- round((obj_prediction-obj_solver)/obj_solver,2)
hist(error)
