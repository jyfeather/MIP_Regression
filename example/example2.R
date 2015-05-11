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

source("./R/PredictionVector.R")
source("./R/BAB.R")
source("./R/EpsilonSimilar.R")
source("./R/MIPSolution.R")
source("./R/PredictionModel.R")
####################################
#  input
####################################
path <- paste(getwd(), "/example/Data/Sample/air05", sep = "")
#path <- paste(getwd(), "/example/Data/Sample/p0033.mps", sep = "")
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
  vs_train <- matrix(NA, ncol = 2*ncol)
  sols <- matrix(NA, ncol = ncol)
  beg_t <- proc.time()
  for (i in 1:numtree) {
    newobjs <- similar_obj(objective)
    newcons <- similar_cons(constraint_3)
    res <- Rglpk_solve_LP(obj = newobjs, mat = constraint_1, dir = constraint_2, rhs = newcons, 
                          bounds = rep("C",nrow), max = maximum)  
    v <- c(v_train(res), as.vector(newobjs))
    sol <- res$solution
    vs_train <- rbind(vs_train, v)
    sols <- rbind(sols, sol)
    cat(i)
  }
  exc_t <- proc.time() - beg_t
  vs_train <- na.omit(vs_train)
  sols <- na.omit(sols)
  write.csv(sols, file = "./example/temp/train_sols.csv")
  write.csv(vs_train, file = "./example/temp/train_vector.csv")
  rm(numtree, beg_t, newobjs, newcons, res, v, sol, i)
} else {
  sols <- read.csv(file = "./example/temp/train_sols.csv", header = T)  
  vs_train <- read.csv(file = "./example/temp/train_vector.csv", header = T)  
  sols <- sols[,-1]
  vs_train <- vs_train[,-1]
}

####################################
#   training process
####################################
# ncol prediction models, time consuming, 0.5 seconds per prediction
# model_list <- list()
# for (i in 1:ncol) {
#   reg <- cbind(vs, sols[,i])
#   reg <- as.data.frame(reg)
#   names(reg)[ncol(reg)] <- "y"
#   model <- TrainModel(reg, "gbm")  #  Error: protect(): protection stack overflow 
#   #model <- TrainModel(reg, "linear")
#   #model <- TrainModel(as.matrix(reg), "xgboost")
#   model_list[[i]] <- model
# }

####################################
#   testing process
####################################
# test solver
num = 100
vs_test <- matrix(NA, ncol = 2*ncol)
beg_t <- proc.time()
while(num > 0) {
  testobj <- similar_obj(objective) 
  testcon <- similar_cons(constraint_3)

  # solution from LP solver
  sol_solver <- Rglpk_solve_LP(obj = testobj, mat = constraint_1, dir = constraint_2, rhs = testcon, 
                               bounds = rep("C",nrow), max = maximum)
  vs_test <- rbind(vs_test, c(sol_solver$sol, as.vector(testobj)))
  num <- num - 1
}
exc_t_solver <- proc.time() - beg_t
vs_test <- na.omit(vs_test)

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
  loc <- which(distance %in% sort(distance)[1:3])[1:3]
  sol <- rbind(sol, colSums(vs_train[loc,1:ncol])/3)
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
hist(error)

cor <- c()
for (i in 1:num) {
  cor <- c(cor, cor(sol_prediction[i,], sol_solver[i,]))  
}
plot(cor)
