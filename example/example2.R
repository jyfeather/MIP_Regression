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
#path <- paste(getwd(), "/example/Data/Sample/air04", sep = "")
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
first = FALSE 
numnode <- 2
if (first) {
  numtree <- 300
  #reg <- array(NA, c(numtree, ncol*2, ncol)) too large
  vs <- matrix(NA, ncol = 2*ncol)
  sols <- matrix(NA, ncol = ncol)
  beg_t <- proc.time()
  for (i in 1:numtree) {
    newobjs <- similar_obj(objective)
    newcons <- similar_cons(constraint_3)
    for (j in 1:numnode) {
      res <- sub_prob(newobjs, constraint_1, constraint_2, newcons, bounds, maximum, nrow, ncol)  
      v <- v_train(res)
      sol <- res$sol
      vs <- rbind(vs, v)
      sols <- rbind(sols, sol)
    }
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
  model <- TrainModel(reg, "gbm")  #  Error: protect(): protection stack overflow 
  #model <- TrainModel(reg, "linear")
  #model <- TrainModel(as.matrix(reg), "xgboost")
  model_list[[i]] <- model
}

####################################
#   testing process
####################################
# test solver
num = 10000
beg_t <- proc.time()
while(num > 0) {
  testobj <- similar_obj(objective) 
  testcon <- similar_cons(constraint_3)

  # solution from LP solver
  sol_solver <- sub_prob(numnode, testobj, constraint_1, constraint_2, testcon, bounds, maximum, nrow, ncol)
  sol_solver <- sol_solver$sol  
  num <- num - 1
}
exc_t_solver <- proc.time() - beg_t

# test prediction
num = 10000
vs <- matrix(NA, ncol = 2*ncol)
while (num > 0) {
  testobj <- similar_obj(objective) 
  testcon <- similar_cons(constraint_3)
  #d <- predictor_vector(matrix(0, nrow = nrow, ncol = ncol), as.matrix(constraint_1), constraint_3, sols[1,])
  v <- v_test(constraint_1, testcon, testobj)
  vs <- rbind(vs, v)
  num <- num - 1
}
vs <- na.omit(vs)

# solution from regression
names(vs) <- names(reg)[-length(reg)]

beg_t <- proc.time()
sol <- c()
for (i in 1:ncol) {
  sol <- c(sol, TestPrediction(model_list[[i]], vs, method = "gbm"))
}
exc_t_prediction <- proc.time() - beg_t