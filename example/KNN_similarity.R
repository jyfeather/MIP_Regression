########################################
#   This example is for Guan's idea,
#   which adopts epsilon similar definition
#   as criteria to measure the similarity
########################################
rm(list=ls())

library("Rglpk")
library("lpSolveAPI")

source("./R/EpsilonSimilar.R")

########################################
# function defined
########################################
findSimilar <- function(d, d_set) {
  n <- nrow(d_set)
  m <- ncol(d_set)
  dist <- c()
  for(i in 1:n) {
    dist <- c(dist, sum((d-d_set[i,])^2)/(sum(d^2)+sum(d_set[i,])^2))
  }
  loc <- which(dist %in% sort(dist)[1:5])[1:5]
  return(loc)
}

########################################
#   Input
########################################
path <- paste(getwd(), "/example/Data/Sample/air05", sep = "")
model <- read.lp(filename = path, type = "mps")
ncol <- length(get.type(model))
set.bounds(model, upper = rep(Inf, ncol), columns = 1:ncol) # remove upper bound, currently limited
set.type(model, columns = 1:ncol, type = "real") # change MIP to LP 
nrow <- length(get.rhs(model))

# get cvec, bvec
model_Rglpk <- Rglpk_read_file(path, type = "MPS_fixed")
cvec <- model_Rglpk$objective
bvec <- model_Rglpk$constraints[[3]]
dir <- as.vector(model_Rglpk$constraints[[2]])

########################################
#   solve training LPs
########################################
numtree <- 300
d_train <- matrix(NA, ncol = (ncol+nrow)) # store vector d = [c, b]
sols <- matrix(NA, ncol = ncol)           # store solution
for (i in 1:numtree) {
  ## generate new LP based on origin
  newobjs <- as.vector(similar_obj(cvec))
  newcons <- as.vector(similar_cons(bvec))
  
  ## solve LP 
  set.objfn(model, obj = newobjs)
  set.rhs(model, b = newcons)
  solve(model)
  
  ## keep useful information
  d_train <- rbind(d_train, c(newobjs, newcons))
  sol <- get.variables(model)
  sols <- rbind(sols, sol)
  cat(i)
}
d_train <- na.omit(d_train)
sols <- na.omit(sols)

#######################################
#  construct testing set and solve
#######################################
numtree <- 100
test_objs <- matrix(NA, ncol = ncol)
test_cons <- matrix(NA, ncol = nrow)
sol_solv <- matrix(NA, ncol = ncol)
sol_pred <- matrix(NA, ncol = ncol) 
for (i in 1:numtree) {
  ## generate new LP based on origin
  newobjs <- as.vector(similar_obj(cvec))
  newcons <- as.vector(similar_cons(bvec))
  
  ## keep useful information
  test_objs <- rbind(newobjs, test_objs)
  test_cons <- rbind(newcons, test_cons)
  
  ## solve via sovler
  set.objfn(model, obj = newobjs)
  set.rhs(model, b = newcons)
  solve(model)
  sol_solv <- rbind(get.variables(model), sol_solv)  
  
  ## solve via KNN
  d <- c(newobjs, newcons)
  loc <- findSimilar(d, d_train)    
  sol_pred <- rbind(sol_pred, colSums(sols[loc,])/5)
  
  cat(i)
}
test_objs <- na.omit(test_objs)
test_cons <- na.omit(test_cons)
sol_solv <- na.omit(sol_solv)
sol_pred <- na.omit(sol_pred)

###########################################
#  prediction precision
###########################################
# compare the objective value
obj_solv <- c()
obj_pred <- c()
for (i in 1:numtree) {
  c <- test_objs[i,]
  obj_pred <- c(obj_pred, c %*% sol_pred[i,]) 
  obj_solv <- c(obj_solv, c %*% sol_solv[i,])
}
error <- round((obj_pred-obj_solv)/obj_solv,2)
hist(error, main = "Objective Value Error of 5-NN in LP problem")
