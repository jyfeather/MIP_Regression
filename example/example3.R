rm(list=ls())

library("Rglpk")
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

res <- Rglpk_solve_LP(obj = objective, mat = constraint_1, dir = constraint_2, 
                      rhs = constraint_3, bounds = rep("C",nrow), max = maximum)

loc <- which(res$solution != 0)

# find a basic feasible solution