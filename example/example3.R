######################################################
#     To find feasible basis A_B
#     Using lpSolveAPI
#     get.basis()
#         return a vector: m(constraints, slack variables) + n(decision variables)
######################################################

library("lpSolveAPI")
library("Rglpk")

######################################################
#  air05 example
######################################################
rm(list=ls())
path <- paste(getwd(), "/example/Data/Sample/air05", sep = "")
lp <- read.lp(path, type = "mps")

# change MIP to LP, and solve
set.type(lp, columns = 1:n, type = "real")
solve(lp)
n <- length(get.type(lp))
m <- length(get.basis(lp))
B_loc <- abs(get.basis(lp))

# read constraint matrix via Rglpk package
model <- Rglpk_read_file(path, type = "MPS_fixed")
constraint <- as.matrix(model$constraints[[1]])
constraint <- cbind(diag(m), constraint)

# construct B matrix
B <- constraint[,B_loc]
rhs <- get.rhs(lp)
d_B <- round(solve(B) %*% rhs,2)

######################################################
#   standard example
######################################################
rm(list=ls())
lprec <- make.lp(0, 4)
set.objfn(lprec, c(1, 3, 6.24, 0.1))
add.constraint(lprec, c(0, 78.26, 0, 2.9), ">=", 92.3)
add.constraint(lprec, c(0.24, 0, 11.31, 0), "<=", 14.8)
add.constraint(lprec, c(12.68, 0, 0.08, 0.9), ">=", 4)
set.bounds(lprec, lower = c(28.6, 18), columns = c(1, 4))
set.bounds(lprec, upper = 48.98, columns = 4)
RowNames <- c("THISROW", "THATROW", "LASTROW")
ColNames <- c("COLONE", "COLTWO", "COLTHREE", "COLFOUR")
dimnames(lprec) <- list(RowNames, ColNames)
solve(lprec)
get.objective(lprec)
get.variables(lprec)
B_loc <- abs(get.basis(lprec, nonbasic = FALSE))
constraint <- matrix(c(1,0,0,0,1,0,0,0,1, 0,0.24,12.68, 78.26,0,0,0,11.31,0.08,2.9,0,0.9), nrow = 3)
B <- constraint[,B_loc]
rhs <- get.rhs(lprec)
d_B <- solve(B) %*% rhs

######################################################
#   small example
######################################################
rm(list=ls())
x <- make.lp(0, 2)
set.objfn(x, c(1,1))
add.constraint(x, c(2, 1), ">=", 5)
add.constraint(x, c(2, 2), "=", 6)
set.bounds(x, lower = c(0,0),columns = c(1,2))
set.bounds(x, upper = c(4,4),columns = c(1,2))
solve(x)
B_loc <- abs(get.basis(x))
constraint <- matrix(c(1,0,0,1,2,2,1,1), nrow = 2)
rhs <- get.rhs(x)
d_B <- solve(constraint[,B_loc]) %*% rhs
d_B[B_loc > 2]
get.variables(x)