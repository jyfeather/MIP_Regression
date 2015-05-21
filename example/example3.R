######################################################
#     To find feasible basis A_B
#     refer to http://cran.r-project.org/web/views/Optimization.html
#
#     Using lpSolveAPI
#     get.basis()
#         return a vector: m(constraints, slack variables) + n(decision variables)
######################################################
library("lpSolveAPI")
library("Rglpk")
library("linprog")
library("clpAPI")

######################################################
#  air05 example
######################################################
rm(list=ls())
path <- paste(getwd(), "/example/Data/Sample/air05", sep = "")
#path <- paste(getwd(), "/example/Data/Sample/p0033.mps", sep = "")
lp <- read.lp(path, type = "mps")
#write.lp(lp, filename = "./example/Data/Sample/p0033.lp", type = "lp")

# change MIP to LP, and solve
n <- length(get.type(lp))
set.bounds(lp, upper = rep(Inf, n), columns = 1:n)
set.type(lp, columns = 1:n, type = "real")
solve(lp)
m <- length(get.basis(lp))
B_loc <- abs(get.basis(lp))
res_lpSolveAPI <- get.variables(lp)

# read constraint matrix via Rglpk package
model <- Rglpk_read_file(path, type = "MPS_fixed")
constraint <- as.matrix(model$constraints[[1]])
constraint <- cbind(diag(m), constraint)

# construct B matrix
B <- constraint[,B_loc]
rhs <- get.rhs(lp)
d_B <- round(solve(B) %*% rhs,2)
B_loc <- sort(B_loc)
d_B[B_loc > m]

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
NB <- constraint[, -B_loc]
rhs <- get.rhs(lprec)
d_NB <- c(0,28.6,0,0)
d_B <- solve(B) %*% (rhs - NB %*% d_NB)

###########################################
# https://homepages.rpi.edu/~mitchj/handouts/upperbounds/
###########################################
lp <- make.lp(0,4)
set.objfn(lp, c(0,0,-2,-3))
add.constraint(lp, c(1,0,1,1), "=", 4)
add.constraint(lp, c(0,1,-1,-2), "=", 1)
set.bounds(lp, upper = c(10,4,5,1), columns = 1:4)
solve(lp)
get.variables(lp)
B_loc <- abs(get.basis(lp))
constraint <- matrix(c(1,0, 0,1, 1,0, 0,1, 1,-1, 1,-2), nrow = 2)
B <- constraint[,B_loc]
NB <- constraint[,-B_loc]
d_NB <- c(0,0,4,0)
rhs <- c(4,1)
d_B <- solve(B) %*% (rhs - NB %*% d_NB)
