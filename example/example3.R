######################################################
#     To find feasible basis A_B
#     refer to http://cran.r-project.org/web/views/Optimization.html
#
#     Using lpSolveAPI
#     get.basis()
#         return a vector: m(constraints, slack variables) + n(decision variables)
#
#     Using linprog
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
res_lpSolveAPI[res_lpSolveAPI!=0]

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
rhs <- get.rhs(lprec)
d_B <- solve(B) %*% rhs

######################################################
#   small example
######################################################
rm(list=ls())
x <- make.lp(0, 2)
set.objfn(x, c(1,1))
add.constraint(x, c(2, 1), "=", 5)
add.constraint(x, c(2, 2), "=", 6)
set.bounds(x, lower = c(0,0),columns = c(1,2))
#set.bounds(x, upper = c(4,4),columns = c(1,2))
solve(x)
B_loc <- abs(get.basis(x))
constraint <- matrix(c(1,0,0,1,2,2,1,2), nrow = 2)
rhs <- get.rhs(x)
d_B <- solve(constraint[,B_loc]) %*% rhs
d_B[B_loc > 2]
get.variables(x)

#####################################################
#  Package 'linprog' exmaple  
#####################################################
cvec <- as.vector(model$objective)
bvec <- as.vector(model$constraints[[3]])
bvec <- c(bvec, rep(1, 33))
Amat <- as.matrix(model$constraints[[1]])
Amat <- rbind(Amat, diag(33))
res_linprog <- solveLP( cvec, bvec, Amat, maximum = FALSE)

Amat <- cbind(diag(49), Amat)
basis <- rownames(res_linprog$basvar)
basis_var <- basis[grep("^[0-9].*", basis)]; basis_var <- as.numeric(basis_var)
basis_slack <- basis[grep("S.*", basis)]; basis_slack <- as.numeric(substring(basis_slack,3)) 
A_b <- Amat[,c(basis_slack, basis_var+49)]
d_b <- round(solve(A_b) %*% bvec, 2)
d_b_solved <- res_linprog$basvar


#####################################################
#####################################################
rm(list=ls())
# preparing the model
lp <- initProbCLP()
nrows <- 5
ncols <- 8
# objective function
obj <- c(1, 0, 0, 0, 2, 0, 0, -1)
# upper and lower bounds of the rows
rlower <- c(2.5, -1000, 4, 1.8, 3)
rupper <- c(1000, 2.1, 4, 5, 15)
# upper and lower bounds of the columns
clower <- c(2.5, 0, 0, 0, 0.5, 0, 0, 0)
cupper <- c(1000, 4.1, 1, 1, 4, 1000, 1000, 4.3)
# constraint matrix
ia <- c(0, 4, 0, 1, 1, 2, 0, 3, 0, 4, 2, 3, 0, 4)
ja <- c(0, 2, 4, 6, 8, 10, 11, 12, 14)
ar <- c(3.0, 5.6, 1.0, 2.0, 1.1, 1.0, -2.0, 2.8,
        -1.0, 1.0, 1.0, -1.2, -1.0, 1.9)
# direction of optimization
setObjDirCLP(lp, 1)
# load problem data
loadProblemCLP(lp, ncols, nrows, ia, ja, ar,
               clower, cupper, obj, rlower, rupper)
# solve lp problem
solveInitialCLP(lp)
# retrieve the results
getSolStatusCLP(lp)
getObjValCLP(lp)
getColPrimCLP(lp)
getRowPrimCLP(lp)