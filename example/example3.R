######################################################
#     To find feasible basis A_B
#     Using lpSolveAPI
#     get.basis()
#         return a vector: m(constraints, slack variables) + n(decision variables)
######################################################
rm(list=ls())

library("lpSolveAPI")
path <- paste(getwd(), "/example/Data/Sample/air05", sep = "")
lp <- read.lp(path, type = "mps")

# change MIP to LP, and solve
n <- length(get.type(lp))
m <- length(get.basis(lp))
set.type(lp, columns = 1:n, type = "real")
solve(lp)
B_loc <- get.basis(lp)
rhs <- get.rhs(lp)
B <- matrix(NA, nrow = m)
for (i in B_loc) {
  if (abs(i) <= m) {
    B_i <- rep(0, m); B_i[abs(i)] = 1
    B <- cbind(B, B_i)
  } else {
    B <- cbind(B, get.column(lp, abs(i)-m)$column[-1])
  }
}

# find a basic feasible solution
constraint_1 <- as.matrix(constraint_1)
B <- constraint_1[,sample(loc, nrow)]

# lpsolve
lprec <- make.lp(0, 4)
set.objfn(lprec, c(1, 3, 6.24, 0.1))
add.constraint(lprec, c(0, 78.26, 0, 2.9), ">=", 92.3)
add.constraint(lprec, c(0.24, 0, 11.31, 0), "<=", 14.8)
add.constraint(lprec, c(12.68, 0, 0.08, 0.9), ">=", 4)
set.bounds(lprec, lower = c(28.6, 18), columns = c(1, 4))
set.bounds(lprec, upper = 48.98, columns = 4)
solve(lprec)
get.objective(lprec)
get.variables(lprec)
get.basis(lprec, nonbasic = FALSE)

x <- make.lp(0, 2)
set.objfn(x, c(1,1))
add.constraint(x, c(2, 1), ">=", 5)
add.constraint(x, c(2, 2), "=", 6)
set.bounds(x, lower = c(0,0),columns = c(1,2))
set.bounds(x, upper = c(4,4),columns = c(1,2))
solve(x)
get.basis(x)
solve(matrix(c(2,2,1,0),nrow = 2)) %*% c(5,6)