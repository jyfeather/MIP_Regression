############################################
# wonder the difference when facing an infeasible LP
# http://lpsolve.sourceforge.net/5.5/
############################################

library(lpSolveAPI)

############################################
# simple problem
############################################
# w/o trick
lprec <- make.lp(0, 2)
set.objfn(lprec, c(1, 1))
add.constraint(lprec, c(1, 0), ">=", 6)
add.constraint(lprec, c(0, 1), ">=", 6)
add.constraint(lprec, c(1, 1), "<=", 11)
beg <- proc.time()
n <- 10000
while(n > 0) {
  solve(lprec)
  n <- n - 1
}
proc.time() - beg

# w/ trick
lprec <- make.lp(0, 5)
set.objfn(lprec, c(1, 1, 1000, 1000, 1000))
add.constraint(lprec, c(1, 0, 1, 0, 0), ">=", 6)
add.constraint(lprec, c(0, 1, 0, 1, 0), ">=", 6)
add.constraint(lprec, c(1, 1, 0, 0, -1), "<=", 11)
beg <- proc.time()
n <- 10000
while(n > 0) {
  solve(lprec)
  n <- n - 1
}
proc.time() - beg

###########################################
# complex problem
###########################################
path <- paste(getwd(), "/example/Data/Sample/air05", sep = "")
model <- read.lp(filename = path, type = "mps")
set.type(model, columns = 1:ncol, type = "real") # change MIP to LP 

# w/o trick
beg <- proc.time()
solve(model)
proc.time() - beg

# w/ trick

