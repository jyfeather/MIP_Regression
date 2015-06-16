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
add.constraint(lprec, c(1, 0), "<=", 1)
add.constraint(lprec, c(0, 1), "<=", 1)
#add.constraint(lprec, c(1, 1), ">", 100)
#add.constraint(lprec, c(0, 1), ">", 2)
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
ncol <- length(get.bounds(model)$lower) 
set.type(model, columns = 1:ncol, type = "real") # change MIP to LP 
add.constraint(model, c(1,2, rep(0, ncol-2)), "=", 4)

# reorganize the LP
nrow <- length(get.constr.type(model))
mat <- matrix(0, nrow = nrow+1, ncol = ncol)
for (i in 1:ncol) {
  tmp <- get.column(model, i)
  tmp.column <- tmp$column
  tmp.nzrow <- tmp$nzrow + 1
  mat[tmp.nzrow,i] <- tmp.column  
}
rhs <- get.rhs(model)
dir <- get.constr.type(model)

# w/o trick
model2 <- make.lp(0, ncol)
set.objfn(model2, mat[1,])
for(i in 1:nrow) {
  add.constraint(model2, mat[i+1,], "=", rhs[i])
}
add.constraint(model2, c(1,2, rep(0, ncol-2)), "=", 4) # to make it infeasible
set.bounds(model2, lower = rep(0,ncol), upper = rep(1,ncol), columns = 1:ncol)

beg <- proc.time()
solve(model2)
proc.time() - beg

# w/ trick
model3 <- make.lp(0, ncol + nrow + 1)
set.objfn(model3, c(mat[1,], rep(100000,nrow+1)))
for(i in 1:nrow) {
  tmp <- rep(0, nrow+1)
  tmp[i] <- 1
  lhs <- c(mat[i+1,], tmp)
  add.constraint(model3, lhs, "=", rhs[i])
}
add.constraint(model3, c(1,2, rep(0, ncol+nrow-2), 1), "=", 4) # to make it infeasible
set.bounds(model3, lower = rep(0,ncol), upper = rep(1,ncol), columns = 1:ncol)

beg <- proc.time()
solve(model3)
proc.time() - beg
