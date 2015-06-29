rm(list=ls())

library(lpSolveAPI)
path <- paste(getwd(), "/example/Data/Sample/air05", sep = "")
lp <- read.lp(path, type = "mps")

n <- length(get.type(lp))
set.bounds(lp, upper = rep(Inf, n), columns = 1:n)
set.type(lp, columns = 1:n, type = "real")
t.start <- proc.time()
solve(lp)
t.end <- proc.time()-t.start


library(quadprog)
# reorganize the LP
nrow <- length(get.constr.type(lp))
ncol <- length(get.bounds(lp)$lower) 
Amat <- matrix(0, nrow = nrow+1, ncol = ncol)
for (i in 1:ncol) {
  tmp <- get.column(lp, i)
  tmp.column <- tmp$column
  tmp.nzrow <- tmp$nzrow + 1
  Amat[tmp.nzrow,i] <- tmp.column  
}
Amat <- Amat[-1,]
bvec <- get.rhs(lp)
dir <- get.constr.type(lp)
Dmat <- diag(ncol)
dvec <- rnorm(ncol)
solve.QP(Dmat, dvec, t(Amat), bvec, meq = 0)
