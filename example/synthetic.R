################################################
# LP & QP comparison
################################################

################################################
# generate synthetic problem w/ box constraints
################################################
rm(list=ls())
library(Matrix)
library(quadprog)
source(file = "./R/ADMM.R")

n <- 100;
P = matrix(rnorm(n*n),n,n)
P = nearPD(P)$mat

q = rnorm(n)
r = 0

l = rnorm(n)
u = rnorm(n)
lb = pmin(l, u)
ub = pmax(l, u)

# quadprog
print("solve QP using quadprog")
Amat <- t(rbind(diag(n),-diag(n)))
bvec <- c(lb, -ub)
t_start <- proc.time()
obj <- solve.QP(P, q, Amat, bvec)$val
print(proc.time()-t_start)

print("-------------------------------")
# ADMM
print("solve QP using ADMM")
t_start <- proc.time()
obj2 <- QP_ADMM(P, q, r, lb, ub, 1, 1)
print(proc.time()-t_start)

print("-------------------------------")
# linprog
print("solve LP using linprog")
cc <- rnorm(n)
t_start <- proc.time()
obj3 <- linprog(cc, lb = lb, ub = ub)
print(proc.time()-t_start)