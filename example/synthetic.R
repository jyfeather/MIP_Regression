################################################
# LP & QP comparison
################################################
rm(list=ls())

################################################
# generate synthetic problem
################################################
n <- 1000 # number of variables
m <- 800  # number of constraints

A <- matrix(rnorm(n*m),m,n)
b <- rnorm()

# for LP
C <- 

library(lpSolveAPI)
