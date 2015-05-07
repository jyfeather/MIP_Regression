# Branch and Bound implementation (BAB)
# http://www.mat.univie.ac.at/~neum/glopt/software_g.html#bb_codes

# solve sub problems
# currently, only for binary programming problem
sub_prob <- function(numnode, obj, con1, con2, con3, bound, max, nrow, ncol) {
  i <- con1$i
  j <- con1$j
  v <- con1$v
  if (numnode > 0) {
    # branch procedure
    for (k in 1:numnode) {
      i <- c(i, nrow+k)
      loc <- sample(ncol, 1)
      j <- c(j, loc)
    }
    v <- c(v, rep(1, numnode))
    con1 <- simple_triplet_matrix(i, j, v, nrow = max(i), ncol = max(j), dimnames = NULL)
    con2 <- c(con2, rep("==", numnode))
    con3 <- c(con3, round(runif(numnode)))
  }
  # change to LP problems from IP/BP 
  type <- rep("C", max(j))
  res <- Rglpk_solve_LP(obj, con1, con2, con3, bound, type, max)
  # return A,b,c,x*  
  return(list(A=con1, b=con3, c=obj, sol=res$sol))
}
