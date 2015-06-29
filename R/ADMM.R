norm_vec <- function(x) sqrt(sum(x^2))
#####################################
# min (1/2)*x'*P*x + q'*x + r
# s.t. lb <= x <= ub
# 
# Return: solution in a vector
# 
# rho is the augmented Lagrangian paramter
# alpha is the over-relaxation parameter (typical values are [1,0, 1.8])
#####################################
QP_ADMM <- function(P, q, r, lb, ub, rho, alpha) {
  library(pracma)
  
  t_start <- proc.time() 
  
  QUIET = 1;
  MAX_ITER = 1000;
  ABSTOL = 1e-4;
  RELTOL = 1e-2;
  
  n = nrow(P)
  X = rep(0,n)
  z = rep(0,n)
  u = rep(0,n)
  
  history.objval = c()
  history.r_norm = c()
  history.s_norm = c()
  history.eps_pri = c()
  history.eps_dual = c()
  
  if(!QUIET) {
    sprintf("%3s    %10s    %10s    %10s    %10s    %10s", 'iter',
            'r norm', 'eps pri', 's norm', 'eps dual', 'objective')
  }
  
  for (k in 1:MAX_ITER) {
    if (k>1) {
      x = mldivide(R, mldivide(t(R), rho*(z-u)-q))
    } else {
      R = as.matrix(chol(P + rho*diag(n)))
      x = mldivide(R, mldivide(t(R), rho*(z-u)-q))
    }
    
    # z-update with relaxation
    zold = z
    x_hat = alpha*x+(1-alpha)*zold
    z = min(ub, max(lb, x_hat+u))
    
    # u-update
    u = u + x_hat - z
    
    # diagnostics, reporting, termination checks
    history.objval = c(history.objval, 1/2*t(x)%*%P%*%x + t(q)%*%x + r)
    
    history.r_norm = c(history.r_norm, norm_vec(x-z))
    history.s_norm = c(history.s_norm, norm_vec(-rho*(z-zold)))
    
    history.eps_pri = c(history.eps_pri, sqrt(n)*ABSTOL+RELTOL*max(norm_vec(x),norm_vec(-z)))
    history.eps_dual = c(history.eps_dual, sqrt(n)*ABSTOL+RELTOL*norm_vec(rho*u))
    
    if (!QUIET) {
      sprintf('%3d    %10.4f    %10.4f    %10.4f    %10.4f    %10.2f', k, 
              history.r_norm[k], history.eps_pri[k], 
              history.s_norm[k], history.eps_dual[k], history.objval[k])
    }
    
    if (history.r_norm[k]<history.eps_pri[k] && history.s_norm[k]<history.eps_dual[k]) {
      break
    }
  }
  
  if (!QUIET) {
    print(proc.time()-t_start)
  }
  
  return(history.objval[k])
}
