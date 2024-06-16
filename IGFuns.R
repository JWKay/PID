## Functions for implementation of the IG PID for MV Gaussians

IG_GaussM_PQR = function(sizes, P, Q, R) {
  
  n0=sizes[1]; n1=sizes[2]; n2=sizes[3]
  
  ind0 = seq(1,n0); ind1=seq(n0+1, n0+n1); ind2= seq(n0+n1+1, n0+n1+n2)
  
  I0 = diag(1, nrow = n0, ncol = n0)
  I1 = diag(1, nrow = n1, ncol = n1)
  I2 = diag(1, nrow = n2, ncol = n2)
  
  P1 = t(P)
  Q1 = t(Q)
  R1 = t(R)
  
  dP = det(I1 - P1%*%P)
  dQ = det(I2 - Q1%*%Q)
  dR = det(I2 - R1%*%R)
  dQR = det(I1 - R%*%Q1%*%Q%*%R1)
  
  
  
  
  r1 = cbind(I0, P, Q)
  r2 = cbind(P1, I1, R)
  r3 = cbind(Q1, R1, I2)
  
  ## Full covariance matrix
  
  Sig = rbind(r1, r2, r3)
  
  # Test for covariance matrix being positive definite
  
  ev = eigen(Sig)$values
  if(!all(ev>0)) stop("covariance matrix is not positive definite")
  
  PD = "yes"
  
  
  ## Compose Sig5
  
  R5 = t(P) %*% Q
  r51 = cbind(I0, P, Q)
  r52 = cbind(P1, I1, R5)
  r53 = cbind(Q1, t(R5), I2)
  
  Sig5 = rbind(r51, r52, r53)
  
  ## Compose Sig6
  
  Q6 = P %*% R
  r61 = cbind(I0, P, Q6)
  r62 = cbind(P1, I1, R)
  r63 = cbind(t(Q6), R1, I2)
  
  Sig6 = rbind(r61, r62, r63)
  
  
  
  sig5_inv = solve(Sig5)
  sig6_inv = solve(Sig6)
  
  ## Find the range of feasible values for t.
  
  feas_test=function(t){m = (1-t)* sig5_inv + t*sig6_inv; x = all(eigen(m)$values>0); 2*x -1}
  root_hi = uniroot(feas_test, c(1, 100))$root
  root_lo = uniroot(feas_test, c(-100,0))$root
  
  feas = c(root_lo, root_hi)
  
  
  KLdiv = function(t){
    mm = solve((1-t)*sig5_inv + t*sig6_inv)
    0.5*(sum(log(eigen(mm)$values))  - sum(log(eigen(Sig)$values)))
    
  }
  
  ## can blow up near the roots, so add a nudge inwards
  
  ## plot of g(t) appear in the plot window in RStudio
  
  xlo = root_lo + 0.1 ; xhi = root_hi - 0.1
  xvals = as.matrix(seq(xlo, xhi, length.out=100), ncol =1)
  fvals0 = apply(xvals, 1, KLdiv)
  
  plot(xvals, fvals0, type ="l", xlab = "t", ylab ="g(t)", col = "blue", lwd =2)
  
  # Initial value for t-star
  
  x0 = (root_lo + root_hi)/2
  
  # Perform the constrained optimisation
  
  out = optim(x0, KLdiv, method = "L-BFGS-B", lower = root_lo, upper = root_hi)
  
  tstar = out$par
  syn = out$value
  
  ### Mutual informations and interaction information
  
  i13 = 0.5*log(1/dQ)
  i23 = 0.5*log(1/dR)
  i13G2 = 0.5*log(dP*dR/det(Sig))
  i23G1 =  0.5*log(dP*dQ/det(Sig))
  jmi = 0.5*log(dP/det(Sig))
  ii = jmi - i13 - i23
  
  
  ## convert output to bits
  
  inf = c(i13, i23, i13G2, i23G1, jmi, ii)/log(2)
  
  unq1 = i13G2 - syn
  unq2 = i23G1 - syn
  
  red = i13 - unq1
   
  ## convert pid components to bits
  
  pid = c(unq1, unq2, red, syn)/log(2)
  
  ## return the numerical results
  
  list(PD = PD, feas = feas, t_star = tstar, inf = inf, pid = pid)
  
}


IG_GaussM_Dat = function(sizes, mat){
  
  # Convert input matrix to correlation form
  
  d = diag(1/sqrt(diag(mat)))
  mat = d%*%mat%*%d
  
  n0=sizes[1]; n1=sizes[2]; n2=sizes[3]
  
  ind0 = seq(1,n0); ind1=seq(n0+1, n0+n1); ind2= seq(n0+n1+1, n0+n1+n2)
  
  I0 = diag(1, nrow = n0, ncol = n0)
  I1 = diag(1, nrow = n1, ncol = n1)
  I2 = diag(1, nrow = n2, ncol = n2)
  
  
  S00 = mat[ind0, ind0]
  S01 = mat[ind0, ind1]
  S02 = mat[ind0, ind2]
  S11 = mat[ind1, ind1]
  S12 = mat[ind1, ind2]
  S22 = mat[ind2, ind2]
  
  InvSq00 = backsolve(chol(S00), I0)
  InvSq11 = backsolve(chol(S11), I1)
  InvSq22 = backsolve(chol(S22), I2)
  
  P = t(InvSq00)%*%S01%*%InvSq11
  
  Q= t(InvSq00)%*%S02%*%InvSq22
  
  R = t(InvSq11)%*%S12%*%InvSq22
  
  
  
  P1 = t(P)
  Q1 = t(Q)
  R1 = t(R)
  
  dP = det(I1 - P1%*%P)
  dQ = det(I2 - Q1%*%Q)
  dR = det(I2 - R1%*%R)
  dQR = det(I1 - R%*%Q1%*%Q%*%R1)
  
  
  
  
  r1 = cbind(I0, P, Q)
  r2 = cbind(P1, I1, R)
  r3 = cbind(Q1, R1, I2)
  
  ## Full covariance matrix
  
  Sig = rbind(r1, r2, r3)
  
  ## Check that Sig is positive definite
  
  ev = eigen(Sig)$values
  if(!all(ev>0)) stop("covariance matrix is not positive definite")
  
  PD = "yes"
  
  ## Compose Sig5
  
  R5 = t(P) %*% Q
  r51 = cbind(I0, P, Q)
  r52 = cbind(P1, I1, R5)
  r53 = cbind(Q1, t(R5), I2)
  
  Sig5 = rbind(r51, r52, r53)
  
  ## Compose Sig6
  
  Q6 = P %*% R
  r61 = cbind(I0, P, Q6)
  r62 = cbind(P1, I1, R)
  r63 = cbind(t(Q6), R1, I2)
  
  Sig6 = rbind(r61, r62, r63)
  
  
  
  sig5_inv = solve(Sig5)
  sig6_inv = solve(Sig6)
  
  ## Find the range of feasible values for t.
  
  feas_test=function(t){m = (1-t)* sig5_inv + t*sig6_inv; x = all(eigen(m)$values>0); 2*x -1}
  root_hi = uniroot(feas_test, c(1, 100))$root
  root_lo = uniroot(feas_test, c(-100,0))$root
  
  feas = c(root_lo, root_hi)
  
  
  KLdiv = function(t){
    mm = solve((1-t)*sig5_inv + t*sig6_inv)
    0.5*(sum(log(eigen(mm)$values))  - sum(log(eigen(Sig)$values)))
    
  }
  
  ## can blow up near the roots, so add a nudge inwards
  
  ## plot of g(t) appear in the plot window in RStudio
  
  xlo = root_lo + 0.1 ; xhi = root_hi - 0.1
  xvals = as.matrix(seq(xlo, xhi, length.out=100), ncol =1)
  fvals0 = apply(xvals, 1, KLdiv)
  
  plot(xvals, fvals0, type ="l", xlab = "t", ylab ="g(t)", col = "blue", lwd =2)
  
  # Initial value for t-star
  
  x0 = (root_lo + root_hi)/2
  
  # Perform the constrained optimisation
  
  out = optim(x0, KLdiv, method = "L-BFGS-B", lower = root_lo, upper = root_hi)
  
  tstar = out$par
  syn = out$value
  
  ### Mutual informations and interaction information
  
  i13 = 0.5*log(1/dQ)
  i23 = 0.5*log(1/dR)
  i13G2 = 0.5*log(dP*dR/det(Sig))
  i23G1 =  0.5*log(dP*dQ/det(Sig))
  jmi = 0.5*log(dP/det(Sig))
  ii = jmi - i13 - i23
  
  
  ## convert output to bits
  
  inf = c(i13, i23, i13G2, i23G1, jmi, ii)/log(2)
  
  unq1 = i13G2 - syn
  unq2 = i23G1 - syn
  
  red = i13 - unq1
  
  ## convert pid components to bits
  
  pid = c(unq1, unq2, red, syn)/log(2)
  
  ## return the numerical results
  
  list(PD = PD, feas = feas, t_star = tstar, inf = inf, pid = pid)
  
}

 

## Code for the trivariate NQ Gaussian

IG_GaussU_pqr = function(p, q, r){
  
  mat = matrix(c(1, p, q, p, 1, r, q, r, 1), ncol =3)
  ev = eigen(mat)$values
  if(!all(ev>0)) stop("covariance matrix is not positive definite")
  
  pp =1-p^2; qq=1-q^2; rr=1-r^2
  dd = 1 -p^2 -q^2 -r^2 + 2*p*q*r
  
  Sig = mat
  ### Mutual informations
  
  i13 = 0.5*log(1/qq)
  i23 = 0.5*log(1/rr)
  i13G2 = 0.5*log(pp*rr/det(Sig))
  i23G1 =  0.5*log(pp*qq/det(Sig))
  jmi = 0.5*log(pp/det(Sig))
  ii = jmi - i13 - i23
  
  B = (q -p*r)^2 + r^2 * (1-p^2)
  A = (q -p*r)^2 + r^2 * (1-p^2)*(1-q^2)
  
  num = (1-p^2)*(1-q^2)*(1-r^2)*B
  den = det(Sig) *A
  syn = 0.5*log(num/den)
  
  t1 = q*(q - p*r) ; t2 = r*(r -p*q)
  tval = t1/(t1 + t2)
  
  ## convert output to bits
  
  inf = c(i13, i23, i13G2, i23G1, jmi, ii)/log(2)
  
  unq1 = i13G2 - syn
  unq2 = i23G1 - syn
  
  red = i13 - unq1
  
  
  pid = c(unq1, unq2, red, syn)/log(2)
  
  list( inf = inf, pid = pid)
  
}


