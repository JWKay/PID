# Computation of Idep and MMI PIDs for
# univariate Gaussian predictors and target.

# Inputs are p, q, r, the Pearson correlations between
# X_0 and X_1, X_0 and Y, X_1 and Y.

idepGU <- function(p, q, r){

mat = matrix(c(1, p, q, p, 1, r, q, r, 1), ncol =3)
ev = eigen(mat)$values
if(!all(ev>0)) stop("covariance matrix is not positive definite")

pp =1-p^2; qq=1-q^2; rr=1-r^2
dd = 1 -p^2 -q^2 -r^2 + 2*p*q*r

i0y = 0.5*log2(1/qq)
i1y= 0.5*log2(1/rr)
mitot = 0.5*log2(pp/dd)

b = i0y
i =0.5*log2((1-q^2*r^2)/qq)
k = 0.5*log2(pp*rr/dd)



# Idep PID

unq0 = min(b, i, k)
unq1 =  i1y -i0y + unq0
red = i0y - unq0
syn = mitot -i1y - unq0


idep=c(unq0, unq1, red, syn)

# MMI PID

redM = min(i0y, i1y)
synM = ifelse(i0y < i1y, mitot - i1y, mitot - i0y) 
unq0M = ifelse(i0y < i1y, 0, i0y - i1y) 
unq1M = ifelse(i1y < i0y, 0, i1y - i0y)
mmi =c(unq0M, unq1M, redM, synM)

list(idep=idep, mmi=mmi)


}









# Computation of  Idep and MMI PIDs for
# multivariate Gaussian predictors and target.

# Inputs are:

#	sizes	a numeric list containing the values of n0, n1, n2
# 	mat		a positive definite covariance or correlation matrix




idepGM <-  function(sizes, mat){

# Test for matrix being positive definite

ev = eigen(mat)$values
if(!all(ev>0)) stop("input matrix is not positive definite")


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
SigZ = rbind(r1, r2, r3)

dSigZ = det(SigZ)

b = 0.5* log2(1/dQ)
i = 0.5* log2(dQR/dQ)
k = 0.5* log2( dP*dR/dSigZ)

mitot = 0.5* log2(dP/dSigZ)
i0y =  b
i1y = 0.5* log2(1/dR)
 
# Idep PID

unq0 = min(b, i, k)
red = i0y - unq0
unq1 = i1y - i0y + unq0
syn = mitot - i1y - unq0
idep = c(unq0, unq1, red, syn)

# MMI PID

redM = min(i0y, i1y)
synM =  ifelse(i0y < i1y, mitot - i1y, mitot - i0y)
unq0M = ifelse(i0y < i1y, 0, i0y - i1y)
unq1M = ifelse(i1y < i0y, 0, i1y - i0y)
mmi = c(unq0M, unq1M, redM, synM)

list(idep = idep, mmi = mmi)

}








# Code to compute p-values for the tests of deviance,
# as described in Appendix D

# The univariate case:


# Inputs are p, q, r, the Pearson correlations between
# X_0 and X_1, X_0 and Y, X_1 and Y, and also n, the 
# number of observations.


DevTestU <- function(p, q, r, n){

mat = matrix(c(1, p, q, p, 1, r, q, r, 1), ncol =3)
ev = eigen(mat)$values
if(!all(ev>0)) stop("covariance matrix is not positive definite")

pp =1-p^2; qq=1-q^2; rr=1-r^2
dd = 1 -p^2 -q^2 -r^2 + 2*p*q*r

dev1 = n*log(1/dd)
dev2= dev1 + n*log(pp)
dev3= dev1 + n*log(qq)
dev4= dev1 + n*log(rr)
dev5 = dev1 + n*log(pp*qq)
dev6= dev1 + n*log(pp*rr)
dev7 = dev1 + n*log(qq*rr)

deviances =c(dev1, dev2, dev3, dev4, dev5, dev6, dev7)
df =c(3, 2, 2, 2, 1, 1, 1)

pvals = 1 - pchisq(deviances, df)

models =c("U1", "U2", "U3", "U4", "U5", "U6", "U7")

list(models = models, pvalues = pvals)

}

# Code to compute p-values for the tests of deviance,
# as described in Appendix D

# The multivariate case:
# Inputs are :

#	sizes	a numeric list containing the values of n0, n1, n2
# 	mat		a positive definite covariance or correlation matrix
#   n		the number of observations


DevTestM <- function(sizes, mat, n){

ev = eigen(mat)$values
if(!all(ev>0)) stop("covariance matrix is not positive definite")


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

r1 = cbind(I0, P, Q)
r2 = cbind(P1, I1, R)
r3 = cbind(Q1, R1, I2)
SigZ = rbind(r1, r2, r3)

dSigZ = det(SigZ)


dev1 = n *log(1/dSigZ)
dev2 = dev1 + n*log(dP)
dev3 = dev1 + n*log(dQ)
dev4 = dev1 + n*log(dR)
dev5 = dev1 + n*log(dP*dQ)
dev6 = dev1 + n*log(dP*dR)
dev7 = dev1 + n*log(dQ*dR)

deviances= c(dev1, dev2, dev3, dev4, dev5, dev6, dev7)

N01=n0*n1; N02=n0*n2; N12=n1*n2

df = c(N01+N02+N12, N02+N12, N01 + N12, N01+N02, N12, N02, N01)

pvals = 1-pchisq(deviances, df)

models =c("M1", "M2", "M3", "M4", "M5", "M6", "M7")

list(models =models, pvalues = pvals)

}











