# A script to run the various functions available to compute
# the unnormalised IG, Idep and Immi PIDs.


# First, we set the working directory

 setwd("/Users/jk/Desktop/IG_IDEP_MMI")
 
# Then, we load the Idep and Immifunction definitions into RStudio

  source("IdepFuns/IdepGauss.R")

# We load the IG function definitions into RStudio
  
  source("IGFuns.R")
  

# Example 1: Calcium in Bones

ca = read.csv("bones.csv", header=T)
caf = ca[ca$G=="f", 2:12]


# First analysis: (H, W, A), (DC, DR), (CF, CH)

ind=c(1,2,3, 7, 9, 10, 11)

sizes=c(3,2,2); nobs =73

caf1 = caf[, ind]

mat = cov(caf1) 

out1 = IG_GaussM_Dat(sizes, mat)
out1

out2 =idepGM(sizes, mat)
out2




# Second analysis: 

ind=c(1,2,3, 4, 5, 6, 8, 10, 11)

sizes=c(3,4,2); nobs=73

caf1 = caf[, ind]

mat = cov(caf1)

out3 = IG_GaussM_Dat(sizes, mat)
out3

out4= idepGM(sizes, mat)
out4


# Example 2:  Inputting the matrices P, Q, R

sizes = c(3,4,3)

P = matrix(rep(-0.15, 12), byrow =T, ncol =4)
Q = matrix(rep(0.15, 9), byrow =T, ncol =3)
R = matrix(rep(0.15, 12), byrow =T, ncol =3)
out5= IG_GaussM_PQR(sizes, P, Q, R)
out5

## Setting up the covariance matrix

n0=sizes[1]; n1=sizes[2]; n2=sizes[3]

ind0 = seq(1,n0); ind1=seq(n0+1, n0+n1); ind2= seq(n0+n1+1, n0+n1+n2)

I0 = diag(1, nrow = n0, ncol = n0)
I1 = diag(1, nrow = n1, ncol = n1)
I2 = diag(1, nrow = n2, ncol = n2)

P1 = t(P)
Q1 = t(Q)
R1 = t(R)

r1 = cbind(I0, P, Q)
r2 = cbind(P1, I1, R)
r3 = cbind(Q1, R1, I2)

## Full covariance matrix

Sig8 = rbind(r1, r2, r3)

out6 = idepGM(sizes, Sig8)
out6

## Additional example: scalar inouts and output

p=0.8484; q = 0.6383; r=0.7168 

out7 = IG_GaussU_pqr(p, q, r)
out7

out8 = idepGU(p, q, r)
out8





