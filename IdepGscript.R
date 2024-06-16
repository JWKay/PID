# A script to run the various functions available to compute
# the Idep and Immi PIDs, and also perform tests of deviance


# First, we set the working directory

 setwd("/Users/jk/Desktop/IdepFuns")
 
# Then, we load the function definitions ito R

  source("IdepGaussR.txt")
  
  
  


# Example 5

 p=0.8484; q = 0.6383; r=0.7168; nobs =84

gsout = idepGU(p, q, r)
gsout

gstests =DevTestU(p, q, r, nobs)
gstests



# Example 6

ca = read.csv("bones.csv", header=T)
caf = ca[ca$G=="f", 2:12]


# First analysis: (H, W, A), (DC, DR), (CF, CH)

ind=c(1,2,3, 7, 9, 10, 11)

sizes=c(3,2,2); nobs =73

caf1 = caf[, ind]

mat = cov(caf1) 


out1 =idepGM(sizes, mat)
out1
tests1 =DevTestM(sizes, mat, 73)
tests1




# Second analysis: 

ind=c(1,2,3, 4, 5, 6, 8, 10, 11)

sizes=c(3,4,2); nobs=73

caf1 = caf[, ind]

mat = cov(caf1)



out2= idepGM(sizes, mat)
out2

tests2 =DevTestM(sizes, mat, nobs)
tests2













