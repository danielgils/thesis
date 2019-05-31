#' @title D, A, G and V optimality criteria by hand
#' @author Daniel Gil
#' November 2018
#' 
#' This script is to replicate the numerical example of the first chapter from 
#' the thesis of Roselinde Kessels. I have it on paper, page 27.
#' 

rm(list = ls())

#' The design has 3 choice sets with 2 alternatives each
#' The alternatives has 2 attributes:
#' Attribute 1 has 3 levels
#' Attribute 2 has 2 levels
#' Effects-type coding is used in both
#' 
#' Design matrix (X_0) in numerical value
(x_0 <- matrix(c(1, 2, 2, 1,
                2, 2, 3, 1,
                3, 2, 1, 1), nrow = 6, ncol = 2, byrow = T))

#' Design matrix (X) in effects-type coding
(x <- matrix(c(1, 0, 1,
               0, 1, -1,
               0, 1, 1,
               -1, -1, -1,
               -1, -1, 1,
               1, 0, -1),nrow = 6, ncol = 3, byrow = T))
#' Remember, the number of columns of this design matrix is the number of 
#' parameters k. Here, k = 3
#' For illustration, we use only three prior parameters, randomly drawn from a
#' muultivariate normal distribution with mean c(-1, 0, -1) and covariance 
#' matrix equals to the identity of size 3 (I_3)

#' Information matrix of the multinomial logit model is
#' \sum_{s=1}^3 X'_s(P_s - p_s p'_s) X_s
#' where p_s and P_s are the probability of each choice set and the matrix 
#' respectively.
#' p_js = exp{X'_s \Beta} / (sum_{j=1}^j exp{X'_s \Beta}), with j the number of
#' profiles for each choice set
#' First divide the design matrix by each choice set
x_1 <- x[1:2,]
x_2 <- x[3:4,]
x_3 <- x[5:6,]

#' The draws of the multinormal distribution are
b <- matrix(c( 0.805, -0.080, -0.603,
              -2.083,  2.238, -1.624,
              -0.486, -0.087, -1.594), 3, 3, byrow = T)

#' Using the first draw the probabilities are
#' First compute the exponential for each alternative. This objects are called
#' num because they are the numerator in the probability formula
num_1_1 <- exp(t(x_1[1,]) %*% b[1,1:3])
num_1_2 <- exp(t(x_1[2,]) %*% b[1,1:3])
den_1 <- num_1_1 + num_1_2
p_1_1 <- num_1_1 / den_1 
p_1_2 <- num_1_2 / den_1 
p_1 <- matrix(c(p_1_1, p_1_2), 2, 1)

#' Or alteratively
num_1 <- exp(x_1 %*% b[1,])
den_1_ <- sum(num_1)
p_1_ <- num_1 / den_1_
all.equal(p_1, p_1_)
# P_s = diag(p_s)
P_1 <- diag(as.vector(p_1_))

#' P_s - p_s p'_s
PP_1 <- P_1 - (p_1 %*% t(p_1))
#' X'_s(P_s - p_s p'_s) X_s. 
I_1 <- t(x_1) %*% PP_1 %*% x_1

#' For the SECOND choice set
num_2 <- exp(x_2 %*% b[1,])
den_2 <- sum(num_2)
p_2 <- num_2 / den_2

# P_s = diag(p_s)
P_2 <- diag(as.vector(p_2))

#' P_s - p_s p'_s
PP_2 <- P_2 - (p_2 %*% t(p_2))
#' X'_s(P_s - p_s p'_s) X_s. 
I_2 <- t(x_2) %*% PP_2 %*% x_2

#' For the THIRD choice set
num_3 <- exp(x_3 %*% b[1,])
den_3 <- sum(num_3)
p_3 <- num_3 / den_3

# P_s = diag(p_s)
P_3 <- diag(as.vector(p_3))

#' P_s - p_s p'_s
PP_3 <- P_3 - (p_3 %*% t(p_3))
#' X'_s(P_s - p_s p'_s) X_s. 
I_3 <- t(x_3) %*% PP_3 %*% x_3

# Information Matrix
(I <- I_1 + I_2 + I_3)

# D-optimal criterion
(D <- det(solve(I))^(1/3))
(DD <- det(I)^(-1/3))
# A-optimal criterion
(A <- sum(diag(solve(I))))


# From the package
(x <- matrix(c(1, 0, 1,
               0, 1, -1,
               0, 1, 1,
               -1, -1, -1,
               -1, -1, 1,
               1, 0, -1),nrow = 6, ncol = 3, byrow = T))
(b <- c( 0.805, -0.080, -0.603))
group <- c(1,1,2,2,3,3)
u <- x %*% diag(b)
u <- .rowSums(u, m = nrow(x), n = length(b))
p <- exp(u) / rep(rowsum(exp(u), group), each = 2)
# information matrix
# t(x) %*% (diag(p) - (p %*% t(p))) %*% x this is just for one choice set
(info <- crossprod(x * p, x) - crossprod(rowsum(x * p, group)))
I

