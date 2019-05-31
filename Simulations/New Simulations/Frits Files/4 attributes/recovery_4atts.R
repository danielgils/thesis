

### Recovery sct ANA R sims 4 attributes 

#true parameters
TRUEtheta1 <- c(0.1, 0.3, 0.5, 0.55, 0.33, 0.01, 0.08, 0.66, 0.2, 0.11, 0.22, 0.33)
TRUEalpha1 <- c(0.70, 0.80, 0.90, 0.99)

TRUEtheta2 <- c(0.1, 0.5, 0.99, 0.65, 0.33, 0.01, 0.08, 0.56, 0.75, 0.15, 0.45, 0.66)
TRUEalpha2 <- c(0.50, 0.70, 0.90, 0.99)

TRUEtheta3 <- c(0.1, 0.5, 0.99, 0.65, 0.53, 0.21, 0.08, 0.56, 0.75, 0.88, 0.5, 0.9)
TRUEalpha3 <- c(0.40, 0.60, 0.90, 0.99)

TRUEtheta4 <- c(0.3, 0.5, 0.99, 0.65, 0.53, 0.71, 0.88, 0.56, 0.75, 0.88, 0.5, 0.9)
TRUEalpha4 <- c(0.20, 0.50, 0.70, 0.99)


## estimated params 
est1.1.theta <- c(0.1034 ,
            0.3129 ,
            0.5522 ,
            0.4773 ,
            0.2965 ,
            0.001 ,
            0.0054 ,
            0.6548,
            0.1797 ,
            0.1101 ,
            0.244 ,
            0.3325) 

est1.1.ana <- c(0.9536, 
                0.8443, 
            0.9669, 
            0.6503 )


plot(TRUEtheta1, est1.1.theta, ylim=c(0,1), xlim = c(0,1), col='blue')
points(TRUEalpha1, est1.1.ana, col='green')
cor(c(TRUEtheta1, TRUEalpha1), c(est1.1.theta, est1.1.ana))


est1.2.theta <- c(0.001, 
                  0.2346, 
                  0.5011 ,
                  0.584 ,
                  0.4401 ,
                  0.2048 ,
                  0.001 ,
                  0.6531 ,
                  0.1759 ,
                  0.001 ,
                  0.1513 ,
                  0.2507 )
est1.2.ana <- c(
                  0.9612, 
                  0.9384 ,
                  0.768 ,
                  0.7308)



plot(TRUEtheta1, est1.2.theta, ylim=c(0,1), xlim = c(0,1), col='blue')
points(TRUEalpha1, est1.2.ana, col='green')
cor(c(TRUEtheta1, TRUEalpha1), c(est1.2.theta, est1.2.ana))



est1.3.theta <- c(
0.0809, 
0.2957, 
0.541, 
0.4773, 
0.2965, 
0.001, 
0.0632, 
0.6748 ,
0.2274, 
0.0783, 
0.217, 
0.3087)

est1.3.ana <- c(
0.784, 
0.762, 
0.947, 
0.8949) 


plot(TRUEtheta1, est1.3.theta, ylim=c(0,1), xlim = c(0,1), col='blue')
points(TRUEalpha1, est1.3.ana, col='green')
cor(c(TRUEtheta1, TRUEalpha1), c(est1.3.theta, est1.3.ana))


################
#attr 1
a=1
b=3
c=1

plot(c(TRUEtheta1[a:b], TRUEalpha1[c]), c(est1.3.theta[a:b], est1.3.ana[c]), col=c('blue', 'blue', 'blue','green'), type = 'l', xlim = c(0,1), ylim=c(0,1))
points(c(TRUEtheta1[a:b], TRUEalpha1[c]), c(est1.2.theta[a:b], est1.2.ana[c]), col=c('blue', 'blue', 'blue','green'), type='l')

plot((TRUEtheta1[a:b] * TRUEalpha1[c]), (est1.3.theta[a:b] * est1.3.ana[c]), col='blue', xlim = c(0,1), ylim=c(0,1))
points((TRUEtheta1[a:b] * TRUEalpha1[c]), (est1.2.theta[a:b] * est1.2.ana[c]), col='blue', xlim = c(0,1), ylim=c(0,1))



#attr 2
a=4
b=6
c=2


plot(c(TRUEtheta1[a:b], TRUEalpha1[c]), c(est1.3.theta[a:b], est1.3.ana[c]), col=c('blue', 'blue', 'blue','green'), type = 'l', xlim = c(0,1), ylim=c(0,1))
points(c(TRUEtheta1[a:b], TRUEalpha1[c]), c(est1.2.theta[a:b], est1.2.ana[c]), col=c('blue', 'blue', 'blue','green'), type='l')

plot((TRUEtheta1[a:b] * TRUEalpha1[c]), (est1.3.theta[a:b] * est1.3.ana[c]), col='blue', xlim = c(0,1), ylim=c(0,1))
points((TRUEtheta1[a:b] * TRUEalpha1[c]), (est1.2.theta[a:b] * est1.2.ana[c]), col='blue', xlim = c(0,1), ylim=c(0,1))


#attr 3
a=7
b=9
c=3


plot(c(TRUEtheta1[a:b], TRUEalpha1[c]), c(est1.3.theta[a:b], est1.3.ana[c]), col=c('blue', 'blue', 'blue','green'), type = 'l', xlim = c(0,1), ylim=c(0,1))
points(c(TRUEtheta1[a:b], TRUEalpha1[c]), c(est1.2.theta[a:b], est1.2.ana[c]), col=c('blue', 'blue', 'blue','green'), type='l')

plot((TRUEtheta1[a:b] * TRUEalpha1[c]), (est1.3.theta[a:b] * est1.3.ana[c]), col='blue', xlim = c(0,1), ylim=c(0,1))
points((TRUEtheta1[a:b] * TRUEalpha1[c]), (est1.2.theta[a:b] * est1.2.ana[c]), col='blue', xlim = c(0,1), ylim=c(0,1))





0.0649 6623.141 0
0.5231 6623.141 0
0.999 6623.141 0
0.6581 6623.141 0
0.3941 6623.141 0
0.0339 6623.141 0
0.001 6623.141 0
0.5346 6623.141 0
0.6808 6623.141 0
0.1574 6623.141 0
0.4698 6623.141 0
0.6357 6623.141 0
0.4769 6623.141 0
0.7351 6623.141 0
0.9464 6623.141 0
0.9674 6623.141 0

# T=3


driet<- c(0.0717, 
0.4836, 
0.9936, 
0.6274, 
0.5424, 
0.216, 
0.0855, 
0.5667 ,
0.7828 ,
0.8552 ,
0.461 ,
0.8686 ,
0.4256 ,
0.5549 ,
0.9185 ,
0.9375 )

plot(c(TRUEtheta3, TRUEalpha3),  driet)
cor(c(TRUEtheta3, TRUEalpha3),  driet)
0.63/0.13
0.13/0.63
