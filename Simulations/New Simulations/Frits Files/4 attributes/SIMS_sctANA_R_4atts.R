###
# DATAGEN 3 attributes 
##


rm(list= ls())


source('C:\\Users\\u0105757\\Desktop\\SUPER COMP\\SCT ANA R sims\\4 attributes\\settings_Rsims_4atts.R')


#DESIGN 
bindes <- sct.alts[sample(nrow(sct.alts), size = M * S, replace = TRUE), ]


#####################
####### DATA ########
#####################

#data.ISN <- datagen(TRUEalpha = TRUEalpha, TRUEtheta = TRUEtheta, I = I, TrueT = TrueT, bindes = bindes)

#procent opt out
#print(sum(data.ISN[ , , 4]) / sum(data.ISN))

# #information dataset 
# data.choice <- vector(mode = 'list') 
# data.choice[["data.ISN"]] <- data.ISN
# data.choice[["TrueTheta"]] <- TRUEtheta
# data.choice[["TrueAlpha"]] <- TRUEalpha
# data.choice[["TrueT"]] <- TrueT

#save dataset 
#name <- paste("C:\\Users\\u0105757\\Desktop\\SUPER COMP\\SCT ANA R sims\\3 attributes\\data.ISN_SCT_R_ANA_3atts_Trule_",TrueT,".RData", sep = "")
#save(data.choice, file = name)

#data.ISNO <- array(rep(data.ISN , O), c(I, S, N, O))

#############             ##################
########   Estimation  ########
############### ######################


#true parameters
TRUEtheta1 <- c(0.1, 0.3, 0.5, 0.55, 0.33, 0.01, 0.08, 0.66, 0.2, 0.11, 0.22, 0.33)
TRUEalpha1 <- c(0.70, 0.80, 0.90, 0.99)

TRUEtheta2 <- c(0.1, 0.5, 0.99, 0.65, 0.33, 0.01, 0.08, 0.56, 0.75, 0.15, 0.45, 0.66)
TRUEalpha2 <- c(0.50, 0.70, 0.90, 0.99)

TRUEtheta3 <- c(0.1, 0.5, 0.99, 0.65, 0.33, 0.01, 0.08, 0.56, 0.75, 0.1, 0.5, 0.9)
TRUEalpha3 <- c(0.40, 0.60, 0.90, 0.99)

TRUEtheta4 <- c(0.3, 0.5, 0.99, 0.65, 0.33, 0.71, 0.88, 0.56, 0.75, 0.2, 0.5, 0.9)
TRUEalpha4 <- c(0.40, 0.60, 0.70, 0.99)

I = 500


parsT <- valuesT <- vector(mode = 'list')


start <- Sys.time()
for(TRUET in 1:P){
  
  
  if (TRUET == 1){
    TRUEtheta <- TRUEtheta1 
    TRUEalpha <- TRUEalpha1 
  }
  
  if (TRUET == 2){
    TRUEtheta <- TRUEtheta2 
    TRUEalpha <- TRUEalpha2 
  }
  
  if (TRUET == 3){
    TRUEtheta <- TRUEtheta3 
    TRUEalpha <- TRUEalpha3
  }
  
  if (TRUET == 4){
    TRUEtheta <- TRUEtheta4 
    TRUEalpha <- TRUEalpha4
  }
  
  
  ### datagen
  data.ISN <- datagen(TRUEalpha = TRUEalpha, TRUEtheta = TRUEtheta, I = I, TrueT = TRUET, bindes = bindes)
  data.ISNO <- array(rep(data.ISN , O), c(I, S, N, O))
  
  
  value <- numeric()
  pars <- vector(mode = 'list')
  
  for (r in 1:10){
    
    #startwaarden
    startparams <- runif(PK+P)
    #startparams <- c(TRUEtheta, TRUEalpha)
    
    # optimizer
    m1 <- optim(startparams, sct.ana.R.3atts, bindes = bindes, data.ISNO = data.ISNO,
                I = I,
                SCmat = SCmat, O.pat.SMO = O.pat.SMO, Trule = TRUET,
                O = O, S = S, M = M, A = A, 
                hessian = FALSE, method="L-BFGS-B", lower = rep(0.001, (PK+P)), upper = rep(0.999, (PK+P)))
    
    print(r)
    
    #STORE
    value[r]<- m1$value
    pars[[r]]<- m1$par
    
  }
  
  ##store
  valuesT[[t]]  <- value
  index <- which.min(value) 
  parsT[[t]] <- pars[[index]]
  print(t)
}
end <- Sys.time()


####
# parsT[[1]] <- m1$par
# plot(parsT[[1]], c(TRUEtheta, TRUEalpha))
# points(parsT[[1]][(1+PK):(PK+P)], TRUEalpha, col='blue')
# cor(parsT[[1]], c(TRUEtheta, TRUEalpha))
# 
# plot(parsT[[2]], c(TRUEtheta, TRUEalpha))
# points(parsT[[2]][(1+PK):(PK+P)], TRUEalpha, col='blue')
# 
# 
# plot(parsT[[3]], c(TRUEtheta, TRUEalpha), main='T3')
# points(parsT[[3]][(1+PK):(PK+P)], TRUEalpha, col='blue')






