
############################
### settings choice data ###
############################

library(idefix)
library(abind)

#DESIGN 

# for (i in 1:2000){
# bindes <- sct.alts[sample(nrow(sct.alts), size = M * S, replace = TRUE), ]
# 
# dif <- max(colSums(bindes)) - min(colSums(bindes))
# 
# if(dif < bestdif){
#   bestdif <- dif
#   bestDES <- bindes 
# }
# }
# bindes <- bestDES
# save(bindes, file = "C:\\Users\\u0105757\\Desktop\\SUPER COMP\\SCT ANA R sims\\4 attributes\\bindes_4atts.RData")




#####################
####### DATA ########
#####################

#I = 300
# 
# #true parameters
# TRUEtheta1 <- c(0.1, 0.3, 0.5, 0.55, 0.33, 0.01, 0.08, 0.66, 0.2, 0.11, 0.22, 0.33)
# TRUEalpha1 <- c(0.70, 0.80, 0.90, 0.99)
# 
# TRUEtheta2 <- c(0.1, 0.5, 0.99, 0.65, 0.33, 0.01, 0.08, 0.56, 0.75, 0.15, 0.45, 0.66)
# TRUEalpha2 <- c(0.50, 0.70, 0.90, 0.99)
# 
# TRUEtheta3 <- c(0.1, 0.5, 0.99, 0.65, 0.53, 0.21, 0.08, 0.56, 0.75, 0.88, 0.5, 0.9)
# TRUEalpha3 <- c(0.40, 0.60, 0.90, 0.99)
# 
# TRUEtheta4 <- c(0.3, 0.5, 0.99, 0.65, 0.53, 0.71, 0.88, 0.56, 0.75, 0.88, 0.5, 0.9)
# TRUEalpha4 <- c(0.20, 0.50, 0.70, 0.99)
# 
# dataT = 4
#   
#   if (dataT == 1){
#     TRUEtheta <- TRUEtheta1 
#     TRUEalpha <- TRUEalpha1 
#   }
#   
#   if (dataT == 2){
#     TRUEtheta <- TRUEtheta2 
#     TRUEalpha <- TRUEalpha2 
#   }
#   
#   if (dataT == 3){
#     TRUEtheta <- TRUEtheta3 
#     TRUEalpha <- TRUEalpha3
#   }
#   
#   if (dataT == 4){
#     TRUEtheta <- TRUEtheta4 
#     TRUEalpha <- TRUEalpha4
#   }
#   
#   
#   ### datagen
#   data.ISN <- datagen(TRUEalpha = TRUEalpha, TRUEtheta = TRUEtheta, I = I, TrueT = dataT, bindes = bindes)
# 
#   #procent opt out
#   print(sum(data.ISN[ , , 4]) / sum(data.ISN))
#   
#   # #information dataset
#    data.choice <- vector(mode = 'list')
#    data.choice[["data.ISN"]] <- data.ISN
#    data.choice[["TrueTheta"]] <- TRUEtheta
#    data.choice[["TrueAlpha"]] <- TRUEalpha
#    data.choice[["dataT"]] <- dataT
# 
#   #save dataset
#   name <- paste("C:\\Users\\u0105757\\Desktop\\SUPER COMP\\SCT ANA R sims\\4 attributes\\data.ISN_SCT_R_ANA_4atts_Trule_",
#                 dataT,".RData", sep = "")
#   save(data.choice, file = name)



S = 20 # keuzesets 
M = 3 # 3 alternatieven
N = 4 # totaal (+ opt out) alternatieven per set
J = M * S #totaal alternatieven
PK = 12 #parameters
P = 4 #attributen
O <- 2^M # aantal mogelijke consideration sets
A <- 2^P
levels <- c(3, 3, 3, 3)

I = 300

#####################
####   Design   #####
#####################

load("bindes_4atts.RData")

# all alternatives
alts <- Profiles(lvls = levels, coding = c('D', 'D', 'D', 'D'))

#dummy coding
st <- seq(1, 7, 2) 
end <- st+1
x <- matrix(nrow = nrow(alts), ncol = P)
for(c in 1:P){
  
  a <- st[c]
  b <- end[c]
  
    for(i in 1:nrow(alts)){
      if(isTRUE(all.equal(c(0, 0), as.numeric(alts[i, a:b])))){
        x[i, c] <- 1
      } else{x[i, c] <- 0}
    }
}
sct.alts <- cbind(x[ ,1], alts[ ,1:2], x[ ,2], alts[ ,3:4], x[ ,3], alts[ ,5:6], x[, 4], alts[ ,7:8])
rownames(sct.alts) <- NULL

### Other settings
# mogelijke consideration sets (combinatie alternatieven)
O.pat <- Profiles(rep(2, M), coding = c(rep("D", each = M)))

O.pat.SMO <- array(data = NA, c(S, M, O))
for(i in 1: nrow(O.pat)){
  O.pat.SMO[, , i]  <- matrix(rep(O.pat[i, ], each = S), nrow = S, ncol = M)
}

# SC regel (oplopend aantal attribuutlevels)
pat <- Profiles(rep(2, P), coding = c(rep("D", each = P)))
SCmat <- cbind(pat, rowSums(pat))
SCmat <- SCmat[order(SCmat[ , P + 1]), ]

# mogelijke ANA sets (combinatie attributen)
ana.pat <- Profiles(rep(2, P), coding = c(rep("D", each = P)))
ana.pat <- cbind(ana.pat, rowSums(ana.pat))
ana.pat <- ana.pat[order(ana.pat[ , P + 1]), ]
ana.order <- as.numeric(ana.pat[ , P + 1])
ana.pat <- ana.pat[ , -(P + 1)]


##########
######### Functions 
#########
logit <- function(x) {log(x / (1 - x))}
logist <- function(x) {exp(x) / (1 + exp(x))}

m_acceptable_R_4atts <- function (ana.set, Trule, bindes, theta, SCmat, S, M, H, r1, r2) {
  
  ana.set <- rep(ana.set, times = levels)
  ana.theta <- theta * ana.set
  ana.theta[ana.theta == 0] <- 1
  

  # 1 - kans dat geen enkel attribuut(level) van alternatief geaccepteerd wordt 
  if (Trule == 1) {
    
      temp <- (1 - ana.theta) %o% rep(1, J)
      # matrix met attribuutlevelkansen per alternatief
      temp2 <- t(matrix(temp[t(bindes) == 1], byrow = FALSE, nrow = P))
      
      #kans dat alternatief m in set s in consideration set zit
      probc.SM <- 1 - matrix(apply(temp2, 1, prod), nrow = S, byrow = TRUE)
    
  }
  
  else if (Trule == 2) {
    
    temp <- ana.theta %o% rep(1 , J)
    # matrix met attribuutlevelkansen per alternatief
    temp2 <- t(matrix(temp[t(bindes) == 1], byrow = FALSE, nrow = P))
    
    temp2.JPH <- temp2 %o% rep(1, H)
    SCpat.JPH <- rep(1 , J) %o% t(SCmat[r1:r2, 1:P])
    
    #kans dat alternatief m in set s in consideration set zit
    probc.SM <- 1- matrix(apply(apply(SCpat.JPH * temp2.JPH + ((1 - temp2.JPH) * (1 - SCpat.JPH)),
                                       c(1, 3), prod), 1, sum),
                           nrow = S, byrow = T)
  }
  
  
  else if (Trule == 3) {
    
    temp <- ana.theta %o% rep(1 , J)
    # matrix met attribuutlevelkansen per alternatief
    temp2 <- t(matrix(temp[t(bindes) == 1], byrow = FALSE, nrow = P))
    
    temp2.JPH <- temp2 %o% rep(1, H)
    SCpat.JPH <- rep(1 , J) %o% t(SCmat[r1:r2, 1:P])
    
    #kans dat alternatief m in set s in consideration set zit
    probc.SM <- matrix(apply(apply(SCpat.JPH * temp2.JPH + ((1 - temp2.JPH) * (1 - SCpat.JPH)),
                                   c(1, 3), prod), 1, sum),
                       nrow = S, byrow = T)
  }
  
  
  else if (Trule == 4) {
    
    return(matrix(exp(bindes %*% log(ana.theta)), nrow = S, byrow = TRUE))
    
  }
  return(probc.SM)
}

#data generation function
datagen <- function(TRUEalpha, TRUEtheta, I, TrueT, bindes){
  
  
  data.ISN <- array(NA, c(I, S, N))
  Trule <- TrueT
  #H = aantal patronen voor subsetconj regel, zie SCmat
  if (Trule == 1){ r1 <- 2; r2 <- 16} #not relevant
  if (Trule == 2){ r1 <- 1; r2 <- 5}
  if (Trule == 3){ r1 <- 12; r2 <- 16}
  if (Trule == 4){ r1 <- 16; r2 <- 16} #not relevant
  H <- r2 - r1 + 1
  
  for (i in 1:I){
    
    # which attributes acceptable (alpha)
    ana.set <- rbinom(n = P, size = 1 , prob = TRUEalpha)
    ana.des <- t(diag(rep(ana.set, times = levels)) %*% t(bindes))
    
    #Probabilities accepting alternative given ANA pattern 
    probs.acc <- m_acceptable_R_4atts(ana.set = ana.set, Trule = TrueT, bindes = bindes, theta = TRUEtheta, SCmat = SCmat, 
                                          S = S, M = M, H = H, r1 = r1, r2 = r2 )
    #Realize acceptability pattern 
    # which levels acceptable (theta)
    theta.m <- matrix(rep(TRUEtheta, each = J), nrow = J)
    a.levels <- t(apply(theta.m, 1, function(x) {rbinom(n = PK, size = 1 , prob = x)}))
    
    #combined
    ana.a.des <- ana.des * a.levels
    #reversed
    ana.a.des[, which(rep(ana.set, levels) == 0)] <- bindes[ , which(rep(ana.set, levels) == 0)]
    
    #consideration sets for each choice set  
    c.sets <- matrix(as.integer(rowSums(ana.a.des) >= TrueT), nrow = S, byrow = TRUE) 
    
    pi_wa <- (probs.acc * c.sets) / rowSums(probs.acc * c.sets) 
    pi_wa[is.na(pi_wa)] <- 0
    
    #Include opt-out
    pi_wa <- cbind(pi_wa, (1 - round(rowSums(pi_wa), digits = 4)))
    
    #Draw from multinomial 
    data.ISN[i, , ] <- t(apply(pi_wa, 1, function (x) {rmultinom(1, 1, x)}))
    
  }
  
  return(data.ISN)
  
}
  
#likelihood function 
sct.ana.R.4atts <- function(params, bindes, data.ISNO, I, SCmat, O.pat.SMO, Trule, O, S, M, A) {
  
  #parameters
  theta <- params[1:PK]  
  theta[theta < 0.001] <- 0.001
  alpha <- params[(PK+1):(PK + P)]
  
  
  #H = aantal patronen voor subsetconj regel, zie SCmat
  if (Trule == 1){ r1 <- 2; r2 <- 16} #not relevant
  if (Trule == 2){ r1 <- 1; r2 <- 5}
  if (Trule == 3){ r1 <- 12; r2 <- 16}
  if (Trule == 4){ r1 <- 16; r2 <- 16} #not relevant
  H <- r2 - r1 + 1
  
  
  #initiate lists 
  probc.ASM <- probc.ASM.vm <- probc.ASMO.vm <- P.sm_oa <- P.sn_oa <- prob.IS <- Pi.OS <- prob.I <- vector(mode = "list")
  P_ia <- matrix(NA, I, A)
  
  ### for every ANA pattern 
  for (a in 1:nrow(ana.pat)){
    
    #ana.set 
    ana.set <- ana.pat[a, ]
    
    if((P - ana.order[a]) >= Trule){
      
      probc.ASM[[a]] <- matrix(1, S, M)
      
    } else {
      # bereken kans dat alternatief m acceptable is volgens SC-T model voor elk ana patroon (probc.SM)
      probc.ASM[[a]] <- m_acceptable_R_4atts(ana.set = ana.set, Trule = Trule, bindes = bindes, 
                                       theta = theta, SCmat = SCmat, S = S, M = M, H = H, r1 = r1, r2 = r2)
    }

    
    ## CONSIDERATIION
    #copy each ana pat slice 32 times for each consideration pattern
    probc.ASMO.vm[[a]] <- array(rep(probc.ASM[[a]] , O), c(S, M, O))
    
    #O pattern times Vm
    probc.ASMO.vm[[a]] <- probc.ASMO.vm[[a]] * O.pat.SMO
    
    #probability for each M | S, A, O
    tot.sum.Vm.SMO <- tot.sum.vm. <- aperm(replicate(M, array(apply(probc.ASMO.vm[[a]], c(1, 3), sum),
                                                              c(S, O)), simplify="array"), c(1, 3, 2)) 
    P.sm_oa[[a]] <- probc.ASMO.vm[[a]] / tot.sum.Vm.SMO
    
    
    #replace NaN by zero (result of 0/0)
    P.sm_oa[[a]][is.nan(P.sm_oa[[a]])] <- 0
    
    #include opt-out 
    coptout <- array(rep(0, O * S), c(O, S, 1))
    coptout[1, , ] <- 1
    P.sn_oa[[a]] <- abind(P.sm_oa[[a]], aperm(coptout, c(2, 3, 1)), along = 2)
    
    
    #Pi.OS kans consideration pattern given ana pattern for each choice set  
    funp <- function(p) {apply(t(p * t(O.pat)) + t((1 - p) * t(1 - O.pat)), 1, prod)}
    Pi.OS[[a]] <- apply(probc.ASM[[a]], 1, funp )
    
    #for every I calculate prob of data for each o
    probfun <- function(data.SNO){apply(t(Pi.OS[[a]]) * apply(data.SNO * P.sn_oa[[a]], c(1,3), sum), 1, sum)}
    prob.IS[[a]] <- t(apply(data.ISNO, 1, probfun))
    
    #multiply over sets
    prob.I[[a]] <- apply(prob.IS[[a]], 1, prod)
    
    #Probs all individuals for every ana pattern 
    P_ia[ , a] <- prob.I[[a]]
    
  }
  
  #probabilities voor elk ana patroon 
  ana.probs <- apply(t(t(ana.pat) * alpha) +  t(t(1 - ana.pat) * (1 - alpha)), 1 , prod)
  
  LL <- sum(log(apply(t(ana.probs * t(P_ia)), 1, sum)))          
  
  #if(is.nan(LL)){print(params)}
  #print(-LL)
  return(-LL)
  
}
  

#############             ##################
########   Estimation  ########
############### ######################


estim_T4 <- function(dataT){
  
  if (dataT == 1){
    load("data.ISN_SCT_R_ANA_4atts_Trule_1.RData")
  }
  
  if (dataT == 2){
    load("data.ISN_SCT_R_ANA_4atts_Trule_2.RData")
  }
  
  if (dataT == 3){
    load("data.ISN_SCT_R_ANA_4atts_Trule_3.RData")
  }
  
  if (dataT == 4){
    load("data.ISN_SCT_R_ANA_4atts_Trule_4.RData")
  }
  
  data.ISNO <- array(rep(data.choice$data.ISN , O), c(I, S, N, O))
  
  ##### 
  for(estT in 1:P){
    
    
    for (r in 1:4){
      
      #startwaarden
      startparams <- runif(PK + P)
      #startparams <- c(TRUEtheta, TRUEalpha)
      
      # optimizer
      m1 <- optim(startparams, sct.ana.R.4atts, bindes = bindes, data.ISNO = data.ISNO,
                  I = I,
                  SCmat = SCmat, O.pat.SMO = O.pat.SMO, Trule = estT,
                  O = O, S = S, M = M, A = A, 
                  hessian = FALSE, method="L-BFGS-B", lower = rep(0.001, (PK+P)), upper = rep(0.999, (PK+P)))
      
      
      #write
      par <- round(m1$par, digits = 4)
      value <- rep(round(m1$value, digits = 4), PK + P)
      conv <- rep(m1$convergence, PK + P)
      
      table <- cbind(par, value, conv)
      rownames(table) <- NULL
      
      filio <- paste("MODEL_sct_ana_R_4att_estT_", estT, "_DATA_T", dataT,".txt", sep = "")
      
      write.table(data.frame(table), file = filio, append = TRUE, col.names = FALSE,  row.names = FALSE)
      cat('\n', file = filio, append = TRUE)  
      
    }
  }
}


