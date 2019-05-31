###
# Simulation study Marcel
###

rm(list=ls())

#setwd
setwd('C:/Users/u0105757/Desktop/Rotterdam study/simulation paper')

#libraries
source(file = "functions.R")


### settings (adjust)
Nresp<-200
Ntotal<-12
Nfixed<-6 
prior<-"inf"  
hetero<-"high"
iterations<-10  


#initialize 
Ntasks<- Ntotal - Nfixed

truePOPmean<-c(-0.1, 0.1, -0.4,	0.4, -0.8, 0.8)
truePOPstd<-diag(length(truePOPmean))
if (hetero=="low"){diag(truePOPstd)<-c(0.025, 0.025, 0.1, 0.1, 0.2, 0.2)
} else if (hetero=="high") {diag(truePOPstd)<-c(0.1, 0.1, 0.4, 0.4, 0.8, 0.8)
} else {print("wrong heterogeneity specification")}

profs<-profiles(lvls= c(3,3,3), coding ="contr.treatment")[[2]] #dummy coding
fcombs<-full_sets(cand= profs, n_alts = 2)
n_alts<-2

#draw samples from true distribution (each sample = respondent)
if (Nresp == 100){
  if (hetero=="low"){
    betas<-read.table(file = "lowH_100.txt", sep = "", dec=",")
  } else if (hetero=="high") {betas<-read.table(file = "highH_100.txt", sep = "", dec=",")
  } else {print("wrong heterogeneity specification")}
  
} else if (Nresp == 200) {
  if (hetero=="low"){
    betas<-read.table(file = "lowH_200.txt", sep = "", dec=",")
  } else if (hetero == "high") {betas<-read.table(file = "highH_200.txt", sep = "", dec=",")}
} else {print("Nresp is wrong")}


#prior distribution
if (prior == "uninf"){
  prior_m<-c(0,0,0,0,0,0)
  prior_std<-diag(length(truePOPmean))
  diag(prior_std)<- 10
} else if (prior== "inf") {
  prior_m<-truePOPmean
  prior_std<-truePOPstd
} 

prior_cov<-(prior_std)^2


##empty vectors
total_des<-list(Nresp)
total_Y<-list(Nresp)
fullX<-list(iterations)
fully<-list(iterations)

#initial fixed designs (2 blocks of 6)
if (hetero =="low"){
fix1<-read.table(file="lowH_fix6_block1.txt", sep="")
fix2<-read.table(file="lowH_fix6_block2.txt", sep="")
} else if (hetero =="high"){
fix1<-read.table(file="highH_fix6_block1.txt", sep="")
fix2<-read.table(file="highH_fix6_block2.txt", sep="")
}

#############
######## SIMULATE
#############

#repeat simulations
for (i in 1: iterations){
  
  #for every respondend 
  for (r in 1: Nresp){
   
    #initial design (2 different blocks of 6 initial sets)
    #even participant --> set 1, odd --> set 2
    ifelse (r %% 2 == 0, design_I<-fix1, design_I<-fix2)
    
    #Individual betas
    beta_I<-as.matrix(betas[r, ])
    
    #generate individual responses for inital set(s)
    Y_I<-numeric()
    
    for (s in 1:(nrow(design_I)/n_alts)){
      setn<-as.matrix(design_I[((s*2)-1) : (s*2), ])
      Y_I<-c(Y_I, respond(par = beta_I, set = setn, n_alts = 2))
    }
    
    #For each respondent get optimal design sequentially (Ntasks)
    design_I<-rbind(design_I, matrix(data = NA, nrow =Ntasks*n_alts, ncol = ncol(profs)))
    
    for (c in 1:Ntasks){
      
      ### update parameters
      #take sets and responses so far
      ind<-apply(design_I, 1, function(x) all(is.na(x)))
      des_i <- as.matrix(design_I[ !ind, ])
      y<-Y_I[1:nrow(des_i)]
      
      #sample from posterior
      samples<-impsampling(prior_mean = prior_m, prior_covar = prior_cov, des = des_i, n_alts = n_alts, Y=y, m=8)
      sam<-samples[[1]]
      w<-samples[[2]]
      
      
      #new set
      new_set<-KL_info(fp= profs, fcomb= fcombs, par_samples = sam, weights = w, n_alts = n_alts)
      
      #update and respond
      design_I[((Nfixed*n_alts) +(2*c)-1):((Nfixed*n_alts) +(2*c)) , ]<- new_set
      resp<-respond(par = beta_I, set = new_set, n_alts = n_alts, bin = TRUE)
      Y_I<-c(Y_I, resp)
      
    }
    
    #save all individual designs together
    total_des[[r]]<-design_I
    total_Y[[r]]<-Y_I  #binary respons
    
  }
  
  
  #bind together
  fullX[[i]]<-do.call(rbind, total_des)
  fully[[i]]<-do.call(rbind, total_Y)
  
}

X<-do.call(rbind, fullX)
Y<-do.call(rbind, fully)

#binary response to discrete
D<-transf(Y)

#write to excel file
#write.table(X, file="X_sim24.txt", quote = FALSE, sep=" ", row.names = FALSE)
#write.table(D, file="Y_sim24.txt", quote = FALSE, sep=" ", row.names = FALSE)


