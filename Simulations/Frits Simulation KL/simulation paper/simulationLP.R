###
# Simulation study Marcel
###

rm(list=ls())

#setwd
setwd('C:/Users/u0105757/Desktop/Rotterdam study/simulation paper')

#libraries
#source(file = "functions.R")
library(IASB)

### settings (adjust)
Nresp<-100
Ntotal<-12
Nfixed<-1 
prior<-"uninf"  
hetero<-"low"
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

#initial fixed designs 
fix1<-read.table(file="init1.txt", sep="")
fix2<-read.table(file="init2.txt", sep="")
fix3<-read.table(file="init3.txt", sep="")
fix4<-read.table(file="init4.txt", sep="")
fix5<-read.table(file="init5.txt", sep="")
fix6<-read.table(file="init6.txt", sep="")
fix7<-read.table(file="init7.txt", sep="")
fix8<-read.table(file="init8.txt", sep="")
fix9<-read.table(file="init9.txt", sep="")
fix10<-read.table(file="init10.txt", sep="")
fix11<-read.table(file="init11.txt", sep="")
fix12<-read.table(file="init12.txt", sep="")

#############
######## SIMULATE
#############

#repeat simulations
for (i in 1: iterations){
  
  #for every respondend 
  for (r in 1: Nresp){
    
    #initial design (12 different initial sets)
    if( 0 <  r & r <  9){design_I<-fix1
    }else if( 9 <= r & r < 17){design_I<-fix2
    }else if(17 <= r & r < 25){design_I<-fix3
    }else if(25 <= r & r < 33){design_I<-fix4
    }else if(33 <= r & r < 41){design_I<-fix5
    }else if(41 <= r & r < 49){design_I<-fix6
    }else if(49 <= r & r < 57){design_I<-fix7
    }else if(57 <= r & r < 65){design_I<-fix8
    }else if(65 <= r & r < 73){design_I<-fix9
    }else if(73 <= r & r < 81){design_I<-fix10
    }else if(81 <= r & r < 89){design_I<-fix11
    }else {design_I<-fix12}
    
    
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
      
      #update parameters
      ind<-apply(design_I, 1, function(x) all(is.na(x)))
      des_i <- as.matrix(design_I[ !ind, ])
      y<-Y_I[1:nrow(des_i)]
      
      samples<-imp_sampling(prior_mean = prior_m, prior_covar = prior_cov, des = des_i, n_alts = n_alts, Y=y, m=9)
      sam<-samples[[1]]
      w<-samples[[2]]
      
      
      #new set
      new_set<-KL_info(fp= profs, fcomb= fcombs, par_samples = sam, weights = w, n_alts = n_alts)
      
      #update
      design_I[((Nfixed*n_alts) +(2*c)-1):((Nfixed*n_alts) +(2*c)) , ]<- new_set
      resp<-respond(par = beta_I, set = new_set, n_alts = n_alts, bin = TRUE)
      Y_I<-c(Y_I, resp)
      
    }
    
    #save all individual designs together
    total_des[[r]]<-design_I
    total_Y[[r]]<-Y_I
    
  }

  
  #bind together
  fullX[[i]]<-do.call(rbind, total_des)
  fully[[i]]<-do.call(rbind, total_Y)
  
}

X<-do.call(rbind, fullX)
Y<-do.call(rbind, fully)

source(file = "functions.R")
D<-transf(Y)


#write to excel file
 write.table(X, file="X_sim1_lattice2_512.txt", quote = FALSE, sep=" ", row.names = FALSE)
 write.table(D, file="Yd_sim1_lattice2_512.txt", quote = FALSE, sep=" ", row.names = FALSE)

