

### Functions
## read in all possible choice sets 

#log posterior
logPost<-function(par, prior_mean, prior_covar, des,  n_alts, Y ){
  
  p<-t(t(des) * par)
  p<-.rowSums(p, m= nrow(des), n=length(par))
  expp<-exp(p)
  p<-expp/rep(rowsum(expp, rep(seq(1, nrow(des)/n_alts, 1), each = n_alts)), each=n_alts)
  
  log_L<-sum(Y*log(p))
  
  logprior2=-0.5*(par -t(prior_mean))%*%solve(prior_covar)%*%(as.matrix(par) - prior_mean)
  
  logpost<-(length(prior_mean)/2*log(2*pi)-0.5*log(det(prior_covar))) + logprior2 + log_L
  
  return(logpost)
  
}

#hessian 
hessian<-function (par, des, covar, n_alts){
  
  des<-as.matrix(des)
  p<- des %*% diag(par)
  p<-.rowSums(p, m= nrow(des), n=length(par))
  expp<-exp(p)
  p<-expp/rep(rowsum(expp, rep(seq(1, nrow(des)/n_alts, 1), each = n_alts)), each=n_alts)
  
  info<-crossprod(des*p, des) - crossprod(rowsum(des*p, rep(seq(1, nrow(des)/n_alts, 1), each=n_alts)))
  
  hess<-(-info - solve(covar))
  
  return(hess)
  
}

#likelihood
Lik <- function(par, des, n_alts, Y){
  
  des<-as.matrix(des)
  
  p<-t(t(des) * par)
  p<-.rowSums(p, m= nrow(des), n=length(par))
  expp<-exp(p)
  p<-expp/rep(rowsum(expp, rep(seq(1, nrow(des)/n_alts, 1), each = n_alts)), each=n_alts)
  
  L<-prod(p^Y)
  
  return(L)
  
}

# density proposal distribution (multivariate t)
g_dens<-function (par, g_mean, g_covar){
  
  df=length(g_mean)
  n<-length(par)
  dif<-g_mean-par
  invcov<-solve(g_covar)
  differ<-as.numeric(t(dif)%*% invcov %*% dif)
  
  iMVSTd=1/(det(g_covar)^(0.5))*(1+((1/df)*differ))^(-(df+length(par))/2)
}

#importance sampling algoritme 
impsampling<- function (prior_mean, prior_covar, des, n_alts, Y, m, ...) {
  prior1 <- (2 * pi)^(-length(prior_mean)/2) * (det(prior_covar))^(-0.5)
  maxest <- maxLik::maxNR(logPost, start = prior_mean, prior_mean = prior_mean, 
                          prior_covar = prior_covar, des = des, Y = Y, n_alts = n_alts)$estimate
  H <- hessian(par = maxest, des = des, covar = prior_covar, n_alts = n_alts)
  g_covar <- -solve(H)
  
  g_draws <- lattice_mvt(mean = maxest, cvar= g_covar, df= length(maxest), m=m)
  
  prior = LK = dens_g = (weights <- numeric(nrow(g_draws)))
  
  for (r in 1:nrow(g_draws)) {
    prior[r] <- prior1 * exp(-0.5 * (g_draws[r, ] - prior_mean) %*% 
                               solve(prior_covar) %*% as.matrix(g_draws[r, ] - prior_mean))
    LK[r] <- Lik(par = g_draws[r, ], des = des, Y = Y, n_alts = n_alts)
    dens_g[r] <- g_dens(par = g_draws[r, ], g_mean = maxest, 
                        g_covar = g_covar)
  }
  w <- LK * prior/dens_g
  w <- w/sum(w)
  
  return(list(g_draws, w, maxest, g_covar))
}

KL <- function (set, par.draws, weights){
  # Probabilities.
  num2 <- tcrossprod(set, par.draws)
  mmat2 <- as.matrix(t(apply(num2, 2, max)))
  numm2 <- exp(sweep(num2, 2, mmat2, FUN = "-"))
  nummax <- exp(sweep(num2, 2, mmat2, FUN = "-"))
  denom <- colSums(numm2)
  ps <- sweep(nummax, 2, denom, FUN = "/")
  lgp <- log(ps)
  wprob <- sweep(ps, 2, weights, FUN="*")
  twp <- rowSums(wprob)
  lgwp <- sweep(lgp, 2, weights, FUN="*")
  tlwp <- rowSums(lgwp)
  #kullback Leibler information
  klinfo <- sum(twp * (log(twp) - tlwp))
  return (as.numeric(klinfo))
}



KLs <- function(full.comb, par.draws, cte.des, cand.set, weights) {
  #take set
  set <- as.matrix(cand.set[as.numeric(full.comb), ])
  #Alternative specific constants 
  set <- cbind(cte.des, set)
  #calculate for all sets the KLinfo.
  kl <- KL(set, par.draws, weights)
  return(kl)
}

Altspec <- function(alt.cte, n.sets) {
  # create matrix
  mat <- diag(length(alt.cte))
  n.zero <- which(alt.cte == 0)
  mat[n.zero, n.zero] <- 0
  # delete zero columns
  del.col <- c(which(apply(mat, 2,   function(x) all(x == 0))))
  mat <- mat[, -del.col]
  #rbind for full design 
  mat <- as.matrix(mat)
  cte.mat <- do.call(rbind, replicate(n.sets, mat, simplify = FALSE)) 
  #return
  return(cte.mat)
}


SeqKL2 <- function(cand.set, full.comb, n.alts, alt.cte, par.draws, weights, reduce = TRUE) {
  # Handling par.draws.
  if (!(is.matrix(par.draws))) {
    par.draws <- matrix(par.draws, nrow = 1)
  }
  # Error alternative specific constants. 
  if (length(alt.cte) != n.alts) {
    stop("n.alts does not match the alt.cte vector")
  }
  # Create alternative specific design.
  cte.des <- Altspec(alt.cte = alt.cte, n.sets = 1)
  # Error handling cte.des
  if (ncol(cand.set) + ncol(cte.des) != ncol(par.draws)) {
    stop("dimension of par.draws does not match the dimension of alt.cte + cand.set.")
  }
  
  # If no weights, equal weights.
  if (is.null(weights)) {
    weights <- rep(1, nrow(par.draws))
  }
  # Calculate KL for each set. 
  kl.infos <- apply(full.comb, 1, KLs, par.draws, cte.des, cand.set, weights)
  # Select maximum.
  comb.nr <- as.numeric(full.comb[which.max(kl.infos), ])
  set <- cand.set[comb.nr, ]
  # Add alternative specific constants if necessary
  if (!is.null(cte.des)) {
    set <- cbind(cte.des, set)
  }
  row.names(set) <- NULL
  # return.
  return(list(set = set, cnr = comb.nr ))
}

















#profiles generation 
profiles<- function (lvls, coding, intercept = FALSE) {
  
  lvls <- lapply(X = as.list(lvls), function(x) (1:x))
  D <- as.data.frame(expand.grid(lvls))
  cn <- names(D)
  D[, cn] <- lapply(D[, cn],  factor)
  con <- list()
  
  for (i in 1:length(lvls)) {
    name <- paste("Var", i, sep = "")
    con[name] <- coding
  }
  
  CD <- as.data.frame(model.matrix(~., D, contrasts = con))
  if (intercept == F) { CD <- CD[, -1]}
  
  return(list(D, as.matrix(CD)))
}

#all posible sets 
full_sets<- function(candset, n_alts){
  
  fun<-function(x){return(1:x)}
  
  s1<-as.list(rep(nrow(candset), n_alts))
  s2<-lapply(X = s1, fun)
  fc<-as.data.frame(expand.grid(s2))
  
  return(fc)
}

#respond function 
respond<-function (par, set, n_alts, bin=TRUE){
  
  par<-as.matrix(par)
  d <- as.matrix(set)
  
  #prob
  U <- d %*% t(par)
  expU <- exp(U)
  p <- expU/sum(expU)
  
  #choice
  choice<-findInterval(x=runif(1), vec=c(0,cumsum(p)))
  
  Y<-rep(0,length(p))
  Y[choice] <- 1
  
  #return
  ifelse(bin, return(Y), return(choice))
  
}

###### convert binary respons to discrete 
transf<- function(Y){
  
  D<-matrix(data = NA, ncol= ncol(Y)/2, nrow= nrow(Y))
  
  for (r in 1:nrow(Y)){
    for (c in 1:ncol(D)){ifelse(all(Y[r ,(c*2-1):(c*2)] == c(0,1)) ,  D[r, c] <- 2  ,  D[r, c] <- 1)}
  }
  return(D)
}

###lattice draws
lattice <- function (K, b, m){
  
  
  base <- function(num){
    
    a1 <- c1 <- rep(0, m)
    a1[m] <- c1[m] <- num %% b
    
    for(j in seq(m - 1, 1, by = -1)){
      tem <- (num - c1[j + 1]) / (b^(m - j))
      a1[j] <- tem %% b
      c1[j] <- c1[j + 1] + a1[j] * (b^(m - j))
    }
    
    return(a1)
  }
  
  a <- 1571
  N <- b^m # number of lattice points
  u <- runif(K)
  av <- rep(1, K)
  
  for(i in 2:K) {av[i] <- (a * av[i - 1]) %% N}
  
  e <- matrix(0, N, K)
  seq <- rep(0, m)
  kk <- -m
  
  for(k in 1:m) {seq[k] <- b^kk; kk <- kk + 1;}
  
  for(i in 1:N){
    ei <- c(crossprod(seq, base(i - 1))) * av + u
    e[i, ] <- ei - floor(ei)
  }
  
  latt <- matrix(0, N, K)
  for(i in 1:K){
    for(j in 1:N){latt[j, i] <- 1 - abs(2 * e[j, i] - 1)}
  }
  
  latt <- qnorm(latt)
  
  return(latt)
}

###lattice mvt 
lattice_mvt<- function (mean, cvar, df, m, b=2){
  
  dim<-length(mean)
  
  #gen lattice from standard normal
  lat<- lattice(K = dim, b, m)
  
  mean <- t(mean)
  X <- matrix(NA, nrow(lat), dim)
  A <- chol(cvar)
  
  for(i in 1:nrow(lat))
  {
    
    inv<-rgamma(1, df/2)
    invw<-inv/(df/2)
    W<-1/invw
    
    Z<-lat[i,]
    
    r<-mean + sqrt(W) * (Z %*% t(A))
    X[i, ]<-t(r)
  }
  
  return (X)
}


