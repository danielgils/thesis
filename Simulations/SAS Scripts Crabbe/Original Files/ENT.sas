/*
1) Construction of individualized choice designs for the mixed logit choice model with the ENT criterion.
2) Gibbs sampling for the mixed logit choice model.
*/

libname mc 'LOCATION OF DATASETS';

%let nalts = 2;
%let niset = 1;
%let tset = 15;

%let nres = 50;
%let df = 6;
%let dfv = 8;

%let ndraws = 512;
%let nrep = 100;

proc plan ordered;
factors x1=3 x2=3 x3=3/noprint;
output out=candidate;
run;
proc transreg design data=candidate;
model class(x1 x2 x3/effects);
output out=candid;
run;

proc iml WORKSIZE=6000000;
options nonotes;

use mc.betatrue1; read all into truebeta;
use mc.mutrue1; read all into pmean;
use mc.sigmatrue1; read all into pcov;

use candid(keep = &_trgind); read all into candidat; 
do i = 1 to nrow(candidat)-1;
	do j = 1 to nrow(candidat)-i;
    	holddesignfull = holddesignfull // candidat[i, ];
    	holddesignfull = holddesignfull // candidat[i + j, ];
	end;
end;

start choice(beta,Xlast,y1,cset,y);
y=j(cset,1,0);
if cset=1 then y=y;
else do;
y[1:cset-1,1]=y1;
end;
p=j(&nalts,1,0);
ran=uniform(0);
utis=exp(Xlast*beta);
p=utis/sum(utis);
  if ran<=p[1] then y[cset,1]=1;
  else if ran<=p[1]+p[2] then y[cset,1]=2;
  else y[cset,1]=3;
finish;

start RANDMVT(N, DF, Mean, Cov, RANZ);	
	/* check parameters */
	if N<1 then do;
	print "The requested number of observations should be at least 1:" N; stop;
	end; 
	if DF<1 then do;
	print "The degrees of freedome should be at least 1:" DF; stop;
	end;
	mMean = colvec(Mean);
	p = nrow(mMean); 
	X = j(N,p,.);
	*Z = j(p,1,.);
	A = root( Cov );
	inv = 0; 
	do i=1 to N;
		call randgen(inv,'GAMMA',DF/2);   /* Generate inv~ Gamma(DF/2,1)*/
		invW = inv/(DF/2);      /* invW ~ Gamma(DF/2, 2/DF ) */
	 	W = 1/invW;    /* W ~ Inverse Gamma(DF/2, DF/2) */
		*call randgen( Z, 'NORMAL'); 
		Z=RANZ[,i];
		r = mMean + sqrt(W)* A` * Z ; /*A` Z ~ the standard univariate normal(if p = 1) or Multivariate normal(if p>1) */
		X[i, ] = r`;
	end;
	return(X);
finish;

start prob(yn,X,b);
row=nrow(X);
cset=row/&nalts;
p=1;
do ip=1 to cset;
xx=X[(ip-1)*&nalts+1:ip*&nalts,];
	m = max(xx*b);
	u=(xx*b)-m;
  	psub=exp(u)/sum(exp(u));
  	p=p*psub[yn[ip]];
end;
return(p);
finish;

start MVST(n,beta,meanbeta,cov);
dif=beta-meanbeta;
invcov=inv(cov);
diff=t(dif)*invcov*dif;
iMVSTd=1/(det(cov)**(0.5))*(1+((1/n)*diff))**(-(n+&df)/2);
return(iMVSTd);
finish MVST;

/************************main loop*******************************/

mu_rep_all    = j(&df,&nrep,.);
sigma_rep_all = j((&nrep*&df),&df,.);
beta_rep_all  = j((&nrep*&df),&nres,.);

RMSE = j(&nrep,3,.);

do rep=1 to &nrep; 

obs=&nres#&tset#&nalts;
Xmat=j(obs,&df,0);
ymat=j(&nres,&tset,0);
use mc.lat; read all into zfull;

logprior1=-&df/2*log(2*3.1415926)-0.5*log(det(pcov));
prior1=1/((2*3.1415926)**(&df/2)*det(pcov)**(0.5));

do res=1 to &nres; 
holddesign = holddesignfull;
nrowcand=nrow(holddesign);
beta=truebeta[,res]; 

points = (res-1)*&ndraws+1:res*&ndraws;
z = zfull[points,];
z=z`;

chol=t(root(pcov)); 
priordraws=pmean+chol*z;  

priordr=j(&ndraws,1,0); 
do ndr=1 to &ndraws;
priordr[ndr,]=prior1*exp(-0.5*t(priordraws[,ndr]-pmean)*inv(pcov)*(priordraws[,ndr]-pmean));
end;
klpriorsmall=-10000;
nsetcand=nrowcand/&nalts;
do u=1 to nsetcand;
	l=(u-1)*&nalts+1:u*&nalts;
	X=holddesign[l,];
		numer = X*priordraws;
		mmat = numer[<>,];
		numermax = exp(numer - mmat);
		denom = numermax[+,];
		probs = numermax/denom;
		logprobs = log(probs);
		elemfin = probs[,:];
		logelem = probs # logprobs;
		logelemfin = logelem[,:];
		likelielem = probs # log(t(priordr));
		likelielemfin = likelielem[,:];
		klprior = 0;
		do klelem = 1 to &nalts;
		klprior = klprior + ( logelemfin[klelem,] - (log(elemfin[klelem,]) * elemfin[klelem,]) + likelielemfin[klelem,] );
		end;
if klprior > klpriorsmall then;
do;
klpriorsmall = klprior;
startset=X;
end;
end; 
Xorg=startset;

y1=0;
do iset=1 to &niset;
as=(iset-1)*&nalts+1:iset*&nalts;
Xlast=Xorg[as,];
call choice(beta,Xlast,y1,iset,y);
y1=y;
end;
cset=&niset;

do j=1 to &niset;
k=(j-1)*&nalts+1:j*&nalts;
ncands = nrowcand/&nalts;
	do h = 1 to ncands until(found);
	l=(h-1)*&nalts+1:h*&nalts;
	found = (Xorg[k,] = holddesign[l,]);
	end;
idx = setdif(1:nrowcand,l);
holddesign = holddesign[idx,];
nrowcand=nrow(holddesign);
end;
print nrowcand;

do jj=&niset+1 to &tset;

start logpos(b) global(Xorg,y,logprior1,pmean,pcov,cset);
logl=0;
	do si = 1 to cset;
    k = (si-1)#&nalts+1 : si#&nalts;  
   	p = exp(Xorg[k,]*t(b));
	p = p/sum(p);
	logl=logl+log(p[y[si,]]);
   	end;
logprior2=-0.5*(b-t(pmean))*inv(pcov)*(t(b)-pmean);
logpost=logl+logprior1+logprior2;
return(logpost);
finish logpos;

start grd(b) global(Xorg,y,pmean,pcov,cset);
nobs=&nalts*cset;
choicematrix=j(nobs,1,0);
do ire=1 to cset;
ses=(ire-1)*&nalts;
ind=ses+y[ire];
 choicematrix[ind,]=1;
end;
   g = j(1, &df, 0);
   utils = Xorg*t(b);
   exputils = exp(utils); 
       	do j = 1 to cset;
	 	som = 0;
     	ssom = 0; 
	 		do k = (j-1)#&nalts+1 to j#&nalts;
	 		som = som + exputils[k];
      		ssom = ssom + exputils[k] * Xorg[k, ];
      		end;
     		do k = (j-1)#&nalts+1 to j#&nalts;  
      		der = choicematrix[k,] * (Xorg[k, ] - ssom/som);
      		g[, ] = g[, ] + der; 
      		end;  
	 	end;
gprior=inv(pcov)*(pmean-t(b));
gpos=(t(g)+gprior)`;
return (gpos); 
finish grd;

start hes(b) global(Xorg,pcov,cset);
info=j(&df,&df,0);
do scc=1 to cset;
   xx=Xorg[(scc-1)#&nalts+1:scc#&nalts,];
   u=xx*t(b);
   p=exp(u)/sum(exp(u));
   pp=diag(p);
   info=info+(t(xx)*(pp-p*t(p))*xx);
end;
hes=-info-inv(pcov);
return(hes);
finish hes;

b=j(1,&df,0);
bmode = j(1, &df, 0.);
optn = {1, 0};
tc=repeat(.,1,12);
tc[1]=600;
tc[7]=1.e-8;
call nlpnrr(rc,bmode,"logpos",b,optn,,,,,"grd","hes") tc=tc;

covst=-inv(hes(bmode));
prbeta=RANDMVT(&ndraws,&dfv,bmode,covst,z); 
prbeta=prbeta`;
denpi0=j(&ndraws,1,0);
denlikeli=j(&ndraws,1,0);
deng=j(&ndraws,1,0);
wsub=j(&ndraws,1,0);
do ndr=1 to &ndraws;
denpi0[ndr,]=prior1*exp(-0.5*t(prbeta[,ndr]-pmean)*inv(pcov)*(prbeta[,ndr]-pmean));
denlikeli[ndr,]=prob(y,Xorg,prbeta[,ndr]);
deng[ndr,]=MVST(&dfv,prbeta[,ndr],t(bmode),covst);
wsub[ndr,]=denlikeli[ndr,]*denpi0[ndr,]/deng[ndr,];
end;
weight=wsub/wsub[+,];

cset=cset+1;
klinfosmall=-10000;

loglikeliprior = log(t(denlikeli) # t(denpi0));
priorlikeli=j(&ndraws,1,0);
do ndr=1 to &ndraws;
priorlikeli[ndr,]=prob(y,Xorg,priordraws[,ndr]);
end;
term = priorlikeli[:,];
logterm = log(term);

nsetcand=nrowcand/&nalts;
do uu=1 to nsetcand;
	ij=(uu-1)*&nalts+1:uu*&nalts;
	X=holddesign[ij,];
		numer = X*prbeta;
		mmat = numer[<>,];
		numermax = exp(numer - mmat);
		denom = numermax[+,];
		probs = numermax/denom;
		logprobs = log(probs);
		elem = probs # t(weight);
		elemfin = elem[,+];
		logelem = probs # logprobs # t(weight);
		logelemfin = logelem[,+];
		likelielem = probs # loglikeliprior # t(weight);
		likelielemfin = likelielem[,+];
		klinfo = 0;
		do klelem = 1 to &nalts;
		klinfo = klinfo + ( logelemfin[klelem,] - (log(elemfin[klelem,]) * elemfin[klelem,]) + likelielemfin[klelem,] -
								(logterm * elemfin[klelem,]) );
		end;
if klinfo > klinfosmall then;
do;
klinfosmall = klinfo;
Xnewsub=X;
index=ij;
end;
end;
Xorg=Xorg//Xnewsub;
y1=y;
call choice(beta,Xnewsub,y1,cset,y);
idx = setdif(1:nrowcand,index);
holddesign = holddesign[idx,];
nrowcand=nrow(holddesign);
end;

print nrowcand;

y=y`;
nxmatsub=&nalts#&tset;
Xmat[(res-1)*nxmatsub+1:res*nxmatsub,]=Xorg;
ymat[res,]=y;
end;

/*MCMC estimation*/

%let nsample=30000;
%let nburnin=70000;
%let nskip=1;

use mc.betatrue1; read all into betatrue;
use mc.mutrue1; read all into mutrue;
use mc.sigmatrue1; read all into sigmatrue;

X=Xmat;
Y=ymat;
Y = Y`;

start loglikeli(betan,y,X);
logl=j(1,&nres,0);
do res=1 to &nres;
tt=(res-1)*&nalts*&tset+1:res*&nalts*&tset;
Xsub=X[tt,];
row=nrow(Xsub);
nset=row/&nalts;
logl[,res] = 0;
	do s=1 to nset;
	xx=Xsub[(s-1)*&nalts+1:s*&nalts,];
	m = max(xx*betan[,res]);
	u=(xx*betan[,res])-m;
  	p=exp(u)/sum(exp(u));
  	logl[,res] = logl[,res] + log(p[y[s,res],]);    
	end;
end;
return(logl);
finish;

start nextmu(betan,Sigma);
lb=t(root(Sigma)); 
zr = j(&df,1,.);
call randgen(zr,'normal');
munew=betan[,:]+lb*zr/sqrt(&nres);    
return(munew);
finish;

start nextbeta(betan,rho,mu,Sigma,y,X,logLold,bnew,logLnew,newrho,meanind);
lb=t(root(Sigma));
zr = j(&df,&nres,.);
call randgen(zr,'normal');
betanew=betan+sqrt(rho)*(lb*zr);
bn=betanew-repeat(mu,1,&nres);
bo=betan-repeat(mu,1,&nres);
logLnew=loglikeli(betanew,y,X);
varmat1=bn#(inv(Sigma)*bn);
varmat2=bo#(inv(Sigma)*bo);
r=exp(logLnew-logLold-0.5*(varmat1[+,]-varmat2[+,]));
ind=(ranuni(repeat(0,1,&nres))<=r);
meanind=ind[,+]/&nres;
newrho=rho-0.001*(meanind<0.3)+0.001*(meanind>0.3);
logLnew=logLnew#ind+logLold#(1-ind);
bnew=betanew#repeat(ind,&df,1)+betan#repeat((1-ind),&df,1);
finish;

start nextSigma(betan,mu);
s=betan-repeat(mu,1,&nres);
Sigmanew=s*s`+&df*I(&df);
lb=t(root(inv(Sigmanew)));
zr = j(&df,&nres+&df,.);
call randgen(zr,'normal');
ss=lb*zr;
Sigmanew=inv(ss*ss`);
return(Sigmanew);
finish;

/*******starting values********/
mu=j(&df,1,0);
betao=repeat(mu,1,&nres);
Sigma=&df*I(&df);
rho=0.01;

mudraw=j(&df,&nsample,0);
Sigmadraw=j(&df*&nsample,&df,0);
betamean=j(&df,&nres,0);

logLold=loglikeli(betao,Y,X);

ndraw=&nburnin+&nsample*&nskip;
do r=1 to ndraw;
call nextbeta(betao,rho,mu,Sigma,Y,X,logLold,bnew,logLnew,newrho,accept);
betao=bnew;
logLold=logLnew;
rho=newrho;
 mu=nextmu(betao,Sigma);
 Sigma=nextSigma(betao,mu);
nk=(r-&nburnin)/&nskip;
if r>&nskip & int(nk)=nk & nk>0 then do;
mudraw[,nk]=mu;
nnk=(nk-1)*&df+1:nk*&df;
Sigmadraw[nnk,]=Sigma;
betamean=betamean+betao;
end;
end;
betamean=betamean/&nsample;
mumean=mudraw[,:];
sigmamean=j(&df,&df,0);
do i = 1 to &nsample;
k=(i-1)*&df+1:i*&df;
	sigmamean = sigmamean + Sigmadraw[k,];
end;
sigmamean = sigmamean/&nsample;

RMSE_mu = sqrt((t(mumean - mutrue)*(mumean - mutrue))/&df);

nlong = (&df*(&df+1))/2;
sigmatruelong = j(nlong,1,.);
sigmameanlong = j(nlong,1,.);
j2 = 0;
do i = 1 to &df;
	j1 = j2 + 1;
	j2 = j1 + (&df - i);
	sigmatruelong[j1:j2,] = sigmatrue[i:&df,i];
	sigmameanlong[j1:j2,] = sigmamean[i:&df,i];
end; 
RMSE_sigma   = sqrt((t(sigmameanlong - sigmatruelong)*(sigmameanlong - sigmatruelong))/nlong);

betasum = 0;
do m = 1 to &nres;
betasum = betasum + (t(betamean[,m] - betatrue[,m]) * (betamean[,m] - betatrue[,m]));
end;
RMSE_beta = sqrt(betasum / (&nres*&df));

RMSE[rep,1] = RMSE_mu; RMSE[rep,2] = RMSE_sigma; RMSE[rep,3] = RMSE_beta;

mu_rep_all[,rep] = mumean;
k = (rep-1)*&df+1:rep*&df;
sigma_rep_all[k,] = sigmamean;
beta_rep_all[k,] = betamean;

end;

meanRMSE = RMSE[:,];
print meanRMSE;
create mc.RMSE from RMSE; append from RMSE;
create mc.mu from mu_rep_all; append from mu_rep_all;
create mc.sigma from sigma_rep_all; append from sigma_rep_all;
create mc.beta from beta_rep_all; append from beta_rep_all;

options notes;
quit;
