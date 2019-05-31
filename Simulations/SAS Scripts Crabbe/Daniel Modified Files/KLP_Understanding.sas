/*Daniel: I modified and comment this script to understand what has been done
	I created an excel file called Objects_created_KL.xlsx to paste some output that I needed to understand and
	easily manipulate
	Dan: All my comments start with "Dan:"	*/
/*
1) Construction of individualized choice designs for the mixed logit choice model with the KLP criterion.
2) Gibbs sampling for the mixed logit choice model.
*/

libname mc 'C:\Users\danie\Documents\Daniel Gil\KULeuven\Stage 2\Thesis\Simulations\SAS Scripts Crabbe\Daniel Modified Files';

%let nalts = 2;
%let niset = 1;
%let tset = 15;

%let nres = 50;
*%let nres = 2;
%let df = 6;
%let dfv = 8;

*%let ndraws = 512;
%let ndraws = 10;
*%let nrep = 100;
%let nrep = 2;

* Daniel: This proc is to create all possible combinations between attributes/factors;
* Daniel: It looks like this is scenario 1 mentioned in the paper;
proc plan ordered;
factors x1=3 x2=3 x3=3/noprint;
output out=candidate;
run;
* Daniel: This proc is to create columns with effects coding for each attribute/factor;
proc transreg design data=candidate;
model class(x1 x2 x3/effects);
output out=candid;
run;

proc iml WORKSIZE=6000000;
options nonotes;

use mc.betatrue1; read all into truebeta; * Dan: (6x50) Betatrue is the true betas for each of the 50 respondents
												There are six betas because there are 6 variables, 2 for each factor;
use mc.mutrue1; read all into pmean; * Dan: (6x1) Mutru1 is the true mean of the prior distribution;
use mc.sigmatrue1; read all into pcov; * Dan: (6x6) Sigmatrue1 is the true variance (or sd, not sure) of the prior dist.;

* Dan: candid just select the columns from the candidate set (candid) that are represented with effects coding;
use candid(keep = &_trgind); read all into candidat; 

* Dan: this loop creates a table with all possible choice sets of two alternatives that can be created from
	candidat. Here, each pair of observations/rows is a choice set.
	holddesignfull is 702x6; 
do i = 1 to nrow(candidat)-1;
	do j = 1 to nrow(candidat)-i;
	* Dan: // is a concatenation operator, vertical;
	    holddesignfull = holddesignfull // candidat[i, ];
		holddesignfull = holddesignfull // candidat[i + j, ];
	end;
end;

* Dan: NOTE: (IMPORTANT) I moved all modules/functions (begin with start) to the beginning. Before, some of 
		them were inside loops, so each function was defined more than once;
* Dan: start statement is to create functions within proc iml;
* Dan: choice is a function that selects the next choice set based on the estimation of the betas,
	the last choice set, the last response (at least that's what I think so far);
start choice(beta,Xlast,y1,cset,y);
* Dan: Parameters:
		- Beta: Matrix 6x1. Betas for the corresponding respondent
		- Xlast: Matrix 2x6. This is the best candidate set for the respondent
		- y1: Scalar. In the beginning is zero and is updated later
		- iset: Scalar. Is the iteration number correspoding to the loop that runs from 1 to niset
		- y: Scalar. Is not defined previously, but i think is the output;
* Dan: The J function creates a matrix with nrow rows and ncol columns with all elements equal to value 
	[J (nrow <, ncol> <, value> )];
	y=j(cset,1,0);
	if cset=1 then y=y;
	else do;
		y[1:cset-1,1]=y1;
	end;
	p=j(&nalts,1,0); * Dan: Matrix 2x1;
	ran=uniform(0); * Dan: Scalar. Random number (I think it's always zero);
	utis=exp(Xlast*beta); * Dan: Matrix 2x1. Utility of each alternative;
	p=utis/sum(utis); * Dan: Matrix 2x1. Probability of each alternative;
  	if ran<=p[1] then y[cset,1]=1; * Dan: I think y is always going to be 1, because ran is always zero;
  	else if ran<=p[1]+p[2] then y[cset,1]=2;
  	else y[cset,1]=3;
finish;

* Dan: RANDMVT is a function that generates samples from a multivariate normal distribution or a 
	multivariate t distribution, because of the degrees of freedom asked.
	Note 2: Actually it seems that in generates samples from an inverse gamma, which are them summed
	to the input mean vector, Mean, multiplied by a normal distribution.(at least that's what I think so far)
	Note 3: According to the article (pag 14) this function generates samples from a MULTIVARIATE T DISTRIBUTION
	Note 4: from R this function are the lattice points;

start RANDMVT(N, DF, Mean, Cov, RANZ);	
	* Dan: RANZ are the corresponding lattice points;
	/* check parameters */
	if N<1 then do;
		print "The requested number of observations should be at least 1:" N; stop;
	end; 
	if DF<1 then do;
		print "The degrees of freedome should be at least 1:" DF; stop;
	end;
	* Dan: colvec converts a matrix into a column vector.;
	mMean = colvec(Mean);
	p = nrow(mMean);  * Number of parameters;
	X = j(N,p,.);
	*Z = j(p,1,.);
	A = root( Cov ); * Cholesky decomposition;
	inv = 0; 
	do i=1 to N;
	* Dan: randgen generates random numbers from a specified distribution.;
		call randgen(inv,'GAMMA',DF/2);   /* Generate inv~ Gamma(DF/2,1)*/
		invW = inv/(DF/2);      /* invW ~ Gamma(DF/2, 2/DF ) */
	 	W = 1/invW;    /* W ~ Inverse Gamma(DF/2, DF/2) */
		*call randgen( Z, 'NORMAL'); 
		Z=RANZ[,i]; * Dan: matrix 6x10;
		r = mMean + sqrt(W)* A` * Z ; /*A` Z ~ the standard univariate normal(if p = 1) or Multivariate normal(if p>1) */
		X[i, ] = r`;
	end;
	return(X);
finish;

* Dan: prob is a function that computes the probability of the alternative chosen. y is the alternative chosen
	and is selected randomly by choice function;
start prob(yn,X,b);
	* Dan: Parameters:
		- yn: response, but i think is always 1 (scalar)
		- X: is the best candidate set (Matrix 2x6)
		- b: I think is a draw from the posterior of beta (6x1);
	row=nrow(X);
	cset=row/&nalts; * Dan: in this case cset = 1;
	p=1;
	* Dan: For each candidate set,;
	do ip=1 to cset;
		xx=X[(ip-1)*&nalts+1:ip*&nalts,]; * Dan: In this case xx = X (matrix 2x6);
		m = max(xx*b);
		u=(xx*b)-m;
  		psub=exp(u)/sum(exp(u));
  		p=p*psub[yn[ip]];
	end;
	return(p);
finish;

* Dan: MVST is a function that; 
start MVST(n,beta,meanbeta,cov);
	* Dan: Parameters:
		- n: is equal to macro variable dfv. It's 8 in this case
		- beta: I think is a draw from the posterior of beta (6x1)
		- meanbeta: I think is the result from the optimization. Matrix 6x1
		- cov = is the inverse of the hessian of meanbeta. Matrix 6x6;
	dif=beta-meanbeta;
	invcov=inv(cov);
	diff=t(dif)*invcov*dif;
	iMVSTd=1/(det(cov)**(0.5))*(1+((1/n)*diff))**(-(n+&df)/2);
	return(iMVSTd);
finish MVST;

* Dan: logpos is the function to optimize;
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

* Dan: grd is the gradient of logpos;
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
	
* Dan: hes is the hessian of logpos;	
start hes(b) global(Xorg,pcov,cset);
	info=j(&df,&df,0);
	do scc=1 to cset; * Dan: this loop is for each candidate set that has been answered;
		xx=Xorg[(scc-1)#&nalts+1:scc#&nalts,];
		u=xx*t(b);
		p=exp(u)/sum(exp(u));
		pp=diag(p);
		info=info+(t(xx)*(pp-p*t(p))*xx);
	end;
	hes=-info-inv(pcov);
	return(hes);
finish hes;

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

/************************main loop*******************************/
* Dan: These four matrices are exported at the end of the proc;
mu_rep_all    = j(&df,&nrep,.); * Dan: Matrix 6x2 (original dimension is 6x100) 
										Correspond mu estimates in each repetition;
sigma_rep_all = j((&nrep*&df),&df,.); * Dan: Matrix 12x6 (original dimension is 600x6) 
											First 6 rows are Cov matrix of rep 1, and so on;
beta_rep_all  = j((&nrep*&df),&nres,.); * Dan: Matrix 12x50 (original dimension is 600x50)
											First 6 rows are betas for rep1, and so on;;

RMSE = j(&nrep,3,.); * Dan: Matrix 2x3 (original dimension is 100x3);

* Dan: Loop for each repetition;
*do rep=1 to &nrep; 
do rep=1 to 1; 
	/*submit rep;
    %put Repetition= &rep;
    endsubmit;*/

	* Dan: # is a multiplication operator, elementwise;
	obs=&nres#&tset#&nalts; * Dan: obs = 50x15x2 = 1500;
	Xmat=j(obs,&df,0);	* Dan: Matrix 1500x6. I think this is the desing matrix;
	ymat=j(&nres,&tset,0); * Dan: Matrix 50x15;
	use mc.lat; read all into zfull; * Dan: Matrix 25600x6. lat are the lattice points used 
									to sample from the prior;

	logprior1=-&df/2*log(2*3.1415926)-0.5*log(det(pcov)); * Dan: logprior1 = -3.43419. Constant value;
	prior1=1/((2*3.1415926)**(&df/2)*det(pcov)**(0.5));	* Dan: prior1 = 0.0322515. Constant value;

	do res=1 to &nres; 
	*do res=1 to 2; 
		/*submit res;
    	%put Respondeeeeeeeeeeeeeeeent = &res;
    	endsubmit;*/
		* Dan: This section until next loop (u loop) is to get the prior draws for the corresponding respondent;
		holddesign = holddesignfull; * Dan: Copy of holddesignfull 702x6;
		nrowcand=nrow(holddesign); * Dan: nrowcand = 702;
		beta=truebeta[,res]; * Dan: Matrix 6x1. Betas for the corresponding respondent;

		points = (res-1)*&ndraws+1:res*&ndraws; * Dan: Indicates what rows to choose from lattice points for each respondent;
		z = zfull[points,]; * Dan: 10x6. These are the lattice points for the respondent. In this case ndraws=10;
		* Dan: ` is the transpose operator;
		z=z`; * Dan: 6x10;
		
		* Dan: t function returns the transpose of the argument;
		* Dan: root function returns cholesky decomposition;
		chol=t(root(pcov)); * Dan: Matrix 6x6 (same dimensions as pcov). In this case is a diagonal matrix;
		priordraws=pmean+chol*z;  * Dan: Matrix 6x10. Each respondent has a different matrix;
		klpriorsmall=-10000; * Dan: Don't know where this number came from;
		nsetcand=nrowcand/&nalts; * Dan: 702/2 = 351 candidate choice sets;
		do u=1 to nsetcand; * Dan: This loop is for each candidate set. To choose the best next candidate set;
		*do u=1 to 10; * Dan: This loop is for each candidate set;
			/*submit u;
    		%put nsetcand_u = &u;
    		endsubmit;*/
			
			l=(u-1)*&nalts+1:u*&nalts; * Dan: Indicates the rows of each candidate set: (1,2),(3,4),(5,6) and so on;
			X=holddesign[l,]; * Dan: Matrix 2x6. Correspoding candidate set;
			numer = X*priordraws; * Dan: Matrix 2x10 (2x6 * 6x10). Each alternative is multiplied by each prior draw;
			* Dan: <> is the element maximum operator which in this case keeps the maximum in each row;
			mmat = numer[<>,]; * Dan: Matrix 1x10. Takes the maximum in each column/Draw;
			numermax = exp(numer - mmat); * Dan: Matrix 2x10;
			* Dan: + is the sum of rows of the matrix for each column;
			denom = numermax[+,]; * Dan: Matrix 1x10. Sum of rows in each column;
			probs = numermax/denom; * Dan: Matrix 2x10;
			logprobs = log(probs); * Dan: Matrix 2x10;
			elemfin = probs[,:]; * Dan: 2x1. It gives the average for each row; 
			logelemfin = logprobs[,:]; * Dan: 2x1. It gives the average for each row; 
			klprior = 0;
			do klelem = 1 to &nalts; * Dan: This loop is just to compute the sum of that operation in both alts;
				klprior = klprior + ( elemfin[klelem,] * (log(elemfin[klelem,]) - logelemfin[klelem,]) );
			end;
			* Each time a new candidate set gives a better value, it's chosen to continue;
			* If it does give a better value then the candidate set is discarded;
			if klprior > klpriorsmall then do; 
				klpriorsmall = klprior; * Dan: updates klpriorsmall number;
				startset=X; * Dan: Matrix 2x6;
			end;
		end; * Dan: end u loop;
		Xorg=startset; * Dan: Matrix 2x6 This is the best candidate set for the respondent;

		y1=0;
		do iset=1 to &niset; * Dan: I think niset is the number of next sets or something like that. Dont know for sure;
			as=(iset-1)*&nalts+1:iset*&nalts; * Dan: Indicates always the same rows (1,2). I dont know what is niset;
			Xlast=Xorg[as,]; * Dan: Matrix 2x6 This is the best candidate set for the respondent. So far it's the same as Xorg;
			call choice(beta,Xlast,y1,iset,y); * Dan: I think y is always 1;
			y1=y;
		end;
		cset=&niset; * Dan: scalar. Candidate set = 1;

		* Dan: this loop is to find the chosen choice set and remove it from the list of all possible choice sets;
		do j=1 to &niset; * Dan: I think niset is the number of next sets or something like that. Dont know for sure;
			k=(j-1)*&nalts+1:j*&nalts; * Dan: Indicates always the same rows (1,2). I dont know what is niset;
			ncands = nrowcand/&nalts; * Dan: 702/2 = 351 candidate choice sets;
			do h = 1 to ncands until(found); * Dan: This loop stops until it find the best candidate set in the whole
													list of possible candidate sets. For documentation:
											https://newonlinecourses.science.psu.edu/stat481/node/44/;
				l=(h-1)*&nalts+1:h*&nalts; * Dan: Indicates the rows of each candidate set: (1,2),(3,4),(5,6) and so on;
				found = (Xorg[k,] = holddesign[l,]); * Dan: I think this a boolean object to stop the loop;
			end;
			* Dan: setdif function returns all element values present in A but not in B SETDIF (A, B);
			idx = setdif(1:nrowcand,l); * Dan: sequence of values from 1 to nrowcand, removing l;
			holddesign = holddesign[idx,]; * Dan: Because the chosen choice set is used, then it is removed from the
													list of all possible choice sets;
			nrowcand=nrow(holddesign); * Dan: Recompute the number of rows in the possible choice set list;
		end;
		*print nrowcand;
		do jj=&niset+1 to &tset; * Dan: This loop is to repeat the whole process for the remaining choice sets in the survey;
			/*submit cset;
    		%put csetttt = &cset;
    		endsubmit;
			submit jj;
    		%put tset_jj = &jj;
    		endsubmit;*/
			b=j(1,&df,0); * Dan: Matrix 1x6;
			bmode = j(1, &df, 0.); * Dan: Matrix 1x6;
			optn = {1, 0}; * Dan: Column Vector 2x1. These are the options for nlpnrr function: 1 is that the problem
								is to minimize or maximize and zero in the second part means that No printed output 
								is produced. 
					http://support.sas.com/documentation/cdl/en/imlug/65547/HTML/default/viewer.htm#imlug_nonlinearoptexpls_sect018.htm;
			tc=repeat(.,1,12); * Dan: Matrix 1x12 filled with dots;
			tc[1]=600; 
			tc[7]=1.e-8;
			call nlpnrr(rc,bmode,"logpos",b,optn,,,,,"grd","hes") tc=tc; * Function to optimize logpos function
									logpos, grd (gradient) and hes (hessian) are defined above. See help to know what is it for.
									bmode is the result of the optimization (I guess it is, not sure) Matrix 1x6.
				From Sas Help: The row vector xr (bmode), which has length $n$, contains the optimal point when rc $>0$.
							bmode is the posterior mode and it is used as mean in the multivariate t distribution (pag14article);
			
			* Dan: inv function computes the inverse of a matrix;
			covst=-inv(hes(bmode)); * Dan: Matrix 6x6. This is the var-cov matrix for the multivariate t distribution (pag14article);

			prbeta=RANDMVT(&ndraws,&dfv,bmode,covst,z); * Dan: returns a matrix of 10x6 (ndraws, df). 
									Draws from a MULTIVARIATE T DISTRIBUTION, I think these are from the posterior of beta
									According to R these are lattice points;
			prbeta=prbeta`; * Dan: Matrix 6x10;
			denpi0=j(&ndraws,1,0); * Dan: Matrix 10x1;
			denlikeli=j(&ndraws,1,0); * Dan: Matrix 10x1;
			deng=j(&ndraws,1,0); * Dan: Matrix 10x1;
			wsub=j(&ndraws,1,0); * Dan: Matrix 10x1;
			do ndr=1 to &ndraws;
				denpi0[ndr,]=prior1*exp(-0.5*t(prbeta[,ndr]-pmean)*inv(pcov)*(prbeta[,ndr]-pmean)); * Dan: looks like a multivariate normal;
				denlikeli[ndr,]=prob(y,Xorg,prbeta[,ndr]); * Dan: returns a probability for each draw from the posterior
											of beta. Prob function returns the prob of the alternative chosen by y, using the design and 
											the draws from the posterior of betas. Matrix 10x1;
				deng[ndr,]=MVST(&dfv,prbeta[,ndr],t(bmode),covst); * Dan: it makes a calculation but dont know yet;
				wsub[ndr,]=denlikeli[ndr,]*denpi0[ndr,]/deng[ndr,]; * Dan: Scalar;
			end;
			weight=wsub/wsub[+,]; * Dan: Matrix 10x1. Sum of rows in each column;
	
			* Dan: from here until the jj loop ends, is the same as with the first candidate set;
			cset=cset+1; 
			klinfosmall=-10000;
			
			nsetcand=nrowcand/&nalts; * Dan: Recompute the number of rows in the possible choice set list;
			do uu=1 to nsetcand; * Dan: This loop is for each candidate set. To choose the best next candidate set;
				ij=(uu-1)*&nalts+1:uu*&nalts; * Dan: Indicates the rows of each candidate set: (1,2),(3,4),(5,6) and so on;
				X=holddesign[ij,]; * Dan: Matrix 2x6. Correspoding candidate set;
				numer = X*prbeta; * Dan: Matrix 2x10 (2x6 * 6x10). Each alternative is multiplied by each posterior draw;
				* Dan: <> is the element maximum operator which in this case keeps the maximum in each row;
				mmat = numer[<>,]; * Dan: Matrix 1x10. Takes the maximum in each column/Draw;
				numermax = exp(numer - mmat); * Dan: Matrix 2x10;
				* Dan: + is the sum of rows of the matrix for each column;
				denom = numermax[+,]; * Dan: Matrix 1x10. Sum of rows in each column;
				probs = numermax/denom; * Dan: Matrix 2x10;
				logprobs = log(probs); * Dan: Matrix 2x10;
				* Dan: # is a multiplication operator, elementwise;
				elem = probs # t(weight); * Dan: Multiply each prob with corresponding weight. Matrix 2x10;
				elemfin = elem[,+]; * Dan: 2x1. It gives the sum for each row across all columns; 
				logelem = logprobs # t(weight); * Dan: Multiply each 
				with corresponding weight. Matrix 2x10;
				logelemfin = logelem[,+]; * Dan: 2x1. It gives the average for each row; 
				klinfo = 0;
				do klelem = 1 to &nalts; * Dan: This loop is just to compute the sum of that operation in both alts;
					klinfo = klinfo + ( elemfin[klelem,] * (log(elemfin[klelem,]) - logelemfin[klelem,]) );
				end;
				* Each time a new candidate set gives a better value, it's chosen to continue;
				* If it does give a better value then the candidate set is discarded;
				if klinfo > klinfosmall then;
					do;
					klinfosmall = klinfo; * Dan: updates klinfosmall number;
					Xnewsub=X; * Dan: Matrix 2x6;
					index=ij; * Dan: save the position of the best candidate set in the list of possible candsets;
				end;
			end; * Dan: end uu loop;
			Xorg=Xorg//Xnewsub; * Dan: append of new choice set to the existing one. 
								The final matrix has to be 30x6 if 15 choice sets are considered;
			y1=y;
			call choice(beta,Xnewsub,y1,cset,y); * Dan: I think y is always 1;
			* Dan: setdif function returns all element values present in A but not in B SETDIF (A, B);
			idx = setdif(1:nrowcand,index); * Dan: sequence of values from 1 to nrowcand, removing index;
			holddesign = holddesign[idx,]; * Dan: Because the chosen choice set is used, then it is removed from the
													list of all possible choice sets;
			nrowcand=nrow(holddesign); * Dan: Recompute the number of rows in the possible choice set list;
		end;  * Dan: end jj loop;
		*print nrowcand;
		y=y`; * Dan: Matrix 1x15 filled with 1s;
		nxmatsub=&nalts#&tset;
		Xmat[(res-1)*nxmatsub+1:res*nxmatsub,]=Xorg;
		ymat[res,]=y;
		*submit res;
		*%put end_respoooooondeeent = &res;
		*endsubmit;
	end; * Dan: end res loop;
	* Dan: The result of this loop is Xmat, which is a huge matrix 1500x6 with the respective 15 choice sets for each
			respondent (50res x 15choicesets x 2alts = 1500);
	*print res;
	*print Xmat;
	*print ymat;
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
end; * Dan: end rep loop;

meanRMSE = RMSE[:,];
print meanRMSE;
create mc.RMSE from RMSE; append from RMSE;
create mc.mu from mu_rep_all; append from mu_rep_all;
create mc.sigma from sigma_rep_all; append from sigma_rep_all;
create mc.beta from beta_rep_all; append from beta_rep_all;

options notes;
quit;
