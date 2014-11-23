// CONSTRAINED LINEAR REGRESSION OF 3 SEGMENTS with 1 INTERVAL
// V14_1a: BUILDS ON 12_1d
// MODEL UNCERTAINTY: INTRINSIC (unknown), YOBS (linear, bounded), XOBS (linear, bounded)
// DISTRIBUTIONS
// CHI ~ MN(chi_theta,chi_mu,chi_sd) params specified
// NU  ~ N(Ymodel,sdI) T[0,inf]
// Yobs~ N(nu, sdy)  	sdy = f(nu) = sdy_min + ay*nu
// Xobs~ N(chi,sdx)   sdx = f(chi)= sdx_min + ax*chi

data {
	int<lower=0> 	N;   					// number of cases
	vector[N] 		xobs;     		// observed predictor (covariate)
	vector[N] 		yobs;     		// observed outcome (variate)
	real	 				s2lim[2];			// bounds for S2 [S2_lower, S2_upper]
	real<lower=0> oy[3];				// yobs uncert [oy_min, oy_mult_min, oy_mult_max]
	real<lower=0> ox[3];				// xobs uncert [ox_min, ox_mult_min, ox_mult_max]
	real<lower=0>	oI[2]; 				// [lower, upper] bounds for intrinsic uncertainty of model
	real 					chi_mu[2];		// mean of chi (true x) distribution
  real<lower=0> chi_sig[2];		// sd of chi distrib
  real<lower=0> chi_theta[2];	// mixing proporiton of chi distrib	
  real 					PETlim; 			// max PET, limits ETo and other vars
}
transformed data{
	real			Xmax; 						// max x value
	real 			Ymax;
	Xmax <- max(xobs);
	Ymax <- max(yobs);
}

parameters { 	// PARAMETERS ARE ASSIGNED UNIFORM PROB ON VALID INTERVAL UNLESS SPECIFIED OTHERWISE IN MODEL
	real<lower=s2lim[1],upper=s2lim[2]>			S2;			// slope of segment 2: constrained by input s2lims
	real<lower=0> 													K1;			// x location of 1st segment break, constrained by [0, max(obs)]
	real<lower=K1,upper=PETlim>							ETo;		// x intercept of 3rd segment	
	real<lower=ETo-K1,upper=(ETo-K1)/(1-((s2lim[1])/S2))>	 	dK;			// span from K1 to 2nd segment break; limit by lower bound S2
	vector<lower=0,upper=1.5*Ymax>[N]  			nu; 		// true, error-free values of y (unknown)
	vector<lower=0,upper=1.5*Xmax>[N]				chi;		// true, error-free values of x (unknown)
	real<lower=oI[1],upper=oI[2]> 					sdI;		// intrinsic uncertainty of model
	real<lower=ox[2],upper=ox[3]>						ax;			// x uncertainty multiplier range - stddev/unit
	real<lower=oy[2],upper=oy[3]>						ay;			// y uncertainty multiplier range
	real<lower=0,upper=10>									A;			// linear segment knot transition parameter, larger = sharper
}

model {
	real 			mu_nu;							// temp value of mu_nu
	real 			pchi[2];						// temp var for log component densities
	real 			S1;									// temp var derived each MCMC iteration
	real 			K2;									// temp var derived each MCMC iteration
	
	// NON-UNIFORM PRIORS
	sdI ~ cauchy(oI[1],100);
	A 	~ lognormal(1,1);
	// broad, uninformative priors on main parameters
	K1  ~ cauchy(0,2.5*Xmax);
	ETo ~ cauchy(K1,2.5*Xmax);
	dK  ~ cauchy(ETo-K1,2.5*Xmax);
	
	K2 <- K1 + dK;
	S1 <- S2*(K2-ETo)/dK;
	// loop through each observation t
	for (t in 1:N){
		// chi ~ MixedNormal(chi_theta | chi_mu,chi_sig)
		for (g in 1:2) {
			pchi[g] <- log(chi_theta[g]) + normal_log(chi[t],chi_mu[g],chi_sig[g]);
		}
		increment_log_prob(log_sum_exp(pchi));
		
		// modeled nu: inv_logit() is analytic approx to Heaviside fcns
		mu_nu <- ( chi[t]-K1 ) *S1 *inv_logit( A*(chi[t]-K1) ) + (chi[t]-K2) *(S2-S1) *inv_logit( A*(chi[t]-K2) );
		nu[t] ~ normal(mu_nu, sdI) T[0,]; 					// add noise to modeled
	}	
	yobs ~ lognormal(log(nu), 	ay); 							// yobs with noise
	xobs ~ lognormal(log(chi), 	ax);							// xobs with noise
}

generated quantities{
	real 				K2;		// derived K2
	real 				S1; 	// derived S1
	real 				mu_nu;
	real 				nu_lp[N];

	K2 <- K1 + dK;
	S1 <- S2*(K2-ETo)/dK;	
	for (t in 1:N){
		mu_nu <- ( chi[t]-K1 ) *S1 *inv_logit( A*(chi[t]-K1) ) + (chi[t]-K2) *(S2-S1) *inv_logit( A*(chi[t]-K2) );
		nu_lp[t] <- normal_log(nu[t],mu_nu,sdI);
	}
}