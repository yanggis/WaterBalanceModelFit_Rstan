// CHANGEPOINT DETECTION: 1D, SINGLE CHANGE, GAUSSIAN SPREAD
// 
// This RStan model detects a step change in the location and/or spread of 
// a 1D Gaussian distribution as a function of a latent change point variable
// Tcp.  The location and scale parameters are assumed to be constant before 
// and after Tcp. This is equivalent to a mixed Gaussian where the source 
// distribution changes at time Tcp. An example use is to detect changes in 
// the residual of a fitted 2D model. 
//
// The only trick here is to use the inv_logit() function as an approximation
// to the Heaviside function, with an appropriate scale parameter A that strikes
// a compromise between the continuity needed for the Stan HMC sampler and 
// a quick transition to capture the change effect. The alternative approach of 
// using a categorical variable to specify the source distributions is not 
// supported in Stan. 
//
// (c) TC Moran, 2014
// Distributed under the MIT License.

data {
	int<lower=1> 	N; 				// number of data points
	real 			y[N]; 				// observations
	real 			minIntvl;			// minimum interval length
	real 			Tidx[N];			// time index vector
}
transformed data {
	real dy;
	dy <- max(y)-min(y); 		// used to bound mu and sig priors 
}
parameters {
	real<lower=minIntvl,upper=N-minIntvl>	Tcp;	// Tcp signifies last time index of first interval
	real<lower=min(y),upper=max(y)>  		mu0; 		// first distribution location
	real<lower=-dy,upper=dy>						dmu1;		// change in location after changepoint Tcp 
	real<lower=0,upper=dy> 							sig0; 	// scale of first distribution
	real<lower=-sig0,upper=dy>					dsig1;	// change in scale after Tcp
}
model {
	real H; 			// Heaviside function approximation
	real A; 			// scale param for inv_logit
	real mu;			// temp var for location param
	real sig; 		// temp var for scale param 
	
	A <- 10; 			// appropriate for unit time step, smooth but full transition
	for (n in 1:N) {
		H <- inv_logit(A*(Tidx[n]-Tcp));
		mu <- mu0 + dmu1*H;
		sig<- sig0+dsig1*H;
		y[n] ~ normal(mu,sig);
	}
}