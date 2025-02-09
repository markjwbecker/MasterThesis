data {
	int<lower=0> N;         // No. of time periods observed
	int<lower=0> K;         // No. of variables in y_t
	int<lower=0> P;         // Lag order VAR
	int<lower=0> Q;         // Lag order VMA
	int<lower=0> PF;        // Number of free VAR
	int<lower=0> QF;        // Number of free VMA
	int<lower=0> KF;        // Number of free VMA
	vector[K] Y[N];         // Data
	matrix[K, K] AR[P];     // Restriction matrix for VAR
	matrix[K, K] BR[Q];     // Restriction matrix for VMA
	matrix[K, K] KSIR;      // Restriction matrix for KSI
	                        // 0 Set to zero 
	                        // 1 Set to one
	                        // 2 Estimate
}
parameters {
	vector<lower=0>[PF] AF;                // Vector of free VAR parameters
	vector<lower=0>[QF] BF;                // Vector of free VMA parameters
	vector<lower=0>[KF] KSIF;              // Vector of free KSE parameters
	cholesky_factor_corr[K] L_corr_noise;  // Cholesky of correlation matrix for errors
	vector<lower=0>[K] sd_noise;           // Sd of errors
	vector[K] intercept;                   // Expected value of y (mu not c)
}
transformed parameters {
	matrix[K,K] L_sigma;                   // Cholesky of covariance matrix
  L_sigma = diag_pre_multiply(sd_noise, L_corr_noise);
}
model {
	vector[K] mus[N];                     // Mu
  vector[K] YM[N];                      // Centered Data
  vector[K] err[N];                     // Residuals 
  matrix[K, K] A[P];                    // VAR matrices 
  matrix[K, K] B[Q];                    // VMA matrices 
  matrix[K, K] KSI;                     // KSI matrix
  int kk;
  
  // Put the free VAR parameters in the corrsponding position in VAR matrices
  kk = 1;
	for (p in 1:P){
	  for (i in 1:K){
  	  for (j in 1:K){
  	    if(AR[p,i,j]!=2){
  	      A[p,i,j] = AR[p,i,j];
  	    } else {
  	      A[p,i,j] = AF[kk];
  	      kk = kk + 1;
  	    }
  	  }
	  }
	}
	
	// Put the free VMA parameters in the corrsponding position in VMA matrices
	kk = 1;
	for (q in 1:Q){
	  for (i in 1:K){
  	  for (j in 1:K){
  	    if(BR[q,i,j]!=2){
  	      B[q,i,j] = BR[q,i,j];
  	    } else {
  	      B[q,i,j] = BF[kk];
  	      kk = kk + 1;
  	    }
  	  }
	  }
	}
	
	// Put the free KSI parameters in the corrsponding position in KSI matrix
  kk = 1;
  for (i in 1:K){
	  for (j in 1:K){
	    if(KSIR[i,j]!=2){
	      KSI[i,j] = KSIR[i,j];
	    } else {
	      KSI[i,j] = KSIF[kk];
	      kk = kk + 1;
	    }
	  }
  }

  // Center Data as intercept is unconditional mean
  for (n in 1:N)
	  YM[n] = Y[n] - intercept;
	// Calculate conditional forecast for centered data
	for (n in 1:N) {
	  for ( k in 1:K)
	    mus[n,k] = 0;
		for (p in 1:min(P,n-P)) 
			mus[n] = mus[n] + A[p] * YM[n-p];
		for (q in 1:min(Q,n-Q)) 
			mus[n] = mus[n] + B[q] * err[n-q];
	  err[n] = YM[n] - 	mus[n];
	}
	// Priors...
	L_corr_noise ~ lkj_corr_cholesky(2.0);
	sd_noise ~ normal(0,1);
	intercept ~ normal(0,1);
	AF ~ normal(0,1);
	BF ~ normal(0,1);
	KSIF ~ normal(0,1);
  // Likelihood
	YM ~ multi_normal_cholesky(mus,KSI*L_sigma);
}
generated quantities {
  // Define bew output to simplify 
	matrix[K,K] Sigma;
	matrix[K,K] Aout[P];
	matrix[K,K] Bout[Q];
	matrix[K,K] KSIout;
	int kk;
	Sigma = L_sigma * L_sigma';
	kk = 1;
	for (p in 1:P){
	  for (i in 1:K){
  	  for (j in 1:K){
  	    if(AR[p,i,j]!=2){
  	      Aout[p,i,j] = AR[p,i,j];
  	    } else {
  	      Aout[p,i,j] = AF[kk];
  	      kk = kk + 1;
  	    }
  	  }
	  }
	}
	
	kk = 1;
	for (q in 1:Q){
	  for (i in 1:K){
  	  for (j in 1:K){
  	    if(BR[q,i,j]!=2){
  	      Bout[q,i,j] = BR[q,i,j];
  	    } else {
  	      Bout[q,i,j] = BF[kk];
  	      kk = kk + 1;
  	    }
  	  }
	  }
	}

	kk = 1;
  for (i in 1:K){
	  for (j in 1:K){
	    if(KSIR[i,j]!=2){
	      KSIout[i,j] = KSIR[i,j];
	    } else {
	      KSIout[i,j] = KSIF[kk];
	      kk = kk + 1;
	    }
	  }
  }
	
}
