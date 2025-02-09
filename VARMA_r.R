cat("\014") 
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(MTS)
da=read.table("q-gdp-ukcaus.txt",header=T)
gdp <- log(da[,3:5])
zt <- diffM(gdp) # Txm datamatris
colnames(zt) <- c("uk","ca","us")

kdx <- Kronid(zt)$index
spec <- Kronspec(kdx)

Phi_spec <- spec$PhiID
Theta_spec <- spec$ThetaID

fit1 <- Kronfit(zt,kdx) #frekventistiska estimaten

N <- nrow(zt)
K <- 3
P <- 1
Q <- 1


KSIR <- Phi_spec[1:K,1:K]
AR <- array(0,dim=c(P,K,K))
if(P>0){
  for(i in 1:P){
    AR[i,,] <- Phi_spec[,((i-1)*K+1):(i*K)+K]
  }
} 

BR <- array(0,dim=c(Q,K,K))
if(Q>0){
  if(Q>0){
    for(i in 1:Q){
      BR[i,,] <- Theta_spec[,((i-1)*K+1):(i*K)+K]
    }
  }
}


PF <- sum(AR==2)
QF <- sum(BR==2)
KF <- sum(KSIR==2)

varma_data <- list(N = N,      # Sample size
                   K = K,      # number of variables
                   P = P,      # Lags for VAR
                   PF=PF,      # Number of free VAR parameters
                   Q = Q,      # Lags for VMA
                   QF=QF,      # Number of free VMA parameters
                   KF=KF,      # Number of free KSI parameters
                   Y = zt,     # Data
                   KSIR=KSIR,  # Restriction matrix KSI
                   AR=AR,      # Restriction matrix VAR
                   BR=BR)      # Restriction matrix VMA

fit <- stan(
  file = "varma2.stan", # Stan program
  data = varma_data,    # named list of data
  chains = 4,           # number of Markov chains
  warmup = 10000,       # number of warmup iterations per chain
  iter = 40000,         # total number of iterations per chain
  cores = 4             # number of cores (could use one per chain)
  # refresh = 0         # no progress shown
)

print(fit)

plot(fit,pars=c("intercept","Aout","Bout","Sigma"))
