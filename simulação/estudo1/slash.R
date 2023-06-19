library(rstan)
options( mc.cores = 4 )
source( 'https://raw.githubusercontent.com/holtz27/rbras/main/source/slash_estudo1.R' )

model1 = stan_model( model_code = 'data{
  int<lower=0> T;
  real y0;
  vector[T] y;
  real mu;                         // mean log volatility
  real<lower=-1,upper=1> phi;      // persistence of volatility
  real<lower=0> sigma;
  vector[T] h;                     // std log volatility time t
  vector<lower=0, upper=1>[T] l;
  real<lower=-1,upper=1> b1;
  real b0;
  real b2;
}
parameters{
  real<lower=0> v;
}
model{
  // --- prioris
  v ~ gamma(2.0, 0.25);
  
  //--- Sampling volatilitys:
    for(t in 1:T) l[t] ~ beta( v, 1 );
  
  //--- Sampling observations:
    y[1] ~ normal( b0 + b1 * y0 + b2 * exp( h[1] ), exp( 0.5 * h[1] ) / sqrt( l[1] ) );
  for(t in 2:T) y[t] ~ normal( b0 + b1 * y[t-1] + b2 * exp(h[t]), exp( 0.5 * h[t] ) / sqrt( l[t] ) );
}' )

model2 = stan_model( model_code = 'data{
  int<lower=0> T;
  real y0;
  vector[T] y;
  real mu;                         // mean log volatility
  real<lower=-1,upper=1> phi;      // persistence of volatility
  real<lower=0> sigma;
  vector[T] h;                     // std log volatility time t
  vector<lower=0, upper=1>[T] l;
  real<lower=-1,upper=1> b1;
  real b0;
  real b2;
}
parameters{
  real<lower=0> v;
}
model{
  // --- prioris
  target += 0.5 * log( v /( v+3 ) ) + 0.5 * log( trigamma( 0.5 * v ) - trigamma( 0.5 * (v+1) ) - 2 * ( v + 3 ) / ( v * square(v+1) ) );
  
  //--- Sampling volatilitys:
    for(t in 1:T) l[t] ~ beta( v, 1 );
  
  //--- Sampling observations:
    y[1] ~ normal( b0 + b1 * y0 + b2 * exp( h[1] ), exp( 0.5 * h[1] ) / sqrt( l[1] ) );
  for(t in 2:T) y[t] ~ normal( b0 + b1 * y[t-1] + b2 * exp(h[t]), exp( 0.5 * h[t] ) / sqrt( l[t] ) );
}' )

model3 = stan_model( model_code = 'data{
  int<lower=0> T;
  real y0;
  vector[T] y;
  real mu;                         // mean log volatility
  real<lower=-1,upper=1> phi;      // persistence of volatility
  real<lower=0> sigma;
  vector[T] h;                     // std log volatility time t
  vector<lower=0, upper = 1>[T] l;
  real<lower=-1,upper=1> b1;
  real b0;
  real b2;
}
parameters{
  real<lower=0.02,upper=0.5> lambda;
  real<lower=0> v; 
}
model{
  // --- prioris
  lambda ~ uniform(0.02, 0.5);
  v ~ exponential( lambda );
  
  //--- Sampling volatilitys:
    for(t in 1:T) l[t] ~ beta( v, 1 );
  
  //--- Sampling observations:
    y[1] ~ normal( b0 + b1 * y0 + b2 * exp( h[1] ), exp( 0.5 * h[1] ) / sqrt( l[1] ) );
  for(t in 2:T) y[t] ~ normal( b0 + b1 * y[t-1] + b2 * exp(h[t]), exp( 0.5 * h[t] ) / sqrt( l[t] ) );
}' )

r = 200
set.seed( 852465 )
seeds = sample(1:1e6, r)
x = slash_estudo1( M = 500,
                   warmup = 500,
                   model1 = model1, model2 = model2, model3 = model3,
                   r = r,
                   tails = c( 1, 1.5, 2, 2.5, 3 ),
                   mu = 1.0,
                   phi = 0.985,
                   sigma = 0.15,
                   b0 = 0.2,
                   b1 = 0.03,
                   b2 = -0.05,
                   y0 = 0,
                   T = 2e3,
                   seeds = seeds )

#save(x, file = 'slash_estudo1.RData')
load(file = 'slash_estudo1.RData')
summary = x$summary
summary[2:4, ] = round( apply( x$summary[2:4, ], MARGIN = 2, as.numeric ), 5 )
summary

par(mfrow=c(3,1))
boxplot(x$ess1[1, ], x$ess1[2, ], x$ess1[3, ], x$ess1[4, ], x$ess1[5, ],
        names = c('v1','v2','v3','v4','v5'), main = 'priori1')

boxplot(x$ess2[1, ], x$ess2[2, ], x$ess2[3, ], x$ess2[4, ], x$ess2[5, ],
        names = c('v1','v2','v3','v4','v5'), main = 'priori2')

boxplot(x$ess3[1, ], x$ess3[2, ], x$ess3[3, ], x$ess3[4, ], x$ess3[5, ],
        names = c('v1','v2','v3','v4','v5'), main = 'priori3')
par(mfrow=c(1,1))

par(mfrow=c(3,1))
boxplot(x$time1[1, ], x$time1[2, ], x$time1[3, ], x$time1[4, ], x$time1[5, ],
        names = c('v1','v2','v3','v4','v5'), main = 'priori1')

boxplot(x$time2[1, ], x$time2[2, ], x$time2[3, ], x$time2[4, ], x$time2[5, ],
        names = c('v1','v2','v3','v4','v5'), main = 'priori2')

boxplot(x$time3[1, ], x$time3[2, ], x$time3[3, ], x$time3[4, ], x$time3[5, ],
        names = c('v1','v2','v3','v4','v5'), main = 'priori3')
par(mfrow=c(1,1))




















