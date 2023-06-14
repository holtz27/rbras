library(rstan)
options( mc.cores = 4 )
source( 'https://raw.githubusercontent.com/holtz27/rbras/main/source/estudo1.R' )

model1 = stan_model( model_code = 'data{
                                   int<lower=0> T;
                                   real y0;
                                   vector[T] y;
                                   real mu;                         // mean log volatility
                                   real<lower=-1,upper=1> phi;      // persistence of volatility
                                   real<lower=0> sigma;
                                   vector[T] h;                     // std log volatility time t
                                   vector<lower=0>[T] l;
                                   real<lower=-1,upper=1> b1;
                                   real b0;
                                   real b2;
                                   }
                                   parameters{
                                     real<lower=2> v;
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
                                   vector<lower=0>[T] l;
                                   real<lower=-1,upper=1> b1;
                                   real b0;
                                   real b2;
                                 }
                                 parameters{
                                   real<lower=2> v;
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
                                   vector<lower=0>[T] l;
                                   real<lower=-1,upper=1> b1;
                                   real b0;
                                   real b2;
                                 }
                                 parameters{
                                   real<lower=0.02,upper=0.5> lambda;
                                   real<lower=2> v; 
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

x = estudo1( M = 500,
             warmup = 500,
             model1 = model1, model2 = model2, model3 = model3,
             r = 10,
             tails = c( 5, 10, 15, 20 ),
             mu = 1.0,
             phi = 0.985,
             sigma = 0.15,
             b0 = 0.2,
             b1 = 0.03,
             b2 = -0.05,
             y0 = 0,
             T = 2e3 )

