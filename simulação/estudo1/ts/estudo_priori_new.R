library(rstan)
options( mc.cores = 4 )

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
                                     v ~ gamma(12, 0.8);
                                     
                                     //--- Sampling volatilitys:
                                     for(t in 1:T) l[t] ~ gamma(0.5 * v, 0.5 * v);
                                    
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
                                   for(t in 1:T) l[t] ~ gamma(0.5 * v, 0.5 * v);
                                  
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
                                   for(t in 1:T) l[t] ~ gamma(0.5 * v, 0.5 * v);
                                  
                                   //--- Sampling observations:
                                   y[1] ~ normal( b0 + b1 * y0 + b2 * exp( h[1] ), exp( 0.5 * h[1] ) / sqrt( l[1] ) );
                                   for(t in 2:T) y[t] ~ normal( b0 + b1 * y[t-1] + b2 * exp(h[t]), exp( 0.5 * h[t] ) / sqrt( l[t] ) );
                                 }' )

estudo1 = function(M = 500,
                   warmup = 500,
                   r = 2,
                   tails = c( 5, 10, 15, 20 ) ){
  
  v1 = v2 = v3 = rep(0, r)
  summary = data.frame()
  
  for( v in tails ){
    if( v == tails[1] ) time.init = Sys.time()
    for( i in 1:r ){
      
      cat( paste0('réplica ', i, ' com o parâmetro v = ', v ) )
      mu = 1.0
      phi = 0.985
      sigma = 0.13
      b0 = 0.01
      b1 = 0.1
      b2 = -0.02
      y0 = 0
      v = v #( log(v) )
      T = 2e3
      y = h = l = rep(0, T)
      
      data_seed = sample( 1:1e6, 1 )
      set.seed( data_seed )
      for(t in 1:T){
        if(t == 1){
          l[t] = rgamma(1, shape = v/2, rate = v/2)
          h[t] = rnorm(1, mean = mu, sd = sigma * 1 / sqrt( (1 - phi * phi) ) )
          y[t] = rnorm(1, b0 + b1 * y0 + b2 * exp( h[t] ), exp(h[t]/2) / sqrt( l[t] ))
        }else{
          l[t] = rgamma(1, shape = v/2, rate = v/2)
          h[t] = rnorm(1, mean = (mu + phi * ( h[t-1] - mu )), sd = sigma)
          y[t] = rnorm(1, b0 + b1 * y[t-1] + b2 * exp(h[t]), exp(h[t]/2) / sqrt( l[t] ))
        }
      }
      
      fit1 = sampling(model1, 
                      list(T = T, 
                           y = y,
                           y0 = 0,
                           mu = 1.0,
                           phi = 0.985,
                           sigma = 0.13,
                           h = h,
                           l = l,
                           b0 = 0.01,
                           b1 = 0.1,
                           b2 = -0.02),
                      iter = warmup + M,
                      warmup = warmup,
                      thin = 1,
                      chains = 4 )
      v1[ i ] = mean( extract( fit1 )$v )
      
      fit2 = sampling(model2, 
                      list(T = T, 
                           y = y,
                           y0 = 0,
                           mu = 1.0,
                           phi = 0.985,
                           sigma = 0.13,
                           h = h,
                           l = l,
                           b0 = 0.01,
                           b1 = 0.1,
                           b2 = -0.02),
                      iter = warmup + M,
                      warmup = warmup,
                      thin = 1,
                      chains = 4 )
      v2[ i ] = mean( extract( fit2 )$v )
      
      fit3 = sampling(model3, 
                      list(T = T, 
                           y = y,
                           y0 = 0,
                           mu = 1.0,
                           phi = 0.985,
                           sigma = 0.13,
                           h = h,
                           l = l,
                           b0 = 0.01,
                           b1 = 0.1,
                           b2 = -0.02),
                      iter = warmup + M,
                      warmup = warmup,
                      thin = 1,
                      chains = 4 )
      v3[ i ] = mean( extract( fit3 )$v )
      
      cat( '\r' )
      
    }
    
    vies = matrix( c(vies1 = mean( v1 - v ),
                     vies2 = mean( v2 - v ),
                     vies3 = mean( v3 - v )), ncol = 1 )
    
    smse = matrix( c(smse1 = mean( (v1 - v)**2 ),
                     smse2 = mean( (v2 - v)**2 ),
                     smse3 = mean( (v3 - v)**2 )), ncol = 1)
    
    data = cbind(vies, smse, v)
    data = data.frame( data )
    data = round( data, digits = 3 )
    summary = rbind(summary, data)
    if( v == tails[ length(tails) ] ) time.final = Sys.time()
    
  } 
  
  #row.names( summary ) = rep(c('priori1', 'priori2','priori3'), length(tails))  
  colnames( summary ) = c( 'vies', 'smse', 'v' )
  return( list(summary = summary, time = time.final - time.init) )
}

estudo1()

