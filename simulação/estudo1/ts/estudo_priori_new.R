library(rstan)
options( mc.cores = 4 )

path1 = 'C:/Users/8936381/Documents/rstan/Estudo_simulação/Estudo1/ts/ts_priori1.stan'
model1 = stan_model( path1 )

path2 = 'C:/Users/8936381/Documents/rstan/Estudo_simulação/Estudo1/ts/ts_priori2.stan'
model2 = stan_model( path2 )

path3 = 'C:/Users/8936381/Documents/rstan/Estudo_simulação/Estudo1/ts/ts_priori3.stan'
model3 = stan_model( path3 )

M = 500
warmup = 500
r = 10
v1 = v2 = v3 = rep(0, r)

for( v in c(5, 10, 15, 20) ){
  for( i in 1:r ){
    
    cat( paste0('réplica ', i, ' com o parâmetro v = ', v ) )
    
    if( i == 1 ) time.init = Sys.time()
    
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
    if( i == r ) time.final = Sys.time()
    
  }
} 
  
  
  
  


# Total time
time.final - time.init

vies = matrix( c(vies1 = mean( v1 - v ),
                 vies2 = mean( v2 - v ),
                 vies3 = mean( v3 - v )), ncol = 1 )

smse = matrix( c(smse1 = mean( (v1 - v)**2 ),
                 smse2 = mean( (v2 - v)**2 ),
                 smse3 = mean( (v3 - v)**2 )), ncol = 1)

summary = cbind(vies, smse)
summary = data.frame( summary )
row.names( summary ) = c('priori1', 'priori2','priori3')
colnames( summary ) = c( 'vies', 'smse')
round( summary, digits = 3 )

