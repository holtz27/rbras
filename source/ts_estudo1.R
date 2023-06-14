ts_estudo1 = function(M,
                   warmup,
                   model1, model2, model3,
                   r,
                   tails,
                   mu,
                   phi,
                   sigma,
                   b0,
                   b1,
                   b2,
                   y0,
                   T){
  
  v1 = v2 = v3 = rep(0, r)
  
  for( v in tails ){
    if( v == tails[1] ) time.init = Sys.time()
    for( i in 1:r ){
      
      cat( paste0('réplica ', i, ' com o parâmetro v = ', v ) )
      
      # Data generation
      y = h = l = rep(0, T)
      data_seed = sample( 1:1e6, 1 )
      set.seed( data_seed )
      for(t in 1:T){
        if(t == 1){
          l[t] = rgamma(1, shape = 0.5 * v, rate = 0.5 * v)
          h[t] = rnorm(1, mean = mu, sd = sigma * 1 / sqrt( (1 - phi * phi) ) )
          y[t] = rnorm(1, b0 + b1 * y0 + b2 * exp( h[t] ), exp(h[t]/2) / sqrt( l[t] ))
        }else{
          l[t] = rgamma(1, shape = 0.5 * v, rate = 0.5 * v)
          h[t] = rnorm(1, mean = (mu + phi * ( h[t-1] - mu )), sd = sigma)
          y[t] = rnorm(1, b0 + b1 * y[t-1] + b2 * exp(h[t]), exp(h[t]/2) / sqrt( l[t] ))
        }
      }
      
      # models fit
      fit1 = sampling(model1, 
                      list(T = T, 
                           y = y,
                           y0 = y0,
                           mu = mu,
                           phi = phi,
                           sigma = sigma,
                           h = h,
                           l = l,
                           b0 = b0,
                           b1 = b1,
                           b2 = b2),
                      iter = warmup + M,
                      warmup = warmup,
                      thin = 1,
                      chains = 4 )
      v1[ i ] = mean( extract( fit1 )$v )
      
      fit2 = sampling(model2, 
                      list(T = T, 
                           y = y,
                           y0 = y0,
                           mu = mu,
                           phi = phi,
                           sigma = sigma,
                           h = h,
                           l = l,
                           b0 = b0,
                           b1 = b1,
                           b2 = b2),
                      iter = warmup + M,
                      warmup = warmup,
                      thin = 1,
                      chains = 4 )
      v2[ i ] = mean( extract( fit2 )$v )
      
      fit3 = sampling(model3, 
                      list(T = T, 
                           y = y,
                           y0 = y0,
                           mu = mu,
                           phi = phi,
                           sigma = sigma,
                           h = h,
                           l = l,
                           b0 = b0,
                           b1 = b1,
                           b2 = b2),
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
    vies = round( vies, digits = 3 )
    smse = matrix( c(smse1 = mean( (v1 - v)**2 ),
                     smse2 = mean( (v2 - v)**2 ),
                     smse3 = mean( (v3 - v)**2 )), ncol = 1)
    smse = round( smse, digits = 3 )
    
    if( v == tails[1] ){
      summary = data.frame(vies, smse)
      summary = rbind(c('v', v), summary)
    }else{
      data = data.frame(vies, smse)
      data = rbind(c('v', v), data)
      summary = cbind(summary, data)
    }
    
    
    if( v == tails[ length(tails) ] ) time.final = Sys.time()
    
  } 
  
  row.names( summary ) = c('','priori1', 'priori2','priori3' )
  
  return( list(summary = summary, time = time.final - time.init) )
}
