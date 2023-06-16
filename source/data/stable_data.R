stable_data = function(mu, phi, sigma,
                       b0, b1, b2,
                       y0,
                       a, #( log(v) )
                       T,
                       seed = NULL){
  #library( stabledist )
  # er
  #a = 1.85    # a e ( 0, 2 ] 
  #T = 2e3
  y = h = rep(0, T)
  if( is.null(seed) ) set.seed( sample(1:1e6,1) ) 
  set.seed( seed )
  for(t in 1:T){
    if(t == 1){
      h[t] = rnorm(1, mean = mu, sd = sigma * 1 / sqrt( (1 - phi * phi) ) )
      y[t] = rstable(n = 1, 
                     alpha = a, 
                     beta = 0, 
                     gamma = exp(h[t]/2),
                     delta = b0 + b1 * y0 + b2 * exp( h[t] ) 
      )
    }else{
      h[t] = rnorm(1, mean = (mu + phi * ( h[t-1] - mu )), sd = sigma)
      y[t] = rstable(n = 1, 
                     alpha = a, 
                     beta = 0, 
                     gamma = exp(h[t]/2),
                     delta = b0 + b1 * y[t-1] + b2 * exp( h[t] ) 
      )
    }
  }
  return( y )
}
