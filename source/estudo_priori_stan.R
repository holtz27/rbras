estudo_priori_stan = function(y0,
                              model,
                              v,
                              M,
                              warmup,
                              nchains,
                              lags,
                              r, 
                              seeds){
  
  source( 'https://raw.githubusercontent.com/holtz27/rbras/main/source/ts_data.R' )
  source( 'https://raw.githubusercontent.com/holtz27/rbras/main/source/svmsmn_stan_fit.R' )
  
  param = matrix(c(1.0, 0.985, 0.13,
                   0.01, 0.1, -0.02,
                   v), ncol = 1)
  
  param_hat = matrix(nrow = 7, ncol = r)
  
  for( i in 1:r ){
    
    if( i == 1 ) time.init = Sys.time()
    
    data = ts_data(mu = 1.0, phi = 0.985, sigma = 0.13,
                b0 = 0.01, b1 = 0.1, b2 = -0.02,
                y0 = 0,
                v = v, 
                T = 2e3,
                seed = seeds[ i ]
    )
    
    cat( paste0('réplica ', i, ' com o parâmetro v = ', v ) )
    param_hat[, i] = svmsmn_stan_fit(data = data$y, 
                                     y0,
                                     mu = param[1, 1], phi = param[2, 1], sigma = param[3, 1],
                                     h = data$h, l = data$l,
                                     b0 = param[4, 1], b1 = param[5, 1], b2 = param[6, 1],
                                     model, 
                                     M,
                                     warmup,
                                     nchains,
                                     lags)
    cat( '\r' )
    if( i == r ) time.final = Sys.time()
  }
  
  #vies
  vies = mean( param_hat - v )
  vies = data.frame( vies )
  row.names(vies) = c( 'v' )
  #vies = apply( param_hat - matrix( rep(param, r), ncol = r ), 1, mean )
  #vies = data.frame( vies )
  #row.names(vies) = c('mu', 'phi', 'sigma',
  #                    'b0', 'b1', 'b2', 'v')
  #smse
  smse = mean( (param_hat - v)**2 )
  smse = sqrt( smse )
  #smse = apply( (param_hat - matrix( rep(param, r), ncol = r ))**2, 1, mean )
  #smse = sqrt( smse )
  
  summary = cbind(vies, smse)
  
  return( list( summary = summary, 
                time = time.final - time.init) )
}
