slash_estudo1 = function(M,
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
                         T,
                         seeds = NULL){
  source( 'https://raw.githubusercontent.com/holtz27/rbras/main/source/data/slash_data.R' )
  v1 = v2 = v3 = rep(0, r)
  ess1 = ess2 = ess3 = time1 = time2 = time3 = matrix(nrow = length( tails ), ncol = r)
  lin = 1
  for( v in tails ){
    if( v == tails[1] ) time.init = Sys.time()
    for( i in 1:r ){
      cat( paste0('réplica ', i, ' com o parâmetro v = ', v ) )
      # Data generation
      data.gen = slash_data(mu = mu, phi = phi, sigma = sigma,
                            b0 = b0, b1 = b1, b2 = b2,
                            y0 = y0, v = v, T = T, seed = seeds[ i ])
      # Data generation
      y = data.gen$y
      h = data.gen$h
      l = data.gen$l
      # models fit
      time.fit.init = Sys.time()
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
      time1[ lin, i ] = Sys.time() - time.fit.init
      v1[ i ] = mean( extract( fit1 )$v )
      ess1[ lin, i ] = as.numeric( summary( fit1 )$summary[1, 'n_eff'] )
      
      time.fit.init = Sys.time()
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
      time2[ lin, i ] = Sys.time() - time.fit.init
      v2[ i ] = mean( extract( fit2 )$v )
      ess2[ lin, i ] = as.numeric( summary( fit2 )$summary[1, 'n_eff'] )
      
      time.fit.init = Sys.time()
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
      time3[ lin, i ] = Sys.time() - time.fit.init
      v3[ i ] = mean( extract( fit3 )$v )
      ess3[ lin, i ] = as.numeric( summary( fit3 )$summary[1, 'n_eff'] )
      cat( '\r' )
      
    }
    lin = lin + 1
    vies = matrix( c(vies1 = mean( v1 - v ),
                     vies2 = mean( v2 - v ),
                     vies3 = mean( v3 - v )), ncol = 1 )
    #vies = round( vies, digits = 3 )
    smse = matrix( c(smse1 = mean( (v1 - v)**2 ),
                     smse2 = mean( (v2 - v)**2 ),
                     smse3 = mean( (v3 - v)**2 )), ncol = 1)
    #smse = round( smse, digits = 3 )
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
  row.names( ess1 ) = paste0( 'v=', tails )
  colnames( ess1 ) = paste0( 'it ', 1:r )
  row.names( ess2 ) = paste0( 'v=', tails )
  colnames( ess2 ) = paste0( 'it ', 1:r )
  row.names( ess3 ) = paste0( 'v=', tails )
  colnames( ess3 ) = paste0( 'it ', 1:r )
  row.names( time1 ) = paste0( 'v=', tails )
  colnames( time1 ) = paste0( 'time it ', 1:r )
  row.names( time2 ) = paste0( 'v=', tails )
  colnames( time2 ) = paste0( 'time it ', 1:r )
  row.names( time3 ) = paste0( 'v=', tails )
  colnames( time3 ) = paste0( 'time it ', 1:r )
  return( list(summary = summary, 
               time = time.final - time.init,
               ess1 = ess1,
               ess2 = ess2,
               ess3 = ess3,
               time1 = time1,
               time2 = time2,
               time3 = time3) )
}
