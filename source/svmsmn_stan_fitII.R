svmsmn_stan_fitII = function(data, 
                             y0, 
                             model, 
                             M,
                             warmup,
                             nchains = 1, 
                             lags = 1,
                             model_name,
                             normal,
                             cores){
  
  model_selection = data.frame()
  T = length( data )
  fit = sampling(model, 
                 list(T = length( data ), 
                      y = data,
                      y0 = y0),
                 iter = warmup + M,
                 warmup = warmup,
                 thin = lags,
                 chains = nchains)
  parameters = extract( fit )
  source( 'https://raw.githubusercontent.com/holtz27/rbras/main/source/model_selection.R' )
  ############################################################################
  ############################## Model Selection
  # p( y | theta ) = p( y | b, theta_h, v, h ) = p( y | b, h )
  # data = ( y0 , y )
  # param_hat = ( b = b_hat, h = h_hat, l = l_hat )
  ############################################################################
  if( normal ){
    # Theta draws
    draws = matrix(c( parameters$mu,
                      parameters$phi,
                      parameters$sigma,
                      parameters$b0, 
                      parameters$b1, 
                      parameters$b2 ), nrow = 6, byrow = TRUE)
    draws = rbind( draws, t( as.matrix( parameters$h ) ) )
    draws = rbind( draws, matrix(1, nrow = T, ncol = length(parameters$mu)) )  
    # convergence check
    ############### Numeric Analysis
    ############################### theta
    mcmc = coda::as.mcmc( t( draws[1:6, ] ) )
    ####### Geweke Statistic
    # |G| > 1.96 evidencia não convergencia
    CD_theta = coda::geweke.diag( mcmc )
    convergence = sum( as.numeric( c( abs( CD_theta$z ) < 1.96 ) ) ) == 6 
    # model selection
    dic = svmsmn_dic(data = y, y0 = 0, 
                     param_draws = draws[4:(nrow(draws)), ])
    ############### waic
    waic = svmsmn_waic(data = y, y0 = 0,
                       draws = draws[4:(nrow(draws)), ])
    ############### loo
    loo = svmsmn_loo(data = y, y0 = 0, 
                     draws = draws[4:(nrow(draws)), ],
                     cores = cores)
  }else{
    # Theta draws
    draws = matrix(c( parameters$mu,
                      parameters$phi,
                      parameters$sigma,
                      parameters$b0, 
                      parameters$b1, 
                      parameters$b2,
                      parameters$v ), nrow = 7, byrow = TRUE)
    draws = rbind( draws, t( as.matrix( parameters$h ) ) )
    draws = rbind( draws, t( as.matrix( parameters$l ) ) )
    # convergence check
    ############### Numeric Analysis
    ############################### theta
    mcmc = coda::as.mcmc( t( draws[1:7, ] ) )
    ####### Geweke Statistic
    # |G| > 1.96 evidencia não convergencia
    CD_theta = coda::geweke.diag( mcmc )
    convergence = sum( as.numeric( c( abs( CD_theta$z ) < 1.96 ) ) ) == 7
    ############### DIC
    dic = svmsmn_dic(data = y, y0 = 0, 
                     param_draws = draws[ c( 4:6, 
                                             8:(T+7),
                                             (T+8):(nrow(draws)) ), ]
    )
    ############### waic
    waic = svmsmn_waic(data = y, y0 = 0,
                       draws = draws[ c( 4:6, 
                                         8:(T+7),
                                         (T+8):(nrow(draws)) ), ]
    )
    ############### loo
    loo = svmsmn_loo(data = y, y0 = 0, 
                     draws = draws[ c( 4:6, 
                                       8:(T+7),
                                       (T+8):(nrow(draws)) ), ],
                     cores = cores
    )
  }
  
  model_selection = rbind( model_selection, 
                           c(dic,
                             waic$estimates[3,1],
                             loo$estimates[3,1],
                             convergence)
  )
  row.names( model_selection ) = model_name
  colnames( model_selection ) = c( 'dic', 'waic', 'loo', 'convergence' )
  
  return( model_selection )
}
