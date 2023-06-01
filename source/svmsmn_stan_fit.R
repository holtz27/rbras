svmsmn_stan_fit = function(data, y0, 
                           model, M, nchains = 3, 
                           lags = 1){
  
  fit = sampling(model, 
                 list(T = length( data ), 
                      y = data,
                      y0 = y0),
                 iter = M,
                 chains = nchains)
  
  parameters = extract( fit )
  draws = matrix(c( parameters$mu,
                    parameters$phi,
                    parameters$sigma,
                    parameters$b0, 
                    parameters$b1, 
                    parameters$b2,
                    parameters$v ), nrow = 7, byrow = TRUE)
  # jumps
  jumps = seq(1, 0.5 * M * nchains, by = lags)
  draws = draws[, jumps ]
  
  return( matrix( apply(draws, 1, mean), ncol = 1 ) )
}
