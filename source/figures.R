trace_plots = function(draws, names){
  n = nrow( draws )
  par(mfrow = c(1,1))
  mat = matrix(seq(1, 3 * n), nrow = 3, byrow = FALSE)
  layout( mat )
  
  for(i in 1:n){
    plot( draws[i, ], type = 'l', main = names[i], xlab = '', ylab = '')
    plot(acf( draws[i, ], lag.max = 100, plot = FALSE)[1:100], main ='', xlab = '', ylab = 'ACF')
    plot(density( draws[i, ] ), main ='', xlab = '', ylab = '')
  }
  
  par(mfrow = c(1,1))
}

abs_plots = function( draws_h, y ){
  library(ggplot2)
  T = length( y )
  e.vol_hat = apply( exp( draws_h ) , MARGIN = 1, FUN = mean )
  e.vol_min = apply( exp( draws_h ) , MARGIN = 1, FUN = quantile, probs = c(0.025) )
  e.vol_max = apply( exp( draws_h ) , MARGIN = 1, FUN = quantile, probs = c(0.975) )
  data = matrix(c(1:T, abs( y ), e.vol_hat, e.vol_min, e.vol_max), ncol = 5)
  data = data.frame(data)
  names(data) = c('obs', 'y.abs', 'e.hat', 'e.min','e.max')
  h = ggplot(data)
  h = h + geom_line(aes(obs, y.abs), color = 'grey70')
  h = h + geom_ribbon(aes(x = 1:T, ymax = e.vol_max, ymin = e.vol_min), 
                      fill = 'blue' ,alpha = 0.2)
  h = h + geom_line(aes(obs, e.hat), linewidth = 0.75)
  h = h + theme_test() + xlab('') + ylab('|Retornos|')
}

tail_plot = function(draws, model_name){
  library(ggplot2)
  library(latex2exp)
  l_hat = apply(draws, MARGIN = 1, mean)
  df = data.frame(obs = 1:(length(l_hat)), l = l_hat)
  h = ggplot(df) + geom_line(aes(x=obs, y = l)) 
  h = h + theme_test() + xlab('') 
  h = h + ylab(TeX(paste0('SVM-', model_name, ': ', '$\\lambda_{t}$')))
}

