freq.models = function( M ){
  models = rep(0, 4)
  convergence = 0
  
  for( i in 1:length(M) ){
    
    if( (sum( df[[i]]$convergence ) == 4) && (!is.na(sum( df[[i]]$convergence ))) ){
      
      convergence = convergence + 1
      
      x = apply( M[[ i ]], MARGIN = 2, which.min )
      y = duplicated( x )
      indice = NULL
      
      for(i in 2:3){
        
        if( y[i] ){ 
          indice = x[ i ]
          models[ indice ] = models[ indice ] + 1
          break
        }
      } 
    }
    
  }
  
  data = data.frame( models )
  row.names( data ) = c( 'normal', 'ts', 'slash', 'vgamma' )
  colnames( data ) = c('frequencia')
  
  return( list(data = data, convergence = convergence) )
}
