freq.models = function( M ){
  models = rep(0, 4)
  
  for( i in 1:length(M) ){
    
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
  
  data = data.frame( models )
  row.names( data ) = c( 'normal', 'ts', 'slash', 'vgamma' )
  colnames( data ) = c('frequencia')
  
  return( data )
}
