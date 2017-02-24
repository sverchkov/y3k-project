library("adjacencymatrixoperations")
library("Matrix")
library("graph")

matrixFromNEM = function( nem.object ) {
  
  edge.list = edges( nem.object$graph )
  
  nodes = nodes( nem.object$graph )
  n = length( nodes )
  
  effects = nem.object$selected
  m = length( effects )
  
  result = Reduce( rbind, Map(
    function ( edge.vector, effect.vector )
      c( nodes %in% edge.vector, effects %in% effect.vector ),
    edge.list, nem.object$mappos[nodes] ) )
  
  colnames( result ) = c( nodes, effects )
  rownames( result ) = nodes
  
  cbind( transitivelyClose( result[1:n,1:n], fill.diagonal = TRUE ), result[1:n, (n+1):m] )
}

effectMatrixFromNEM = function( nem.object ) {
  
  nem.matrix = matrixFromNEM( nem.object )
  n = nrow( nem.matrix )
  m = ncol( nem.matrix ) - n
  nem.matrix[1:n,1:n] %&% nem.matrix[1:n,(n+1):m]
}