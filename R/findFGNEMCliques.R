#' Find FGNEM cliques
#' 
#' Finds Cliques of +ve interactions in and FGNEM
#' 
#' @param acc the accessibility matrix
#' @return a list of character arrays, each array is a clique
findFGNEMCliques <- function ( acc ) {
  #node.names <- rownames( acc )
  bin.acc <- (acc == 1)
  clique.acc <- bin.acc & t( bin.acc )
  clique.ptrs <- apply( clique.acc, 1, function ( col ) which( col )[1] )
  
  Map( function ( clique.id ) names( which( clique.ptrs == clique.id ) ), unique( clique.ptrs ) )
}