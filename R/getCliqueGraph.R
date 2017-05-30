#' Get a graph of cliques
#' 
#' @param cliques a list of character arrays representing cliques of nodes
#' @param adj an adjacency matrix among nodes, with +ve and -ve edges
#' @return an adjacenct matrix among cliques, with +ve and -ve edges
getCliqueGraph <- function ( cliques, adj ) {
  clique.strs <- Reduce( c, Map( paste, cliques, collapse = ";" ) )
  n <- length( clique.strs )
  clique.graph <- matrix( data = 0, nrow = n, ncol = n, dimnames = list( clique.strs, clique.strs ) )

  for ( a in rownames( adj ) ){
    clique.a <- which( Reduce( c, Map( function ( c ) a %in% c, cliques ) ) )
    
    for ( b in colnames( adj ) ){
      edge <- adj[a,b]
      if ( edge != 0 ){
        clique.b <- which( Reduce( c, Map( function ( c ) b %in% c, cliques ) ) )
        old.edge = clique.graph[ clique.a, clique.b ]
        
        if ( edge != old.edge ){
          if ( old.edge != 0 ){
            warning( paste( "Prior edge between", clique.strs[ clique.a ], "and", clique.strs[ clique.b ], "was", old.edge, "but new edge between", a, "and", b, "is", edge ) )
          }else{
            clique.graph[ clique.a, clique.b ] <- edge
          }
        }
      }
    }
  }
  
  clique.graph
}