## Tools for producing output to cytoscape

#' Output nem to cytoscape sif and node attribute tsv
#' @param nem.object The NEM object
#' @param sif.file The SIF file for output
#' @param atr.file The attribute file for output
nemToFiles = function( nem.object, sif.file, atr.file ){
  
  # Write edges
  
  G = nem.object$graph

  sif.lines = character()
  edge.list = edges( G )
  for ( source.node in names( edge.list ) ){
    target.nodes = paste0( edge.list[[ source.node ]], collapse = " " )
    sif.lines = c( sif.lines, paste0( source.node, " aa ", target.nodes ) )
  }
  
  effect.list = nem.object$mappos
  for ( source.node in names( effect.list ) ) if ( source.node != "null" ){
    effects = paste0( effect.list[[ source.node ]], collapse = " " )
    sif.lines = c( sif.lines, paste0( source.node, " ae ", effects ) )
  }
  
  writeLines( sif.lines, con = sif.file )

  # Make attribute table
  actions = nem.object$control$Sgenes
  effects = names( nem.object$LLperGene )
  attr.tab = data.frame(
    node = c( actions, effects ),
    type = c( rep( "action", length( actions ) ), rep( "effect", length( effects ) ) ),
    group = "none",
    stringsAsFactors = FALSE )
  # Assign effects to groups
  for ( action in names( effect.list ) ){
    attr.tab[ which( attr.tab$node %in% effect.list[[ action ]] ), "group"] = action
  }
  attr.tab[ which( attr.tab$node %in% actions ), "group" ] = paste0( actions, ".group" )
  for ( node in actions )
    for ( child in edge.list[[ node ]] )
      if ( node %in% edge.list[[ child ]] )
        attr.tab[ which( attr.tab$node == child ), "group" ] = attr.tab[ which( attr.tab$node == node ), "group" ]

  write.table( attr.tab, file = atr.file, quote = FALSE, sep = "\t", row.names = FALSE ) 
}