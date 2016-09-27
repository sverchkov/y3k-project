# Output nem to cytoscape

load("results/ferm.nem.RData")

G = ferm.nem$graph

file.conn = file("results/ferm.nem.sif")

lines = character()
for ( a in G@nodes ){
  children = paste0( G@nodes[G@edgeL[[a]]$edges], collapse = " " )
  lines = c( lines, paste0( a, " aa ", children ) )
}

for ( a in G@nodes ){
  effects = paste0( ferm.nem$mappos[[a]], collapse = " " )
  lines = c( lines, paste0( a, " ae ", effects ) )
}

writeLines( lines, file.conn )
close( file.conn )

file.conn = file("results/ferm.nem.orphans.txt")
writeLines( ferm.nem$mappos$null, file.conn )
close( file.conn )