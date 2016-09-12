# NEM Analysis for fermentation data

load( "clean-data/ferm-for-nem.RData" )

D = getDensityMatrix( ferm.p, dirname="tmp" )

# getDensityMatrix forgets row names for some reason
dimnames( D ) = dimnames( ferm.p )

save( D, file = "intermediates/ferm.RData" )

# NEM code can't handle infinities
D[ D == -Inf ] = -.Machine$double.xmax

library( nem )
nem.control = set.default.parameters( colnames( D ) )
nem.control$type = "CONTmLLBayes"
  
ferm.nem = nem( D, inference = "triples", control = nem.control )

save( ferm.nem, file = "results/ferm.RData" )