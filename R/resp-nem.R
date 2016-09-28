# NEM Analysis for fermentation data
library( nem )

load( "clean-data/resp-for-nem.RData" )

D = getDensityMatrix( resp.p, dirname="tmp" )

# getDensityMatrix forgets row names for some reason
dimnames( D ) = dimnames( resp.p )

save( D, file = "intermediates/resp.RData" )

# NEM code can't handle infinities
D[ D == -Inf ] = -.Machine$double.xmax

nem.control = set.default.parameters( colnames( D ) )
nem.control$type = "CONTmLLBayes"
  
resp.nem = nem( D, inference = "triples", control = nem.control )

save( resp.nem, file = "results/resp.nem.RData" )