# NEM Analysis for fermentation data
library( nem )

load( "clean-data/resp-for-nem.RData" )

D = getDensityMatrix( resp.p.proteins, dirname="tmp" )

# getDensityMatrix forgets row names for some reason
dimnames( D ) = dimnames( resp.p.proteins )

save( D, file = "intermediates/resp-proteins.RData" )

# NEM code can't handle infinities
D[ D == -Inf ] = -.Machine$double.xmax

nem.control = set.default.parameters( colnames( D ) )
nem.control$type = "CONTmLLBayes"
  
resp.proteins.nem = nem( D, inference = "triples", control = nem.control )

save( resp.proteins.nem, file = "results/resp-proteins.nem.RData" )