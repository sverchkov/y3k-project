# NEM Analysis for fermentation data
library( nem )

load( "clean-data/resp-rdr-for-nem.RData" )

D = getDensityMatrix( resp.rdr.p, dirname="tmp" )

# getDensityMatrix forgets row names for some reason
dimnames( D ) = dimnames( resp.rdr.p )

save( D, file = "intermediates/resp-rdr.RData" )

# NEM code can't handle infinities
D[ D == -Inf ] = -.Machine$double.xmax

nem.control = set.default.parameters( colnames( D ) )
nem.control$type = "CONTmLLBayes"
  
resp.rdr.nem = nem( D, inference = "triples", control = nem.control )

save( resp.rdr.nem, file = "results/resp-rdr.nem.RData" )