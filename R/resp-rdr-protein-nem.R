# NEM Analysis for fermentation data
library( nem )

load( "clean-data/resp-rdr-for-nem.RData" )

D = getDensityMatrix( resp.rdr.p.proteins, dirname="tmp" )

# getDensityMatrix forgets row names for some reason
dimnames( D ) = dimnames( resp.rdr.p.proteins )

save( D, file = "intermediates/resp-rdr-proteins.RData" )

# NEM code can't handle infinities
D[ D == -Inf ] = -.Machine$double.xmax

nem.control = set.default.parameters( colnames( D ) )
nem.control$type = "CONTmLLBayes"
  
resp.rdr.proteins.nem = nem( D, inference = "triples", control = nem.control )

save( resp.rdr.proteins.nem, file = "results/resp-rdr-proteins.nem.RData" )