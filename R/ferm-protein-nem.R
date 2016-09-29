# NEM Analysis for fermentation data
library( nem )

load( "clean-data/ferm-for-nem.RData" )

D = getDensityMatrix( ferm.p.proteins, dirname="tmp" )

# getDensityMatrix forgets row names for some reason
dimnames( D ) = dimnames( ferm.p.proteins )

save( D, file = "intermediates/ferm-proteins.RData" )

# NEM code can't handle infinities
D[ D == -Inf ] = -.Machine$double.xmax

nem.control = set.default.parameters( colnames( D ) )
nem.control$type = "CONTmLLBayes"
  
ferm.proteins.nem = nem( D, inference = "triples", control = nem.control )

save( ferm.proteins.nem, file = "results/ferm-proteins.nem.RData" )