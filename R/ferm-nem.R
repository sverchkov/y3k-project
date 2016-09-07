# NEM Analysis for fermentation data

load( "clean-data/ferm-for-nem.RData" )

D = getDensityMatrix( ferm.p, dirname="tmp" )

dimnames( D ) = dimnames( ferm.p )

save( D, file = "intermediates/ferm.RData" )

ferm.nem = nem( D )

save( ferm.nem, file = "results/ferm.RData" )