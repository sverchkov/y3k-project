## Run FGNEM on RESP data
library("KnockoutNets")
library("dplyr")

# Load respiratory data
respiration.data = read.delim( "raw-data/AvgKO_Resp_KO.txt", stringsAsFactors = FALSE ) %>% select( -X )

colnames( respiration.data ) = sub( "X.", "", colnames( respiration.data ) )

# Save to clean data
saveRDS( respiration.data, file = "clean-data/respiration.v2.rds" )

expr = as.matrix( select( respiration.data, -Molecule.Type, -Molecule.Name ) )
rownames( expr ) = respiration.data$Molecule.Name

# Limiter for experimentation
# expr = expr[,1:10]

# Build object for input to FGNEM
eg = list( egenes = expr
         , knockdown.cols = colnames( expr )
         , lof = colnames( expr )
         , stddev = apply( expr, 2, sd, na.rm = TRUE ) )

# FGNEM settings
params = paramGen( 1.5 , 1 ) # Defaults, maybe worth changing

# Run FGNEM
results <- scoreBestModelEstimate( eg
                                 , params = params
                                 , doTransitivity = TRUE
                                 , summarization = max # or logsum
           )

# Save results
saveRDS( results, file = "results/respiration-fgnem.rds" )