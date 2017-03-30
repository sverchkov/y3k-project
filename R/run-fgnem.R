## Run FGNEM on RESP data
library("KnockoutNets")
library("dplyr")

# Load respiratory data
resp.data = read.delim( "raw-data/AvgKO_Resp_KO.txt", stringsAsFactors = FALSE ) %>%
  select( -X ) %>%
  rowwise() %>% mutate( Standard.Name = strsplit( Molecule.Name, " " )[[1]][1] )

colnames( resp.data ) = sub( "X.", "", colnames( resp.data ) )

# Save to clean data
saveRDS( resp.data, file = "clean-data/respiration.rds" )

expr = as.matrix( select( resp.data, -Molecule.Type, -Molecule.Name, -Standard.Name ) )
rownames( expr ) = resp.data$Standard.Name

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
                                 , doTransitivity = FALSE
                                 , summarization = max # or logsum
           )

# Save results
saveRDS( results, file = "results/respiration-fgnem.rds" )

## Unused:

# Write data to a table in the format that KnockoutNets likes
#
# Particularly:
# knockdown.cols Vector of knockdown names (one per expression matrix column)
# lof Loss-of-function genes (vector)
# stddev Vector of standard deviations (one per expression matrix column)
# expr Matrix of gene expression levels
# file Output file
# append Whether to append to the file
write.egene.tab( knockdown.cols = colnames( expr ),
                 lof = colnames( expr ),
                 stddev = apply( expr, 2, sd, na.rm = TRUE ),
                 expr,
                 file = "results/respiration.fgnem.tab",
                 append = FALSE)
