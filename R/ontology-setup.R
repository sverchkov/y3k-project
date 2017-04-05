# Set up ontology files
# This creates files in data_clean
library("dplyr")
# DOES NOT WORK library("ontoCAT") # For parsing the ontology file, downloaded from Bioconductor

# Function to extract systematic name
#extractSystematicName = function (v) sapply( strsplit( v, '[|]' ), function (x) x[1] )

# We will take the non-destructive approach and keep the junk columns
associations =
  read.delim("raw-data/gene_association.sgd", header = FALSE, comment.char = "!", stringsAsFactors = FALSE ) %>%
  # We care about the GO term, which is column 5. Name it
  rename( go.term = V5 ) %>%
  # Column 3 is the Standard Gene Name
  rename( standard.name = V3 ) #%>%
  # Extract systematic name
  #mutate( systematic.name = extract.systematic.name( V11 ) )

# Save the annotations.
saveRDS( associations, "clean-data/associations.rds" )

# Get a convenient name map view
name.map = distinct( associations, standard.name )

# Ideally, we'd use a proper ontology parser
# DOES NOT WORK ontology = getOntology( "raw-data/goslim_yeast.txt" )

# Extract go term names
parseOntology = function( file ){
  lines = readLines( file )
  idLines = grep( "^id: ", lines )
  list( ids = sub( "^id: ", "", lines[ idLines ] ),
        names = sub( "^name: ", "", lines[ idLines+1 ] ),
        namespaces = sub("^namespace: ", "", lines[ idLines+2 ] ) )
} 

id.name.pairs = parseOntology( "raw-data/goslim_yeast.obo" )

go.term.names = data.frame( id = id.name.pairs$ids, name = id.name.pairs$names, space = id.name.pairs$namespaces, stringsAsFactors = FALSE )

saveRDS( go.term.names, "clean-data/go-term-names.rds" )