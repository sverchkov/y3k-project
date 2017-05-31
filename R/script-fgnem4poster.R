# Script: prepare FGNEM for poster
library("dplyr")
library("gsubfn")

# Assuming that results is the fgnem, respiration.data is the data
stopifnot( exists( 'results' ) )
stopifnot( exists( 'respiration.data' ) )

# Get cliques
cliques <- findFGNEMCliques( results$acc )

# Get clique graph
clique.graph <- getCliqueGraph( cliques, results$adj )
# As a dataframe
graph.df = data.frame( source = character(0), sign = character(0), dest = character(0), stringsAsFactors = FALSE )
for ( s in rownames( clique.graph ) ) for ( d in colnames( clique.graph ) )
  if ( clique.graph[s,d] == 1 ) {
    graph.df[ nrow( graph.df ) + 1 ,] <- c( s, 'pos', d )
  }else if ( clique.graph[s,d] == -1 ) {
    graph.df[ nrow( graph.df ) + 1 ,] <- c( s, 'neg', d )
  }

# Get attachments of effects to KOs
attachment.array <- getMAPFGNEMAttachments( results$egenes.logprobs )

# Get attachment signs
attachment.signs <- getMAPFGNEMAttachmentSigns( results$egenes.logprobs )
attachment.signs.df <- data.frame(
  effect = rownames( attachment.signs ),
  sign = attachment.signs[,'list.name'],
  action = attachment.signs[,'gene'],
  stringsAsFactors = FALSE )

# Get attachments of effects to KO cliques
clique.attachment.array = getCliqueAttachments( cliques, attachment.array )

# Convert clique attachments to a dataframe
clique.attachment.df <-
  data.frame( action.clique = clique.attachment.array, effect = names( clique.attachment.array ), stringsAsFactors = FALSE ) %>%
  filter( action.clique != "null" ) %>% # Filter out null-attachments
  left_join( attachment.signs.df, by = c( 'effect' = 'effect' ) ) # Add up/down regulation

# Filter to proteins only
clique.attachment.proteins <- clique.attachment.df %>%
  inner_join( respiration.data %>%
               filter( Molecule.Type == 'Protein' ) %>%
               select( Molecule.Name ),
             by = c( 'effect' = 'Molecule.Name' ) ) %>%
  mutate( effect.std = toupper( strapplyc( effect, '^[^ ]+ \\((.*)\\)$', simplify = TRUE ) ) ) # Standardze name

# Make attachment groups
clique.attachment.groups <- clique.attachment.proteins %>% group_by( action.clique, sign ) %>%
  summarize( effect.group = paste( effect.std, collapse = ";" ) ) %>% ungroup()

## Make dataframe for sif
##
# First, we get a mapping of cliques names to IDs
action.mapping.df <- data.frame( set = rownames( clique.graph ), id = paste0( 'A', 1:nrow( clique.graph ) ), stringsAsFactors = FALSE )
effect.mapping.df <- clique.attachment.groups %>% distinct( effect.group ) %>% mutate( id = paste0( 'E', 1:length( effect.group ) ) )
all.mapping.df <- rbind.data.frame( action.mapping.df, effect.mapping.df %>% rename( set = effect.group ) )

for.sif.df <- rbind(
  graph.df,
  clique.attachment.groups %>% ungroup() %>% select( source = action.clique, sign = sign, dest = effect.group ) ) %>%
  left_join( all.mapping.df, by = c( 'source' = 'set' ) ) %>% rename( source.id = id ) %>%
  left_join( all.mapping.df, by = c( 'dest' = 'set' ) ) %>% rename( dest.id = id ) %>%
  filter( source.id != dest.id ) %>%
  select( source.id, sign, dest.id )

# Out to files
write.csv( all.mapping.df, file = "results/resp-fgnem-nodes.csv", row.names = FALSE, quote = FALSE )
write.table( for.sif.df, file = "results/resp-fgnem-edges.sif", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t" )

# Get lists of proteins from sets
effect.groups <- strsplit( effect.mapping.df$effect.group, ';' )

# Get enrichments!
alpha = 0.1
gamma = 0.9

gs <- gs2edge( character(), n.upp = 50, idtype = "ENSEMBL", lib = "org.Sc.sgd" )
I <- gs$I

for ( i in 2:length(effect.groups) ){
  gene.group <- effect.groups[[i]]
  
  y <- gs$y
  y[ intersect( names( gs$y ), gene.group ) ] <- 1
  
  p <- estP( I, y, alpha, gamma )
  
  res <- NULL
  try( res <- sequentialRM( I, y, nupstart = 10, by=1, alpha, gamma, p ) )
  if( !is.null( res ) ){
    saveRDS( res, file = paste0( "results/go-plots/result-",i,".rds" ) )
    
    pdf( paste0( "results/go-plots/effects-",i,".pdf" ), width=15, height=5)
    rmPlot( I, y, res$onwholes )
    dev.off()
  }
}

### The test zone
# Get map enrichment groups for effect.groups[[1]]
gene.group = effect.groups[[1]]

alpha = 0.01
gamma = 0.02

# For some reason the ribosomal proteins are also off-limits?
gene.group <- gene.group[ !grepl( "^R", gene.group ) ]

gs <- gs2edge( gene.group, n.upp = 50, idtype = "ENSEMBL", lib = "org.Sc.sgd" )
I <- gs$I
y <- gs$y

p <- estP( I, y, alpha, gamma )

res <- sequentialRM( I, y, nupstart = 10, by=1, alpha, gamma, p )

pdf("results/go-plots/effects-1.pdf", width=15, height=5)
rmPlot( I, y, res$onwholes )
dev.off()

### Using individual attachments below

# Convert attachments to a dataframe
attachment.df <- data.frame( action = attachment.array, effect = names( attachment.array ), stringsAsFactors = FALSE )

deduplication.table <-
  data.frame(
    effect.std1 = c( "YDR385wYOR133w", "YDR385wYOR133w",
                     "YHR215w;YAR071w", "YHR215w;YAR071w",
                     "YHR021C;YKL156W", "YHR021C;YKL156W",
                     "YHR043C;YHR044C", "YHR043C;YHR044C",
                     "YJL221C;YIL172C;YOL157c", "YJL221C;YIL172C;YOL157c", "YJL221C;YIL172C;YOL157c",
                     "YER117W;YBL087C", "YER117W;YBL087C" ),
    effect.std = c( "YDR385w", "YOR133w",
                    "YHR215w", "YAR071w",
                    "YHR021C", "YKL156W",
                    "YHR043C", "YHR044C",
                    "YJL221C", "YIL172C", "YOL157c",
                    "YER117W", "YBL087C" ),
    stringsAsFactors = FALSE
  )

protein.attachments <-
  attachment.df %>%
  left_join( respiration.data %>% select( Molecule.Name, Molecule.Type ), by = c( "effect" = "Molecule.Name" ) ) %>%
  filter( Molecule.Type == 'Protein' ) %>%
  mutate( effect.std1 = strapplyc( effect, '^[^ ]+ \\((.*)\\)$', simplify = TRUE ) ) %>%
  left_join( deduplication.table, by = c( "effect.std1" = "effect.std1" ) ) %>%
  mutate( effect.std = toupper( ifelse( is.na( effect.std ), effect.std1, effect.std ) ) )

# Get sets of effects
effect.sets <- list()
for ( a in unique( protein.attachments$action ) ) if( a != "null" ) {
  effect.sets[[ length(effect.sets) + 1 ]] <- ( protein.attachments %>% filter( action == a ) )$effect.std
}

# OK now we can get enrichments.

# Using Rolemodel for gene ontology
library("Rolemodel")
# The yeast database is org.Sc.sgd.db
library("org.Sc.sgd.db")
# Hack to get Rolemodel to play nice with org.Sc.sgd
org.Sc.sgdGO2ALLEGS <- org.Sc.sgdGO2ALLORFS

enrichments <- Map( function ( gene.set ){
  res <- NULL
  try( res <- rmTable( gene.set, idtype = "ENSEMBL", lib = "org.Sc.sgd", n.upp = 20, n.low = 1, nupstart = 10, by = 1 ) )
}, effect.sets )
