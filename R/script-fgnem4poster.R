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

# Get attachments of effects to KOs
attachment.array <- getMAPFGNEMAttachments( results$egenes.logprobs )

# Get attachment signs
attachment.signs <- getMAPFGNEMAttachmentSigns( results$egenes.logprobs )

# Get attachments of effects to KO cliques
clique.attachment.array = getCliqueAttachments( cliques, attachment.array )

# Convert clique attachments to a dataframe
clique.attachment.df <-
  data.frame( action.clique = clique.attachment.array, effect = names( clique.attachment.array ) ) %>%
  filter( action.clique != "null" ) # Filter out null-attachments
  
# Attachments with up/down regulation


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
  result <- NULL
  try( result <- rmTable( gene.set, idtype = "ENSEMBL", lib = "org.Sc.sgd", n.upp = 20, n.low = 1, nupstart = 10, by = 1 ) )
}, effect.sets )
