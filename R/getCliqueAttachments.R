#' Group attachments by cliques
#' 
#' @param cliques list of character arrays of genes, representing cliques
#' @param attachment.array array of attachments, where names are effects and values are genes ( that correspond to cliques )

getCliqueAttachments <- function ( cliques, attachment.array ){
  result <- attachment.array
  for ( clique in cliques )
    result[ attachment.array %in% clique ] <- paste( clique, collapse = ";" )
  
  result
}