#' Get MAP FGNEM attachments
#' 
#' @param egenes.logprobs E gene log probability tables (named list where each element is an S gene with a table of probabilities of pos/null/neg attachments)
#' @return a character vector whose names are E genes and values are the MAP S-genes
getMAPFGNEMAttachments <- function ( egenes.logprobs ) {
  max.log.ratio <- Map( function ( egene.table ){
    non.nulls <- egene.table[,c('pos','neg')] - egene.table[,'null']
    apply( non.nulls, 1, max, na.rm = TRUE )
  }, egenes.logprobs )
  max.log.ratio.table <- cbind( do.call( cbind, max.log.ratio ), null = 0 )
  
  apply( max.log.ratio.table, 1, function ( scores ) names( which.max( scores ) ) )
}