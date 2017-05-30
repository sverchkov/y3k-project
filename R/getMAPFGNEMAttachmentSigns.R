#' Helper
#' 
#' @param scores.list list of score arrays
#' @return matrix with highest score, name of highest score
maxAcrossList <- function ( scores.list ) {
  table <- do.call( cbind, scores.list )
  cbind( value = apply( table, 1, max, na.rm = TRUE ), list.name = apply( table, 1, function ( scores ) names( which.max( scores ) ) ) )
}

#' Get MAP FGNEM attachments
#' 
#' @param egenes.logprobs E gene log probability tables (named list where each element is an S gene with a table of probabilities of pos/null/neg attachments)
#' @return a numeric vector of -1, 0, and 1's whose names are E genes and values the sign of attachment to the MAP S-genes
getMAPFGNEMAttachmentSigns <- function ( egenes.logprobs ) {
  
  pos.log.ratio <- Map( function ( egene.table ){
    egene.table[,'pos'] - egene.table[,'null']
  }, egenes.logprobs )

  neg.log.ratio <- Map( function ( egene.table ){
    egene.table[,'neg'] - egene.table[,'null']
  }, egenes.logprobs )

  pos.maxima <- maxAcrossList( pos.log.ratio )
  neg.maxima <- maxAcrossList( neg.log.ratio )
  
  pos.or.neg <- maxAcrossList( list( pos = pos.maxima[,'value'], neg = neg.maxima[,'value'] ) )
  
  cbind( pos.or.neg, gene = ifelse( pos.or.neg[,'list.name'] == "pos", pos.maxima[,'list.name'], neg.maxima[,'list.name'] ) )
}