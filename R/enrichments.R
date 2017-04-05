library("dplyr")

#' Get category hits
#' 
#' @param category The category whose enrichment we're computing
#' @param objects The group of objects
#' @param mapping The mapping defining object memberships in categories
getCategoryHits = function ( category, objects, mapping )
  mapping %>% filter( .[[2]] %in% category, .[[1]] %in% objects ) %>% distinct( hits = .[[1]] )

#' Count category hits
#' 
#' @param category The category whose enrichment we're computing
#' @param objects The group of objects
#' @param mapping The mapping defining object memberships in categories
countCategoryHits = function ( category, objects, mapping )
  nrow( getCategoryHits( category, objects, mapping ) )
  
#' Get hypergeometric enrichment of a category in a group of objects
#' 
#' @param category The category whose enrichment we're computing
#' @param objects The group of objects
#' @param mapping The mapping defining object memberships in categories
hyperGeometricEnrichment = function( category, objects, mapping ){
  bg.size = as.numeric( mapping %>% distinct( .[[1]] ) %>% summarize( n() ) )
  q = countCategoryHits( category, objects, mapping ) # Success in sample
  m = as.numeric( mapping %>% filter( .[[2]] == category ) %>% summarize( n() ) ) # Success in bg
  n = bg.size - m # Failure in bg
  k = length( objects ) # Sample size
  phyper( max( 0, q-1 ), m, n, k, lower.tail = FALSE )
}

#' Get a ranked list of enrichments
#'
#' @param objects The group of objects
#' @param mapping The mapping defining object memberships in categories
hyperGeometricEnrichmentRanking = function( objects, mapping ){
  mapping %>% distinct( category = .[[2]] ) %>%
    rowwise() %>%
    mutate( pVal = hyperGeometricEnrichment( category, objects, mapping ) ) %>%
    arrange( pVal )
}

#' Benjamini-Hochberg FDR filter
#' 
#' @param rankedList: data frame with ascending ranked column "pVal"
#' @param alpha: the false discovery rate (default is 0.05)
filterFDR = function ( rankedList, alpha = 0.05 ){
  # for a goven alpha, find largest k such that P_k <= k/m * alpha where
  # k is the position down the list, P_k is the pVal, and m is the length of the list
  m = nrow( rankedList )
  k = 0
  for ( i in 1:m ) {
    if ( rankedList[i,"pVal"] > (i/m)*alpha ) break;
    k = i
  }
  rankedList[ seq( by=1, length.out=k ), ]
}