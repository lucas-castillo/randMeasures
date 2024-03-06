.manage_alternatives <- function(a, s){
  if (!is.null(a) & length(a) == 1){
    stop("a should be a vector of alternatives, not the number of alternatives")
  }
  if (is.null(a)){
    a <- min(s):max(s)
  }
  return(a)
}


#' Redundancy
#'
#' @param s Sequence vector.
#' @param a Optional. Vector of alternatives (e.g. 1:6 in a mental dice task).
#'
#' @return Redundancy value from 0 to 1.
#'
#' @export
#'
#' @examples
#' redundancy(sample(1:10, 100, T))
redundancy <- function(s, a=NULL){
  a <- .manage_alternatives(a, s)

  n <- length(s)
  H_max <- log2(length(a))
  counts <- sapply(a, \(x){sum(s == x)})
  H_single <- log2(n) -
    sum(counts * log2(counts), na.rm = T) / n
  return(1 - H_single/H_max)
}
response_frequencies <- function(s, a=NULL){
  a <- .manage_alternatives(a,s)
  n <- length(s)
  frequencies <- sapply(a, \(x){sum(s == x)}) / n
  return(data.frame(alternative = a, frequency=frequencies))
}
coupon <- function(s,a=NULL){
  # Ginsburg and Karpiuk (1994)
  a <- .manage_alternatives(a,s)
  n <- length(s)
  alternatives_found <- rep(0, length(a))
  accumulator <- c()
  for (i in 1:length(s)){
    alternatives_found[which(a == s[i])] <- alternatives_found[which(a == s[i])] + 1
    if (!any(alternatives_found == 0)){
      accumulator <- c(accumulator, sum(alternatives_found))
      alternatives_found <- rep(0, length(a))
    }
  }
  mean(accumulator)
}
gap <- function(s, a=NULL){
  a <- .manage_alternatives(a,s)
  diffs <- c()
  for (alt in a){
    diffs <- c(diffs, diff(which(s == alt)))
  }
  return(median(diffs))
}
