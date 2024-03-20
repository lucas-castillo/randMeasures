.manage_alternatives <- function(a, s){
  if (!is.null(a) & length(a) == 1){
    stop("a should be a vector of alternatives, not the number of alternatives")
  }
  if (is.null(a)){
    a <- min(s):max(s)
  }
  return(a)
}
.trim <- function(s){
  r <- repetitions(s, T)
  s[c(T, !r[2:length(r)])]
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

#' Repetitions
#'
#' @param s Sequence vector.
#' @param full Optional. Boolean.
#'
#' @return If full returns a vector of length = length(s) with whether the ith item is a repetition of the last. Otherwise returns the mean repetition rate.
#' @export
#'
#' @examples
#' s <- sample(1:10, 100, T)
#' repetitions(s)
#' repetitions(s, full=T)
repetitions <- function(s, full=F){
  v <- s[2:length(s)] == s[1:(length(s) - 1)]
  if (full) return (c(NA, v)) else return (mean(v, na.rm=T))
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
  if (is.null(accumulator)){return(NA)}else{return(mean(accumulator))}
}
repetition_gap <- function(s, a=NULL, measure="median"){
  a <- .manage_alternatives(a,s)
  diffs <- c()
  for (alt in a){
    diffs <- c(diffs, diff(which(s == alt)))
  }
  if (measure == "median"){
    return(median(diffs))
  } else if (measure == "mean"){
    return(mean(diffs))
  } else if (measure == "mode"){
    as.numeric(names(which(table(diffs) == max(table(diffs)))))
  } else if (measure == "distances"){
    return(table(diffs, dnn = "distances"))
  }
}

poker <- function(s,a=NULL){
  #GK94
}
.response_matrix <- function(s, a, wrap=T){
  M <- matrix(0, length(a), length(a))
  for (i in 2:length(s)){
    last_index <- which(a == s[i-1])
    this_index <- which(a == s[i])
    M[last_index, this_index] <- M[last_index, this_index] + 1
  }
  if (wrap){
    last_index <- which(a == s[i])
    this_index <- which(a == s[1])
    M[last_index, this_index] <- M[last_index, this_index] + 1
  }
  return(M)
}
RNG <- function(s, a=NULL){
  # Evans' RNG. Note T+Neil has a typo here
  a <- .manage_alternatives(a,s)
  M <- .response_matrix(s,a, wrap = T)

  row_sums <- rowSums(M)
  top <- sum(M * log10(M), na.rm = T)
  if (top == 0) return(0)
  bottom <- sum(row_sums * log10(row_sums), na.rm = T)
  top / bottom
}

NSQ <- function(s, a=NULL){
  a <- .manage_alternatives(a,s)
  M <- .response_matrix(s,a, wrap = T)
  sum(M == 0) / (length(a) ** 2 - 1)
}

adjacency <- function(s, full=F, unpack=F, drop.r=NULL){
  if (is.null(drop.r)){
    drop.r <- T
    message("Calculating adjacency percentage after dropping repetitions")
  }
  if (drop.r){
    r <- repetitions(s, T)
  }
  # note that t&n use count / n, not count / (n-1)
  diffs <- diff(s)
  descending <- diffs == -1
  ascending <- diffs == 1
  combined <- ascending | descending

  if (drop.r){
    r <- repetitions(s, T)
    descending[r] <- NA
    ascending[r] <- NA
    combined[r] <- NA
  }
  if (full & unpack){
    return(data.frame(
      descending = c(NA, descending),
      ascending = c(NA, ascending),
      combined = c(NA, combined)
    ))
  } else if (full){
    return(c(NA, combined))
  } else if (unpack){
    return(list(
      "descending" = mean(descending, na.rm=T),
      "ascending" = mean(ascending, na.rm=T),
      "combined" = mean(combined, na.rm=T)
    ))
  } else{
    return(mean(combined, na.rm=T))
  }
}

turning_points <- function(s, full=F, drop.r=NULL){
  if (is.null(drop.r)){
    drop.r <- T
    message("Calculating turning point percentage after dropping repetitions")
  }
  if (drop.r){
    n <- length(s)
    r <- repetitions(s, T)
    one <- s[1:(n-2)]
    two <- s[2:(n-1)]
    three <- s[3:(n)]
    TP <- (one > two & two < three) | (one < two & two > three)
    TP[r] <- NA
    if (full){return(TP)} else{return(mean(TP, na.rm=T))}
  }
  d <- diff(s)
  d[d > 0] <- 1
  d[d < 0] <- -1
  zero_idx <- which(d == 0)
  d <- d[d!=0]
  tp <- abs(diff(d)) == 2
  if (length(zero_idx) > 0){
    for (i in 1:length(zero_idx)){
      if (i == 1){
        tp2 <- append(tp, F, after=zero_idx[i]-1)
      } else{
        tp2 <- append(tp2, F, after=zero_idx[i]-1)
      }
    }
  } else {
    tp2 <- tp
  }
  if (full){return(c(NA, tp2, NA))} else return(mean(tp2, na.rm = T))
}
distances <- function(s, full=F, drop.r=NULL){
  n <- length(s)
  if (is.null(drop.r)){
    drop.r <- T
    message("Calculating distances after dropping repetitions")
  }
  if (drop.r){
    r <- repetitions(s, T)
    D <- abs(s[2:n] - s[1:(n-1)])
    D[r] <- NA
  } else{
    D <- abs(s[2:n] - s[1:(n-1)])
  }

  if (full){
    return(c(NA, D))
  } else{
    return(mean(D, na.rm = T))
  }
}

phase_length <- function(s){
  phases <- diff(which(turning_points(s, T, F)))
  table(phases)
}

cluster_ratio <- function(s, a=NULL){
  a <- .manage_alternatives(a,s)
  M <- .response_matrix(s,a, wrap = T)
  var(as.vector(M))
}
runs <- function(s){
  ended <- F
  i <- 2
  k <- 1
  run_counter <- c()
  trend <- s[2] - s[1]
  while (!ended){
    if (trend > 0){
      start <- i-1
      while(trend > 0){
        i <- i + 1
        if (i > length(s)) break
        trend <- s[i] - s[i-1]
      }
      end <- i-1
      run_counter <- c(run_counter, end-start+1)
    } else if (trend <=0){
      start <- i-1
      while(trend <=0){
        i <- i + 1
        if (i > length(s)) break
        trend <- s[i] - s[i-1]
      }
      end <- i-1
    }

    if (i > length(s)){
      ended <- T
    }
  }
  unaccounted_items <- length(s) - sum(run_counter)
  return(var(c(run_counter, rep(1, unaccounted_items))))
}

fod <- function(s, a=NULL){
  TT <- table(diff(s))
  data.frame(difference=as.integer(names(TT)), frequency=as.vector(TT) / (length(s) - 1))
}
.interleaved_digrams <- function(s, a){
  M <- matrix(0, length(a), length(a))
  for (i in 3:length(s)){
    last_index <- which(a == s[i-2])
    this_index <- which(a == s[i])
    M[last_index, this_index] <- M[last_index, this_index] + 1
  }
  return(M)
}
RNG2 <- function(s, a=NULL){
  a <- .manage_alternatives(a,s)
  M <- .interleaved_digrams(s,a)
  row_sums <- rowSums(M)
  top <- sum(M * log10(M), na.rm = T)
  if (top == 0) return(0)
  bottom <- sum(row_sums * log10(row_sums), na.rm = T)
  top / bottom
}

#' Get all measures
#'
#' @param s Sequence vector.
#' @param a Optional. Vector of alternatives (e.g. 1:6 in a mental dice task).
#'
#' @return A dataframe with all measures in the package
#' @export
#'
#' @examples
#' all_measures(sample(1:5, 100, T))
all_measures <- function(s, a=NULL, drop.r=NULL){
  if (is.null(drop.r)){
    drop.r <- T
    message("Calculating Adj, TPs and distances")
  }
  a <- .manage_alternatives(a,s)
  response_frequencies <- response_frequencies(s, a)
  phase_length <- phase_length(s)
  phase_length <- phase_length / sum(phase_length)

  fod <- fod(s, a)



  data.frame(
    redundancy = redundancy(s, a),
    repetitions = repetitions(s),
    adjacency=adjacency(s, drop.r = drop.r),
    turning_points=turning_points(s, drop.r = drop.r),
    distances=distances(s, drop.r = drop.r),
    coupon = coupon(s, a),
    repetition_gap_med =repetition_gap(s,a, "median"),
    repetition_gap_mean= repetition_gap(s,a, "mean"),
    RNG = RNG(s,a),
    RNG2 = RNG2(s, a),
    NSQ = NSQ(s,a),
    cluster_ratio = cluster_ratio(s,a),
    phase_length_1=if ("1" %in% names(phase_length)){phase_length["1"]}else{0},
    phase_length_2=if ("2" %in% names(phase_length)){phase_length["2"]}else{0},
    phase_length_3=if ("3" %in% names(phase_length)){phase_length["3"]}else{0},
    phase_length_4=if ("4" %in% names(phase_length)){phase_length["4"]}else{0},
    phase_length_5=if ("5" %in% names(phase_length)){phase_length["5"]}else{0},
    fod_m5 = if (-5 %in% fod$difference) {fod$frequency[fod$difference == -5]}else{0},
    fod_m4 = if (-4 %in% fod$difference) {fod$frequency[fod$difference == -4]}else{0},
    fod_m3 = if (-3 %in% fod$difference) {fod$frequency[fod$difference == -3]}else{0},
    fod_m2 = if (-2 %in% fod$difference) {fod$frequency[fod$difference == -2]}else{0},
    fod_m1 = if (-1 %in% fod$difference) {fod$frequency[fod$difference == -1]}else{0},
    fod_0  = if ( 0 %in% fod$difference) {fod$frequency[fod$difference ==  0]}else{0},
    fod_p1 = if (+1 %in% fod$difference) {fod$frequency[fod$difference == +1]}else{0},
    fod_p2 = if (+2 %in% fod$difference) {fod$frequency[fod$difference == +2]}else{0},
    fod_p3 = if (+3 %in% fod$difference) {fod$frequency[fod$difference == +3]}else{0},
    fod_p4 = if (+4 %in% fod$difference) {fod$frequency[fod$difference == +4]}else{0},
    fod_p5 = if (+5 %in% fod$difference) {fod$frequency[fod$difference == +5]}else{0}
  )
}
