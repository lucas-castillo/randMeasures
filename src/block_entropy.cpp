#include <Rcpp.h>
using namespace Rcpp;
bool strings_equal(String a, String b){
  return (a == b);
}
LogicalVector vector_string_equal(StringVector a, String b){
  LogicalVector result(a.size());
  for (int i = 0; i < a.size(); i++){
    result(i) = strings_equal(a(i), b);
  }
  return(result);
}
String paste0(NumericVector x, String collapse="."){
  String result("");
  for (int i=0; i < x.size(); i++){
    result += (String)((x(i)));
    if (i != (x.size() - 1)){
      result += collapse;
    }
  }
  return(result);
}
NumericVector subset(NumericVector x, int start, int end){
  NumericVector result(end - start + 1);
  for (int i=start; i <= end; i++){
    result(i-start) = x(i);
  }
  return(result);
}

NumericVector relative_frequencies(NumericVector v){
  NumericVector unique_values = sort_unique(v);
  NumericVector counts(unique_values.size());
  for (int i = 0; i < unique_values.size(); i++){
    int uv = unique_values(i);
    counts(i) = sum(uv == v);
  }
  return(counts / v.size());
}
NumericVector relative_frequencies_s(CharacterVector v){
  StringVector unique_values = sort_unique(v);
  NumericVector counts(unique_values.size());
  for (int i = 0; i < unique_values.size(); i++){
    String uv = unique_values(i);
    counts(i) = sum(vector_string_equal(v, uv));
  }
  return(counts / v.size());
}
NumericVector vlog2(NumericVector x){
  NumericVector res(x.size());
  for (int i = 0; i < x.size(); i++){
    res(i) = std::log2(x(i));
  }
  return(res);
}

//' Calculate block entropy (as in randfindR package)
//' @param x sequence
//' @param block_size length of blocks in which the original sequence should be divided for analysis
//' @export
// [[Rcpp::export]]
double block_entropy(NumericVector x, int block_size = 1){
  // compute entropy directly from observed frequencies for block_size = 1
  if (block_size == 1) {
    NumericVector p_xi = relative_frequencies(x);
    double result = sum(p_xi * vlog2(p_xi)) * (-1);
    return(result);
  }
  int max_index = x.size() - block_size;
  CharacterVector blocks(max_index + 1);
  for (int i = 0; i <= max_index; i++){
    NumericVector block = subset(x, i, i+block_size-1);
    String block_string = paste0(block);
    blocks(i) = block_string;
  }

  NumericVector p_xi = relative_frequencies_s(blocks);
  double result = sum(p_xi * vlog2(p_xi)) * (-1);
  return(result);
}
