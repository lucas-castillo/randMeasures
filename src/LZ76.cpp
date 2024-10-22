#include <Rcpp.h>
using namespace Rcpp;


//' Calculate Lempel-Ziv complexity as implemented in Angelike & Musch, 2024
//'
//' @param v sequence
//' @export
// [[Rcpp::export]]
double lz76_complexity(NumericVector v)
{
  int c = 1;
  int l = 0;
  int i = -1;
  int k = 1;
  int k_max = 1;
  int n = v.size();
  bool continue_loop = true;

  while (continue_loop){
    if (v(i + k) == v(l + k)){
      k += 1;
      if ((l+k) >= n){
        c += 1;
        continue_loop = false;
      }
    } else{
      if (k > k_max){
        k_max = k;
      }
      i += 1;

      if (i == l){
        c += 1;
        l = l + k_max;
        if ((l+1) >= n){
          continue_loop = false;
        } else {
          i = 0;
          k = 1;
          k_max = 1;
        }
      } else {
        k = 1;
      }
    }
  }
  return(c);
}
