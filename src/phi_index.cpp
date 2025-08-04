#include <Rcpp.h>
using namespace Rcpp;


NumericVector sliceVector(NumericVector x, int start, int end){
  NumericVector result(end - start + 1);

  for (int i = start; i <= end; i++){
    result(i-start) = x(i);
  }
  return(result);
}

double contingMatrixToPhi(NumericMatrix M, int setSize, int Nresp) {
  double a1 = M(0,0);
  double a2 = M(0,1);
  double a3 = M(1,0);
  double a4 = M(1,1);

  double Den = sum(M);

  // if(Den > 0){
  double AA1 = ((a1 + a2) * (a1 + a3)) / Den;
  double AA2 = ((a2 + a4) * (a1 + a2)) / Den;
  double AA3 = ((a3 + a4) * (a1 + a3)) / Den;
  double AA4 = ((a2 + a4) * (a3 + a4)) / Den;
  // }

  double B1 =(pow(a1-AA1, 2) / AA1);
  double B2 =(pow(a2-AA2, 2) / AA2);
  double B3 =(pow(a3-AA3, 2) / AA3);
  double B4 =(pow(a4-AA4, 2) / AA4);

  double s  = 1;
  if (a2>a4){
    s = -1;
  }
  return s*(sqrt((B1 + B2 + B3 + B4)/(setSize*Nresp))*100);
}

NumericVector classify_window(NumericVector y, int A, int B, int C){
  NumericVector O(A);
  int Pos = A;
  for (int D = 0; D <= (A-1); D++){
    if (y(C+D-1) == B){
      O(Pos-1) = 1;
    }
    Pos -= 1;
  }
  return rev(O);
}

LogicalVector whichRowEqualsVector(NumericVector x, NumericMatrix M){
  LogicalVector result(M.nrow());
  for (int i = 0; i < M.nrow(); i++){
    NumericVector temp = M(i, _);
    result(i) = is_true(all(temp == x));
  }
  return result;
}

NumericVector tallyOccurrences(
    NumericVector y,
    NumericMatrix options,
    int A, int B, int Nresp
){
  NumericVector frequencies(options.nrow());

  for (int C=1; C <= ((Nresp-A)+1); C++){
    NumericVector O = classify_window(y, A, B, C);
    LogicalVector index = whichRowEqualsVector(O, options);

    for (int i = 0; i < index.size(); i++){
      if (index(i)){
        frequencies(i) += 1;
      }
    }
  }
  return frequencies;
}

//[[Rcpp::export]]
NumericVector phi_index_cpp(
    NumericVector seq,
    int minScale, int maxScale,
    int maxOrder,
    List list_options
){
  NumericVector phi_indices(maxOrder - 1);
  int n = seq.size();
  int setSize = maxScale - minScale + 1;

  List phiArrayList(maxOrder);
  List phiObsPredList(maxOrder);
  List phiObsFreqList(maxOrder);
  List phiList(maxOrder);

  // Compute base frequencies (for order A=1)
  IntegerVector alternatives = seq_len(setSize) + minScale - 1;
  IntegerVector baseFreq(alternatives.size());
  for (int i=0; i<alternatives.size(); i++){
    int opt = alternatives(i);
    LogicalVector temp = seq == opt;
    baseFreq(i) = sum(temp);
  }
  IntegerVector baseInvFreq = n - baseFreq;
  // Compute frequencies and phi for other orders
  for (int A = 2; A <= maxOrder; A++){
    // Matrix
    //        repeat alternate
    // obs        0          0
    // pred       0          0
    NumericMatrix phiObsPred(2,2);
    List phiObsFreq(setSize);
    NumericMatrix options = list_options[A-2];
    NumericVector phiArrayFreq(options.nrow());

    LogicalVector sameInd = options(_, 0) == options(_, A-1);

    for (int B = 1; B<=setSize; B++){
      NumericVector freq = tallyOccurrences(seq, options, A, B, n);
      phiObsFreq[B-1] = freq;
      phiArrayFreq += freq;

      for (int r = 0; r < options.nrow(); r++){
        NumericVector BinString = options( r, _ );

        int freqMid, freq1, freq2;
        double result;
        NumericVector tmp;
        if (A == 2){
          freqMid = n;
          if (BinString(0) == 1){
            freq1 = baseFreq(B-1);
          } else{
            freq1 = baseInvFreq(B-1);
          }
          if (BinString(1) == 1){
            freq2 = baseFreq(B-1);
          } else{
            freq2 = baseInvFreq(B-1);
          }
        } else {
          List oldPhiObsFreq = phiObsFreqList[A - 2];
          tmp = oldPhiObsFreq[B - 1];
          NumericMatrix oldOptions = list_options[A-3];

          NumericVector string1 = sliceVector(BinString, 0, A-2);
          NumericVector string2 = sliceVector(BinString, 1, A-1);

          NumericVector tmp2 = tmp[whichRowEqualsVector(string1, oldOptions)];
          freq1 = tmp2(0);
          tmp2 = tmp[whichRowEqualsVector(string2, oldOptions)];
          freq2 = tmp2(0);

          if (A == 3) {
            if (BinString(1) == 1){
              freqMid = baseFreq(B - 1);
            } else{
              freqMid = baseInvFreq(B - 1);
            }
          } else {
            oldPhiObsFreq = phiObsFreqList[A - 3];
            tmp = oldPhiObsFreq[B - 1];
            NumericMatrix oldOptions = list_options[A-4];

            NumericVector stringMid = sliceVector(BinString, 1, A-2);
            tmp2 = tmp[whichRowEqualsVector(stringMid, oldOptions)];

            freqMid = tmp2(0);

          }
        }
        if (freqMid > 0){
          result = freq1 * freq2 / (double) freqMid;
        } else {
          result = 0;
        }



        if (BinString(0) == BinString(A-1)){
          phiObsPred(1,0) += result;
        } else {
          phiObsPred(1,1) += result;
        }
      }
    }

    phiArrayList[A - 1] = phiArrayFreq;
    phiObsFreqList[A - 1] = phiObsFreq;

    phiObsPred(0,0) = sum((NumericVector) phiArrayFreq[sameInd]);
    phiObsPred(0, 1) = sum((NumericVector) phiArrayFreq[!sameInd]);
    phiObsPredList[A - 1] = phiObsPred;

    double phi = contingMatrixToPhi(phiObsPred, setSize, n);
    phi_indices(A-2) = phi;
  }

  return phi_indices;
}
