#include <Rcpp.h>
using namespace Rcpp;


inline double MAX2(double x, double y) {
  if (x > y){
    return x;
  }
  else{
    return y;
  }
}

inline double MAX3(double x, double y, double z) {
  return MAX2(MAX2(x,y), MAX2(y,z));
}

inline double ARGMAXI2(double x, double y) {
  if (x > y){
    return 1.0;
  }
  else{
    return 2.0;
  }
}

inline double ARGMAXD2(double x, double y) {
  if (x > y){
    return 1.0;
  }
  else{
    return 3.0;
  }
}

inline double ARGMAXID2(double x, double y) {
  if (x > y){
    return 2.0;
  }
  else{
    return 3.0;
  }
}

inline double ARGMAX3(double x, double y, double z) {
  double k = ARGMAXI2(x,y);
  double l = ARGMAXD2(x,z);
  if (k == l){
    return 1.0;
  } else if (k == 2.0 && l == 1.0){
    return 2.0;
  } else if (k == 1.0 && l == 3.0){
    return 3.0;
  } else{
    return ARGMAXID2(y,z);  
  }
}
  
inline NumericVector nucleotideToInt(std::string x) {
  
  double n = x.size();  
  CharacterVector read(n);
  for( int i=0; i<n; i++ ) {
      read[i] = x.substr( i, 1);
  }
  
  NumericVector readInt(n);
  for (int i = 0; i < n; i++) {
    if (read[i] == "A"){
      readInt[i] = 0;
    } 
    else if (read[i] == "C"){
      readInt[i] = 1;
    } 
    else if (read[i] == "G"){
      readInt[i] = 2;
    }
    else if (read[i] == "T"){
      readInt[i] = 3;
    }
    else{
      readInt[i] = 4;
    }
  }
  return readInt;
}


inline CharacterVector IntToState(NumericVector x, int t) {
  CharacterVector path(t);
  for (int i = 0; i < t; i++) {
    if (x[i] == 1){
      path[t-i-1] = "M";
    } else if (x[i] == 2) {
      path[t-i-1] = "I";
    } else {
      path[t-i-1] = "D";
    }
  }
  return path;
}

//' Generate trace back.
//' 
//' This function returns the  trace back, i.e. the sequence of most probable 
//' states in the alignment between read x and reference sequence y.
//'
//' @param q double, index of the state for the last two nucleotides in the read and reference sequence.
//' @param i double, length of the read.
//' @param j double, length of the reference sequence.
//' @param Tm matrix, GPHMM probabilities at each location in the read and reference sequence in the M state.
//' @param Tx matrix, GPHMM probabilities at each location in the read and reference sequence in the insertion state.
//' @param Ty matrix, GPHMM probabilities at each location in the read and reference sequence in the deletion state.
//'
CharacterVector traceBack(double q, double i, double j, 
NumericMatrix Tm, NumericMatrix Tx, NumericMatrix Ty) {
  
  Rcpp::NumericVector path(i + j);
  int t = 0;
  path[t] = q;
  while (i > 1 && j > 1) {
    if (path[t] == 1) {
      path[t+1] = Tm(i,j);
      i = i - 1;
      j = j - 1;
    } else if (path[t]==2) {
      path[t+1] = Tx(i,j);
      i = i - 1;
    } else {
      path[t+1] = Ty(i,j);
      j = j - 1;
    }
    t = t + 1;
  }
  CharacterVector pathChar = IntToState(path, t+1); 
  return pathChar;
}




//' Calculate GPHMM probability.
//' 
//' This function returns the GPHMM probability that a read x could have been
//' sequenced from a reference sequence y.
//'
//' @param x string, with the sequence of the read.
//' @param y string, with the sequence of the reference.
//' @param tau double, probability to transition from any state to end state.
//' @param pp matrix, emission probabilities in the M state.
//' @param qX vector, emission probabilities in the insertion state.
//' @param qY vector, emission probabilities in the deletion state.
//' @param dX double, transition probability from the M to the insertion state.
//' @param dY double, transition probability from the M to the deletion state.
//' @param eX double, transition probability from the insertion to the insertion state.
//' @param eY double, transition probability from the deletion to the deletion state.
//' 
//' @export
//' @examples
//' param <- initializeGphmm()
//' tau <- param[['tau']]
//' pp <- param[['pp']]
//' qX <- param[['qX']]
//' qY <- param[['qY']]
//' dX <- 1/(1+exp(-sum(param[['deltaX']] * c(1, 20))))
//' dY <- 1/(1+exp(-sum(param[['deltaY']] * c(1, 20))))
//' eX <- param[['epsX']]
//' eY <- param[['epsY']]
//' calculategphmm('ATCG', 'ATGG', tau, pp, qX, qY, dX, dY, eX, eY)
// [[Rcpp::export]]
List calculategphmm(std::string x, std::string y,
double tau, NumericMatrix pp, NumericVector qX, NumericVector qY,
double dX, double dY, double eX, double eY) {
  
  double n = x.size(); //length of read x
  double m = y.size(); //length of template y
  
  // convert nucleotides "A", "C", "G", "T" to integer 0, 1, 2, 3
  // to be easily called in matrices s, dX, dY
  NumericVector readInt = nucleotideToInt(x);
  NumericVector templInt = nucleotideToInt(y);
  
  //initialize matrices
  NumericMatrix Vm(n+1,m+1);
  NumericMatrix Vx(n+1,m+1);
  NumericMatrix Vy(n+1,m+1);
  NumericMatrix Tm(n+1,m+1);
  NumericMatrix Tx(n+1,m+1);
  NumericMatrix Ty(n+1,m+1);
  const double minusInf = -10000.0;
  const double zero = 0.0;
  
  Vx(1, 1) = minusInf;
  Vy(1, 1) = minusInf;
  Tx(1, 1) = zero;
  Ty(1, 1) = zero;
  NumericVector initCol = rep( minusInf, m+1);
  Vm(0,_) = initCol;
  Vx(0,_) = initCol;
  Vy(0,_) = initCol;
  NumericVector initColZero = rep( zero, m+1);
  Tm(0,_) = initColZero;
  Tx(0,_) = initColZero;
  Ty(0,_) = initColZero;
  NumericVector initRow = rep( minusInf, n+1);
  Vm(_,0) = initRow;
  Vx(_,0) = initRow;
  Vy(_,0) = initRow;
  NumericVector initRowZero = rep( zero, n+1);
  Tm(_,0) = initRowZero;
  Tx(_,0) = initRowZero;
  Ty(_,0) = initRowZero;
  
  //viterbi algo
  for (int i = 1; i < n+1; i++){
    for (int j = 1; j < m+1; j++){
      if (i == 1 & j == 1){
        Vm(i, j) = zero;
        Tm(i, j) = 1.0;
      }
      else{
        Vm(i, j) = pp(readInt[i-1], templInt[j-1]) + MAX3((Vm(i-1, j-1) + log(1- dX - dY - tau)), (Vx(i-1, j-1) + log(1- eX - tau)), (Vy(i-1, j-1) + log(1- eY - tau)));
        Vx(i, j) = qX[readInt[i-1]] + MAX2((Vm(i-1, j) + log(dX)), (Vx(i-1, j) + log(eX)));
        Vy(i, j) = qY[templInt[j-1]] + MAX2((Vm(i, j-1) + log(dY)), (Vy(i, j-1) + log(eY)));
        Tm(i, j) = ARGMAX3((Vm(i-1, j-1) + log(1- dX - dY - tau)), (Vx(i-1, j-1) + log(1- eX - tau)), (Vy(i-1, j-1) + log(1- eY - tau)));
        Tx(i, j) = ARGMAXI2((Vm(i-1, j) + log(dX)), (Vx(i-1, j) + log(eX)));
        Ty(i, j) = ARGMAXD2((Vm(i, j-1) + log(dY)), (Vy(i, j-1) + log(eY)));
      }
    }
  }
  //termination
  double V = log(tau) + MAX3(Vm(n, m), Vx(n, m), Vy(n, m));
  double Q = ARGMAX3(Vm(n, m), Vx(n, m), Vy(n, m));
  //traceback
  CharacterVector trace = traceBack(Q, n, m, Tm, Tx, Ty);
  
  return Rcpp::List::create(Rcpp::Named("V", V), Rcpp::Named("path", trace));
}





