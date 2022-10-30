#include <RcppArmadillo.h>

using namespace Rcpp;



RcppExport SEXP linearcppaopt(SEXP z, SEXP v) {
  
  Rcpp::NumericMatrix zr(z);
  Rcpp::NumericMatrix vr(v);

  
  int nb = zr.ncol();
  int n = zr.nrow();
  
  arma::mat Z(zr.begin(), n, nb, false);
  arma::mat V(vr.begin(), nb, nb, false);

    
  #define ARMA_DONT_PRINT_ERRORS
  #include<exception>
  #include <armadillo>
  
  std::ostream nullstream(0);
  arma::set_cerr_stream(nullstream);
  
  arma::mat ZWZ = arma::zeros(nb,nb);
  arma::mat ZWZv = arma::zeros(nb,nb);
  arma::mat ZWZi = arma::zeros(nb,nb);
  
  double out;
    
    ZWZ.zeros();
    for (int ii=0; ii<nb; ii++) {
      for (int jj=ii; jj<nb; jj++) {
        for (int kk=0; kk<n; kk++) {
          ZWZ(ii,jj) += Z(kk,ii)*Z(kk,jj);}
        ZWZ(jj,ii) = ZWZ(ii,jj);}}

    ZWZv = ZWZ + V;
    
    bool flag = inv_sympd(ZWZi,ZWZv);
    
    if(flag==1){
    out = trace(ZWZi);
	} else {
    out = 100000;}

  return as<NumericVector>(wrap(out));
  
}





RcppExport SEXP linearcppdopt(SEXP z, SEXP v) {
  
  Rcpp::NumericMatrix zr(z);
  Rcpp::NumericMatrix vr(v);
  
  int nb = zr.ncol();
  int n = zr.nrow();
  
  arma::mat Z(zr.begin(), n, nb, false);
  arma::mat V(vr.begin(), nb, nb, false);

  #define ARMA_DONT_PRINT_ERRORS
  #include<exception>
  #include <armadillo>
  
  std::ostream nullstream(0);
  arma::set_cerr_stream(nullstream);
  
  arma::mat ZWZ = arma::zeros(nb,nb);
  arma::mat ZWZv = arma::zeros(nb,nb);
  arma::mat ZWZi = arma::zeros(nb,nb);
  
  double out;
  double eout;
  double outlogdet;
    
    ZWZ.zeros();
    for (int ii=0; ii<nb; ii++) {
      for (int jj=ii; jj<nb; jj++) {
        for (int kk=0; kk<n; kk++) {
          ZWZ(ii,jj) += Z(kk,ii)*Z(kk,jj);}
        ZWZ(jj,ii) = ZWZ(ii,jj);}}

  ZWZv = ZWZ + V;
    
  bool flag = log_det_sympd(outlogdet,ZWZv);

if(flag==1){
    out = - outlogdet / nb;
} else {
    out = log(100000);}

eout = exp(out);

  return as<NumericVector>(wrap(eout));
  
}





RcppExport SEXP logisticcppdopt(SEXP z, SEXP v, SEXP a, SEXP w) {
  
  Rcpp::NumericMatrix zr(z);
  Rcpp::NumericMatrix vr(v);
  Rcpp::NumericMatrix ar(a);
  Rcpp::NumericVector wr(w);
  
  int nb = zr.ncol();
  int n = zr.nrow();
  int M = ar.nrow();
  int p = ar.ncol();
  
  arma::mat A(ar.begin(), M, p, false);
  arma::mat Z(zr.begin(), n, nb, false);
  arma::mat V(vr.begin(), nb, nb, false);
  arma::vec weights(wr.begin(), M, false);
  
  #define ARMA_DONT_PRINT_ERRORS
  #include<exception>
  #include <armadillo>
  
  std::ostream nullstream(0);
  arma::set_cerr_stream(nullstream);


  arma::mat At = A.t();

  arma::vec LogDet = arma::zeros(M);
  arma::vec ExpLogDet = arma::zeros(M);
  arma::vec eta = arma::zeros(n);
  arma::vec prob = arma::zeros(n);
  arma::vec W = arma::zeros(n);

  arma::mat ZWZ = arma::zeros(nb,nb);
  arma::mat ZWZv = arma::zeros(nb,nb);

  double out;
  double outlogdet;

  for (int k=0; k<M; k++) {

   eta = Z*At.unsafe_col(k);
   prob = 1/(1+exp(-eta));
   W = prob%(1-prob);

   ZWZ.zeros();
   for (int ii=0; ii<nb; ii++) {
    for (int jj=ii; jj<nb; jj++) {
     for (int kk=0; kk<n; kk++) {
      ZWZ(ii,jj) += W(kk)*Z(kk,ii)*Z(kk,jj);}
      ZWZ(jj,ii) = ZWZ(ii,jj);}}

    ZWZv = ZWZ + V;

    bool flag = log_det_sympd(outlogdet,ZWZv);

    if(flag==1){
     LogDet(k) = - outlogdet / nb;
    } else {
      return wrap(100000);}

  }

  ExpLogDet = exp(LogDet);
  out = dot(ExpLogDet, weights);

  return as<NumericVector>(wrap(out));

}





RcppExport SEXP logisticcppaopt(SEXP z, SEXP v, SEXP a, SEXP w) {
  
  Rcpp::NumericMatrix zr(z);
  Rcpp::NumericMatrix vr(v);
  Rcpp::NumericMatrix ar(a);
  Rcpp::NumericVector wr(w);
  
  int nb = zr.ncol();
  int n = zr.nrow();
  int M = ar.nrow();
  int p = ar.ncol();
  
  arma::mat A(ar.begin(), M, p, false);
  arma::mat Z(zr.begin(), n, nb, false);
  arma::mat V(vr.begin(), nb, nb, false);
  arma::vec weights(wr.begin(), M, false);
  
  #define ARMA_DONT_PRINT_ERRORS
  #include<exception>
  #include <armadillo>
  
  std::ostream nullstream(0);
  arma::set_cerr_stream(nullstream);

  arma::mat At = A.t();

  arma::vec Trace = arma::zeros(M);
  arma::vec eta = arma::zeros(n);
  arma::vec prob = arma::zeros(n);
  arma::vec W = arma::zeros(n);

  arma::mat ZWZ = arma::zeros(nb,nb);
  arma::mat ZWZv = arma::zeros(nb,nb);
  arma::mat ZWZi = arma::zeros(nb,nb);

  double out;

  for (int k=0; k<M; k++) {

   eta = Z*At.unsafe_col(k);
   prob = 1/(1+exp(-eta));
   W = prob%(1-prob);

   ZWZ.zeros();
   for (int ii=0; ii<nb; ii++) {
    for (int jj=ii; jj<nb; jj++) {
     for (int kk=0; kk<n; kk++) {
      ZWZ(ii,jj) += W(kk)*Z(kk,ii)*Z(kk,jj);}
      ZWZ(jj,ii) = ZWZ(ii,jj);}}

    ZWZv = ZWZ + V;

    bool flag = inv_sympd(ZWZi,ZWZv);

    if(flag==1){
     Trace(k) = trace(ZWZi);
    } else {
    return wrap(100000);}

  }

  out = dot(Trace, weights);

  return as<NumericVector>(wrap(out));

}





RcppExport SEXP poissoncppdopt(SEXP z, SEXP v, SEXP a, SEXP w) {
  
  Rcpp::NumericMatrix zr(z);
  Rcpp::NumericMatrix vr(v);
  Rcpp::NumericMatrix ar(a);
  Rcpp::NumericVector wr(w);
  
  int nb = zr.ncol();
  int n = zr.nrow();
  int M = ar.nrow();
  int p = ar.ncol();
  
  arma::mat A(ar.begin(), M, p, false);
  arma::mat Z(zr.begin(), n, nb, false);
  arma::mat V(vr.begin(), nb, nb, false);
  arma::vec weights(wr.begin(), M, false);
  
  #define ARMA_DONT_PRINT_ERRORS
  #include<exception>
  #include <armadillo>
  
  std::ostream nullstream(0);
  arma::set_cerr_stream(nullstream);

  arma::mat At = A.t();

  arma::vec LogDet = arma::zeros(M);
  arma::vec ExpLogDet = arma::zeros(M);
  arma::vec eta = arma::zeros(n);
  arma::vec W = arma::zeros(n);

  arma::mat ZWZ = arma::zeros(nb,nb);
  arma::mat ZWZv = arma::zeros(nb,nb);

  double out;
  double outlogdet;

  for (int k=0; k<M; k++) {
  
   eta = Z*At.unsafe_col(k);
   W = exp(eta);

   ZWZ.zeros();
   for (int ii=0; ii<nb; ii++) {
    for (int jj=ii; jj<nb; jj++) {
     for (int kk=0; kk<n; kk++) {
       ZWZ(ii,jj) += W(kk)*Z(kk,ii)*Z(kk,jj);}
       ZWZ(jj,ii) = ZWZ(ii,jj);}}

    ZWZv = ZWZ + V;

    bool flag = log_det_sympd(outlogdet,ZWZv);
  
    if(flag==1){
      LogDet(k) = - outlogdet / nb;
    } else {
      return wrap(100000);}
    
  }

  ExpLogDet = exp(LogDet);
  out = dot(ExpLogDet, weights);

  return as<NumericVector>(wrap(out));

}





RcppExport SEXP poissoncppaopt(SEXP z, SEXP v, SEXP a, SEXP w) {
  
  Rcpp::NumericMatrix zr(z);
  Rcpp::NumericMatrix vr(v);
  Rcpp::NumericMatrix ar(a);
  Rcpp::NumericVector wr(w);
  
  int nb = zr.ncol();
  int n = zr.nrow();
  int M = ar.nrow();
  int p = ar.ncol();
  
  arma::mat A(ar.begin(), M, p, false);
  arma::mat Z(zr.begin(), n, nb, false);
  arma::mat V(vr.begin(), nb, nb, false);
  arma::vec weights(wr.begin(), M, false);
  
  #define ARMA_DONT_PRINT_ERRORS
  #include<exception>
  #include <armadillo>
  
  std::ostream nullstream(0);
  arma::set_cerr_stream(nullstream);

  arma::mat At = A.t();

  arma::vec Trace = arma::zeros(M);
  arma::vec eta = arma::zeros(n);
  arma::vec W = arma::zeros(n);

  arma::mat ZWZ = arma::zeros(nb,nb);
  arma::mat ZWZv = arma::zeros(nb,nb);
  arma::mat ZWZi = arma::zeros(nb,nb);

  double out;

  for (int k=0; k<M; k++) {

   eta = Z*At.unsafe_col(k);
   W = exp(eta);

   ZWZ.zeros();
   for (int ii=0; ii<nb; ii++) {
    for (int jj=ii; jj<nb; jj++) {
     for (int kk=0; kk<n; kk++) {
      ZWZ(ii,jj) += W(kk)*Z(kk,ii)*Z(kk,jj);}
      ZWZ(jj,ii) = ZWZ(ii,jj);}}

   ZWZv = ZWZ + V;

   bool flag = inv_sympd(ZWZi,ZWZv);

   if(flag==1){
     Trace(k) = trace(ZWZi);
   } else {
    return wrap(100000);}

  }

  out = dot(Trace, weights);

  return as<NumericVector>(wrap(out));

}



