// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
//' @import Rcpp
//' @import RcppArmadillo

SEXP spmat_to_dgCMatrix(arma::sp_mat& spmat) {
  Rcpp::S4 dgCMatrix("dgCMatrix");
  
  // Convert the sp_mat to dgCMatrix
  int n = spmat.n_rows;
  int nc = spmat.n_cols;
  int nnz = spmat.n_nonzero;
  Rcpp::IntegerVector i(nnz);
  Rcpp::IntegerVector p(nc + 1);
  Rcpp::NumericVector x(nnz);
  // Fill in the 'i' and 'x' slots
  int pos = 0;
  for (int j = 0; j < nc; j++) {
    p[j] = pos;
    for (arma::sp_mat::const_iterator it = spmat.begin_col(j); it != spmat.end_col(j); ++it) {
      i[pos] = it.row();
      x[pos] = *it;
      pos++;
    }
  }
  p[nc] = pos;
  // Set the slots of the dgCMatrix object
  dgCMatrix.slot("i") = i;
  dgCMatrix.slot("p") = p;
  dgCMatrix.slot("x") = x;
  dgCMatrix.slot("Dim") = Rcpp::IntegerVector::create(n, nc);
  return dgCMatrix;
}

// [[Rcpp::export]]
Rcpp::S4 GenerateMatrixC( int K_B, int p ) {

  int n = (K_B+1)*p ;
  arma::umat locations (2, (2*K_B+1)*p) ;
  arma::vec values ((2*K_B+1)*p) ;
  int i = 0 ;
  
  for (int u=0; u<(K_B+1)*p; u++){
    locations(0,i) = u ;
    locations(1,i) = u ;
    if(u<p){
      values(i) = 1 ;
    }
    else{
      values(i) = -1 ;
    }
    i++ ;
  }
  for (int j=0; j<p; j++){
    for (int k=1; k<=K_B; k++){
      locations(0,i) = k*p+j ;
      locations(1,i) = j ;
      values(i) = 1 ;
      i++ ;
    }
  }
  arma::sp_mat C(locations, values, n, n) ;
  return spmat_to_dgCMatrix(C) ;
}