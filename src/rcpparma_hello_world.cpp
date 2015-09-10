// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//



//Function to generate a single list of two matrices (train or test), given an index
//  [[Rcpp::export]]
Rcpp::List MatSplit( const arma::mat & A, const arma::mat & B, const arma::uvec & splitind){
  arma::mat tA = A.rows(splitind);
  arma::mat tB = B.rows(splitind);
  return Rcpp::List::create(Rcpp::Named("A")=tA,
                            Rcpp::Named("B")=tB);

}

//Function to estimate RMSE from point estimate, bootstrap estimate, and heldout estimate (across all K, for one feature-pair)
// [[Rcpp::export]]
arma::mat RMSE(const arma::Mat<double> & TrainCor, const arma::Mat<double> & TestCor, const Rcpp::NumericVector & dim){
  arma::mat RMSE = square(TrainCor-TestCor);
  RMSE = sum(RMSE,1);
  RMSE = sqrt(RMSE/TestCor.n_cols);
  RMSE.reshape(dim[0],dim[1]);
  return RMSE;
}


//Function to estimate Median and Mean Absolute Deviation from point estimate, bootstrap estiamte, and heldout estimate
// [[Rcpp::export]]
arma::vec MAD(const arma::mat & TrainCor, const arma::mat & TestCor,bool isMedian, const Rcpp::NumericVector & dim){
  //First generate the abs error
  arma::mat errmat = abs(TrainCor-TestCor);
  arma::vec mad(TrainCor.n_rows);
  if(isMedian){
    mad = median(errmat,1);
  }
  else{
    mad = mean(errmat,1);
  }
  mad.reshape(dim[0],dim[1]);
  return(mad);
}

//Function to perform point estimate of correlation
// [[Rcpp::export]]
arma::mat PointCor(const arma::mat & A, const arma::mat & B){
  arma::mat C= cor(A,B);
  C.elem(find_nonfinite(C)).zeros();
  return C;
}

//Function to generate matrix of bootstrap replicates.
// [[Rcpp::export]]
arma::umat GenBoot(const int samplesize, const int bootstrapnumber){
  arma::umat samplemat = arma::randi<arma::umat>(bootstrapnumber,samplesize,arma::distr_param(0,samplesize-1));
  return(samplemat);
}

//Script that takes two matrices, performs bootstrapped correlation, and returns the median
// [[Rcpp::export]]
arma::mat BeQTL(const arma::mat & A, const arma::mat & B, const arma::umat & Bootmat){
  int bsi= Bootmat.n_rows;
  arma::mat C(A.n_cols*B.n_cols,Bootmat.n_rows);

  arma::mat tA(A.n_rows,A.n_cols);
  arma::mat tB(B.n_rows,B.n_cols);
  arma::mat tC(A.n_rows,B.n_rows);
  for(int i=0; i<bsi; i++){
    tA = A.rows(Bootmat.row(i));
    tB = B.rows(Bootmat.row(i));
    tC = cor(tA,tB);
    tC.elem(find_nonfinite(tC)).zeros();
    C.col(i) = vectorise(tC,0);
  }

 return reshape(median(C,1),A.n_cols,B.n_cols);
}

template<typename T>
T mod(T a, int n)
{
  return a- arma::floor(a/n)*n;
}

//Function to return row and column indices given a matrix and vectorized indices
// [[Rcpp::export]]
arma::mat Ind(const int nrow, const int ncol, const arma::uvec & ind){
  arma::mat rc(ind.n_elem,2);
  arma::uvec::const_iterator b = ind.begin();
  arma::uvec::const_iterator e = ind.end();
  int it=0;
  rc.col(0)= floor(ind/nrow);
  rc.col(1)= mod(ind,nrow)

  rc.col(0) = ind
}

//Function to summarize results from BeQTL and return a dataframe
// [[Rcpp::export]]
Rcpp::DataFrame SumRes(const arma::mat  & cormat, const arma::mat & errmat, const Rcpp::CharacterVector cortype, const Rcpp::CharacterVector errtype, const Rcpp::DataFrame SnpDF, const Rcpp::DataFrame Genedf, const int samplesize, const double tcutoff){

  arma::mat tmat = sqrt(samplesize-2)*(cormat/(sqrt(1-square(cormat))));
  arma::uvec goods = find(abs(tmat)>tcutoff);


}


