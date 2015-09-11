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
arma::mat MAD(const arma::mat & TrainCor, const arma::mat & TestCor,bool isMedian, const arma::ivec  dim){
  //First generate the abs error

  arma::mat errmat = arma::abs(TrainCor-TestCor);

  arma::vec madvec(TrainCor.n_rows);
  if(isMedian){
    madvec = arma::median(errmat,1);
  }
  else{
    madvec = arma::mean(errmat,1);
  }

  int rownum = arma::as_scalar(dim[0]);
  int colnum = arma::as_scalar(dim[1]);
  arma::mat mad = arma::reshape(madvec,rownum,colnum);
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
  Rcpp::Rcout<<"Starting Bootstrap!"<<std::endl;
  arma::mat C(A.n_cols*B.n_cols,Bootmat.n_rows);
  arma::mat tA(A.n_rows,A.n_cols);
  arma::mat tB(B.n_rows,B.n_cols);
  arma::mat tC(A.n_rows,B.n_rows);
  for(int i=0; i<bsi; i++){
    tA = A.rows(Bootmat.row(i));
    tB = B.rows(Bootmat.row(i));
    tC = cor(tA,tB);
    C.col(i) = vectorise(tC,0);
  }
  C.elem(find_nonfinite(C)).zeros();

 return reshape(median(C,1),A.n_cols,B.n_cols);
}

// [[Rcpp::export]]
arma::mat BeQTL2(const arma::mat & A, const arma::mat & B, const arma::umat & Bootmat){
  int bsi= Bootmat.n_rows;
  arma::mat C(A.n_cols*B.n_cols,Bootmat.n_rows);
  arma::mat tC(A.n_rows,B.n_rows);
  for(int i=0; i<bsi; i++){
    tC = cor(A.rows(Bootmat.row(i)),B.rows(Bootmat.row(i)));
    C.col(i) = vectorise(tC,0);
  }
  C.elem(find_nonfinite(C)).zeros();

  return reshape(median(C,1),A.n_cols,B.n_cols);
}

template<typename T>
T mod(T a, int n)
{
  return a - arma::floor(a/n)*n;
}


//Function to return row and column indices given a matrix and vectorized indices
// [[Rcpp::export]]
arma::umat Ind(const int nrow, const arma::uvec & ind){
  arma::umat rc(ind.n_elem,2);
  rc.col(0)= mod(ind,nrow);
  rc.col(1)= arma::floor(ind/nrow);


  return(rc);

}

// [[Rcpp::export]]
arma::mat Nind(const arma::mat inpmat,const arma::uvec rowindex,const arma::uvec colindex){
  arma::mat copmat = inpmat(rowindex,colindex);
  return copmat;

}



//Function to summarize results from BeQTL and return a dataframe
// [[Rcpp::export]]
Rcpp::DataFrame SumRes(const arma::mat  & cormat, const arma::mat & errmat, const Rcpp::DataFrame SnpDF, const Rcpp::DataFrame Genedf, const int samplesize, const double tcutoff){
  Rcpp::Rcout<<"Generating Tmat"<<std::endl;
  arma::mat tmat = sqrt(samplesize-2)*(cormat/(sqrt(1-square(cormat))));
  Rcpp::Rcout<<"Finding strong t"<<std::endl;
  arma::uvec goods = find(abs(tmat)>tcutoff);
  Rcpp::Rcout<<"Generating matrix index"<<std::endl;
  arma::umat goodmat = Ind(tmat.n_rows,goods);
  Rcpp::Rcout<<"Subsetting tmat"<<std::endl;
  Rcpp::Rcout<<"This many good results"<<goodmat.n_rows<<std::endl;
  arma::vec tvec = tmat(goods);
  Rcpp::Rcout<<"Subsetting Errmat"<<std::endl;
  arma::vec errvec = errmat(goods);
  Rcpp::Rcout<<"Generating SNP and Gene lists"<<std::endl;
  Rcpp::IntegerVector GoodGenes = Rcpp::wrap(arma::conv_to<arma::ivec>::from(goodmat.col(0)));
  Rcpp::IntegerVector GoodSNPs = Rcpp::wrap(arma::conv_to<arma::ivec>::from(goodmat.col(1)));

//Subset SNP anno

  Rcpp::CharacterVector SNPnames = SnpDF["rsid"];
  SNPnames = SNPnames[GoodSNPs];
  arma::ivec SNPchr = Rcpp::as<arma::ivec>(SnpDF["Chrom"]);
  SNPchr = SNPchr(goodmat.col(1));
  arma::ivec SNPpos = Rcpp::as<arma::ivec>(SnpDF["Pos"]);
  SNPpos = SNPpos(goodmat.col(1));
//Subset Geneanno
  Rcpp::CharacterVector GeneNames = Genedf["Symbol"];
  GeneNames = GeneNames[GoodGenes];
  arma::ivec Genechr = Rcpp::as<arma::ivec>(Genedf["Chrom"]);
  Genechr = Genechr(goodmat.col(0));
  arma::ivec Genestart = Rcpp::as<arma::ivec>(Genedf["Start"]);
  Genestart = Genestart(goodmat.col(0));
  arma::ivec Genestop = Rcpp::as<arma::ivec>(Genedf["Stop"]);
  Genestop = Genestop(goodmat.col(0));

  arma::ivec CisDist(GoodGenes.length());
  Rcpp::Rcout<<"Calculating Cisdist"<<std::endl;
  CisDist = arma::min(arma::join_cols(abs(Genestop-SNPpos),abs(Genestart-SNPpos)),1);
  Rcpp::Rcout<<"CisDist Calculated"<<std::endl;
  CisDist.elem(find(Genechr!=SNPchr)).fill(-1);
  return  Rcpp::DataFrame::create(Rcpp::Named("SNP")=SNPnames,
                                  Rcpp::Named("Gene")=GeneNames,
                                  Rcpp::Named("t-stat")=Rcpp::wrap(tvec),
                                  Rcpp::Named("err")=Rcpp::wrap(errvec),
                                  Rcpp::Named("CisDist")=Rcpp::wrap(CisDist));

}


