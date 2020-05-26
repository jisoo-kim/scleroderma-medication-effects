//#include <Rcpp.h>
#include <RcppArmadillo.h> // might cause problem, installing gfortran-6.1.pkg from the R website solves the problem: https://stackoverflow.com/questions/35999874/mac-os-x-r-error-ld-warning-directory-not-found-for-option
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double meanC(Rcpp::NumericVector x) {
  int n = x.size();
  double total = 0;
  
  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  return total / n;
}

// [[Rcpp::export]]
Rcpp::NumericVector rowSumsC(Rcpp::NumericMatrix x) {
  int nrow = x.nrow(), ncol = x.ncol();
  Rcpp::NumericVector out(nrow);
  
  for (int i = 0; i < nrow; i++) {
    double total = 0;
    for (int j = 0; j < ncol; j++) {
      total += x(i, j);
    }
    out[i] = total;
  }
  return out;
}

// [[Rcpp::export]]
arma::mat sweep(arma::mat x, int k) {
  //James Goodnight pg 154
  int nrow = x.n_rows, ncol = x.n_cols;
  
  double D = x(k-1, k-1), B = 0;
  
  for (int j = 0; j < ncol; j++) {
    x(k-1, j) = x(k-1, j)/D;
  }
  
  for (int i = 0; i < nrow; i++) {
    if(i != k-1){
      B = x(i, k-1);
      for (int j = 0; j < ncol; j++) {
        x(i,j) = x(i,j) - B*x(k-1,j);
      }
      x(i, k-1) = -B/D;
    }
  }
  x(k-1, k-1) = 1/D;
  return x;
}

// [[Rcpp::export]]
arma::mat adjust(arma::mat x, int k) {
  int nrow = x.n_rows, ncol = x.n_cols;
  
  double D = x(k-1, k-1), B = 0;
  
  for (int j = 0; j < ncol; j++) {
    x(k-1, j) = x(k-1, j)/D;
  }
  
  for (int i = 0; i < nrow; i++) {
    if(i != k-1){
      B = x(i, k-1);
      for (int j = 0; j < ncol; j++) {
        x(i,j) = x(i,j) - B*x(k-1,j);
      }
      
    }
  }
  
  return x;
}

// [[Rcpp::export]]
arma::mat sweep_web(arma::mat x, int k) {
  //https://terrytao.wordpress.com/2015/10/07/sweeping-a-matrix-rotates-its-graph/
  int nrow = x.n_rows, ncol = x.n_cols;
  
  double D = x(k-1, k-1), B = 0;
  
  for (int j = 0; j < ncol; j++) {
    x(k-1, j) = x(k-1, j)/D;
  }
  
  for (int i = 0; i < nrow; i++) {
    if(i != k-1){
      B = x(i, k-1);
      for (int j = 0; j < ncol; j++) {
        x(i,j) = x(i,j) - B*x(k-1,j);
      }
      x(i, k-1) = B/D;
    }
  }
  x(k-1, k-1) = -1/D;
  return x;
}


// [[Rcpp::export]]
arma::mat multsweep(arma::mat x, int startc, int endc) {
  
  // sweep A11 for A = (A11, A12, A21, A22)
  // obtain (inv(A11), inv(A11)*A12, -A21*inv(A11), A22 - A21*inv(A11)*A12)
  int nrow = x.n_rows, ncol = x.n_cols;
  
  for(int k = startc ; k <= endc  ; k++){
    double D = x(k-1, k-1), B = 0;
    
    for (int j = 0; j < ncol; j++) {
      x(k-1, j) = x(k-1, j)/D;
    }
    
    for (int i = 0; i < nrow; i++) {
      if(i != k-1){
        B = x(i, k-1);
        for (int j = 0; j < ncol; j++) {
          x(i,j) = x(i,j) - B*x(k-1,j);
        }
        x(i, k-1) = -B/D;
      }
    }
    x(k-1, k-1) = 1/D;
  }
  
  return x;
}

// [[Rcpp::export]]
arma::mat Blockmat (arma::mat A11, arma::mat A12, arma::mat A21, arma::mat A22){
  return(join_cols(join_rows(A11,A12), join_rows(A21,A22)));
}

// [[Rcpp::export]]
arma::vec vec_ex(arma::vec x, int k) {
  // exclude kth element of vector x
  arma::vec y(x.size()-1);
  int ct = 0;
  for(int i = 0; i < x.size(); i++){
    
    if(i != k){
      y(ct) = x(i);
      ct = ct+1;
    }
  }
  return(y);
}

// [[Rcpp::export]]
arma::mat infl_logistic(arma::mat X, arma::mat dlt, arma::mat Y) {
  //returns matrix with column k = dlt - dlt(k), dlt(k) is coef excluding kth data row
  int nrow = X.n_rows, ncol = X.n_cols; 
  
  arma::mat ph = 1/(1+exp(-X*dlt));
  arma::vec p1 = ph % (1-ph);
  arma::mat sqV = diagmat(sqrt(p1));
  arma::mat Xt = sqV *X;
  arma::mat toswp = Blockmat(Xt.t()*Xt, Xt.t(), Xt, arma::zeros(nrow,nrow)); // A22 = zero matrix
  arma::mat H = multsweep(toswp, 1, Xt.n_cols ); // - A21*inv(A11)*A12
  
  arma::mat H1 = -H.submat(ncol, ncol, H.n_cols-1, H.n_cols-1);//A21*inv(A11)*A12
  
  arma::mat resl = arma::zeros(ncol*nrow,ncol);//012 345 678 start at 3*k
  arma::mat res  = arma::zeros(ncol,ncol);
  for(int k = 0; k < nrow; k++){
    resl.submat(ncol*k, 0,ncol*(k+1)-1, ncol-1) =  p1(k) * X.row(k).t() *X.row(k);
    res = res + resl.submat(ncol*k, 0,ncol*(k+1)-1, ncol-1);
  }
  
  arma::mat resk  = arma::zeros(ncol,ncol);
  arma::mat pert = arma::zeros(ncol,nrow);
  for(int k = 0; k < nrow; k++){
    resk = res - resl.submat(ncol*k, 0,ncol*(k+1)-1, ncol-1) ;
    pert.col(k) =  (multsweep(resk, 1, ncol) * X.row(k).t() ) * (Y(k) - ph(k))/(1-H1(k,k));
  }
  
  return(pert);
}

// [[Rcpp::export]]
arma::mat XtX_k(arma::mat X){
  // X(k) is X taking out the kth row, Xk =  X.row(k)
  // X(k).t() % X(k) = sumj Xj.t() % Xj - Xk.t() % Xk
  //XtX_k returns rbind(X1.t() % X1, ...,Xn.t() % Xn, sumj Xj.t() % Xj)
  
  //list 1:(ncol*nrow) , sum (ncol*nrow+1):(ncol*(nrow+1))
  //arma::mat list = arma::zeros(ncol*nrow,ncol);
  //arma::mat sum  = arma::zeros(ncol,ncol);
  //resl = rbind(list,sum) 
  int nrow = X.n_rows, ncol = X.n_cols; 
  arma::mat resl = arma::zeros(ncol*(nrow+1),ncol);
  
  for(int k = 0; k < nrow; k++){
    resl.submat(ncol*k, 0,ncol*(k+1)-1, ncol-1) =  X.row(k).t() *X.row(k);
    resl.submat(ncol*nrow, 0, ncol*(nrow+1)-1, ncol-1) = resl.submat(ncol*nrow, 0, ncol*(nrow+1)-1, ncol-1) + resl.submat(ncol*k, 0,ncol*(k+1)-1, ncol-1);
  }
  return(resl);
}

// [[Rcpp::export]]
arma::mat inflY0(arma::mat V0, arma::mat V0_l,arma::mat Y0, 
               arma::mat dlt, arma::mat dlt_inf,
               arma::mat V, arma::mat V_l){
  // must guarantee that the first nrow0 columns of dlt_inf correspond to V0
  // and the order of samples match, i.e.
  // the ith col of dlt_inf corresponds to the ith row in V0 when i<nrow0
  //                                 (i - nrow0)th row in V1 when i>=nrow0
  int nrow0 = V0.n_rows, ncol0 = V0.n_cols; 
  //output
  arma::mat Ya0(V.n_rows, dlt_inf.n_cols); //dlt_inf.n_cols = nrow0 + nrow1
  //used within loops
  arma::mat dltk(dlt.n_rows, 1);
  arma::mat imat;
  arma::vec y_ex, l, lk, lk_ex; // lk_ex is taking the kth elem out of lk
  arma::mat denomk ;
  arma::mat denomkvec(V.n_rows,1);
  //used across loops
  arma::mat res0 = XtX_k(V0); // !!!!!!! call function XtX_k
  arma::mat res0k  = arma::zeros(ncol0,ncol0);
  arma::mat V0c_k(nrow0, ncol0);
  //coefficients
  arma::mat alpha0k = arma::zeros(ncol0,nrow0);
  arma::mat tau0k = arma::zeros(ncol0,nrow0);
  
  for(int k = 0; k < nrow0; k++){// Influence of D0 on Y(0)
    y_ex = vec_ex(Y0,k);
    
    dltk = dlt - dlt_inf.col(k);
    lk = V0_l * dltk;
    lk_ex = vec_ex(lk,k);
    l = V_l* dltk;
    
    V0c_k = V0;
    V0c_k.shed_row(k); // V(k)
    
    res0k = res0.submat(ncol0*nrow0, 0, ncol0*(nrow0+1)-1, ncol0-1) - res0.submat(ncol0*k, 0,ncol0*(k+1)-1, ncol0-1); 
    imat = multsweep(res0k, 1, ncol0) * V0c_k.t(); // inv(V(k).t()* V(k)) * V(k).t()
    
    alpha0k.col(k) = imat * y_ex;
    tau0k.col(k)   = imat * lk_ex;
    
    denomk = lk_ex.t() * (lk_ex - V0c_k * tau0k.col(k));
    denomkvec.fill(denomk(0,0));//only takes element-wise division...
    
    Ya0.col(k) = V * alpha0k.col(k) + ( (l - V * tau0k.col(k)) * ( lk_ex.t() * (y_ex - V0c_k* alpha0k.col(k)) )  ) / denomkvec;
    
  }
  
  arma::mat alpha0, tau0;
  
  //using the entire V0
  res0k = res0.submat(ncol0*nrow0, 0, ncol0*(nrow0+1)-1, ncol0-1);
  imat = multsweep(res0k, 1, ncol0) * V0.t(); // inv(V(k).t()* V(k)) * V(k).t()
  
  
  for(int k = nrow0; k < dlt_inf.n_cols  ; k++){// Influence of D1 on Y(0)
   
  dltk = dlt - dlt_inf.col(k);
  lk = V0_l * dltk;
  l = V_l* dltk;
  
  alpha0 = imat * Y0;
  tau0   = imat * lk;
  
  denomk = lk.t() * (lk - V0 * tau0);
  denomkvec.fill(denomk(0,0));//only takes element-wise division...
  
  Ya0.col(k) = V * alpha0 + ( (l - V * tau0) * ( lk.t() * (Y0 - V0 * alpha0) )  ) / denomkvec;
  
  }
  return(Ya0);
}

// [[Rcpp::export]]
arma::mat inflY1(arma::mat V1, arma::mat V1_l,arma::mat Y1, 
                 arma::mat dlt, arma::mat dlt_inf,
                 arma::mat V, arma::mat V_l){
  // must guarantee that the last nrow1 columns of dlt_inf correspond to V1
  // and the order of samples match, i.e.
  // the ith col of dlt_inf corresponds to the ith row in V0 when i<nrow0
  //                                 (i - nrow0)th row in V1 when i>=nrow0
  int nrow1 = V1.n_rows, ncol1 = V1.n_cols, nrow0 = dlt_inf.n_cols - nrow1; 
  int n = dlt_inf.n_cols;
  //output
  arma::mat Ya1(V.n_rows, dlt_inf.n_cols);//dlt_inf.n_cols = nrow0 + nrow1
  //used within loops
  arma::mat dltk(dlt.n_rows, 1);
  arma::mat imat;
  arma::vec y_ex, l, lk, lk_ex; // lk_ex is taking the kth elem out of lk
  arma::mat denomk ;
  arma::mat denomkvec(V.n_rows,1);
  //used across loops
  arma::mat res1 = XtX_k(V1); // !!!!!!! call function XtX_k
  arma::mat res1k  = arma::zeros(ncol1,ncol1);
  arma::mat V1c_k(nrow1, ncol1);
  //coefficients
  arma::mat alpha1k = arma::zeros(ncol1,nrow1);
  arma::mat tau1k = arma::zeros(ncol1,nrow1);
  
  for(int k = nrow0; k < n; k++){// Influence of D1 on Y(1)
    y_ex = vec_ex(Y1,k-nrow0);
    
    dltk = dlt - dlt_inf.col(k);
    lk = V1_l * dltk;
    lk_ex = vec_ex(lk,k-nrow0);
    l = V_l* dltk;
    
    V1c_k = V1;
    V1c_k.shed_row(k-nrow0); // V(k)
    
    res1k = res1.submat(ncol1*nrow1, 0, ncol1*(nrow1+1)-1, ncol1-1) - res1.submat(ncol1*(k-nrow0), 0,ncol1*(k-nrow0+1)-1, ncol1-1); 
    imat = multsweep(res1k, 1, ncol1) * V1c_k.t(); // inv(V(k).t()* V(k)) * V(k).t()
    
    alpha1k.col(k-nrow0) = imat * y_ex;
    tau1k.col(k-nrow0)   = imat * lk_ex;
    
    denomk = lk_ex.t() * (lk_ex - V1c_k * tau1k.col(k-nrow0));
    denomkvec.fill(denomk(0,0));//only takes element-wise division...
    
    Ya1.col(k) = V * alpha1k.col(k-nrow0) + ( (l - V * tau1k.col(k-nrow0)) * ( lk_ex.t() * (y_ex - V1c_k* alpha1k.col(k-nrow0)) )  ) / denomkvec;
    
  }
  
  arma::mat alpha1, tau1;
  
  //using the entire V1
  res1k = res1.submat(ncol1*nrow1, 0, ncol1*(nrow1+1)-1, ncol1-1);
  imat = multsweep(res1k, 1, ncol1) * V1.t(); // inv(V(k).t()* V(k)) * V(k).t()
  //Rcpp:: Rcout << "h6 " << std::endl ;
  for(int k = 0; k < dlt_inf.n_cols - nrow1 ; k++){// Influence of D0 on Y(1)
    
    dltk = dlt - dlt_inf.col(k);
    lk = V1_l * dltk;
    l = V_l* dltk;
    
    alpha1 = imat * Y1;
    tau1   = imat * lk;
    
    denomk = lk.t() * (lk - V1 * tau1);
    denomkvec.fill(denomk(0,0));//only takes element-wise division...
    
    Ya1.col(k) = V * alpha1 + ( (l - V * tau1) * ( lk.t() * (Y1 - V1 * alpha1) )  ) / denomkvec;
    
  }
  return(Ya1);
}

// [[Rcpp::export]]      
NumericVector a5 (arma::mat dlt) {
  Rcpp::Environment env=Rcpp::Environment::global_env();
  Rcpp::Function invlgt("invlgt"); // fun_in_glob() is a function defined existing in global env
  
  //Rcpp::Function myS = myEnv["sigma2"];
  NumericVector y = invlgt(dlt);
  return(y) ;
}                  

// [[Rcpp::export]]      
Rcpp::NumericMatrix a6 (Rcpp::NumericMatrix a, Rcpp::NumericMatrix b) {
  arma::mat d = Rcpp::as<arma::mat>(a);
  arma::mat e = Rcpp::as<arma::mat>(b);
  arma::mat f = d*e;
  Rcpp::NumericMatrix res = Rcpp::wrap(f);
  return(res) ;
} 

// [[Rcpp::export]]      
arma::mat a7 (arma::vec a) {
  return(a*a.t()) ; // matrix;  a.t()*a is scalar
}       
       

// [[Rcpp::export]]             
double timesTwo(double a){
  return(a*2);
}              
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


/*** R
timesTwo(42)
*/
