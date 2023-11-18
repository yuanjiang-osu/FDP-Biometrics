#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

List CPP_Eval(int p, NumericMatrix T_cov_1, NumericMatrix A, NumericVector mu_1, double alpha, double R_var_1, double V_var_1){
  
  double RV_cov_1 = V_var_1;
  
  Rcpp::Function pnorm("pnorm");
  Rcpp::Function qnorm("qnorm");
  Rcpp::Function pbivnorm("pbivnorm"); 
  
  double c = Rcpp::as<double>(qnorm(1-alpha/2));
  
  int len_AA = (p*p - p)/2;
  
  NumericVector AA(len_AA);
  
  NumericVector Upp_B_1(len_AA);
  
  NumericVector Upp_B_2(len_AA);
  
  NumericVector T_1(len_AA);
  
  NumericVector t_1(p);
  
  NumericVector t_2(p);
  
  NumericVector t_3(p);
  
  int count_d = 0;
  
  for (int m1 = 1; m1<p; m1++){
    
    for(int m2 = 0; m2<m1; m2++){
      
      AA[count_d] = A(m1,m2);
      
      Upp_B_1[count_d] = mu_1[m1];
      
      Upp_B_2[count_d] = mu_1[m2];
      
      count_d++;
      
    }
  }
  
//  Rprintf("count_d = %d, len_AA = %d \n", count_d, len_AA);
  
  T_1 = Rcpp::as<NumericVector>(pbivnorm(c-Upp_B_1, c-Upp_B_2, AA)) - 
    Rcpp::as<NumericVector>(pbivnorm(-c-Upp_B_1, c-Upp_B_2, AA)) -
    Rcpp::as<NumericVector>(pbivnorm(c-Upp_B_1, -c-Upp_B_2, AA)) +
    Rcpp::as<NumericVector>(pbivnorm(-c-Upp_B_1, -c-Upp_B_2, AA));
  
  t_1 = Rcpp::as<NumericVector>(pnorm(c-mu_1));
  t_2 = Rcpp::as<NumericVector>(pnorm(-c-mu_1));
  t_3 = Rcpp::as<NumericVector>(pnorm(-c+mu_1));

  count_d = 0;
  
  for (int m1 = 1; m1<p; m1++){
    
    for(int m2 = 0; m2<m1; m2++){
      
      T_cov_1(m1,m2) = 1-(t_1[m1]-t_2[m1])-(t_1[m2]-t_2[m2]) +
        T_1[count_d]-(t_3[m1]+t_2[m1])*(t_3[m2] + t_2[m2]);
      
      R_var_1 = R_var_1 + 2 * T_cov_1(m1, m2);
      
      if ((mu_1[m1] == 0) && (mu_1[m2] == 0)){
        
        V_var_1 = V_var_1 + 2 * T_cov_1(m1, m2);
        RV_cov_1 = RV_cov_1 + 2 * T_cov_1(m1, m2);
      }
      
      if ((mu_1[m1] * mu_1[m2] == 0) && (mu_1[m1]+mu_1[m2] != 0)){
        RV_cov_1 = RV_cov_1 + T_cov_1(m1, m2);
      }
      
      count_d++;
      
    }
  }
  
  return List::create(R_var_1, V_var_1, RV_cov_1);
  
}

