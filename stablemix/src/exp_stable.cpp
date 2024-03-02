// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <Rcpp.h>
#include <numeric>
#include <math.h>
#include <algorithm>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_spline.h>
using namespace Rcpp;

/*
Helper functions for integration and Hougaard density
*/

struct integrand_params { double alpha; double s; };


// [[Rcpp::export]]
double afun(double x, double alpha){
  return pow(sin(alpha*x)/sin(x),(1/(1-alpha)))*(sin((1-alpha)*x)/sin(alpha*x));
}


/*
Calculating Hougaard density by numeric integration
*/


// [[Rcpp::export]]
double dhint(double u, double s, double alpha){
  return (alpha/(1-alpha))*pow(1/s,1/(1-alpha))*afun(M_PI*u, alpha)*
  			exp(-pow(1/s, alpha/(1-alpha))*afun(M_PI*u, alpha));
}


double integrand(double x, void *p) {
  struct integrand_params *params = (struct integrand_params *) p;
  double alpha = (params->alpha);
  double s = (params->s);
  return dhint(x, s, alpha);
}

// [[Rcpp::export]]
double dps_one(double s, double alpha) {
  gsl_set_error_handler_off();
  gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(1000);
  gsl_function F;
  struct integrand_params params = { alpha, s };
  F.function = &integrand;
  F.params = &params;
  double result, error;

  int intresult = gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000,
                        workspace, &result, &error); 
  gsl_integration_workspace_free (workspace);
  if(ISNAN(intresult)) {
    Rcpp::stop("Integration error");
  }

  return result;
}


// [[Rcpp::export]]
NumericVector dps(NumericVector x, double alpha, bool return_log = true){
	int n = x.size();
	NumericVector result(n);
	for(int i = 0; i < n; i++){
		result[i] = dps_one(x[i], alpha);
	}
  if(return_log){
    return log(result);
  }
	return result;
}

// [[Rcpp::export]]
NumericVector dhexpstab(NumericVector x, double alpha, double delta, 
                        double theta, bool return_log = true, double tol = 1e-2) {
  int n = x.size();
  NumericVector result(n);
  if ((alpha > 1-tol) || (alpha < tol) || (delta <= 0) || (theta < 0)) {
    result = R_NegInf;
    return result;
  }
  double beta = pow(delta / alpha, 1 / alpha);
  result = dps(x / beta, alpha, true) - (theta * x) - 
            (log(beta) - pow(beta * theta, alpha));
  if(return_log){
  	return result;
  }
  else{
  	return exp(result);
  }
}



/*
Calculating exp-Hougaard density by numeric integration
(i.e. If Y ~ Hougaard, this calculates the density of log(Y))
*/


// [[Rcpp::export]]
double dlhint(double u, double s, double alpha){
  return exp(log(alpha) - log(1-alpha) - s*alpha/(1-alpha) + 
          log(afun(M_PI*u, alpha)) - afun(M_PI*u, alpha)*exp(-alpha*s/(1-alpha)));
}


double lintegrand(double x, void *p) {
  struct integrand_params *params = (struct integrand_params *) p;
  double alpha = (params->alpha);
  double s = (params->s);
  return dlhint(x, s, alpha);
}


// [[Rcpp::export]]
double dlps_one(double s, double alpha) {
  gsl_set_error_handler_off();
  gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(1000);
  gsl_function F;
  struct integrand_params params = { alpha, s };
  F.function = &lintegrand;
  F.params = &params;
  double result, error;

  int intresult = gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000,
                        workspace, &result, &error); 
  gsl_integration_workspace_free (workspace);
  if(ISNAN(intresult)) {
    Rcpp::stop("Integration error");
  }

  return result;
}

// [[Rcpp::export]]
NumericVector dlps(NumericVector x, double alpha, bool return_log = true){
	int n = x.size();
	NumericVector result(n);
	for(int i = 0; i < n; i++){
		result[i] = std::max(0.0, dlps_one(x[i], alpha));
	}
  if(return_log){
    return log(result);
  }
	return result;
}

// [[Rcpp::export]]
NumericVector dlhexpstab(NumericVector x, double alpha,  
                        double theta, bool return_log = true, double tol = 1e-2) {
  int n = x.size();
  NumericVector result(n);
  if ((alpha > 1-tol) || (alpha < tol) || (theta < 0)) {
    result = R_NegInf;
    return result;
  }
  result = dlps(x, alpha, true) - (theta * exp(x)) +
                pow(theta, alpha);
   for(int i = 0; i <n; i++){
    if(abs(x[i]) > 150){
      result[i] = R_NegInf;
    }
   }

	if(return_log == true){
 	return result;
 	}
  else{
  	return exp(result);
  }
}

// [[Rcpp::export]]
double llhexpstab(NumericVector x, double alpha,  
                        double theta, bool return_log = true, double tol = 1e-2) {
  return sum(dlhexpstab(x, alpha, theta));
}



// [[Rcpp::export]]
NumericMatrix mklA(NumericMatrix lK, NumericMatrix lB, double alpha){
	int nloc = lK.nrow();
	int nbasis = lK.ncol();
	int nyr = lB.ncol();

	NumericMatrix A(nyr, nloc);
	NumericVector sclKB(nbasis);
	double sclmax = 0;
	for(int i = 0; i < nyr; i ++){
		for(int j = 0; j < nloc; j++){
			sclKB = lK(j,_)/alpha + lB(_,i);
			sclmax = max(sclKB);
			A(i,j) = sclmax + log(sum(exp(sclKB-sclmax)));
		}
	}
	return A;
}


// [[Rcpp::export]]
NumericVector plstabmix(NumericVector lz, double alpha, double theta, 
                  NumericVector lK, bool logp = false){
    int n = lz.size();
    NumericVector result(n);
    if(theta ==0){
        result = - exp(- lz);
    }
    else{
      for(int i = 0; i < n; i++){
        result[i] = - sum(pow(theta + exp((lK - lz[i])/alpha), alpha) - 
                       pow(theta, alpha));
      }
    }
    if (logp) {
      return result;
    }
    else{
      return exp(result);
    }
}

// [[Rcpp::export]]
NumericMatrix plstabmixM(NumericMatrix lz, double alpha, double theta, 
                  NumericMatrix lK, bool logp = false){
    int nloc = lz.ncol();
    int nyear = lz.nrow();
    NumericMatrix result(nyear, nloc);
    if (logp) {
        for(int i = 0; i < nloc; i++){
          result(_,i) = plstabmix(lz(_,i), alpha, theta, lK(i,_), true);
      }
    }
    else{
      for(int i = 0; i < nloc; i++){
          result(_,i) = exp(plstabmix(lz(_,i), alpha, theta, lK(i,_), true));
      }
    }
    return result;
}



// [[Rcpp::export]]
double find_upper(double p, double alpha, double theta, 
                  NumericVector lK){
  NumericVector result(1);
  result(0) = 1;
  double lpcur = Rcpp::as<double>(plstabmix(result, alpha, theta, lK, true));
  double lp = log(p);
  while(lpcur < lp){
    result(0) = result(0) + 5;
    lpcur = Rcpp::as<double>(plstabmix(result, alpha, theta, lK, true));
  }
  return Rcpp::as<double>(result);
}

// [[Rcpp::export]]
double find_lower(double p, double alpha, double theta, 
                  NumericVector lK){
  NumericVector result(1);
  result(0) = 1;
  double lpcur = Rcpp::as<double>(plstabmix(result, alpha, theta, lK, true));
  double lp = log(p);
  while(lpcur > lp){
    result(0) = result(0) - 5;
    lpcur = Rcpp::as<double>(plstabmix(result, alpha, theta, lK, true));
  }
  return Rcpp::as<double>(result);
}



// [[Rcpp::export]]
NumericVector dlstabmix(NumericVector lz, double alpha, double theta, 
                  NumericVector lK, bool logp = true){
    int n = lz.size();
    NumericVector result(n);
    if(theta==0){
      for(int i = 0; i < n; i++){
        result[i] = -(exp(-lz[i]) + lz[i]);
      }
    }
    else{
      for(int i = 0; i < n; i++){
        result[i] =  log(sum(pow(theta + exp((lK - lz[i])/alpha), alpha-1)* 
                           exp((lK - lz[i])/alpha)));
      }
  	  result += plstabmix(lz, alpha, theta, lK, true);
    }
    if (logp) {
      return result;
    }
    else{
      return exp(result);
    }
}


// Returns a matrix the same size as lz
// columns of lz are different spatial locations
// rows of lK are different spatial locations
// [[Rcpp::export]]
NumericMatrix dlstabmixM(NumericMatrix lz, double alpha, double theta, 
                  NumericMatrix lK, bool logp = true){
    int nloc = lz.ncol();
    int nyear = lz.nrow();
    NumericMatrix result(nyear, nloc);
    if (logp) {
        for(int i = 0; i < nloc; i++){
          result(_,i) = dlstabmix(lz(_,i), alpha, theta, lK(i,_), true);
      }
    }
    else{
      for(int i = 0; i < nloc; i++){
          result(_,i) = exp(dlstabmix(lz(_,i), alpha, theta, lK(i,_), true));
      }
    }
    return result;
}

// [[Rcpp::export]]
NumericMatrix pevdC(NumericMatrix q, NumericMatrix loc, NumericMatrix scale, 
				NumericMatrix shape, bool lower_tail = true, bool log_p = false){
	int n = q.nrow();
	int m = q.ncol();
	NumericMatrix result(n,m);
	if(!log_p){
		if(lower_tail){
			for(int i = 0; i < n; i++){
				for(int j = 0; j < m; j++){
					if(shape(i,j)!=0){
						result(i,j) = exp(-pow(std::max(1+shape(i,j)*q(i,j), 0.0),(-1/shape(i,j))));
					}
					else{
						result(i,j) = exp(-exp(-q(i,j)));
					}
				}
			}
		}
		else{
			for(int i = 0; i < n; i++){
				for(int j = 0; j < m; j++){
					if(shape(i,j)!=0){
						result(i,j) = 1 - exp(-pow(std::max(1+shape(i,j)*q(i,j), 0.0),(-1/shape(i,j))));
					}
					else{
						result(i,j) = 1 - exp(-exp(-q(i,j)));
					}
				}
			}
		}
	}
	else{
		if(lower_tail){
			for(int i = 0; i < n; i++){
				for(int j = 0; j < m; j++){
					if(shape(i,j)!=0){
						result(i,j) = -pow(std::max(1 + shape(i,j)*q(i,j), 0.0),(-1/shape(i,j)));
					}
					else{
						result(i,j) = -exp(-(q(i,j)));
					}
				}
			}
		}
		else{
			for(int i = 0; i < n; i++){
				for(int j = 0; j < m; j++){
					if(shape(i,j)!=0){
						result(i,j) = log(1-exp(-pow(std::max(1 + shape(i,j)*q(i,j), 0.0),(-1/shape(i,j)))));
					}
					else{
						result(i,j) = log(1-exp(-exp(-(q(i,j)))));
					}
				}
			}
		}
	}

	return(result);
}

// [[Rcpp::export]]
NumericMatrix qevdC(NumericMatrix p, NumericMatrix loc, NumericMatrix scale, NumericMatrix shape){
	int n = p.nrow();
	int m = p.ncol();
	NumericMatrix result(n,m);
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			if(shape(i,j)!=0){
				result(i,j) = loc(i,j) + scale(i,j) * (pow(-log(p(i,j)),-shape(i,j)) - 1)/shape(i,j);
			}
			else{
				result(i,j) =  loc(i,j) - scale(i,j) * log(-log(p(i,j)));
			}
		}
	}
	return(result);
}




// [[Rcpp::export]]
std::vector<size_t> rankC(const std::vector<float>& v_temp) {
    std::vector<std::pair<float, size_t> > v_sort(v_temp.size());

    for (size_t i = 0U; i < v_sort.size(); ++i) {
        v_sort[i] = std::make_pair(v_temp[i], i);
    }

    sort(v_sort.begin(), v_sort.end());

    std::pair<double, size_t> rank;
    std::vector<size_t> result(v_temp.size());

    for (size_t i = 0U; i < v_sort.size(); ++i) {
        if (v_sort[i].first != rank.first) {
            rank = std::make_pair(v_sort[i].first, i);
        }
        result[v_sort[i].second] = rank.second;
    }
    return result;
}

// [[Rcpp::export]]
double calc_jpexcd_one_pairC(NumericVector x, NumericVector y, double p){
  int nx = x.size();
  NumericVector r1(nx) ;
  NumericVector r2(nx) ;
  NumericVector l1(nx), l2(nx);
  LogicalVector lall(nx);
  NumericVector  result(nx) ;

  std::vector<float> x2 = as< std::vector<float> >(x);
  std::vector<float> y2 = as< std::vector<float> >(y);

  r1 =  wrap(rankC(x2));
  r2 =  wrap(rankC(y2));
  l1 = (r1 + 1)/nx;
  l2 = (r2 + 1)/nx;
  lall = (l1 > p)&(l2>p);
  return sum(as<NumericVector>(lall))/nx;
}

// [[Rcpp::export]]
LogicalVector isNA(NumericVector x) {
  int n = x.size();
  LogicalVector out(n);
  
  for (int i = 0; i < n; ++i) {
    out[i] = NumericVector::is_na(x[i]);
  }
  return out;
}






















