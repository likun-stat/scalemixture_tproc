//
//  main.cpp
//  integration
//
//  Created by ZhangLikun on 9/17/18.
//  Copyright © 2018 ZhangLikun. All rights reserved.
//



#include <RcppGSL.h>
#include <math.h>
#include <gsl/gsl_integration.h>

// [[Rcpp::depends(RcppGSL)]]

using namespace Rcpp;


struct my_f_params {double x; double tau_sqd; double gamma; double sigma};

// xval > 1
double mix_me_distn_integrand(double x, void *p) {
    struct my_f_params *params = (struct my_f_params *)p;
    double xval   = (params->x);
    double tau_sqd  = (params->tau_sqd);
    double gamma = (params->gamma);
    double sigma = (params->sigma);
    
    double half_result = R::pt((xval-x)/sigma, gamma);
    
    return R::dnorm(x, 0.0, sqrt(tau_sqd), 0) * half_result;
}


double mix_me_dens_integrand(double x, void *p) {
    struct my_f_params *params = (struct my_f_params *)p;
    double xval   = (params->x);
    double tau_sqd  = (params->tau_sqd);
    double gamma = (params->gamma);
    double sigma = (params->sigma);
    
    double half_result = R::dt((xval-x)/sigma, gamma);
    
    return R::dnorm(x, 0.0, sqrt(tau_sqd), 0) * half_result;
}


// [[Rcpp::export]]
double pmixture_me_uni(double x, double tau_sqd, double gamma, double sigma, double relerr = 1e-10) {
    
    double result = 0.0;
    
    gsl_function F;
    F.function = &mix_me_distn_integrand;
    
    gsl_integration_workspace *work = gsl_integration_workspace_alloc(1e5);
    
    gsl_set_error_handler_off();
    
    struct my_f_params params = { x, tau_sqd, gamma, sigma };
    F.params = &params;
    
    double abserr = 0.0;
    
    // QAGI adaptive integration on infinite intervals
    double err = gsl_integration_qagi(&F, 1e-12, relerr, 1e5, work, &result, &abserr);
    
    if (ISNAN(err)){
        Rcpp::Rcout << "Error in integration. Returning -1" << std::endl;
        Rcpp::Rcout << "Err = " << err << std::endl;
        result = -1.0;
    }
    
    gsl_integration_workspace_free(work);
    
    return result;
}



// [[Rcpp::export]]
NumericVector pmixture_me(NumericVector x, double tau_sqd, double gamma, double sigma, double relerr = 1e-10) {
    
    int n = x.size();
    NumericVector resultVec(n);
    
    gsl_function F;
    F.function = &mix_me_distn_integrand;
    
    gsl_integration_workspace *work = gsl_integration_workspace_alloc(1e5);
    
    gsl_set_error_handler_off();
    
    
    
    
    // Calculate CDF fucntion at each x values
    for(int i = 0; i < n; i++) {
        struct my_f_params params = { x[i], tau_sqd, gamma, sigma };
        F.params = &params;
        
        double result = 0.0;
        double abserr = 0.0;
        
        
        
        // QAGI adaptive integration on infinite intervals
        double err = gsl_integration_qagi(&F, 1e-12, relerr, 1e5, work, &result, &abserr);
        
        if (ISNAN(err)){
            Rcpp::Rcout << "Error in integration. Returning -1" << std::endl;
            Rcpp::Rcout << "Err = " << err << std::endl;
            result = -1.0;
        }
        
        resultVec[i] = result;
    }
    
    gsl_integration_workspace_free(work);
    
    return resultVec;
}


// [[Rcpp::export]]
NumericVector dmixture_me(NumericVector x, double tau_sqd, double gamma, double sigma, double relerr = 1e-6) {
    
    int n = x.size();
    NumericVector resultVec(n);
    
    gsl_function F;
    F.function = &mix_me_dens_integrand;
    
    gsl_integration_workspace *work = gsl_integration_workspace_alloc(1e5);
    
    gsl_set_error_handler_off();
    
    
    // Calculate PDF fucntion at each x values
    for(int i = 0; i < n; i++) {
        struct my_f_params params = { x[i], tau_sqd, gamma, sigma };
        F.params = &params;
        
        double result = 0.0;
        double abserr = 0.0;
        
        
        // QAGI adaptive integration on infinite intervals
        double err = gsl_integration_qagi(&F, 1e-12, relerr, 1e5, work, &result, &abserr);
        
        if (!ISNAN(err)){
            result =  result;
        }
        else {
            Rcpp::Rcout << "Error in integration. Returning -1" << std::endl;
            Rcpp::Rcout << "Err = " << err << std::endl;
            result = -1.0;
        }
        
        
        if (result > 0) {
            resultVec[i] = result;
        } else {
            resultVec[i] = 0;
        }
    }
    
    gsl_integration_workspace_free(work);
    
    return resultVec;
}



/*  Test in R
 # parameter settings
 delta <- 0.3
 tau <- 1   # Always variance = std^2
 
 x_vals <- seq(1.001, 20, length = 10000)
 system.time(pdf_vals<-dmixture_me(x_vals, tau, delta))
 system.time(cdf_vals<-pmixture_me(x_vals, tau, delta))
 
 par(mfrow=c(2,1))
 plot(x_vals, pdf_vals, type="l")
 plot(x_vals, cdf_vals, type="l")
 */






// -------------------------------------------------------------------------- //
// This only makes sense if we want to search over x > 0.  Since we are interested
// in extremes, this seems okay -- things are not so extreme if they are in the
// lower half of the support of the (copula) distribution.
//                                                                            //

// [[Rcpp::export]]
NumericVector find_xrange_pmixture_me(double min_p, double max_p,
                                      NumericVector x_init,
                                      double tau_sqd, double gamma, double sigma, double relerr = 1e-10) {
    NumericVector min_x(1);
    NumericVector max_x(1);
    double p_min_x;
    double p_max_x;
    NumericVector x_range(2);
    
    min_x[0] = x_init[0];
    max_x[0] = x_init[1];
    
    if((min_x[0] <= 0) || (min_p <= 0.5)) {
        Rcpp::stop("This will only work for x > 0, which corresponds to p > 1/2.");
    }
    if(min_x[0] >= max_x[0]) {
        Rcpp::stop("x_init[1] must be smaller than x_init[2].");
    }
    
    
    // First the min
    p_min_x = pmixture_me(min_x, tau_sqd, gamma, sigma, relerr)[0];
    while (p_min_x > min_p) {
        // Rcpp::Rcout << "left x is " << min_x[0] << std::endl;
        // Rcpp::Rcout << "F(" << min_x[0] << ") = " << p_min_x << std::endl;
        min_x[0] = min_x[0]/2;
        p_min_x = pmixture_me(min_x, tau_sqd, gamma, sigma, relerr)[0];
    }
    x_range[0] = min_x[0];
    
    // Now the max
    p_max_x = pmixture_me(max_x, tau_sqd, gamma, sigma, relerr)[0];
    while (p_max_x < max_p) {
        // Rcpp::Rcout << "right x is " << max_x << std::endl;
        // Rcpp::Rcout << "F(" << max_x << ") = " << p_max_x<< std::endl;
        max_x[0] = max_x[0]*2;
        p_max_x = pmixture_me(max_x, tau_sqd, gamma, sigma, relerr)[0];
    }
    x_range[1] = max_x[0];
    
    return x_range;
}
//                                                                            //
// -------------------------------------------------------------------------- //



// -------------------------------------------------------------------------- //
// Update X.s one by one using full conditional without the information of X
// [[Rcpp::export]]
double d_gpd_uni(double x, double loc, double scale, double shape, bool islog = true){
    double tm1 = 1/scale, tml1 = -log(scale);
    double tm2 = shape/scale;
    double tm3 = -1/shape-1;
    double res;
    
    if(abs(shape)>1e-09){
        if(!islog){
            res = tm1*pow(1+tm2*(x-loc), tm3);
        }
        else{
            res = tml1+tm3*log(1+tm2*(x-loc));
        }
    } else{
        if(!islog){
            res = exp((loc-x)/scale);
        }
        else{
            res = (loc-x)/scale;
        }
    }
    
    return res;
}

// [[Rcpp::export]]
NumericVector d_gpd(NumericVector x, double loc, NumericVector scale, double shape, bool islog = true){
    int n = x.size();
    NumericVector resultVec(n);
    NumericVector tm1 = 1/scale, tml1 = -log(scale);
    NumericVector tm2 = shape/scale;
    double tm3 = -1/shape-1;
    
    if(abs(shape)>1e-09){
        if(!islog){
            for(int i=0; i<n; i++){
                resultVec[i] = tm1[i]*pow(1+tm2[i]*(x[i]-loc), tm3);
            }
        }
        else{
            for(int i=0; i<n; i++){
                resultVec[i] = tml1[i]+tm3*log(1+tm2[i]*(x[i]-loc));
            }
        }
    } else{
        if(!islog){
            for(int i=0; i<n; i++){
                resultVec[i] = exp((loc-x[i])/scale[i]);
            }
        }
        else{
            for(int i=0; i<n; i++){
                resultVec[i] = (loc-x[i])/scale[i];
            }
        }
    }
    
    return resultVec;
}


// NOTE: thresh_X is required
// [[Rcpp::export]]
double marg_transform_data_mixture_me_likelihood(NumericVector Y, NumericVector X, NumericVector X_s,
                                                 LogicalVector cen, double prob_below,
                                                 NumericVector theta_gpd, double shape, NumericVector lon, double gamma, double sigma,
                                                 double tau_sqd, double thresh_X){     //time consuming
    NumericVector ll(Y.size());
    int size_Y = Y.size();
    double loc = theta_gpd[0], a = theta_gpd[1], b = theta_gpd[2];
    NumericVector Scale = a+b*lon;
    NumericVector tmp(1), res(1);
    double temp;
    
    for(int i=0; i<size_Y; i++){
        if(cen[i] == TRUE){
            ll[i] = R::pnorm(thresh_X, X_s[i], sqrt(tau_sqd), 1, 1);
        }
        else{
            tmp[0] = X[i]; temp = Y[i];
            res = dmixture_me(tmp, tau_sqd, gamma, sigma);
            ll(i) = R::dnorm(X[i], X_s[i], sqrt(tau_sqd), 1)+ d_gpd_uni(temp, loc, Scale[i], shape, true) -log(res[0]);
            if(!std::isfinite(ll(i))) ll(i) = -std::numeric_limits<double>::infinity();
        }
    }
    
    return sum(ll);
}

// NOTE: thresh_X is required
// [[Rcpp::export]]
double marg_transform_data_mixture_me_likelihood_uni(double Y, double X, double X_s,
                                                 bool cen, double prob_below,
                                                 NumericVector theta_gpd, double shape, double lon, double gamma, double sigma,
                                                 double tau_sqd, double thresh_X, bool missing){     //time consuming
    
    double loc = theta_gpd[0], a = theta_gpd[1], b = theta_gpd[2];
    double scale = a + b*lon;
    NumericVector tmp(1), res(1);
    double ll;
    
    if(cen == TRUE){
        ll = R::pnorm(thresh_X, X_s, sqrt(tau_sqd), 1, 1);
    }
    else{
        if(missing == FALSE){
            tmp[0] = X;
            res = dmixture_me(tmp, tau_sqd, gamma, sigma);
            ll = R::dnorm(X, X_s, sqrt(tau_sqd), 1)+ d_gpd_uni(Y, loc, scale, shape, true) -log(res[0]);
        }
        else{
            ll = 1 - R::pnorm(thresh_X, X_s, sqrt(tau_sqd), 1, 1);
        }
        if(!std::isfinite(ll)) ll = -std::numeric_limits<double>::infinity();
    }
    
    return ll;
}

// [[Rcpp::export()]]
double eig2inv_quadform_vector_cpp (NumericMatrix V, NumericVector d_inv, NumericVector x) {
    int nrow = V.nrow();
    
    NumericVector VtX(nrow);
    for(int c=0; c<nrow; c++){
        for(int r=0; r<nrow; r++){
            VtX[c] = VtX[c] + V(r,c)*x(r);
        }
    }
    
    double res = 0;
    for(int i =0; i<nrow; i++){
        res += d_inv(i)*VtX(i)*VtX(i);
    }
    
    return res;
}

// [[Rcpp::export()]]
double X_s_likelihood_conditional_cpp(NumericVector X_s, double R, NumericMatrix V, NumericVector d){
    NumericVector tmp = X_s/R;
    
    double part1 = -0.5*eig2inv_quadform_vector_cpp(V, 1/d, tmp);
    
    return part1;
}


// NOTE: thresh_X is required
// [[Rcpp::export()]]
List update_X_s_onetime (NumericVector Y, NumericVector X, NumericVector X_s,
                         LogicalVector cen, double prob_below,
                         NumericVector theta_gpd, double shape, NumericVector lon, double gamma, double sigma,
                         double tau_sqd, double thresh_X,
                         NumericVector v_q,  double R, NumericMatrix V, NumericVector d, LogicalVector missing){
    
    int n_s = X.size();
    NumericVector prop_X_s(n_s);
    NumericVector accept(n_s);
    
    double log_num=0, log_denom=0, r=0, temp;
    for(int iter=0; iter<n_s; iter++){
        // tripped : NumericVector X=Y, changing X will change Y as well.
        std::copy(X_s.begin(), X_s.end(), prop_X_s.begin());
        temp = X_s(iter)+v_q(iter)*R::rnorm(0,1);
        prop_X_s(iter) = temp;
        log_num = marg_transform_data_mixture_me_likelihood_uni(Y(iter), X(iter), prop_X_s(iter), cen(iter), prob_below, theta_gpd, shape, lon(iter), gamma, sigma, tau_sqd, thresh_X, missing(iter)) + X_s_likelihood_conditional_cpp(prop_X_s, R, V, d);
        log_denom = marg_transform_data_mixture_me_likelihood_uni(Y(iter), X(iter), X_s(iter), cen(iter), prob_below, theta_gpd, shape, lon(iter), gamma, sigma, tau_sqd, thresh_X, missing(iter)) + X_s_likelihood_conditional_cpp(X_s, R, V, d);
        
        r = exp(log_num - log_denom);
        if(!std::isfinite(r)) r = 0;
        
        if(R::runif(0,1)<r){
            X_s(iter) = temp;
            accept(iter) = accept(iter) + 1;
        }
    }
    
    List result;
    result["X.s"] = X_s;
    result["accept"] = accept;
    
    return result;
}

// -------------------------------------------------------------------------- //


