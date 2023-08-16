#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

#define PI 3.14159265358979323846

// [[Rcpp::export]]
double calculateMVDensity(NumericVector point, NumericMatrix data, double bw) {
    double sum = 0;
    int k = data.ncol();
    for(int row = 0; row < data.nrow(); row++){
        double exponent = 0;
        for(int col = 0; col < k; col++){
            exponent += -0.5 * std::pow(point(col) - data(row,col),2.0) / bw;
        }
        sum += std::exp(exponent);
    }
    //multiply sum by the constant term 2π^(-k/2)*det(Σ)^-1/2
    //note: as sigma is the diagonal matrix of all bw, det(Σ)^-1/2 = Σ^-k/2 
    //      so we can combine with the 2π term
    sum *= std::pow(2.0*PI*bw,-double(k)/2.0);
    return sum;
}
