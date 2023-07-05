#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>
#include <cerrno>
#include "Distributions.hpp"
#include "qpois.h"
#include <limits>
#include <stdexcept>

namespace stats {

    DistributionType str2DistributionType(std::string str){
        DistributionType type;
        if(str == "Beta"){
            type = Beta;
        } else if(str == "Exponential"){
            type = Exponential;
        } else if(str == "Fixed"){
            type = Fixed;
        } else if(str == "HalfNormal"){
            type = HalfNormal;
        } else if(str == "Normal"){
            type = Normal;
        } else if(str == "LogNormal"){
            type = LogNormal;
        } else if(str == "Laplace"){
            type = Laplace;
        } else if(str == "Poisson"){
            type = Poisson;
        } else if(str == "Skellam"){
            type = Skellam;
        } else if(str == "Uniform"){
            type = Uniform;
        } else if(str == "T"){
            type = T;
        } else if(str == "TruncatedNormal"){
            type = TruncatedNormal;
        } else if(str == "DiscreteTruncatedNormal"){
            type = DiscreteTruncatedNormal;
        } else if(str == "ChiSq"){
            type = ChiSq;
        } else {
            throw std::invalid_argument("Attempt to convert unrecognized string to a DistributionType");
        }
        return type;
    }

    //Chi-Squared
    double ChiSqPDF(double x, double df){
        if(df < 1){
            throw std::invalid_argument("Attempt to get a chi-squared density with non-natural degreees of freedom");
        }
        if(x == 0 && df == 1){
            throw std::invalid_argument("Attempt to get chi-squared density at zero with one degree of freedom");
        } else if(x < 0){
            throw std::invalid_argument("Attempt to get chi-squared density for negative value");
        }
        int k = df;
        double denom = std::pow(2.0,k/2.0) * std::tgamma(k/2.0);
        double d = std::pow(x,k/2.0-1.0) * exp(-x / 2.0);
        d /= denom;
        return d;
    }

    double ChiSqCDF(double x, double df){
        if(df < 1){
            throw std::invalid_argument("Attempt to get a chi-squared cumulative probability with non-natural degreees of freedom");
        }
        if(x <= 0){
            return 0;
        }
        double p = 0;
        int k = df;
        try {
            p = boost::math::gamma_p(k/2.0,x/2.0);
        } catch (std::exception & e) {
            throw;
        }
        return p;
    }

    double ChiSqQuantile(double p, double df){
        if(df < 1){
            throw std::invalid_argument("Attempt to get a chi-squared cumulative probability with non-natural degreees of freedom");
        }
        double q = 0;
        int k = df;
        try {
            q = boost::math::gamma_p_inv(k /2.0, p);
        } catch (std::exception & e) {
            throw;
        }
        q *= 2.0;
        return q;
    }

    //Normal
    double NormalPDF(double x, double mu, double sigma){
        if(sigma < 0){
            throw std::invalid_argument("Attempt to get normal density with negative variance");
        }
        if(sigma == 0){ //No variance
            if(x - mu  == 0){ //Infinite density at mu
                return std::numeric_limits<double>::infinity();
            } //No density elsewhere
            return 0.0;
        }
        double result = x - mu;
        result /= sigma;
        result *= result;
        result /= -2;
        result = std::exp(result);
        result /= sigma * root2 * pi;
        return result;
    }

   long double NormalCDF(double x, double mu, double sigma, bool bLogP){
        if(sigma < 0){
            throw std::invalid_argument("Attempt to get normal cumulative density with negative variance");
        }
       if(sigma == 0){ //No variance
           if(x < mu){ //No density prior to the mean
               return 0.0;
           } //All density after the mean
           return 1.0;
       }
       long double z = (x - mu) / sigma;
       long double p = 0.0;
       if(bLogP){
           long double zz = z / root2;
           p = std::log(0.5)+std::log(boost::math::erf(zz) + (long double)1.0);

       } else {
           p = 0.5*(boost::math::erf(z / root2) + 1.0);
       }
       return p;
    } 

    double NormalQuantile(double p, double mu, double sigma){
        if(sigma < 0){
            throw std::invalid_argument("Attempt to get normal quantile with negative variance");
        }
        double erf_invRes;
        try{
            erf_invRes = boost::math::erf_inv(2.0*p-1.0);
        } catch (std::domain_error & e){
            throw std::domain_error("Attempt to get Normal Quantile of non proportion");
        } catch (std::overflow_error & e){
            throw std::overflow_error("Result of Normal Quantile is infinite");
        }
        return mu + sigma * root2 * erf_invRes;
    }

    double NormalCIWidth(double alpha, double sigma){
        if(sigma < 0){
            throw std::invalid_argument("Attempt to get normal ci width with negative variance");
        }
        double res = 0;
        try{
            res = stats::NormalQuantile(1-(alpha/2),0,sigma) * 2;
        } catch(std::domain_error & e){
            throw std::domain_error("Attempt to get Normal CI width non proportional alpha");
        } catch(std::overflow_error & e){
            throw;
        }
        return res;
    }

    //LogNormal
    double LogNormalPDF(double x, double logmu, double logsigma){
        if(logsigma < 0){
            throw std::invalid_argument("Attempt to get lognormal density with negative variance");
        }
        if(x < 0){
            throw std::domain_error("Attempt to get logNormal density of a negative value");
        }
        return NormalPDF(std::log(x), logmu, logsigma)/x;
    }

    double LogNormalCDF(double x, double logmu, double logsigma){
        if(logsigma < 0){
            throw std::invalid_argument("Attempt to get lognormal cumulative density with negative variance");
        }
        if(x < 0){
            throw std::domain_error("Attempt to get cumulative logNormal density of a negative value");
        }
        return NormalCDF(std::log(x),logmu,logsigma);
    }

    double LogNormalQuantile(double p, double logmu, double logsigma){
        if(logsigma < 0){
            throw std::invalid_argument("Attempt to get lognormal quantile with negative variance");
        }
        double normQ = 0;
        try {
            normQ = NormalQuantile(p,logmu,logsigma);
        } catch (std::domain_error & e){
            throw std::domain_error("Attempt to get logNormal Quantile of non proportion");
        } catch (std::overflow_error & e){
            if(p < 0.5){
                return 0;
            }
            throw std::overflow_error("Result of LogNormal Quantile is infinite");
        }
        return std::exp(normQ);
    }

    double HalfNormalPDF(double x, double sigma){
        if(sigma < 0){
            throw std::invalid_argument("Attempt to get half normal density with negative variance");
        }
        if(x < 0){
            return 0;
        }
        return 2*NormalPDF(x,0,sigma);
    }

    double HalfNormalCDF(double x, double sigma){
        if(sigma < 0){
            throw std::invalid_argument("Attempt to get cumulativ half normal density with negative variance");
        }
        if(x <= 0){
            return 0;
        }
        double z = x / sigma;
        return boost::math::erf(z / root2);
    }

    double HalfNormalQuantile(double p, double sigma){
        if(sigma < 0){
            throw std::invalid_argument("Attempt to get half normal quantile with negative variance");
        }
        double q = 0;
        double erf_invRes;
        try{
            erf_invRes = boost::math::erf_inv(p);
        } catch (std::domain_error & e){
            throw std::domain_error("Attempt to get Normal Quantile of non proportion");
        } catch (std::overflow_error & e){
            throw std::overflow_error("Result of Normal Quantile is infinite");
        }
        q = erf_invRes * sigma * root2;
        return q;
    }


    //TruncatedNormal
    double TruncatedNormalPDF(double x, double mu, double sigma, double a, double b){
        if(a > b){
            throw std::invalid_argument("Truncated Normal Distribtion bounds are inverted");
        }
        if(x < a || x > b){
            return 0.0;
        }
        if(sigma < 0){
            throw std::invalid_argument("Attempt to get truncnormal density with negative variance");
        }
        double epsilon = (x - mu) / sigma;
        double alpha = (a - mu) / sigma;
        double beta = (b - mu) / sigma;
        double Zeta = NormalCDF(beta) - NormalCDF(alpha);
        double d = NormalPDF(epsilon) / sigma / Zeta;
        return d;
    }

    double TruncatedNormalCDF(double x, double mu, double sigma, double a, double b){
        if(sigma < 0){
            throw std::invalid_argument("Attempt to get truncnormal cumulative density with negative variance");
        }
        double epsilon = (x - mu) / sigma;
        double alpha = (a - mu) / sigma;
        double beta = (b - mu) / sigma;
        double PhiAlpha = NormalCDF(alpha);
        double Zeta = NormalCDF(beta) - PhiAlpha;
        double p = (NormalCDF(epsilon) - PhiAlpha) / Zeta;
        return p;
    }

    double TruncatedNormalQuantile(double p, double mu, double sigma, double a, double b){
        if(sigma < 0){
            throw std::invalid_argument("Attempt to get truncnormal quantile with negative variance");
        }
        if(a > b){
            throw std::invalid_argument("Truncated Normal Distribtion bounds are inverted");
        }
        double alpha = (a - mu) / sigma;
        double beta = (b - mu) / sigma;
        double PhiAlpha = NormalCDF(alpha);
        double Zeta = NormalCDF(beta) - PhiAlpha;
        double normQ = 0;
        try{
            normQ = NormalQuantile(p * Zeta + PhiAlpha);
        } catch (std::domain_error & e){
            throw std::domain_error("Attempt to get Truncated Normal Quantile of non proportion");
        } catch (std::overflow_error & e){
            if(p > 0.5){
                return b;
            }
            return a;
        }
        double x = mu + sigma * normQ;
        //fprintf(stderr,"α:%0.05f β:%0.05f Φ:%0.05f Ζ:%0.05f x:%0.05f\n",alpha,beta,PhiAlpha,Zeta,x);
        return x;
    }

    //Discrete Truncated Normal
    long double DiscreteTruncatedNormalPMF(double x, double mu, double sigma, double a, double b, bool bLogP){
        if(a > b){
            throw std::invalid_argument("Discrete Truncated Normal Distribtion bounds are inverted");
        }
        if(x < a || x > b){
            return 0.0;
        }
        if(sigma < 0){
            throw std::invalid_argument("Attempt to get discretetruncnormal density with negative variance");
        }
        if(sigma == 0){ //no variance but have to discretize
            return (int(std::round(x)) == int(std::round(mu)));
        }
        long double epsilon1 = (x + 0.5 - mu) / sigma;
        long double epsilon2 = (x - 0.5 - mu) / sigma;
        long double alpha = (a - mu) / sigma;
        long double beta = (b - mu) / sigma;
        long double Zeta = NormalCDF(beta,0.0,1.0,bLogP) - NormalCDF(alpha,0.0,1.0,bLogP);
        long double d1 = NormalCDF(epsilon1,0.0,1.0,bLogP) / Zeta;
        long double d2 = NormalCDF(epsilon2,0.0,1.0,bLogP) / Zeta;
        long double d = d1-d2;
        if(bLogP){
            d = d1 + std::log(1-std::exp(d2-d1));
        }
        return d;
    }

    double DiscreteTruncatedNormalCDF(double x, double mu, double sigma, double a, double b){
        return TruncatedNormalCDF(x,mu,sigma,a,b);
    }

    double DiscreteTruncatedNormalQuantile(double p, double mu, double sigma, double a, double b){
        return std::round(TruncatedNormalQuantile(p,mu,sigma,a,b));
    }

    //Fixed
    double FixedPMF(double x, double mu){
        return ((x == mu) ? 1.0 : 0.0);
    }
    double FixedCDF(double x, double mu){
        return ((x >= mu) ? 1.0 : 0.0);
    }
    double FixedQuantile(double p, double mu){
        return mu;
    }

    //Laplace
    double LaplacePDF(double x, double mu, double b){
        if(b < 0){
            throw std::invalid_argument("Attempt to get laplace density with negative variance");
        }
        double result = std::abs(x - mu);
        if(b == 0){ //Zero variance
            if(result == 0){ //Density at mu with zero variance
                return std::numeric_limits<double>::infinity();
            } //Density not at mu with zero variance
            return 0.0;
        }
        result /= -b;
        result = std::exp(result);
        result /= 2 * b;
        return result;
    }

    double LaplaceCDF(double x, double mu, double b){
        if(b < 0){
            throw std::invalid_argument("Attempt to get laplace cumulative density with negative variance");
        }
        if(b == 0){ //No variance
            if(x < mu){ //No cumulative density before mu
                return 0.0;
            } // All the cumulative density after mu
            return 1.0;
        }
        double z = (x-mu)/b;
        if(x < mu){
            return 0.5 * std::exp(z);
        }
        return 1.0 - 0.5 * std::exp(-z);
    }

    double LaplaceQuantile(double p, double mu, double b){
        if(b < 0){
            throw std::invalid_argument("Attempt to get laplace quantile with negative variance");
        }
        if(p < 0 || p > 1){
            throw std::domain_error("Attempt to get laplace quantile of non proportion");
        }
        if(p == 0 || p == 0){
            throw std::overflow_error("Result of Laplace Quantile is infinite");
        }
        if(p < 0.5){
            return mu + b* std::log(2*p);
        }
        return mu -b*std::log(2 - 2*p);
    }

    double LaplaceCIWidth(double alpha, double sigma){
        if(sigma < 0){
            throw std::invalid_argument("Attempt to get laplace ci width with negative variance");
        }
        if(alpha < 0 || alpha > 1){
            throw std::domain_error("Attempt to get Laplace CI width for non proportion");
        }
        double ci = 0;
        try{
            ci = stats::LaplaceQuantile(1-(alpha/2),0,sigma) * 2;
        } catch (std::overflow_error & e){
            throw;
        }
        return ci;
    }

    //Poisson
    double PoissonPMF(double k, double lambda){
        if(k < 0){
            throw std::invalid_argument("Attempt to get poisson density of negative value");
        }
        if(lambda < 0){
            throw std::invalid_argument("Attempt to get poisson density with negative value");
        }
        if(lambda == 0){
            return (k == 0) ? 1.0 : 0.0;
        }
        double gammaRes = 0;
        try{
            gammaRes = boost::math::tgamma(k + 1);
        } catch (std::overflow_error & e){
            throw std::overflow_error("PoissonPMF result unrepresentable; large gamma");
        } catch (std::underflow_error & e){
            throw std::overflow_error("PoissonPMF result unrepresentable; small gamma");
        } catch (boost::math::evaluation_error & e){
            throw std::logic_error("PoissonPMF result took to long to evaluate gamma function");
        }
        return pow(lambda,k) * exp(-lambda) / gammaRes;
    }

    double PoissonCDF(double k, double lambda){
        if(k < 0){
            throw std::invalid_argument("Attempt to get poisson density of negative value");
        }
        if(lambda < 0){
            throw std::invalid_argument("Attempt to get poisson density with negative lambda");
        }
        if(lambda == 0.0){
            return 1.0;
        }
        if(std::isinf(lambda)){ //With an infinite rate, an infinite number of events will occur
            if(std::isinf(k)){ //Infinite density at +infinity
                return std::numeric_limits<double>::infinity();
            } //Zero density elsewhere
            return 0;
        }
        double p = 0;
        for(double j = 0; j < k; j += 1.0){
            double val = j * std::log(lambda);
            val -= std::lgamma(std::floor(j) + 1);
            p += std::exp(val);
        }
        p *= std::exp(-lambda);
        return p;
    }

    double PoissonQuantile(double p, double lambda){
        if(p < 0 || p > 1){
            throw std::domain_error("Attempt to get poisson quantile of non proportion");
        }
        if(lambda == 0){
            return 0;
        }
        if(lambda < QPois_PreCalc_MaxLambda && p < QPois_PreCalc_MaxP){
            double L = lambda / double(QPois_PreCalc_LambdaIncrement);
            double P = p / double(QPois_PreCalc_PIncrement);
            int L1 = std::floor(L);
            int L2 = std::ceil(L);
            int P1 = std::floor(P);
            int P2 = std::ceil(P);
            double q = 0;
            q += (L2 - L) * (P2 - P) * QPoisson[L1][P1];
            q += (L2 - L) * (P - P1) * QPoisson[L1][P2];
            q += (L - L1) * (P2 - P) * QPoisson[L2][P1];
            q += (L - L1) * (P - P1) * QPoisson[L2][P2];
            q = (q <= lambda) ? std::floor(q) : std::ceil(q);
            return q;
        }
        errno = 0;
        double invQresult = boost::math::gamma_q_inva(lambda,p);
        return std::ceil(invQresult) - 1;
    }

    //Exponential
    double ExponentialPDF(double x, double lambda){
        if(x < 0){
            throw std::domain_error("Attempt to get exponential density of a negative value");
        }
        if(lambda <= 0){
            throw std::invalid_argument("Attempt to get exponential quantile with non-positive lambda");
        }
        if(std::isinf(lambda)){ //With an infinite rate, wait times are always zero
            if(x == 0){ // Infinite density at 0
                return std::numeric_limits<double>::infinity();
            } // No density
            return 0;
        }
        return lambda * std::exp(-x * lambda);
    }

    double ExponentialCDF(double x, double lambda){
        if(x < 0){
            throw std::domain_error("Attempt to get cumulative exponential density of a negative value");
        }
        if(lambda <= 0){
            throw std::invalid_argument("Attempt to get exponential quantile with non-positive lambda");
        }
        return 1.0 - std::exp(-lambda * x);
    }

    double ExponentialQuantile(double p, double lambda){
        if(p < 0 || p > 1){
            throw std::domain_error("Attempt to get quantile of non proportion");
        }
        if(lambda <= 0){
            throw std::invalid_argument("Attempt to get exponential quantile with non-positive lambda");
        }
        if(std::isinf(lambda)){
            return std::numeric_limits<double>::infinity();
        }
        return -1 * std::log(1 - p) / lambda;
    }

    //Beta
    double BetaPDF(double x, double alpha, double beta){
        double betaRes = 0;
        try{
            betaRes = boost::math::beta(alpha,beta);
        } catch (std::domain_error & e){
            throw std::domain_error(std::string("In Beta PDF:") + e.what());
        }
        return pow(x, alpha - 1.0) * pow(x, beta - 1.0) / betaRes;
    }

    double BetaCDF(double x, double alpha, double beta){
        double ibetaRes = 0;
        try{
            ibetaRes = boost::math::ibeta(alpha,beta,x);
        } catch (std::domain_error & e){
            throw std::domain_error(std::string("In Beta CDF:") + e.what());
        }
        return ibetaRes;
    }

    double BetaQuantile(double p, double alpha, double beta){
        double ibetaInvRes = 0;
        try{
            ibetaInvRes = boost::math::ibeta_inv(alpha,beta,p);
        } catch (std::domain_error & e){
            throw std::domain_error(std::string("In Beta Quantile:") + e.what());
        }
        return ibetaInvRes;
    }

    double BetaMaxMonotonicVariance(double x){
        if(x < 0.0 || x > 1.0){
            throw std::domain_error("Attempt to get beta max monotonic variance of non proprotion");
        }
        if(x > 0.5){
            x = 1.0 - x;
        }
        return x * (1.0 - x) * (1.0 - x) / (2.0 - x);
    }

    //Skellam
    double SkellamPMF(double k, double lambda, double kappa){
        //Handle cases where it degenerates to poisson
        if(lambda * kappa == 0){ //Either is zero
            int sign = (kappa == 0) ? 1 : -1;
            if(k * sign < 0){ //If sign of the poisson dist doesn't match that of k, k is impossible
                return 0.0;
            }
            return PoissonPMF(std::abs(k),lambda+kappa);
        }
        k = std::round(k);
        double besselArg = 2*std::sqrt(lambda*kappa);
        double p = std::exp(-(lambda+kappa)) * std::pow(lambda/kappa,k/2.0);
        //For integers, I_k(x) = I_|k|(x)
        p *= boost::math::cyl_bessel_i(std::abs(k),besselArg);
        return p;
    }
    double SkellamCDF(double k, double lambda, double kappa){
        throw std::logic_error("No Closed Form of SkellamCDF");
    }

    double SkellamQuantile(double p, double lambda, double kappa){
        throw std::logic_error("No Closed Form of SkellamQuantile");
    }


    double TPDF(double x, double v){
        if(v < 1){
            throw std::invalid_argument("Attempt to get t pdf with non-positive degrees of freedom");
        }
        double gammaRes1 = 0.0;
        double gammaRes2 = 0.0;
        try{
            gammaRes1 = boost::math::tgamma((v + 1.0)/2.0);
            gammaRes2 = boost::math::tgamma(v/2.0);
        } catch (std::overflow_error & e){
            throw std::overflow_error("TPDF result unrepresentable; large gamma");
        } catch (std::underflow_error & e){
            throw std::overflow_error("TPDF result unrepresentable; small gamma");
        } catch (boost::math::evaluation_error & e){
            throw std::logic_error("TPDF result took to long to evaluate gamma function");
        }
        double coef = gammaRes1 / gammaRes2 / std::sqrt(v*pi);
        double base = 1.0 + std::pow(x,2.0) / v;
        double exponent = -(v+1.0)/2.0;
        return coef * std::pow(base,exponent);
    }

    double TCDF(double x, double v){
        if(v < 1){
            throw std::invalid_argument("Attempt to get t pdf with non-positive degrees of freedom");
        }
        double h = boost::math::hypergeometric_pFq({0.5, (v + 1.0)/2.0},{1.5},-std::pow(x,2.0)/v);
        double gammaRes1 = 0.0;
        double gammaRes2 = 0.0;
        try{
            gammaRes1 = boost::math::tgamma((v + 1.0)/2.0);
            gammaRes2 = boost::math::tgamma(v/2.0);
        } catch (std::overflow_error & e){
            throw std::overflow_error("TCDF result unrepresentable; large gamma");
        } catch (std::underflow_error & e){
            throw std::overflow_error("TCDF result unrepresentable; small gamma");
        } catch (boost::math::evaluation_error & e){
            throw std::logic_error("TCDF result took to long to evaluate gamma function");
        }
        return 0.5 + x * gammaRes1 * h / gammaRes2 / std::sqrt(v*pi);
    }

    double TQuantile(double p, double v){
        throw std::logic_error("T Quantile not implemented");
    }

    //Uniform
    double UniformPDF(double a, double b){
        if(a > b){
            throw std::invalid_argument("Attempt to get uniform pdf with inverted bounds");
        }
        if(std::isinf(a) || std::isinf(b)){
            throw std::invalid_argument("Attempt to get uniform pdf with infinite bounds");
        }
        return 1/(b-a);
    }

    double UniformCDF(double x, double a, double b){
        if(a > b){
            throw std::invalid_argument("Attempt to get uniform cdf with inverted bounds");
        }
        if(std::isinf(a) || std::isinf(b)){
            throw std::invalid_argument("Attempt to get uniform cdf with infinite bounds");
        }
        if(x < a || x > b){
            return 0;
        }
        return (x - a) / (b - a);
    }

    double UniformQuantile(double p, double a, double b){
        if(a > b){
            throw std::invalid_argument("Attempt to get uniform quantile with inverted bounds");
        }
        if(std::isinf(a) || std::isinf(b)){
            throw std::invalid_argument("Attempt to get uniform quantile with infinite bounds");
        }
        return p * (b - a) + a;
    }

    //Generic
    
    SDomain GetDomain(DistributionType distribution, ParamMap parameters){
        SDomain domain;
        domain.min = -std::numeric_limits<double>::infinity();
        domain.max = std::numeric_limits<double>::infinity();
        try{
            switch(distribution){
                case(Beta):
                    domain.min = 0;
                    domain.max = 1;
                    break;
                case(Exponential):
                case(LogNormal):
                case(Poisson):
                case(HalfNormal):
                    domain.min = 0;
                    break;
                case(Fixed):
                    domain.min = parameters.at("mu");
                    domain.max = parameters.at("mu");
                    break;
                case(TruncatedNormal):
                case(DiscreteTruncatedNormal):
                case(Uniform):
                    domain.min = parameters.at("a");
                    domain.max = parameters.at("b");
                    break;
                case(Normal):
                case(Laplace):
                case(Skellam):
                case(T):
                    break;
            }
        } catch (std::out_of_range & e){
            throw std::invalid_argument("Attempt to GetDomain with missing parameters");
        }
        return domain;
    }

    double GetPDF(double X, DistributionType distribution, ParamMap parameters){
        double d = 0;
        try{
            switch(distribution){
                case(Beta):
                    d = BetaPDF(X,parameters.at("alpha"),parameters.at("beta"));
                    break;
                case(Exponential):
                    d = ExponentialPDF(X,parameters.at("lambda"));
                    break;
                case(Fixed):
                    d = FixedPMF(X,parameters.at("mu"));
                    break;
                case(HalfNormal):
                    d = HalfNormalPDF(X,parameters.at("sigma"));
                    break;
                case(Normal):
                    d = NormalPDF(X,parameters.at("mu"),parameters.at("sigma"));
                    break;
                case(LogNormal):
                    d = LogNormalPDF(X,parameters.at("logmu"),parameters.at("logsigma"));
                    break;
                case(TruncatedNormal):
                    d = TruncatedNormalPDF(X,parameters.at("mu"),parameters.at("sigma"),
                                             parameters.at("a"),parameters.at("b"));
                    break;
                case(DiscreteTruncatedNormal):
                    d = DiscreteTruncatedNormalPMF(X,parameters.at("mu"),parameters.at("sigma"),
                                             parameters.at("a"),parameters.at("b"));
                    break;
                case(Laplace):
                    d = LaplacePDF(X,parameters.at("mu"),parameters.at("b"));
                    break;
                case(Poisson):
                    d = PoissonPMF(X,parameters.at("lambda"));
                    break;
                case(Skellam):
                    d = SkellamPMF(X,parameters.at("lambda"),parameters.at("kappa"));
                    break;
                case(T):
                    d = TPDF(X,parameters.at("v"));
                    break;
                case(Uniform):
                    d = UniformPDF(parameters.at("a"),parameters.at("b"));
                    break;
            }
        } catch (std::out_of_range & e){
            throw std::invalid_argument("Attempt to GetPDF with missing parameters");
        } catch (std::domain_error & e){ // LogNormal, Exponential, Beta
            throw;
        } catch (std::invalid_argument & e){ // All dists with scale parameters
            throw;
        } catch (std::overflow_error & e){ //Poisson
            throw;
        } catch (std::underflow_error & e){ //Poisson
            throw;
        } catch (boost::math::evaluation_error & e){ //Poisson or T
            throw;
        }
        return d;
    }

    double GetCDF(double X, DistributionType distribution, ParamMap parameters){
        double q = 0;
        try{
            switch(distribution){
                case(Beta):
                    q = BetaCDF(X,parameters.at("alpha"),parameters.at("beta"));
                    break;
                case(Exponential):
                    q = ExponentialCDF(X,parameters.at("lambda"));
                    break;
                case(Fixed):
                    q = FixedCDF(X,parameters.at("mu"));
                    break;
                case(HalfNormal):
                    q = HalfNormalCDF(X,parameters.at("sigma"));
                    break;
                case(Normal):
                    q = NormalCDF(X,parameters.at("mu"),parameters.at("sigma"));
                    break;
                case(LogNormal):
                    q = LogNormalCDF(X,parameters.at("logmu"),parameters.at("logsigma"));
                    break;
                case(TruncatedNormal):
                    q = TruncatedNormalCDF(X,parameters.at("mu"),parameters.at("sigma"),
                                             parameters.at("a"),parameters.at("b"));
                    break;
                case(DiscreteTruncatedNormal):
                    q = DiscreteTruncatedNormalCDF(X,parameters.at("mu"),parameters.at("sigma"),
                                             parameters.at("a"),parameters.at("b"));
                    break;
                case(Laplace):
                    q = LaplaceCDF(X,parameters.at("mu"),parameters.at("b"));
                    break;
                case(Poisson):
                    q = PoissonCDF(X,parameters.at("lambda"));
                    break;
                case(Skellam):
                    q = SkellamCDF(X,parameters.at("lambda"),parameters.at("kappa"));
                    break;
                case(T):
                    q = TCDF(X,parameters.at("v"));
                    break;
                case(Uniform):
                    q = UniformCDF(X,parameters.at("a"),parameters.at("b"));
                    break;
            }
        } catch (std::out_of_range & e){
            throw std::invalid_argument("Attempt to GetCDF with missing parameters");
        } catch (std::domain_error & e){ //LogNormal, Exponential, Beta
            throw;
        } catch (std::invalid_argument & e){ //Truncated Normal
            throw;
        } catch (std::logic_error & e){ //Skellam
            throw;
        } catch (boost::math::evaluation_error & e){ //T
            throw;
        }
        return q;
    }

    double GetQuantile(double p, DistributionType distribution, ParamMap parameters){
        double x = 0;
        try{
            switch(distribution){
                case(Beta):
                    x = BetaQuantile(p,parameters.at("alpha"),parameters.at("beta"));
                    break;
                case(Exponential):
                    x = ExponentialQuantile(p,parameters.at("lambda"));
                    break;
                case(Fixed):
                    x = FixedQuantile(p,parameters.at("mu"));
                    break;
                case(HalfNormal):
                    x = HalfNormalQuantile(p,parameters.at("sigma"));
                    break;
                case(Normal):
                    x = NormalQuantile(p,parameters.at("mu"),parameters.at("sigma"));
                    break;
                case(LogNormal):
                    x = LogNormalQuantile(p,parameters.at("logmu"),parameters.at("logsigma"));
                    break;
                case(TruncatedNormal):
                    x = TruncatedNormalQuantile(p,parameters.at("mu"),parameters.at("sigma"),
                                                  parameters.at("a"),parameters.at("b"));
                    break;
                case(DiscreteTruncatedNormal):
                    x = DiscreteTruncatedNormalQuantile(p,parameters.at("mu"),
                                    parameters.at("sigma"), parameters.at("a"),
                                    parameters.at("b"));
                    break;
                case(Laplace):
                    x = LaplaceQuantile(p,parameters.at("mu"),parameters.at("b"));
                    break;
                case(Poisson):
                    x = PoissonQuantile(p,parameters.at("lambda"));
                    break;
                case(Skellam):
                    x = SkellamCDF(p,parameters.at("lambda"),parameters.at("kappa"));
                    break;
                case(T):
                    x = TQuantile(p,parameters.at("v"));
                    break;
                case(Uniform):
                    x = UniformQuantile(p,parameters.at("a"),parameters.at("b"));
                    break;
            }
        } catch (std::out_of_range & e){
            throw std::invalid_argument("Attempt to GetQuantile with missing parameters");
        } catch (std::domain_error & e){ // All
            throw;
        } catch (std::overflow_error & e){ // Any Dist with support to any infinity
            throw;
        } catch (std::invalid_argument & e){ //Any dist with a scale parameter
            throw;
        } catch (std::logic_error & e){ //Skellam
            throw;
        }
        return x;
    }


    double GetCIWidth(double alpha, DistributionType distribution, ParamMap parameters){
        double width = 0;
        try{
            switch(distribution){
                case(Normal):
                    width = NormalCIWidth(alpha,parameters.at("sigma"));
                    break;
                case(Laplace):
                    width = LaplaceCIWidth(alpha,parameters.at("b"));
                    break;
                default:
                   throw std::logic_error("CIWidth not yet implemented for chosen distribution");
            }
        } catch (std::out_of_range & e){
            throw std::invalid_argument("Attempt to GetCIWidth with missing parameters");
        } catch (std::domain_error & e){ // All
            throw;
        } catch (std::overflow_error & e){ // Any Dist with support to any infinity
            throw;
        }
        return width;
    }


}


