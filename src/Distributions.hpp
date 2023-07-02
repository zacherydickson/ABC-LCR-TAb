#ifndef _DISTRIBUTIONS_HPP
#define _DISTRIBUTIONS_HPP

#include <algorithm>
#include <cmath>
#include <map>
#include <random>
#include <string>

namespace stats {

    const long double pi = std::atan(1)*4;
    const long double root2 = std::sqrt((long double)2.0);

    enum DistributionType {Beta, ChiSq, Exponential,DiscreteTruncatedNormal, Fixed, HalfNormal, Normal, LogNormal, Laplace, Poisson, Skellam, T, TruncatedNormal, Uniform };

    struct SDomain{
        double min;
        double max;
    };

    typedef std::map<std::string,double> ParamMap; 

    DistributionType str2DistributionType(std::string str);

    double FixedPMF(double x, double mu);
    double FixedCDF(double x, double mu);
    double FixedQuantile(double p, double mu);

    double ChiSqPDF(double x, double df);
    double ChiSqCDF(double x, double df);
    double ChiSqQuantile(double p, double df);

    double NormalPDF(double x, double mu = 0.0, double sigma = 1.0);
    long double NormalCDF(double x, double mu = 0.0, double sigma = 1.0, bool bLogP = false); 
    double NormalQuantile(double p, double mu = 0.0, double sigma = 1.0);
    double NormalCIWidth(double alpha, double sigma);

    double LogNormalPDF(double x, double logmu, double logsigma);
    double LogNormalCDF(double x, double logmu, double logsigma);
    double LogNormalQuantile(double p, double logmu, double logsigma);

    double TruncatedNormalPDF(double x, double mu, double sigma, double a, double b);
    double TruncatedNormalCDF(double x, double mu, double sigma, double a, double b);
    double TruncatedNormalQuantile(double p, double mu, double sigma, double a, double b);

    long double DiscreteTruncatedNormalPMF(double x, double mu, double sigma, double a, double b, bool bLogP = false);
    double DiscreteTruncatedNormalCDF(double x, double mu, double sigma, double a, double b);
    double DiscreteTruncatedNormalQuantile(double p, double mu, double sigma, double a, double b);

    double HalfNormalPDF(double x, double sigma);
    double HalfNormalCDF(double x, double sigma);
    double HalfNormalQuantile(double p, double sigma);

    double LaplacePDF(double x, double mu, double b);
    double LaplaceCDF(double x, double mu, double b);
    double LaplaceQuantile(double p, double mu, double b);
    double LaplaceCIWidth(double alpha, double b);

    double PoissonPMF(double k, double lambda);
    double PoissonCDF(double k, double lambda);
    double PoissonQuantile(double p, double lambda);

    double ExponentialPDF(double x, double lambda);
    double ExponentialCDF(double x, double lambda);
    double ExponentialQuantile(double p, double lambda);

    double BetaPDF(double x, double alpha, double beta);
    double BetaCDF(double x, double alpha, double beta);
    double BetaQuantile(double p, double alpha, double beta);
    double BetaMaxMonotonicVariance(double x);

    double SkellamPMF(double x, double lambda, double kappa);
    double SkellamCDF(double x, double lambda, double kappa);
    double SkellamQuantile(double p, double lambda, double kappa);

    double TPDF(double x, double v);
    double TCDF(double x, double v);
    double TQuantile(double p, double v);

    double UniformPDF(double a, double b);
    double UniformCDF(double x, double a, double b);
    double UniformQuantile(double p, double a, double b);

    template<class Generator>
    double generate_open_canonical(Generator & g){
        double x;
        while((x = std::generate_canonical<double,10>(g)) == 0);
        return x;
    }
    SDomain GetDomain(DistributionType distribution, ParamMap parameters);
    double GetPDF(double X, DistributionType distribution, ParamMap parameters);
    double GetCDF(double X, DistributionType distribution, ParamMap parameters);
    double GetQuantile(double p, DistributionType distribution, ParamMap parameters);
    double GetCIWidth(double alpha, DistributionType distribution, ParamMap parameters);

    ///////Template Functions

    template<typename numericType>
    double GetMedian(std::vector<numericType> vector){
        if(vector.size() == 0)
            return numericType(0);
        int m = std::floor((vector.size()-1) / 2);
        std::nth_element(vector.begin(), vector.begin() + m, vector.end());
        double median = vector[m];
        if(vector.size() % 2 == 0){
            std::nth_element(vector.begin() + m, vector.begin() + m + 1, vector.end());
            median += vector[m+1];
            median /= 2.0;
        }
        return median;
    }

    template<typename numericType>
    double GetMeanAbsoluteDeviation(const std::vector<numericType> & vector, double median){
        double n = vector.size();
        if(vector.size() == 0)
            return numericType(1);
        double mad = 0;
        for(const numericType & x : vector){
            double res = (double) x - median;
            if(res < 0)
                res *= -1;
            mad += res;
        }
        mad /= n;
        return mad;
    }
    
    template<typename numericType>
    double GetMeanAbsoluteDeviation(const std::vector<numericType> & vector){
        numericType median = GetMedian<numericType>(vector);
        return GetMeanAbsoluteDeviation<numericType>(vector,median);
    }
    
    template<typename numericType>
    double GetMean(const std::vector<numericType> & vector){
        double n = vector.size();
        if(vector.size() == 0)
            return numericType(0);
        double mean = 0;
        for(const auto & x : vector){
            mean += x;
        }
        mean /= n;
        return mean;
    }
    
    template<typename numericType>
    double GetStandardDeviation(const std::vector<numericType> & vector, double mean){
        if(vector.size() < 2)
            return numericType(1);
        double stdev = 0;
        for(const numericType & x : vector){
            stdev += (x - mean) * (x - mean);
        }
        return std::sqrt(stdev / (double)(vector.size() - 1));
    }
    
    template<typename numericType>
    double GetStandardDeviation(const std::vector<numericType> & vector){
        double mean = GetMean<numericType>(vector);
        return GetStandardDeviation<numericType>(vector,mean);
    }
}

#endif



