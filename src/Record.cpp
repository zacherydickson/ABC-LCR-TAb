#include <cmath>
#include "Record.hpp"
#include "Distributions.hpp"
#include <stdexcept>

namespace record {

    // Con/Destruction, CRecord

    CRecord::CRecord(const std::vector<std::string> & vParamNames) : 
        nSamples(0), nParams(vParamNames.size()),
        sParamNames(vParamNames.begin(),vParamNames.end()),
        W(CRecord::CalculateThreshold(vParamNames.size())) 
    {
    }

    CRecord::CRecord(const model::CModel & model) :
        CRecord(model.getParamNames)
    {
        this->addSample(model);
    }

    // Static Members, CRecord, protected

    double CRecord::alpha = 0.05;
    double CRecord::epsilon = 0.05;

    // Static Methods, CRecord, protected

    int CRecord::CalculateThreshold(size_t nParam){
        double exponent = 2.0 / double(nParam);
        double num1 = std::pow(2.0,exponent) * stats::pi;
        double denom1 = std::pow(nParam * std:tgamma(nParam / 2.0),exponent);
        double num2 = stats::ChiSqQuantile(1-CRecord::alpha,nParam);
        double denom2 = std::pow(CRecord::epsilon,2.0);
        return num1 /denom1 * num2 / denom2;
    }

    void CRecord::TuneThreshold(double alpha, double epsilon){
        if(alpha <= 0 || alpha > 1){
            throw std::invalid_argument("Attempt to tune CRecord threshold with an alpha outside (0,1]");
        }
        if(espsilon <= 0 || epsilon > 1){
            throw std::invalid_argument("Attempt to tune CRecord threshold with an epsilon outside (0,1]");
        }
        CRecord::alpha = alpha;
        CRecord::epsilon = epsilon;
    }
    

    // Methods, CRecord, public

    void CRecord::addSample(const model::CMode & model) {
        this->vLikelihoods.push_back(model->getNLogP());
        const model::ParamMap & params = model.getParamMap();
        bool vValid = (params.size() == this->nParamss);
        for(const auto & pair : params){
            if(!bValid || !this->sParamNames.contains(pair.first)){
                bValid = false;
                break;
            }
        }
        if(!bValid){
            throw std_invalid_argument("Attempt to add a sample to CRecord with different model than that with which the record  was initialized");
        }
        for(const auto & pair : params){
            this->vSamples[pair.first].push_back(pair.second.value);
        }
        this->nSamples++;
    }

    bool CRecord::isComplete() const {
        double burnin = std::max(std::pow(this->nParams,2.0),this->W);
        if(this->nSamples <= burnin){
            return false;
        }
        double ess = this->estimateEffectiveSampleSize();
        if(ess < this->W){
            return false;
        }
        return true;
    }
    // Methods, CRecord, protected

    double CRecord::calculateSampleCovariance() const {
        Eigen::Matrix<double, this->nParam, this->nParam > sampleCoV();
        //TODO:
        return sampleCov.determinant();
    }

    double CRecord::calculateMultivariateBatchMeans() const {
        //TODO:
        Eigen::Matrix<double, Dynamic, Dynamic> sampleCoV();
        return mBM.determinant();
    }

    double CRecord::estimateEffectiveSampleSizes(){
        double sampleCovDet = this->calculateSampleCovariance();
        double mBMDet = this->calculateMultivariateBatchMeans();
        double ratio = sampleCovDet / mBMDet;
        double ess = this->nSamples * std::pow(ratio,1.0/double(this->nParams));
        return ess;
    }
}
