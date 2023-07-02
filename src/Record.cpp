#include <cmath>
#include <Eigen/LU>
#include "Distributions.hpp"
#include "Logging.h"
#include "Record.hpp"
#include <stdexcept>

#include <sstream>

namespace record {

    // Con/Destruction, CRecord

    CRecord::CRecord(const std::vector<std::string> & vParamNames) : 
        nSamples(0), nParams(vParamNames.size()),
        sParamNames(vParamNames.begin(),vParamNames.end()),
        W(CRecord::CalculateThreshold(vParamNames.size())) 
    {
        this->nStar = std::max(std::pow(this->nParams,2.0),this->W);
        sampleMat.conservativeResize(this->nParams,0);
    }

    CRecord::CRecord(const model::CModel & model) :
        CRecord(model.getParamNames())
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
        double denom1 = std::pow(nParam * std::tgamma(nParam / 2.0),exponent);
        double num2 = stats::ChiSqQuantile(1-CRecord::alpha,nParam);
        double denom2 = std::pow(CRecord::epsilon,2.0);
        return num1 /denom1 * num2 / denom2;
    }

    void CRecord::TuneThreshold(double alpha, double epsilon){
        if(alpha <= 0 || alpha > 1){
            throw std::invalid_argument("Attempt to tune CRecord threshold with an alpha outside (0,1]");
        }
        if(epsilon <= 0 || epsilon > 1){
            throw std::invalid_argument("Attempt to tune CRecord threshold with an epsilon outside (0,1]");
        }
        CRecord::alpha = alpha;
        CRecord::epsilon = epsilon;
        //fprintf(stderr,"alpha: %0.04f\tepsilon: %0.04f\n",CRecord::alpha, CRecord::epsilon);
    }
    

    // Methods, CRecord, public

    void CRecord::addSample(const model::CModel & model) {
        this->vLikelihoods.push_back(model.getNLogP());
        const model::ParamMap & params = model.getParamMap();
        this->sampleMat.conservativeResize(Eigen::NoChange,++(this->nSamples));
        int row = 0;
        for(const std::string & name : this->sParamNames){
            try {
                this->sampleMat(row++,this->nSamples-1) = params.at(name).value;
            } catch (std::out_of_range & e){
                throw std::invalid_argument("Attempt to add a sample to CRecord which is missing paramaters");
            }
        }
    }

    bool CRecord::isComplete() {
        std::stringstream stream;
        //stream << this->sampleMat;
        //fprintf(stderr,"W:%0.04f\ndata:\n%s\n",this->W,stream.str().c_str());
        if(this->nSamples <= this->nStar){
            return false;
        }
        double ess = this->estimateEffectiveSampleSize();
        if(ess < this->W){
            if(ess > 0){
                this->nStar *= this->W / ess;
            } else {
                this->nStar += this->W;
            }
            return false;
        }
        return true;
    }
    // Methods, CRecord, protected

    double CRecord::calculateSampleCovarianceDeterminant() const {
        auto mean = this->sampleMat.rowwise().mean();
        auto unitRow = Eigen::MatrixXd::Constant(1,this->sampleMat.cols(),1.0);
        auto fullmean = mean * unitRow;
        auto demeaned = this->sampleMat - fullmean;
        auto sampleCoV =  demeaned * demeaned.transpose() / (this->nSamples - 1.0);
        //std::stringstream stream;
        //stream << sampleCoV;
        //fprintf(stderr,"sampleCoV:\n%s\n",stream.str().c_str());
        return sampleCoV.determinant();
    }

    double CRecord::calculateMultivariateBatchMeansDeterminant() const {
        size_t batchSize = std::sqrt(this->nSamples); 
        size_t nBatches = this->nSamples / batchSize;
        if(this->nSamples % batchSize){
            nBatches++;
        }
        Eigen::MatrixXd mBM;
        mBM.conservativeResize(this->nParams,nBatches);
        for(int batch = 0; batch < nBatches; batch++){
            int first = batch * batchSize;
            int last = first + batchSize - 1;
            if(last >= this->nSamples){
                last = this->nSamples-1;
            }
            auto subset = this->sampleMat(Eigen::all,Eigen::seq(first,last));
            auto batchMean = subset.rowwise().mean();
            mBM.col(batch) = batchMean;
        }
        //std::stringstream stream;
        //stream << "mBM:\n" << mBM << "\n";
        auto mean = this->sampleMat.rowwise().mean();
        auto unitRow = Eigen::MatrixXd::Constant(1,mBM.cols(),1.0);
        auto fullmean = mean * unitRow;
        auto demeaned = mBM - fullmean;
        auto batchMeanCoV =  demeaned * demeaned.transpose() * batchSize  / (nBatches - 1.0);
        //stream << "batchCov:\n" << batchMeanCoV << "\n";
        //fprintf(stderr,"%s",stream.str().c_str());
        return batchMeanCoV.determinant();
    }

    double CRecord::estimateEffectiveSampleSize() const{
        logger::Log("Calculating multivariate effective sample size ...",logger::INFO);
        double sampleCovDet = this->calculateSampleCovarianceDeterminant();
        double mBMDet = this->calculateMultivariateBatchMeansDeterminant();
        double ratio = sampleCovDet / mBMDet;
        double ess = this->nSamples * std::pow(ratio,1.0/double(this->nParams));
        logger::Log("ESS Calculated at %0.1f",logger::INFO,ess);
        return ess;
    }
}
