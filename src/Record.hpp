#ifndef __RECORD__HPP
#define __RECORD__HPP

#include <Eigen/Core>
#include "Model.hpp"
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace record {

    typedef std::unordered_map<std::string,double> ValueVectorMap;
    typedef std::set<std::string> StringSet;


    //TODO: Change implementation to judge burnin just off of the nlogP
    class CRecord{
        //Cons/Destruction
        CRecord() = delete;
        CRecord(const std::vector<std::string> & vParamNames);
        //Static Members
        protected:
            static double alpha;
            static double epsilon;
        //Members
        protected:
            size_t nSamples;
            Eigen::MatrixXd sampleMat;
            std::vector<double> vLikelihoods;
        //Const Members
        public:
            const double W;
        protected:
            const size_t nParams;
            const StringSet sParamNames;
        //Static Methods
        protected:
            static int CalculateThreshold(size_t nParam);
            Eigen::MatrixXXd CalculateSampleCovariance(const Eigen::MatrixXXd & data);
            static double GetAlpha(){return CRecord::alpha;}
            static double GetEpsilon(){return CRecord::epsilon;}
            static void TuneThreshold(double alpha, double epsilon);
        //Methods
        public:
            void addSample(const model::CModel & model);
            bool isComplete() const;
            size_t size() const {return this->nSamples;}
        protected:
            double calculateSampleCovarianceDeterminant() const;
            double calculateMultivariateBatchMeansDeterminant() const;
            double estimateEffectiveSampleSize() const;
    };

}


#endif //__RECORD__HPP


