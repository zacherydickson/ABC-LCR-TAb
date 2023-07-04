#ifndef __RECORD__HPP
#define __RECORD__HPP

#include <Eigen/Core>
#include "Model.hpp"
#include <string>
#include <set>


//This implementation of a multivariate Effective Sample Size stopping rule is adapted from
//Vats Flegal and Jones 2017 (doi: 10.48550/arXiv.1512.07713)

namespace record {

    typedef std::set<std::string> StringSet;

    class CRecord{
        //Cons/Destruction
        public:
            CRecord() = delete;
            CRecord(const std::vector<std::string> & vParamNames);
            CRecord(const model::CModel & model);
        //Static Members
        protected:
            static double Alpha;
            static double Epsilon;
        //Members
        protected:
            size_t nSamples;
            Eigen::MatrixXd sampleMat;
            std::vector<double> vLikelihoods;
            size_t nStar;
        //Const Members
        protected:
            const size_t nParams;
            const StringSet sParamNames;
            const double W;
        //Static Methods
        public:
            static double GetAlpha(){return CRecord::Alpha;}
            static double GetEpsilon(){return CRecord::Epsilon;}
            static void TuneThreshold(double alpha, double epsilon);
        protected:
            static int CalculateThreshold(size_t nParam);
        //Methods
        public:
            void addSample(const model::CModel & model);
            bool isComplete();
            size_t nextMilestone() const {return this->nStar;}
            size_t size() const {return this->nSamples;}
        protected:
            double calculateSampleCovarianceDeterminant() const;
            double calculateMultivariateBatchMeansDeterminant() const;
            double estimateEffectiveSampleSize() const;
    };

}


#endif //__RECORD__HPP


