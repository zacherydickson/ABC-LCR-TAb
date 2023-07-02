#ifndef __RECORD__HPP
#define __RECORD__HPP

#include <Eigen/Core>
#include "Model.hpp"
#include <string>
#include <set>

namespace record {

    typedef std::set<std::string> StringSet;

    //TODO: Change implementation to judge burnin just off of the nlogP
    class CRecord{
        //Cons/Destruction
        public:
            CRecord() = delete;
            CRecord(const std::vector<std::string> & vParamNames);
            CRecord(const model::CModel & model);
        //Static Members
        protected:
            static double alpha;
            static double epsilon;
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
            static double GetAlpha(){return CRecord::alpha;}
            static double GetEpsilon(){return CRecord::epsilon;}
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


