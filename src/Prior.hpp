#ifndef _PRIOR_HPP
#define _PRIOR_HPP

#include "Distributions.hpp"
#include <iostream>
#include "Model.hpp"
#include <memory>
#include <string>
#include <random>
#include <vector>

namespace prior {

    struct SParameterPriorSpecification {
        stats::DistributionType type;
        stats::ParamMap hyperparameters;
        std::string dependency;
        SParameterPriorSpecification & operator=(const SParameterPriorSpecification & rhs){
            this->type = rhs.type;
            this->hyperparameters = rhs.hyperparameters;
            this->dependency = rhs.dependency;
            return *this;
        }
    };

    typedef std::map<std::string,SParameterPriorSpecification> ParameterPriorMap;

    class CPrior {
        //Con-/Destruction
        public:
            CPrior() = delete;
            CPrior(model::ModelType type, const ParameterPriorMap & mParamPriors);
        //Members
        protected:
            model::vInitialModelState vInitStates;
            model::ModelType modelType;
            ParameterPriorMap mParamPriors;
        //Static Methods
        protected:
            static double CalculatePriorDensity(const std::string & name, double value, stats::DistributionType type, const stats::ParamMap & hyperparameters);
            static CDiscreteFiniteRandomVariable ConstructDFRV(SParameterPriorSpecification paramPrior);
        //Methods
        public:
            double calculateJointPriorDensity(const model::ParamMap & parameters) const;
            model::ModelType getModelType() {return this->modelType;}
            std::vector<std::string> getModelParameterNames(bool bIncludeFixed = true) const;
            size_t getNProt() const {return vInitStates.size();}
            size_t getNParamPriors() const {return mParamPriors.size();}
            std::unique_ptr<model::CModel> GenerateModel(std::mt19937 & gen) const;
            std::ostream& output(std::ostream & os) const;
        protected:
            void collapseDuplicateInitialStates();
            void validateParamPriors() const;
        
    };

    inline std::ostream& operator<<(std::ostream& os, const CPrior & obj){return obj.output(os);}

}


#endif



