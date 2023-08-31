#include <cmath>
#include <limits>
#include "Logging.h"
#include <numeric>
#include "Prior.hpp"
#include <set>
#include <stdexcept>
#include <sstream>


namespace prior {
    
    CPrior::CPrior(model::ModelType type, const ParameterPriorMap & mParamPriors):
        modelType(type),mParamPriors(mParamPriors)
    {
        //Validate the model type is valid
        switch(type){
            case model::StepwiseOU:
                break;
            default:
                throw std::invalid_argument("Attempt to construct CPrior for an unknown model");
        }
        //Determine the number of Initial States, ensure there is at least one
        do {
            std::string A("InitA_" + std::to_string(vInitStates.size()));
            std::string L("InitL_" + std::to_string(vInitStates.size()));
            bool bA = this->mParamPriors.count(A);
            bool bL = this->mParamPriors.count(L);
            if(!bA || !bL){ //If either are not present
                if(vInitStates.size() == 0){ //If no Inits have been seen
                    throw std::invalid_argument("Attempt to construct Prior without an initial abundance and length state for at least one proteins");
                }
                if(bA || bL){ //If one is present and the other is not
                    std::string msg("Attempt to construct Prior with an initial protein missing ");
                    msg += ((bA) ? "Length" : "Abundance");
                    throw std::invalid_argument(msg.c_str());
                }
                break; //Stop processing
            }
            //Construct the initial state
            model::InitialModelState initState;
            initState.abundance = CPrior::ConstructDFRV(this->mParamPriors[A]);
            initState.length = CPrior::ConstructDFRV(this->mParamPriors[L]);
            this->mParamPriors.erase(A);
            this->mParamPriors.erase(L);
            this->vInitStates.push_back(initState);
        } while(true);
        this->validateParamPriors();
        this->collapseDuplicateInitialStates();
    }

    //STATIC METHODS
    
    //protected

    double CPrior::CalculatePriorDensity(const std::string & name, double value, stats::DistributionType type, const stats::ParamMap & hyperparameters){
        double p = 1.0;
        try{
            stats::GetPDF(value,type,hyperparameters);
        } catch (std::exception & e){
            std::stringstream msg;
            msg << "An error occured calculating the prior density for ";
            msg << name << " @ " << value << "; ";
            msg << e.what();
            throw std::runtime_error(msg.str().c_str());
        }
        return p;
    }

    CDiscreteFiniteRandomVariable CPrior::ConstructDFRV(SParameterPriorSpecification paramPrior){
        std::vector<int> vVals;
        std::vector<double> vProbs;
        stats::SDomain domain = stats::GetDomain(paramPrior.type,paramPrior.hyperparameters);
        if(std::isinf(domain.min) || std::isinf(domain.max)){
            throw std::invalid_argument("Attempt to call CPrior::ConstructDFRV for a parameter with an infinite domain");
        }
        for(int i = domain.min; i <= domain.max; i++){
            vVals.push_back(i);
            double p = CPrior::CalculatePriorDensity("InitialState",(double)i,paramPrior.type,paramPrior.hyperparameters);
            vProbs.push_back(p);
        }
        return CDiscreteFiniteRandomVariable(vVals,vProbs);
    }

    //METHODS

    // concrete, public
    
    double CPrior::calculateJointPriorDensity(const model::ParamMap & parameters) const{
        //Ass all model parameters are independent, the joint density is product of the
        //marginal densities. This value might be small so we'll use the log density
        double logDensity = 0;
        for(const auto & pair : parameters){ //Iterate over model parameters
            //Note: pair.first is the parameter name which matches the PriorParameter names
            //      pair.second is a model::SParameterSpecification Object
            const SParameterPriorSpecification & priorSpec = this->mParamPriors.at(pair.first);
            double p = CPrior::CalculatePriorDensity(pair.first,pair.second.value,priorSpec.type,priorSpec.hyperparameters);
            if(p == 0){
                return -std::numeric_limits<double>::infinity();
            }
            logDensity += std::log(p);
        }
        return logDensity;
    }

    std::unique_ptr<model::CModel> CPrior::GenerateModel(std::mt19937 & gen) const{
        logger::Log("Generating a model from prior...",logger::DEBUG+1);
        model::ParamMap parameters;
        //Generate random parameter values for each parameter
        double logDensity = 0;
        for(const auto & pair : this->mParamPriors){
            stats::SDomain domain = stats::GetDomain(pair.second.type,pair.second.hyperparameters);
            model::SParameterSpecification param;
            param.lowerBound = domain.min;
            param.upperBound = domain.max;
            param.dependency = pair.second.dependency;
            double p = 1.0;
            //Initialize values for independent parameters
            if(param.dependency.empty()){
                param.value = stats::GetQuantile(stats::generate_open_canonical(gen),pair.second.type,pair.second.hyperparameters);
                p = CPrior::CalculatePriorDensity(pair.first,param.value,pair.second.type,pair.second.hyperparameters);
            }
            logDensity += std::log(p);
            parameters.insert(std::make_pair(pair.first,param));
        }
        //Initialize values for dependent parameters
        for(auto & pair : parameters){
            if(!pair.second.dependency.empty()){
                pair.second.value = parameters.at(pair.second.dependency).value;
            }
        }
        //Determine which derived model to output
        switch(this->modelType){
            case model::StepwiseOU:
                logger::Log("Returning StepwiseOUModel",logger::DEBUG+1);
                return std::unique_ptr<model::CModel>(new model::CStepwiseOUModel(this->vInitStates,parameters,0,logDensity));
            default:
                logger::Log("Returning Null Base Model",logger::DEBUG+1);
                return std::unique_ptr<model::CModel>(nullptr);
        }
    }

    std::vector<std::string> CPrior::getModelParameterNames(bool bIncludeFixed) const{
        std::mt19937 gen;
        auto ptr = this->GenerateModel(gen);
        std::vector<std::string> vParamNames;
        for(const auto & pair : ptr->getParamMap()){
            if(bIncludeFixed || !pair.second.isFixed()){
                vParamNames.push_back(pair.first);
            }
        }
        return vParamNames;
    }

    std::ostream & CPrior::output(std::ostream & os) const {
        os << "mType{" << this->modelType << "}"; 
        for(int i = 0; i < this->vInitStates.size(); i++){
            os << "\t[A:" << vInitStates[i].abundance.getMode() << ","
               << "L:" << vInitStates[i].length.getMode() << "]";
        }
        for(const auto & pair : mParamPriors){
            os << "\t" << pair.first << "(dType{" << pair.second.type << "}";
            for(const auto & pair2 : pair.second.hyperparameters){
                os << "," << pair2.first << ":" << pair2.second;
            }
            os << ")";
        }
        return os;
    }

    //CPrior: concrete, protected

    void CPrior::collapseDuplicateInitialStates(){
        std::vector<int> vIdxs(this->vInitStates.size());
        std::iota(vIdxs.begin(),vIdxs.end(),0);
        //Put the indexes in ascending order of their state (length then abundance)
        std::sort(vIdxs.begin(),vIdxs.end(),
                [&](int a, int b) -> bool {
                    return this->vInitStates[a] < this->vInitStates[b];
                });
        this->vInitStates[vIdxs[0]].vProtIdxs.push_back(vIdxs[0]);
        for(int i = 0; i < vIdxs.size()-1; i++){
            auto it1 = vIdxs.begin()+i;
            auto it2 = vIdxs.begin()+(i+1);
            //Check if prot[idx[i]]'s state is equal to the subsequent state
            if(this->vInitStates[*it1] == this->vInitStates[*it2]){
                this->vInitStates[*it1].vProtIdxs.push_back(*it2); //Add the next prot to the processing list
                vIdxs.erase(it2); //Stop considering this protein in equality comparisons
                i--; //Take  step back so that the next iteration compares the current state again
            } else { //Indicate the next protein will process itself
                this->vInitStates[*it2].vProtIdxs.push_back(*it2);
            }
        } //End of loop all weight from duplicates is reassigned to a single entry
        //The number of proteins remains constant, just some computation can be skipped later
    }

    void CPrior::validateParamPriors() const {
        for(const auto & pair : this->mParamPriors){
            std::string curParam = pair.first;
            std::set<std::string> sDependents;
            //If it has a dependency make sure it isn't circular or chained
            while(!mParamPriors.at(curParam).dependency.empty()){
                sDependents.insert(curParam);
                //Follow the dependency line until termination, or until encountering a
                //previously encountered parameter
                curParam = mParamPriors.at(curParam).dependency;
                if(mParamPriors.count(curParam) == 0){
                    logger::Log("Identified unknown dependency (%s) for parameter %s",
                            logger::ERROR,curParam.c_str(),pair.first.c_str());
                    throw std::invalid_argument("Encountered unknown parameter dependency");
                }
                if(sDependents.count(curParam) > 0){
                    logger::Log("Identified circular dependency for parameter %s",
                            logger::ERROR,pair.first.c_str());
                    throw std::invalid_argument("Encountered circular parameter dependency");
                }
                if(!mParamPriors.at(curParam).dependency.empty()){
                    logger::Log("Identified chained dependency for parameter %s",
                            logger::ERROR,pair.first.c_str());
                    throw std::invalid_argument("Chained dependencies are currently not allowed");
                }
            }
        }
    }

}
