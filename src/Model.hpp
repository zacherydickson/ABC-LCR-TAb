#ifndef MODEL_HPP
#define MODEL_HPP

#include "ctpl_stl.h"
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include "RandomVariable.hpp"
#include <string>
#include "Tree.hpp"
#include <vector>

namespace model {

    enum ModelType {StepwiseOU};
    ModelType str2ModelType(std::string str);

    struct SParameterSpecification{
        std::string dependency;
        double value;
        double lowerBound;
        double upperBound;
        SParameterSpecification & operator=(const SParameterSpecification & rhs){
            this->value = rhs.value;
            this->lowerBound = rhs.lowerBound;
            this->upperBound = rhs.upperBound;
            this->dependency = rhs.dependency;
            return *this;
        }
        bool isFixed(double at = std::numeric_limits<double>::quiet_NaN()) const {
            if(lowerBound != upperBound && dependency.empty()){
                return false;
            }
            if(!std::isnan(at) && value != at){
                return false;
            }
            return true;
        }
    };

    struct ModelState {
        CDiscreteFiniteRandomVariable abundance;
        CDiscreteFiniteRandomVariable length;
        ModelState & operator=(const ModelState & rhs){
            //fprintf(stderr,"ModelState Assignment\n");
           this->abundance = rhs.abundance;
           this->length = rhs.length; 
           return *this;
        }
        bool operator==(const ModelState &rhs) const{
            if(this->abundance == rhs.abundance && this->length == rhs.length){
                return true;
            }
            return false;
        }
        bool operator<(const ModelState & rhs) const{
            if(this->length < rhs.length){
                return true;
            } else if(this->length == rhs.length && this->abundance < rhs.abundance){
                return true;
            }
            return false;
        }
    };

    //Used for evaluation
    struct SVectorizedModelState {
        std::vector<int> vLength;
        std::vector<int> vAbundance;
        std::ostream& output(std::ostream& os) const {
            os << "[L:{";
            int counter = 0;
            for(int l : vLength){
                os << l;
                if(++counter < vLength.size()){
                    os << ",";
                }
            }
            os << "} A:{";
            counter = 0;
            for(int a : vAbundance){
                os << a;
                if(++counter < vAbundance.size()){
                    os << ",";
                }
            }
            os << "}]";
            return os;
        }
    };


    inline std::ostream& operator<<(std::ostream& os, const SVectorizedModelState & obj){return obj.output(os);}

    struct InitialModelState : ModelState {
        std::vector<int> vProtIdxs;
        InitialModelState & operator=(const InitialModelState & rhs){
            //fprintf(stderr,"ModelState Assignment\n");
           this->abundance = rhs.abundance;
           this->length = rhs.length; 
           this->vProtIdxs = rhs.vProtIdxs;
           return *this;
        }
    };

    struct SEvalBlockIdxPartition {
        SEvalBlockIdxPartition() : nThreads(0) {}
        size_t nThreads;
        std::vector<std::vector<int>> vThreadBlocks;
        SEvalBlockIdxPartition & operator=(const SEvalBlockIdxPartition & rhs){
            this->nThreads = rhs.nThreads;
            this->vThreadBlocks = rhs.vThreadBlocks;
            return *this;
        }
    };

    typedef std::vector<ModelState> vModelState;
    typedef std::vector<InitialModelState> vInitialModelState;
    typedef std::vector<int> EvaluationBlock;
    typedef std::vector<EvaluationBlock> vEvaluationBlock;
    typedef std::map<std::string,vModelState> StateMap;
    typedef std::map<std::string,SParameterSpecification> ParamMap;
    typedef std::map<std::string,double> GradientMap;
    typedef SDepthFirstAccessNodeVector<SVectorizedModelState> SVModelStateNodeVector;
    typedef SBasicNode<SVectorizedModelState> SVModelStateNode;
    typedef std::map<std::string,double> ProposalScaleMap;
    
    //Abstract Base Model Class
    class CModel {
        //Con-/Destruction
        public:
            CModel() = delete;
            CModel(vInitialModelState vInitStates, ParamMap params, double logHastingsRatio = 0, double logJointPriorDensity= 1, const vEvaluationBlock * initEvalBlocks = nullptr, const SEvalBlockIdxPartition * initPartition = nullptr);
        //Members
        public: 
            double logHastingsRatio;
            double logJointPriorDensity;
        protected:
            bool bFixed;
            const vInitialModelState vInitStates;
            ParamMap parameters;
            double nlogJointProbability;
            double minEvalValue;
            vEvaluationBlock vEvalBlocks;
            SEvalBlockIdxPartition evalBlockIdxPartition;
            double maxR2TDivL;
            double maxR2TDivA;
        protected:
        //Static Members
        protected:
            static double MaximumRelativeError;
        private:
            static const std::string modelName;
            static const std::vector<std::string> parameterNames;
        protected:
        //Static Methods  
        public:
            static void TuneEvaluation(double error);
            static double GetEvalRelError(){return CModel::MaximumRelativeError;}
        protected:
            double ExtractValueFromModelString(const std::string & modelStr, const std::string  & keyStr, size_t pos);
            static double ProposeParameterSet(ParamMap & parameters, const ProposalScaleMap & scaleMap, std::mt19937 & gen);
            static double SampleProposalDistribution(SParameterSpecification & param, double scale, std::mt19937 & gen);

        //Methods
        public:
            void evaluate(const Tree & tree, const StateMap & obs, ctpl::thread_pool & threadPool,std::mt19937 & gen, size_t nSim);
            size_t getNEvalBlocks() const {return this->vEvalBlocks.size();}
            double getNLogP() const;
            double getMinEval() const;
            const ParamMap & getParamMap() const {return this->parameters;}
            bool isFixed() const {return this->bFixed;}
            std::ostream& output(std::ostream&) const;
            std::unique_ptr<CModel> proposeJump(const ProposalScaleMap & scaleMap, std::mt19937 & gen) const;
            void setToStr(const std::string & modelStr);
        protected:
            virtual void determineEvaluationBlocks();
            double evaluateBlock(const EvaluationBlock & evalBlock, const Tree & tree, const StateMap & obs,std::mt19937 & gen, size_t nSim) const;
            double evaluateBlocks(int, const std::vector<int> & vBlockIdxs, const Tree & tree, const StateMap & obs,size_t seed, size_t nSim) const;
            virtual size_t initializeSimulationRootNode(const EvaluationBlock & evalBlock, SVModelStateNode & rootNode) const;
            void partitionEvalBlocks(size_t nThreads);
            void validateParameters() const ;
        //Virtual Methods
        public:
            virtual const std::string & getName() const {return modelName;}
            virtual const size_t getNParams() const {return parameterNames.size();}
            virtual const std::vector<std::string> & getParamNames() const {return parameterNames;}
        //Pure virtual functions
        public:
            virtual std::unique_ptr<CModel> constructAdjacentModel(ParamMap & newParams, double logHastingsRatio) const = 0;
        protected:
            virtual bool evalIsLess(int a, int b) const = 0;
            virtual bool evalIsEqual(int a, int b) const = 0;
            virtual void sampleSimulationNode(const EvaluationBlock & evalBlock, SVModelStateNode & node, const SVModelStateNode & parent, const SVModelStateNode & root, double time, std::mt19937 & gen) const = 0;
    };

    inline std::ostream& operator<<(std::ostream& os, const CModel & obj){return obj.output(os);}

    class CStepwiseOUModel : public CModel {
        //Description:
        //  Abundance Evol is modeled as an OU process where the mean value at a node is a
        //  weighted average of the parent value and the selective optimum
        //  The Selective Optimum is influenced by the length of the minimum entropy region
        //  at the parent node
        //  Length Evol Is modeled as a stepwise process where the rate of indels increases
        //  linearly with length, over a time period t there are poisson distribtued events
        //  each with a probability of being an insertion or deletion
        //  The mean number of indels is influenced by the fold change in abundance at the
        //  parent node
        public:
            CStepwiseOUModel() = delete;
            CStepwiseOUModel(vInitialModelState vInitStates, ParamMap params, double logHastingsRatio = 0, double logJointPriorDensity= 1, const vEvaluationBlock * initEvalBlocks = nullptr, const SEvalBlockIdxPartition * initPartition = nullptr);
        //Static Members
        private:
            static const std::string modelName;
            static const std::vector<std::string> parameterNames;
        //Members:
        public:
            bool bFixedZeroTau;
            bool bFixedZeroUpsilon;
        //Methods
        private:
        //Overridden Methods
        public:
            const std::string & getName() const override {return modelName;}
            const size_t getNParams() const override {return parameterNames.size();}
            const std::vector<std::string> & getParamNames() const override {return parameterNames;}
        private:
            std::unique_ptr<CModel> constructAdjacentModel(ParamMap & newParams, double logHastingsRatio) const override;
            bool evalIsLess(int a, int b) const override;
            bool evalIsEqual(int a, int b) const override;
            int sampleAbundance(int protIdx, SVModelStateNode & node, const SVModelStateNode & parent, const SVModelStateNode & root, double time, std::mt19937 & gen) const;
            int sampleLength(int protIdx, SVModelStateNode & node, const SVModelStateNode & parent, const SVModelStateNode & root, double time, std::mt19937 & gen) const;
            void sampleSimulationNode(const EvaluationBlock & evalBlock, SVModelStateNode & node, const SVModelStateNode & parent, const SVModelStateNode & root, double time, std::mt19937 & gen) const override;
    };


      
} //namespace model

#endif



