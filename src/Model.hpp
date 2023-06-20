#ifndef MODEL_HPP
#define MODEL_HPP

#include "ctpl_stl.h"
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include "RandomVariable.hpp"
#include <string>
#include "Tree.hpp"
#include <vector>

namespace model {

    enum ModelType {StepwiseOU,OUStepwise,UnifiedStepwiseOU};
    ModelType str2ModelType(std::string str);

    struct SParameterSpecification{
        double value;
        double lowerBound;
        double upperBound;
        SParameterSpecification & operator=(const SParameterSpecification & rhs){
            this->value = rhs.value;
            this->lowerBound = rhs.lowerBound;
            this->upperBound = rhs.upperBound;
            return *this;
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
    };

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

    typedef std::vector<ModelState> vModelState;
    typedef std::vector<InitialModelState> vInitialModelState;
    typedef std::vector<std::reference_wrapper<const InitialModelState>> EvaluationBlock;
    typedef std::vector<EvaluationBlock> vEvaluationBlock;
    typedef std::map<std::string,vModelState> StateMap;
    typedef std::map<std::string,SParameterSpecification> ParamMap;
    typedef std::map<std::string,double> GradientMap;
    typedef std::pair<ProbabilityMap,ProbabilityMap> SStateNodeValue;
    typedef SDepthFirstAccessNodeVector<SStateNodeValue> SStateNodeVector;
    typedef SDepthFirstAccessNodeVector<double> ModeNodeVector;
    typedef SDepthFirstAccessNodeVector<SVectorizedModelState> SVModelStateNodeVector;
    typedef SBasicNode<SVectorizedModelState> SVModelStateNode;
    typedef std::map<std::string,double> ProposalScaleMap;
    
    //Abstract Base Model Class
    class CModel {
        //Con-/Destruction
        public:
            CModel() = delete;
            CModel(vInitialModelState vInitStates, ParamMap params, double logHastingsRatio = 0, double logJointPriorDensity= 1);
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
            double maxR2TDivL;
            double maxR2TDivA;
        protected:
        //Static Members
        protected:
            static double GoldenSearchBoundMultiple;
            static double GoldenSearchToleranceProportion;
            static double GradientEstimateStepProportion;
            static double MaximumRelativeError;
        private:
            static const std::string modelName;
            static const std::vector<std::string> parameterNames;
        protected:
        //Static Methods  
        public:
            static void TuneEvaluation(double error);
            static void TuneGoldenSearch(double mult, double tol);
            static void TuneGradientEstimation(double prop);
            static double GetGSMult(){return CModel::GoldenSearchBoundMultiple;}
            static double GetGSTolProp(){return CModel::GoldenSearchToleranceProportion;}
            static double GetGEStepProp(){return CModel::GradientEstimateStepProportion;}
            static double GetEvalRelError(){return CModel::MaximumRelativeError;}
        protected:
            double ExtractValueFromModelString(const std::string & modelStr, const std::string  & keyStr, size_t pos);
            static double ProposeParameterSet(ParamMap & parameters, const ProposalScaleMap & scaleMap, std::mt19937 & gen);
            static double SampleProposalDistribution(SParameterSpecification & param, double scale, std::mt19937 & gen);

        //Methods
        public:
            virtual void determineEvaluationBlocks();
            GradientMap estimateGradient(const Tree & tree, const StateMap & obs, ctpl::thread_pool & threadPool,std::mt19937 & gen, size_t nSim) const;
            void evaluate(const Tree & tree, const StateMap & obs, ctpl::thread_pool & threadPool,std::mt19937 & gen, size_t nSim);
            double evaluateBlock(const EvaluationBlock & evalBlock, const Tree & tree, const StateMap & obs,std::mt19937 & gen, size_t nSim) const;
            double evaluateBlocks(int, const std::vector<int> & vBlockIdxs, const Tree & tree, const StateMap & obs,size_t seed, size_t nSim) const;
            size_t getNEvalBlocks() const {return this->vEvalBlocks.size();}
            double getNLogP() const;
            double getMaxVectorDist(const GradientMap & gradient) const;
            double getMinEval() const;
            const ParamMap & getParamMap() {return this->parameters;}
            std::unique_ptr<CModel> goldenSearch(const Tree & tree, const StateMap & obs, ctpl::thread_pool & threadPool,std::mt19937 & gen, size_t nSim) const;
            bool isFixed() const {return this->bFixed;}
            std::ostream& output(std::ostream&) const;
            std::unique_ptr<CModel> proposeJump(const ProposalScaleMap & scaleMap, std::mt19937 & gen) const;
            void setToStr(const std::string & modelStr);
        protected:
            std::vector<std::vector<int>> partitionEvalBlocks(size_t nThreads) const;
        //Virtual Methods
        public:
            virtual const std::string & getName() const {return modelName;}
            virtual const size_t getNParams() const {return parameterNames.size();}
            virtual const std::vector<std::string> & getParamNames() const {return parameterNames;}
        //Pure virtual functions
        public:
            virtual std::unique_ptr<CModel> constructAdjacentModel(ParamMap & newParams, double logHastingsRatio) const = 0;
        protected:
            virtual bool evalIsLess(int a, int b) = 0;
            virtual bool evalIsEqual(int a, int b) = 0;
            virtual size_t initializeSimulationRootNode(const EvaluationBlock & evalBlock, SVModelStateNode & rootNode) const = 0;
            virtual void sampleSimulationNode(const EvaluationBlock & evalBlock, SVModelStateNode & node, const SVModelStateNode & parent, const SVModelStateNode & root, double time, std::mt19937 & gen) const = 0;
    };

    inline std::ostream& operator<<(std::ostream& os, const CModel & obj){return obj.output(os);}

    class CStepwiseOUModel : public CModel {
        //Description:
        //  Length Evol Is modeled as a stepwise process where the rate of indels increases
        //  linearly with length, over a time period t there are poisson distribtued events
        //  each with a probability of being an insertion or deletion
        //  Abundance Evol is modeled as an OU process where the mena value at a node is a
        //  weighted average of the parent value and the selective optimum
        //  The Selective Optimum is influenced by the length of the minimum entropy region
        //  at the current node
        //Con-/Destruction
        public:
            CStepwiseOUModel() = delete;
            CStepwiseOUModel(vInitialModelState vInitStates, ParamMap params, double logHastingsRatio = 0, double logJointPriorDensity= 1);
        //Static Members
        private:
            static const std::string modelName;
            static const std::vector<std::string> parameterNames;
        //Methods
        private:
        //Overridden Methods
        public:
            const std::string & getName() const override {return modelName;}
            const size_t getNParams() const override {return parameterNames.size();}
            const std::vector<std::string> & getParamNames() const override {return parameterNames;}
        private:
            std::unique_ptr<CModel> constructAdjacentModel(ParamMap & newParams, double logHastingsRatio) const override;
            size_t initializeSimulationRootNode(const EvaluationBlock & evalBlock, SVModelStateNode& rootNode) const override;
            bool evalIsLess(int a, int b) override {return this->vInitStates[a].length < this->vInitStates[b].length;}
            bool evalIsEqual(int a, int b) override {return this->vInitStates[a].length == this->vInitStates[b].length;}
            void sampleSimulationNode(const EvaluationBlock & evalBlock, SVModelStateNode & node, const SVModelStateNode & parent, const SVModelStateNode & root, double time, std::mt19937 & gen) const override;
    };

    class COUStepwiseModel : public CModel {
        //Description:
        //  Abundance Evol is modeled as an OU process where the mena value at a node is a
        //  weighted average of the parent value and the selective optimum
        //  Length Evol Is modeled as a stepwise process where the rate of indels increases
        //  linearly with length, over a time period t there are poisson distribtued events
        //  each with a probability of being an insertion or deletion
        //  The mean number of indels is influenced by the fold change in abundance at the
        //  current node relative to the root
        //Con-/Destruction
        public:
            COUStepwiseModel() = delete;
            COUStepwiseModel(vInitialModelState vInitStates, ParamMap params, double logHastingsRatio = 0, double logJointPriorDensity= 1);
        //Static Members
        private:
            static const std::string modelName;
            static const std::vector<std::string> parameterNames;
        //Methods
        private:
        //Overridden Methods
        public:
            const std::string & getName() const override {return this->modelName;}
            const size_t getNParams() const override {return this->parameterNames.size();}
            const std::vector<std::string> & getParamNames() const override {return this->parameterNames;}
        private:
            std::unique_ptr<CModel> constructAdjacentModel(ParamMap & newParams, double logHastingsRatio) const override;
            virtual bool evalIsLess(int a, int b) override {return this->vInitStates[a].abundance < this->vInitStates[b].abundance;}
            virtual bool evalIsEqual(int a, int b) override {return this->vInitStates[a].abundance == this->vInitStates[b].abundance;}
            size_t initializeSimulationRootNode(const EvaluationBlock & evalBlock, SVModelStateNode& rootNode) const override;
            void sampleSimulationNode(const EvaluationBlock & evalBlock, SVModelStateNode & node, const SVModelStateNode & parent, const SVModelStateNode & root, double time, std::mt19937 & gen) const override;
    };

    class CUnifiedStepwiseOUModel : public CModel {
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
            CUnifiedStepwiseOUModel() = delete;
            CUnifiedStepwiseOUModel(vInitialModelState vInitStates, ParamMap params, double logHastingsRatio = 0, double logJointPriorDensity= 1);
        //Static Members
        private:
            static const std::string modelName;
            static const std::vector<std::string> parameterNames;
        //Methods
        private:
        //Overridden Methods
        public:
            const std::string & getName() const override {return this->modelName;}
            const size_t getNParams() const override {return this->parameterNames.size();}
            const std::vector<std::string> & getParamNames() const override {return this->parameterNames;}
        private:
            std::unique_ptr<CModel> constructAdjacentModel(ParamMap & newParams, double logHastingsRatio) const override;
            virtual bool evalIsLess(int a, int b) override {return this->vInitStates[a] < this->vInitStates[b];}
            virtual bool evalIsEqual(int a, int b) override {return this->vInitStates[a] == this->vInitStates[b];}
            size_t initializeSimulationRootNode(const EvaluationBlock & evalBlock, SVModelStateNode& rootNode) const override;
            void sampleSimulationNode(const EvaluationBlock & evalBlock, SVModelStateNode & node, const SVModelStateNode & parent, const SVModelStateNode & root, double time, std::mt19937 & gen) const override;
    };

} //namespace model

#endif



