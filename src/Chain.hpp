#ifndef CHAIN_HPP
#define CHAIN_HPP

#include "AdaptiveParameter.hpp"
#include "ctpl_stl.h"
#include <functional>
#include "Logging.h"
#include <memory>
#include <random>
#include "Prior.hpp"
#include <string>
#include "Tree.hpp"
#include <unordered_map>
#include <unordered_set>

namespace chain {

    class CChain{
        //Con-/Destuction
        public:
            CChain(int id, const prior::CPrior & prior, const Tree & tree, const model::StateMap & obs, unsigned long int seed, size_t nThreads, size_t nSim);
            CChain(int id, const prior::CPrior & prior, const Tree & tree, const model::StateMap & obs, unsigned long int seed, size_t nThreads, size_t nSim,std::string modelStr);
            ~CChain();
        //Members
        public:
            const int id;
        protected:
            std::unordered_map<std::string,aparam::CAdaptiveParameter *> adaptiveScaleMap;
            aparam::CAdaptiveParameter adaptiveSimAlpha;
            bool bFixed;
            double evaluationSD;
            std::mt19937 gen;
            double lastEval;
            size_t nSim;
            size_t iteration;
            std::unique_ptr<model::CModel> model;
            std::reference_wrapper<const model::StateMap> obs;
            std::reference_wrapper<const prior::CPrior> prior;
            ctpl::thread_pool threadPool;
            std::reference_wrapper<const Tree> tree;
        //Static Members
        protected:
            static size_t AcceptScaleHorizon;
            static size_t AcceptAlphaHorizon;
            static bool BGradientDescent;
            static double InitialProposalScale;
            static double InitialSimulationAlpha;
            static double MaximumProposalScaleOoM;
            static size_t SimulationVarianceEstimationN;
            static size_t SimulationVarianceEstimateHorizon;
            static double TargetAcceptRate;
        //Static Methods
        public:
            static bool GetGradientDescentState() {return CChain::BGradientDescent;}
            static size_t GetPSHorizon() {return CChain::AcceptScaleHorizon;}
            static double GetPSInit() {return CChain::InitialProposalScale;}
            static double GetPSOoM() {return CChain::MaximumProposalScaleOoM;}
            static double GetPSRate() {return CChain::TargetAcceptRate;}
            static double GetSimVarAlpha() {return CChain::InitialSimulationAlpha;}
            static size_t GetSimVarAlphaHorizon() {return CChain::AcceptAlphaHorizon;}
            static size_t GetSimVarN() {return CChain::SimulationVarianceEstimationN;}
            static size_t GetSimVarReEvalHorizon() {return CChain::SimulationVarianceEstimateHorizon;}
            static void ToggleGradientDescent() {CChain::BGradientDescent = !CChain::BGradientDescent;}
            static void TunePS(size_t horizon, double init, double oom, double rate);
            static void TuneSimVar(double alpha, size_t alphaHorizon, size_t n, size_t reEvalHorizon);
        //Methods
        public:
            const model::CModel & getModel() const {return *(this->model);}
            bool iterate(int threadId, std::unordered_set<std::string> paramNameSet, double temperature);
            double getEvaluationSD() const {return this->evaluationSD;}
            double getLastEval() const {return this->lastEval;}
            void swapModel(CChain & other);
        protected:
            double calcAcceptanceRatio(const model::CModel & m2, double temperature);
            void doEvaluation(model::CModel * model,bool bQuiet = false);
            double estimateSimulationVariance();
    };

};

#endif //CHAIN_HPP


