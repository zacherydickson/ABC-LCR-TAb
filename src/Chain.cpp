#include "Chain.hpp"
#include "Distributions.hpp"
#include <stdexcept>
#include <sstream>

namespace chain {

//#### Static Members #####

    size_t CChain::AcceptHorizon = 10;
    bool CChain::BGradientDescent = true;
    double CChain::InitialProposalScale = 1.0;
    double CChain::MaximumProposalScaleOoM = 9.0;
    size_t CChain::SimulationVarianceEstimationN = 10;
    double CChain::SimulationAlpha = 0.05;
    //See Mylène Bédard (2008) 10.1016/j.spa.2007.12.005
    //  for aguments why we shouldn't pick this and the real reference for it
    double CChain::TargetAcceptRate = 0.243;

//#### CON/-Destruction ###

    CChain::CChain(int id, const prior::CPrior & prior, const Tree & tree, const model::StateMap & obs, unsigned long int seed, size_t nThreads, size_t nSim):
        id(id), prior(std::cref(prior)), tree(std::cref(tree)),
        obs(std::cref(obs)), gen(), threadPool(nThreads),
        nSim(nSim)
    {

        logger::Log("Chain %d) seed %ld\n",logger::DEBUG,id,seed);
        gen.seed(seed);
        this->model = prior.GenerateModel(this->gen);
        for(const std::string & name : this->model->getParamNames()){
            adaptiveScaleMap[name] = new aparam::CAdaptiveParameter(CChain::TargetAcceptRate,CChain::AcceptHorizon,CChain::TargetAcceptRate,pow(10.0,-CChain::MaximumProposalScaleOoM),pow(10.0,+CChain::MaximumProposalScaleOoM));
        }
        if(id == 0){
            std::stringstream stream;
            stream << *(this->model) << "\n";
            logger::Log("Initial Model\n%s",logger::DEBUG,stream.str().c_str());
        }
        std::vector<double> vEvaluationResults;
        for(int i = 0; i < CChain::SimulationVarianceEstimationN; i++){
            this->doEvaluation(this->model.get());
            vEvaluationResults.push_back(model->getNLogP());
        }
        this->lastEval = stats::GetMean(vEvaluationResults);
        this->evaluationSD = stats::GetStandardDeviation(vEvaluationResults,this->lastEval);
        logger::Log("Chain %d) Estimated std evaluation error: %0.04f",logger::INFO,id,this->evaluationSD);
    }

    CChain::CChain(int id, const prior::CPrior & prior, const Tree & tree, const model::StateMap & obs, unsigned long int seed, size_t nThreads, size_t nSim, std::string modelStr):
        id(id), prior(std::cref(prior)), tree(std::cref(tree)),
        obs(std::cref(obs)), gen(), threadPool(nThreads),
        nSim(nSim)
    {
        gen.seed(seed);
        logger::Log("Chain %d) seed %ld\n",logger::DEBUG,id,seed);
        this->model = prior.GenerateModel(this->gen);
        for(const std::string & name : this->model->getParamNames()){
            adaptiveScaleMap[name] = new aparam::CAdaptiveParameter(CChain::TargetAcceptRate,CChain::AcceptHorizon,CChain::TargetAcceptRate,pow(10.0,-CChain::MaximumProposalScaleOoM),pow(10.0,+CChain::MaximumProposalScaleOoM));
        }
        this->model->setToStr(modelStr);
        if(id == 0){
            std::stringstream stream;
            stream << *(this->model) << "\n";
            logger::Log("Initial Model\n%s",logger::DEBUG,stream.str().c_str());
        }
        std::vector<double> vEvaluationResults;
        for(int i = 0; i < CChain::SimulationVarianceEstimationN; i++){
            this->doEvaluation(this->model.get());
            vEvaluationResults.push_back(model->getNLogP());
        }
        this->lastEval = stats::GetMean(vEvaluationResults);
        this->evaluationSD = stats::GetStandardDeviation(vEvaluationResults,this->lastEval);
        logger::Log("Chain %d) Estimated std evaluation error: %0.04f",logger::INFO,id,this->evaluationSD);
    }

    CChain::~CChain(){
        for(auto & pair : this->adaptiveScaleMap){
            delete pair.second;
            pair.second = nullptr;
        }
    }

//#### Static Methods #####

    //CChain, public

    void CChain::TunePS(size_t horizon, double init, double oom, double rate){
        if(horizon < 1){
            throw std::invalid_argument("Attempt to tune CChain proposal scaling with a non-Natural horizon");
        }
        double initOoM = std::abs(std::log10(init));
        if(initOoM > oom){
            throw std::invalid_argument("Attempt to tune CChain proposal scaling with an initial scale outside of the provided order of magnitude bounds");
        }
        if(oom <= 0.0){
            throw std::invalid_argument("Attempt to tune CChain proposal scaling with non-positive order of magnitude boundaries");
        }
        if(rate <= 0.0 || rate > 1.0){
            throw std::invalid_argument("Attempt to tune CChain proposal scaling with a target acceptance rate outside the bounds of (0,1]");
        }
        CChain::AcceptHorizon = horizon;
        CChain::InitialProposalScale = init;
        CChain::MaximumProposalScaleOoM = oom;
        CChain::TargetAcceptRate = rate;
    }

    void CChain::TuneSimVar(double alpha, size_t n){
        if(alpha <= 0 || alpha > 1.0){
            throw std::invalid_argument("Attempt to tune CChain simulation variance handling with an alpha paramter outside of (0,1]");
        }
        if(n < 2){
            throw std::invalid_argument("Attempt to tune CChain simulation variance handling with insufficient runs to estimate variance");
        }
        CChain::SimulationAlpha = alpha;
        CChain::SimulationVarianceEstimationN = n;
    }

//#### Methods ############

    double CChain::calcAcceptanceRatio(const model::CModel & m2, double temperature){
        logger::Log("Chain %d) Calculating Acceptance Probability ...",logger::DEBUG,this->id);
        if(m2.logJointPriorDensity > 0){
            throw std::invalid_argument("Attempt to calculate acceptance probability for a model which has not had its joint prior density calculated");
        }
        double priorCorrection = m2.logJointPriorDensity - this->model->logJointPriorDensity;
        double sizeCorrection = 0;//m2.getMinEval() - this->model->getMinEval();
        double accept = 0;
        //Shi & Rabosky 2015 (BAMM)
        double logLikelihoodRatio = this->model->getNLogP() - m2.getNLogP();
        //The probability we would see an absolute difference of this size if the true difference was
        //zero
        //double noiseProbability = std::erf(std::abs(logLikelihoodRatio)/2/this->evaluationSD);
        //We shrink the logRatio by the the noise probability
        //logLikelihoodRatio *= (noiseProbability);
        //
        //Consider proposals generously, ie. that the model is inflated and the proposal
        // is deflated, so the difference would be larger
        logLikelihoodRatio += stats::NormalQuantile(1-CChain::SimulationAlpha/2.0) * std::sqrt(2) * this->evaluationSD;
        double logAccept = logLikelihoodRatio + priorCorrection + sizeCorrection;
        logAccept *= temperature;
        if(std::isinf(m2.getNLogP())){
            accept = 0.0;
        } else {
            logAccept += m2.logHastingsRatio;
            if(logAccept < 0.0){
                accept = std::exp(logAccept);
            } else {
                accept = 1.0;
            }
        }
        logger::Log("Chain %d) α = exp(%0.04f(L) + %0.04f(P) + %0.04f(S))^(%0.04f(T))*exp(%0.04f(H)) = %0.02f%%",
                logger::DEBUG,this->id,logLikelihoodRatio,priorCorrection,sizeCorrection,temperature,m2.logHastingsRatio,accept * 100);
        return accept;
    }

    void CChain::doEvaluation(model::CModel * model){
        model->evaluate(this->tree,this->obs,this->threadPool,this->gen,this->nSim);
        logger::Log("Chain %d) evaluated at %0.04f |OoM|",logger::INFO,this->id,model->getNLogP());
    }

    bool CChain::iterate(int threadId, std::unordered_set<std::string> paramNameSet, double temperature){
        model::ProposalScaleMap scaleMap;
        logger::Log("Chain %d) proposing from %0.04f |OoM| with %d sims",logger::INFO,this->id,this->model->getNLogP(),this->nSim);
        //Pull out the relevant scales for proposals
        for(const std::string & name : paramNameSet){
            try{
                scaleMap[name] = this->adaptiveScaleMap.at(name)->getValue();
            } catch(std::out_of_range & e){
                logger::Log("Chain %d) Attempt to iterate chain while proposing an unknown parameter (%s)",
                        logger::ERROR,this->id,name.c_str());
            }
        }
        //Make and evaluate the proposal
        std::unique_ptr<model::CModel> proposal = this->model->proposeJump(scaleMap,this->gen);
        bool bAccepted = false;
        int bOnce= (CChain::BGradientDescent) ? 0 : 1;
        do{
            proposal->logJointPriorDensity = this->prior.get().calculateJointPriorDensity(proposal->getParamMap());
            if(id == 0){
                std::stringstream stream;
                stream << *proposal << "\n";
                logger::Log("Proposed Model\n%s",logger::DEBUG,stream.str().c_str());
            }
            this->doEvaluation(proposal.get());
            //Try accepting the proposal
            double accept = this->calcAcceptanceRatio(*proposal,temperature);
            std::string result = "Rejected";
            if(accept == 1.0 || std::generate_canonical<double,10>(this->gen) < accept){
                this->lastEval = this->model->getNLogP();
                this->model = std::move(proposal);
                bAccepted = true;
                result = "Accepted";
            }
            logger::Log("Chain %d) %s proposal with probability %0.2f%%",logger::INFO,this->id,result.c_str(),accept*100.0);
            if(bAccepted || bOnce){ //Update based on the whole assessment
                //Update the proposal scales as necessary
                for(const std::string & name : paramNameSet){
                    double scaleUpdate = (bAccepted) ? 1.0 : 0.0;
                    if(this->adaptiveScaleMap[name]->update(scaleUpdate,this->gen)){
                        logger::Log("Chain %d) scale for %s proposals updated to %0.04f natural OoMs",logger::INFO,this->id,name.c_str(),std::log(this->adaptiveScaleMap[name]->getValue()));
                    }
                }
            }
            if(!bAccepted && !bOnce){
                logger::Log("Chain %d) Attempting Gradient descent",logger::INFO,this->id);
                size_t quickNSim = (this->nSim > 4) ? this->nSim / 4 : 1; 
                std::unique_ptr<model::CModel> gradientProposal = proposal->goldenSearch(this->tree,this->obs,this->threadPool,this->gen,quickNSim);
                if(gradientProposal){ //Check if a valid model was returned
                    proposal = std::move(gradientProposal);
                } else { //Otherwise Give up
                    bOnce = true;
                }
            }
        } while(!bAccepted && !bOnce++);
        return bAccepted;
    }


    void CChain::swapModel(CChain & other){
        std::swap(this->model,other.model);
    }

}



