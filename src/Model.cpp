#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include "Distributions.hpp"
#include <future>
#include <iterator>
#include "Logging.h"
#include <numeric>
#include "Model.hpp"
#include <queue>
#include <set>
//#include "precalc.h"

namespace model{

    //Modifiable Static Variables
    double CModel::MaximumRelativeError = 0.10; 
    //Constant Static Variables
    const std::string CModel::modelName = "AbstractModel";
    //Currently limited to 64 parameters: proposals for additional parameters will never be
    //proposed
    const std::vector<std::string> CModel::parameterNames = std::vector<std::string>();
    //StepwiseOU
    const std::string CStepwiseOUModel::modelName = "StepwiseOU";
    const std::vector<std::string> CStepwiseOUModel::parameterNames = {"delta","kappa","lambda","muOoM","sigma","tau","upsilon"};

    ModelType str2ModelType(std::string str){
        if(str == "StepwiseOU"){
            return StepwiseOU;
        } else {
            throw std::invalid_argument("Attempt to enumerate unrecognized model type");
        }
    }

    //#### CON-/DESTRUCTION ######################
    
    //CModel
    CModel::CModel(vInitialModelState vinitstates, ParamMap params, double loghastingsratio, double logJointPriorDensity, const vEvaluationBlock * initEvalBlocks, const SEvalBlockIdxPartition * initPartition):
        vInitStates(vinitstates), parameters(params), logHastingsRatio(loghastingsratio),
        nlogJointProbability(-1.0), logJointPriorDensity(logJointPriorDensity),
        bFixed(true)
    {
        for(const auto & pair: params){
            if(!pair.second.isFixed()){
                this->bFixed = false;
                break;
            }
        }
        this->validateParameters();
        if(initEvalBlocks){
            this->vEvalBlocks = *initEvalBlocks;
        }
        if(initPartition){
            this->evalBlockIdxPartition = *initPartition;
        }
    }


    //CStepwiseOUModel
    CStepwiseOUModel::CStepwiseOUModel(vInitialModelState vInitStates, ParamMap params, double logHastingsRatio, double logJointPriorDensity, const vEvaluationBlock * initEvalBlocks, const SEvalBlockIdxPartition * initPartition):
        CModel(vInitStates,params,logHastingsRatio,logJointPriorDensity,initEvalBlocks,initPartition)
    {
        //Validate params
        for(const std::string & name : CStepwiseOUModel::parameterNames){
            if(this->parameters.count(name) == 0){
                throw std::invalid_argument("Attempt to construct CStepwiseOUModel without all required parameters");
            }
        }
        if(this->parameters.size() > CStepwiseOUModel::parameterNames.size()){
            throw std::invalid_argument("Attempt to construct CStepwiseOUModel with unrecognized parameters");
        }
        this->bFixedZeroTau = this->parameters.at("tau").isFixed(0);
        this->bFixedZeroUpsilon = this->parameters.at("upsilon").isFixed(0);
        if(this->vEvalBlocks.size() == 0){
            this->determineEvaluationBlocks();
        }
    }

//#### STATIC METHODS ######################

    //CModel -- public
    
    void CModel::TuneEvaluation(double relError){
        if(relError < 0){
            throw std::invalid_argument("Attempt to tune CModel evaluation with a negative relative error");
        }
        CModel::MaximumRelativeError = relError;
    }


    //CModel -- protected
    
    double CModel::ProposeParameterSet(ParamMap & parameters,const ProposalScaleMap & scaleMap, std::mt19937 & gen) {
        double logPropRatio = 0;
        for(auto & pair : scaleMap){
            SParameterSpecification & param = parameters[pair.first];
            logPropRatio += CModel::SampleProposalDistribution(param,pair.second,gen);
            //Scan through other parameters to see if any depend of the changed parameter and
            //update them as necessary
            for(auto & dependentPair : parameters){
                if(dependentPair.second.dependency == pair.first){
                    dependentPair.second.value = param.value;
                }
            }
        }
        return logPropRatio;
    }

    double CModel::SampleProposalDistribution(SParameterSpecification & param, double scale, std::mt19937 & gen){
        double a = param.lowerBound;
        double b = param.upperBound;
        if(a > b){
            throw std::invalid_argument("Proposal Distribution Bounds are inverted");
        }
        if(a == b){ //The bounds are the same, ie. the value is fixed
            return 0;
        }
        //Adjust the scale to be appropriate for the domain in question 
        //This is done by shinking to a proportion of itself equal to the amount of the normal distribution within the bounds a and b
        double zeta = stats::NormalCDF(b) - stats::NormalCDF(a);
        scale *= zeta;
        double p = stats::generate_open_canonical(gen);
        double q = stats::TruncatedNormalQuantile(p,param.value,scale,a,b);
        double logProposalCorrection = std::log(stats::TruncatedNormalPDF(param.value,q,scale,a,b));
        logProposalCorrection -= std::log(stats::TruncatedNormalPDF(q,param.value,scale,a,b));
        param.value = q;        
        return logProposalCorrection;
    }

    double CModel::ExtractValueFromModelString(const std::string & modelStr, const std::string  & keyStr, size_t pos){
        size_t keyPos = modelStr.find(keyStr,pos);
        if(keyPos + keyStr.size() >= modelStr.size()){
            throw std::invalid_argument(("Provided Model string is missing " + keyStr).c_str());
        }
        //Find the next pair separator in the string
        size_t sepPos = modelStr.find("\t",keyPos);
        //If no tab is found, go to the end of the string
        sepPos = (sepPos == std::string::npos) ? modelStr.size() + 1 : sepPos;
        size_t length = sepPos - keyPos - keyStr.size();
        return atof(modelStr.substr(keyPos+keyStr.size(),length).c_str());
    }

//#### METHODS ######################
    
    //CModel -- concrete, public

    

    void CModel::evaluate(const Tree & tree, const StateMap & obs, ctpl::thread_pool & threadPool,std::mt19937 & gen, size_t nSim){
        using namespace std::placeholders;
        size_t nProt = obs.begin()->second.size();
        size_t nTips = tree.getNLeaves();
        size_t nThreads = threadPool.size();
        logger::Log("Evaluating model with %ld threads ...",logger::DEBUG,nThreads);
        double nLogP = 0;
        if(nThreads > 1){
            //Check if the requested thread paritioning matches the current partitioning
            if(this->evalBlockIdxPartition.nThreads != nThreads){
                this->partitionEvalBlocks(nThreads);
            }
            std::vector<std::vector<int>> & vThreadBlocks = this->evalBlockIdxPartition.vThreadBlocks;
            std::vector<std::future<double>> vFutures;
            for(const std::vector<int> & threadBlock : vThreadBlocks){
                std::future<double> future = threadPool.push(std::bind(&CModel::evaluateBlocks,this,_1,_2,_3,_4,_5,_6),std::cref(threadBlock),std::cref(tree),std::cref(obs),gen(),nSim);
                vFutures.push_back(std::move(future));
            }
            for(auto & future : vFutures){
                nLogP += future.get();
            }
        } else {
            std::vector<int> vIdxs(this->vEvalBlocks.size());
            std::iota(vIdxs.begin(),vIdxs.end(),0);
            nLogP = this->evaluateBlocks(0,vIdxs,tree,obs,gen(),nSim);
        }
        this->nlogJointProbability = nLogP;
        this->minEvalValue = 2*this->vInitStates.size()*obs.size()*std::log(nSim + 2);
        logger::Log("Model evaluated at %0.4f |OoM|",logger::DEBUG,nLogP);
    }

    double CModel::getNLogP() const{
        if(this->nlogJointProbability < 0){
            throw std::logic_error("Attempt to get log Probability of unevaluated model");
        }
        return this->nlogJointProbability;
    }

    double CModel::getMinEval() const{
        if(this->nlogJointProbability < 0){
            throw std::logic_error("Attempt to get minimum log Probability of unevaluated model");
        }
        return this->minEvalValue;
    }

    std::unique_ptr<CModel> CModel::proposeJump(const ProposalScaleMap & scaleMap, std::mt19937 & gen) const {
        //Copy the current Parameter Set
        ParamMap newParams = this->parameters;
        //Update the values of the Copy
        double logHastingsRatio = CModel::ProposeParameterSet(newParams,scaleMap,gen);
        //Return the appropriate Object
        return this->constructAdjacentModel(newParams,logHastingsRatio);
    }

    std::ostream& CModel::output(std::ostream& os) const{
        for(int i = 0; i < this->vInitStates.size(); i++){
            os << "[";
            for(int protIdx : vInitStates[i].vProtIdxs){
                os << protIdx << ";";
            }
            os << "A:" << vInitStates[i].abundance.getMode() << ",L:"
                      << vInitStates[i].length.getMode() << "]\t";
        }
        os << "logPriorDensity: " << this->logJointPriorDensity;
        os << "\tlogHastingsRatio: " << this->logHastingsRatio;
        os << "\tnlogP: " << this->nlogJointProbability;
        for(const std::string & name : this->getParamNames()){
            os << "\t" << name << ": " << this->parameters.at(name).value;
        }
        return os;
    }

    void CModel::setToStr(const std::string & modelStr){
        size_t startPos = modelStr.find("logPriorDensity");
        this->logJointPriorDensity = CModel::ExtractValueFromModelString(modelStr,"logPriorDensity: ",startPos);
        this->logHastingsRatio =  CModel::ExtractValueFromModelString(modelStr,"logHastingsRatio: ",startPos);
        this->nlogJointProbability = CModel::ExtractValueFromModelString(modelStr,"nlogP: ",startPos);
        if(this->nlogJointProbability < 0){
            throw std::invalid_argument("Provided model string contains nLogP less than 0");
        }
        for(auto & pair : this->parameters){
            double value = CModel::ExtractValueFromModelString(modelStr,pair.first + ": ",startPos);
            if(value > pair.second.upperBound || value < pair.second.lowerBound){
                throw std::invalid_argument(("Provided model string contains " + pair.first + "value outside of prior defined domain").c_str());
            }
            pair.second.value = value;
        }
        this->validateParameters();
    }


    //CModel -- concrete, protected
    
    void CModel::determineEvaluationBlocks(){
        logger::Log("Determining evaluation Blocks",logger::DEBUG+1);
        using namespace std::placeholders;
        std::vector<int> vIdxs(this->vInitStates.size());
        std::iota(vIdxs.begin(),vIdxs.end(),0);
        //Put the indexes in ascending order of their length
        std::sort(vIdxs.begin(),vIdxs.end(),std::bind(&CModel::evalIsLess,this,_1,_2));
        EvaluationBlock curBlock;
        curBlock.push_back(vIdxs[0]);
        for(int i = 0; i < vIdxs.size()-1; i++){
            auto it1 = vIdxs.begin()+i;
            auto it2 = vIdxs.begin()+(i+1);
            //Check if prot[idx[i]]'s length is equal to the subsequent state's length
            if(this->evalIsEqual(*it1,*it2)){
            //if(this->vInitStates[*it1].abundance == this->vInitStates[*it2].abundance)
                //Add the next protein to the current block
                curBlock.push_back(*it2);
                vIdxs.erase(it2); //Stop considering this protein
                i--; //Take  step back so that the next iteration compares the current state again
            } else { //The current Block is complete, start a new one
                this->vEvalBlocks.push_back(curBlock);
                curBlock.clear();
                //Start the next block with the next protein
                curBlock.push_back(*it2);
            }
        }
        this->vEvalBlocks.push_back(curBlock);
        //Sort the evaluation blocks in descending order by size
        std::sort(this->vEvalBlocks.begin(),this->vEvalBlocks.end(),
            [&](const EvaluationBlock & a, const EvaluationBlock & b) -> bool {
                return b.size() < a.size();
            });
        logger::Log("Initial states grouped into %zu blocks",logger::DEBUG+1,this->vEvalBlocks.size());
    }

    double CModel::evaluateBlocks(int id, const std::vector<int> & vBlockIdxs, const Tree & tree, const StateMap & obs,size_t seed, size_t nSim) const{
        std::mt19937 gen(seed);
        logger::Log("Processing block of %ld evaluation blocks",logger::DEBUG+1,vBlockIdxs.size());
        double nLogP = 0;
        for(int i : vBlockIdxs){
            nLogP += this->evaluateBlock(this->vEvalBlocks[i],tree,obs,gen,nSim);
        }
        return nLogP;
    }

    void CModel::partitionEvalBlocks(size_t nThreads) {
        this->evalBlockIdxPartition.nThreads = nThreads;
        logger::Log("Partitioning Evaluation Blocks ...",logger::DEBUG+1);
        //Group Evaluation Blocks into a block of Blocks for each thread
        //  Uses a greedy algorithm of adding each evalBlock (in descending order of size)
        //  to the threadBlock which has the smallest total
        std::vector<std::vector<int>> vvEvalBlockIdxs(nThreads);
        std::vector<size_t> threadBlockWeight(nThreads,0);
        for(int i = 0; i < this->vEvalBlocks.size(); i++){
            //Find the smallest threadBlock (could be clever, but the number of elements should be reasonablly small
            auto it = std::min_element(threadBlockWeight.begin(),threadBlockWeight.end());
            size_t pos = it - threadBlockWeight.begin();
            //Add the ith evaluationBlock to the minimally weighted threadBlock
            vvEvalBlockIdxs[pos].push_back(i);
            threadBlockWeight[pos] += this->vEvalBlocks[i].size();
        }
        size_t min = *std::min_element(threadBlockWeight.begin(),threadBlockWeight.end());
        size_t max = *std::max_element(threadBlockWeight.begin(),threadBlockWeight.end());
        this->evalBlockIdxPartition.vThreadBlocks = vvEvalBlockIdxs;
        logger::Log("Partitioned into %ld threadBlocks with %ld to %ld proteins each",logger::DEBUG+1,nThreads,min, max);
    }

    double CModel::evaluateBlock(const EvaluationBlock & evalBlock, const Tree & tree, const StateMap & obs, std::mt19937 & gen, size_t nSim) const{
        logger::Log("Evaluating block of %ld proteins",logger::DEBUG+2,evalBlock.size());
        //Using pseudocounts every probability across data types, proteins, and tips will be a count
        //divided by nSim + 2, so we can pull that term out of the summation and get a baseline
        //We use instead maxSim so that all evaluations have a common lower bound
        double nLogP = 2*evalBlock.size()*obs.size()*std::log(nSim + 2);
        //Set up the Root state for the node vector
        SVModelStateNodeVector vNodes(tree); 
        SVModelStateNode & rootNode = vNodes.atIndex(0);
        size_t nProt = this->initializeSimulationRootNode(evalBlock,rootNode);
        //Initialize Counts across tips, proteins, and data types to 1
        std::vector<std::vector<std::array<double,2>>> vPositiveCounts(obs.size(),std::vector<std::array<double,2>>(nProt,{1.0,1.0}));
        for(int i = 0; i < nSim; i++){ //Perform each simulation
            //Iterate over the nodes after the root, simulating as one goes
            for(int nodeIdx = 1; nodeIdx < vNodes.size(); nodeIdx++){
                SVModelStateNode & node = vNodes.atIndex(nodeIdx);
                //Reset the node for each simulation
                node.value.vLength.clear();
                node.value.vAbundance.clear();
                int parentIdx = node.parent;
                SVModelStateNode & parent = vNodes.atIndex(parentIdx);
                double time = node.height;
                this->sampleSimulationNode(evalBlock,node,parent,rootNode,time,gen);
            }
            if(logger::Verbosity >= logger::DEBUG+3){
                std::stringstream stream;
                stream << vNodes;
                logger::Log("Simulation Result %d\n%s",logger::DEBUG+3,i,stream.str().c_str());
            }
            //Iterate over tips to update counts
            int tipIdx = 0;
            for(const auto & tipPair : obs){
                //First is the label of a tip
                //Second is a vector of all proteins
                auto & node = vNodes.atLabel(tipPair.first);
                int unitIdx = 0;
                for(int j : evalBlock){
                    const InitialModelState & iState = this->vInitStates[j];
                    for(const int protIdx : iState.vProtIdxs){
                        int nodeVal[2] = {node.value.vLength[unitIdx],node.value.vAbundance[unitIdx]};
                        double trueVal[2] = {double(tipPair.second[protIdx].length.getMode()),double(tipPair.second[protIdx].abundance.getMode())};
                        for(int dataIdx = 0; dataIdx < 2; dataIdx++){
                            double relError = (nodeVal[dataIdx] == trueVal[dataIdx]) ? 0.0 : 1.0;
                            if(trueVal[dataIdx] > 0 && relError != 0.0){
                                relError = std::abs(nodeVal[dataIdx]/trueVal[dataIdx] - 1.0);
                                relError /= CModel::MaximumRelativeError;
                                //With an max rel error of 0, for any non-match the error
                                //will be infinite, and therefore 1.0
                                relError = (relError > 1.0) ? 1.0 : relError;
                            }
                            vPositiveCounts[tipIdx][unitIdx][dataIdx] += (1.0 - relError);
                        }
                    }
                    unitIdx++;
                }
                tipIdx++;
            }
        }
        //Iterate over the counts to update the likelihood
        for(const auto & tipVec : vPositiveCounts){
            for(const auto & protVec : tipVec){
                for(int i = 0; i < 2; i++){
                    nLogP -= std::log(protVec[i]);
                }
            }
        }
        logger::Log("Eval block evaluated at %0.4f |OoM|",logger::DEBUG,nLogP);
        return nLogP;
    }

    size_t CModel::initializeSimulationRootNode(const EvaluationBlock & evalBlock, SVModelStateNode& rootNode) const {
        logger::Log("Initializing Sim Root Node ...",logger::DEBUG);
        //Object to keep track of the counts across tips, proteins, and data types
        size_t nProt = 0;
        for(int j : evalBlock){
            const InitialModelState & state = this->vInitStates[j];
            int length = state.length.getMode();
            int abundance = state.abundance.getMode();
            rootNode.value.vLength.push_back(length);
            rootNode.value.vAbundance.push_back(abundance);
            nProt += state.vProtIdxs.size();
        }
        logger::Log("Root Node initialized for block of %zd proteins",logger::DEBUG,nProt);
        return nProt;
    }

    void CModel::validateParameters() const{
        for(const auto & pair : this->parameters){
            std::string curParam = pair.first;
            std::set<std::string> sDependents;
            //If it has a dependency make sure it isn't circular or chained
            while(!parameters.at(curParam).dependency.empty()){
                sDependents.insert(curParam);
                //Follow the dependency line until termination, or until encountering a
                //previously encountered parameter
                curParam = parameters.at(curParam).dependency;
                if(parameters.count(curParam) == 0){
                    logger::Log("Identified unknown dependency (%s) for parameter %s",
                            logger::ERROR,curParam.c_str(),pair.first.c_str());
                    throw std::invalid_argument("Encountered unknown parameter dependency");
                }
                if(sDependents.count(curParam) > 0){
                    logger::Log("Identified circular dependency for parameter %s",
                            logger::ERROR,pair.first.c_str());
                    throw std::invalid_argument("Encountered circular parameter dependency");
                }
                if(!parameters.at(curParam).dependency.empty()){
                    logger::Log("Identified chained dependency for parameter %s",
                            logger::ERROR,pair.first.c_str());
                    throw std::invalid_argument("Chained dependencies are currently not allowed");
                }
            }
            if(pair.second.lowerBound > pair.second.upperBound){
                logger::Log("Identified inverted domain bounds for %s",
                        logger::ERROR,pair.first.c_str());
                throw std::invalid_argument("Proposal Distribution Bounds are inverted");
            }
            if(pair.second.value < pair.second.lowerBound ||
                    pair.second.value > pair.second.upperBound){
                logger::Log("Parameter value (%0.4f) is outside of prior domain for %s",
                        logger::ERROR,pair.second.value,pair.first.c_str());
                throw std::invalid_argument("Parameter value outside of Domain");
            }
        }
    }

    //CModel -- Virtual, public

    //CUnifiedStepwiseOUModel -- concrete, private
    
    int CStepwiseOUModel::sampleAbundance(int protIdx, SVModelStateNode & node, const SVModelStateNode & parent, const SVModelStateNode & root, double time, std::mt19937 & gen) const {
        logger::Log("Sampling Abundance at node %s ...",logger::DEBUG+3,node.label.c_str());
        double selCoef = std::exp(-this->parameters.at("delta").value*time);
        double pTerm = parent.value.vAbundance[protIdx]*selCoef;
        double rTerm = root.value.vAbundance[protIdx]*(1.0-selCoef);
        double lfc = (parent.value.vLength[protIdx]+1.0) / (root.value.vLength[protIdx] + 1.0);
        double tauTerm = std::pow(lfc,this->parameters.at("tau").value);
        double meanAb = pTerm + rTerm * tauTerm;
        double driftCoef = this->parameters.at("sigma").value*time;
        double foldDrift = stats::LogNormalQuantile(stats::generate_open_canonical(gen),0,driftCoef);
        double abn = std::round((meanAb + 1) * foldDrift);
        int abundance = (abn < std::numeric_limits<int>::max()) ? int(abn) : std::numeric_limits<int>::max();

        //int abundance = stats::DiscreteTruncatedNormalQuantile(
        //        stats::generate_open_canonical(gen),
        //        meanAb, driftCoef,0.0,std::numeric_limits<double>::infinity());
        logger::Log("Node %s abundance sampled at %d",logger::DEBUG+3,node.label.c_str(),abundance);
        return abundance;
    }

    int CStepwiseOUModel::sampleLength(int protIdx, SVModelStateNode & node, const SVModelStateNode & parent, const SVModelStateNode & root, double time, std::mt19937 & gen) const {
        logger::Log("Sampling Length at node %s ...",logger::DEBUG+3,node.label.c_str());
        double afc = (1.0+parent.value.vAbundance[protIdx]) / (1.0+root.value.vAbundance[protIdx]);
        double upsilonTerm = std::pow(afc,this->parameters.at("upsilon").value);
        double lambdaTerm = (parent.value.vLength[protIdx]+1)*this->parameters.at("lambda").value*time;
        lambdaTerm *= upsilonTerm;
        double kappaTerm = (parent.value.vLength[protIdx])*this->parameters.at("kappa").value*time;
        kappaTerm *= upsilonTerm;
        double p = stats::generate_open_canonical(gen);
        int nIns = stats::PoissonQuantile(p,lambdaTerm);
        if(errno == EDOM){
            logger::Log("Boost couldn't properly evaluate qPois(p=%0.04f,λ=%0.04f), continuing with %0.04f",logger::WARNING,p,lambdaTerm,nIns);
            errno = 0;
        }
        p = stats::generate_open_canonical(gen);
        int nDel = stats::PoissonQuantile(p,kappaTerm);
        if(errno == EDOM){
            logger::Log("Boost couldn't properly evaluate qPois(p=%0.04f,λ=%0.04f), continuing with %0.04f",logger::WARNING,p,kappaTerm,nDel);
            errno = 0;
        }
        int length = parent.value.vLength[protIdx] + nIns - nDel;
        length = (length > 0) ? length : 0;
        double mutRate = std::exp(this->parameters.at("muOoM").value) * time * length;
        if(mutRate > 0.0){
            int nMut = stats::PoissonQuantile(stats::generate_open_canonical(gen),mutRate);
            if(nMut > 0){
                double sum = 0;
                double max = -1;
                //Generate n+1 random exponentials with mean 1, the length of the longest chunk uf lcr 
                //after n mutations, is the proportion of the max exponentila to the sum
                for(int i = 0; i < nMut+1; i++){
                    double x = stats::ExponentialQuantile(stats::generate_open_canonical(gen),1.0);
                    sum += x;
                    if(x > max){
                        max = x;
                    }
                }
                double prop = max / sum; 
                length = int(std::round(prop * length));
            }
        }
        logger::Log("Node %s length sampled at %d",logger::DEBUG+3,node.label.c_str(),length);
        return length;
    }

    //CStepwiseOUModel -- overridden, private
    
    std::unique_ptr<CModel> CStepwiseOUModel::constructAdjacentModel(ParamMap & newParams, double logHastingsRatio) const{
        return std::unique_ptr<CModel>(new CStepwiseOUModel(this->vInitStates,newParams,logHastingsRatio,1,&(this->vEvalBlocks), &(this->evalBlockIdxPartition)));
    }

    void CStepwiseOUModel::sampleSimulationNode(const EvaluationBlock & evalBlock, SVModelStateNode & node, const SVModelStateNode & parent, const SVModelStateNode & root, double time, std::mt19937 & gen) const{
        logger::Log("Sampling at node %s ...",logger::DEBUG+3,node.label.c_str());
        //Sample for the first protein
        int abundance = this->sampleAbundance(0,node,parent,root,time,gen);
        int length = this->sampleLength(0,node,parent,root,time,gen);
        node.value.vAbundance.push_back(abundance);
        node.value.vLength.push_back(length);
        //check if we have any special cases
        for(int protIdx = 1; protIdx < evalBlock.size(); protIdx++){
            //Built in here is the assumption that all proteins in the eval block have
            //either the same length or abundance
            if(!this->bFixedZeroTau || this->bFixedZeroUpsilon){
                //Update abundance if tau isn't fixed at zero, or if upsilon is also fixed
                //at zero
                abundance = this->sampleAbundance(protIdx,node,parent,root,time,gen);
            }
            if(!this->bFixedZeroUpsilon){
                //Update abundance if upsilon isn't fixed at zero
                length = this->sampleLength(protIdx,node,parent,root,time,gen);
            }
            node.value.vAbundance.push_back(abundance);
            node.value.vLength.push_back(length);
        }
        logger::Log("Node %s sampling complete",logger::DEBUG+3,node.label.c_str());
    } //sampleSimulationNode

    bool CStepwiseOUModel::evalIsLess(int a, int b) const {
        if(!this->bFixedZeroUpsilon && !this->bFixedZeroTau){
            return false;
            //There is no usefully defined sorting if both interactions exist
        }
        if(this->bFixedZeroTau && !this->bFixedZeroUpsilon){
            //Abundance is independant of length, sort by abundance
            if(this->vInitStates[a].abundance < this->vInitStates[b].abundance){
                return true;
            }
            return false;
        }
        //In the remaining two cases, both are independent and length is independent of
        //abundance, we can sort by length; In the former case that is the blocking most
        //likely to group more proteins
        if(this->vInitStates[a].length < this->vInitStates[b].length){
            return true;
        }
        return false;
    }

    bool CStepwiseOUModel::evalIsEqual(int a, int b) const {
        if(!this->bFixedZeroUpsilon && !this->bFixedZeroTau){
            //It is assumed that all initial states provided during construction were
            //unique
            return false;
        }
        if(this->bFixedZeroTau && !this->bFixedZeroUpsilon){
            //Abundance is independant of length test if the abundances are equal
            if(this->vInitStates[a].abundance == this->vInitStates[b].abundance){
                return true;
            }
            return false;
        }
        //In the remaining two cases, both are independent and length is independent of
        //abundance, we can sort by length; In the former case that is the blocking most
        //likely to group more proteins
        if(this->vInitStates[a].length == this->vInitStates[b].length){
            return true;
        }
        return false;
    }

} // model namespace
