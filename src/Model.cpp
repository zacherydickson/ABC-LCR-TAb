#include <algorithm>
#include <array>
#include <cmath>
#include "Distributions.hpp"
#include <future>
#include <iterator>
#include "Logging.h"
#include <numeric>
#include "Model.hpp"
//#include "precalc.h"

namespace model{

    //Modifiable Static Variables
    double CModel::GoldenSearchBoundMultiple = 3;
    double CModel::GoldenSearchToleranceProportion = 0.001;
    double CModel::GradientEstimateStepProportion = 0.001;
    double CModel::MaximumRelativeError = 0.10; 
    //Constant Static Variables
    const std::string CModel::modelName = "AbstractModel";
    //Currently limited to 64 parameters: proposals for additional parameters will never be
    //proposed
    const std::vector<std::string> CModel::parameterNames = std::vector<std::string>();
    //StepwiseOU
    const std::string CStepwiseOUModel::modelName = "StepwiseOU";
    const std::vector<std::string> CStepwiseOUModel::parameterNames = {"delta","kappa","lambda","sigma","tau","muOoM"};
    //OUStepwise
    const std::string COUStepwiseModel::modelName = "OUStepwise";
    const std::vector<std::string> COUStepwiseModel::parameterNames = {"delta","kappa","lambda","sigma","upsilon","muOoM"};

    ModelType str2ModelType(std::string str){
        if(str == "StepwiseOU"){
            return StepwiseOU;
        } else if(str == "OUStepwise"){
            return OUStepwise;
        } else {
            throw std::invalid_argument("Attempt to enumerate unrecognized model type");
        }
    }

    //#### CON-/DESTRUCTION ######################
    
    //CModel
    CModel::CModel(vInitialModelState vinitstates, ParamMap params, double loghastingsratio, double logJointPriorDensity):
        vInitStates(vinitstates), parameters(params), logHastingsRatio(loghastingsratio),
        nlogJointProbability(-1.0), logJointPriorDensity(logJointPriorDensity),
        bFixed(true)
    {
        for(const auto & pair: params){
            if(pair.second.lowerBound != pair.second.upperBound){
                this->bFixed = false;
                break;
            }
        }
    }

    //CStepwiseOUModel
    CStepwiseOUModel::CStepwiseOUModel(vInitialModelState vInitStates, ParamMap params, double logHastingsRatio, double logJointPriorDensity):
        CModel(vInitStates,params,logHastingsRatio,logJointPriorDensity)
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
        this->determineEvaluationBlocks();
    }

    //COUStepwiseModel
    COUStepwiseModel::COUStepwiseModel(vInitialModelState vInitStates, ParamMap params, double logHastingsRatio, double logJointPriorDensity):
        CModel(vInitStates,params,logHastingsRatio,logJointPriorDensity)
    {
        //Validate params
        for(const std::string & name : COUStepwiseModel::parameterNames){
            if(this->parameters.count(name) == 0){
                throw std::invalid_argument("Attempt to construct COUStepwiseModel without all required parameters");
            }
        }
        if(this->parameters.size() > COUStepwiseModel::parameterNames.size()){
            throw std::invalid_argument("Attempt to construct COUStepwiseModel with unrecognized parameters");
        }
        this->determineEvaluationBlocks();
    }

//#### STATIC METHODS ######################

    //CModel -- public
    
    void CModel::TuneEvaluation(double relError){
        if(relError < 0){
            throw std::invalid_argument("Attempt to tune CModel evaluation with a negative relative error");
        }
        CModel::MaximumRelativeError = relError;
    }

    void CModel::TuneGoldenSearch(double mult, double tol){
        if(mult <= 0){
            throw std::invalid_argument("Attempt to tune CModel golden search with boundary multiple of zero");
        }
        if(tol <= 0){
            throw std::invalid_argument("Attempt to tune CModel golden search with non-positive tolerance");
        }
        CModel::GoldenSearchBoundMultiple = mult;
        CModel::GoldenSearchToleranceProportion = tol;
    }

    void CModel::TuneGradientEstimation(double prop){
        if(prop <= 0){
            throw std::invalid_argument("Attempt to tune CModel gradient estimation with non-positive step proportion");
        }
        CModel::GradientEstimateStepProportion = prop;
    }

    //CModel -- protected
    
    double CModel::ProposeParameterSet(ParamMap & parameters,const ProposalScaleMap & scaleMap, std::mt19937 & gen) {
        double logPropRatio = 0;
        for(auto & pair : scaleMap){
            SParameterSpecification & param = parameters[pair.first];
            logPropRatio += CModel::SampleProposalDistribution(param,pair.second,gen);
        }
        return logPropRatio;
    }

    double CModel::SampleProposalDistribution(SParameterSpecification & param, double scale, std::mt19937 & gen){
        double a = param.lowerBound;
        double b = param.upperBound;
        if(a > b){
            //fprintf(stderr, "value: %0.05f, a: %0.05f, b: %0.05f\n",param.value,a,b);
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
//        //fprintf(stderr,"X:%0.05f scale: %0.05f a:%0.05f b:%0.05f zeta:%0.05f p: %0.05f q: %0.05f logP:%0.05f\n",X,scale, a, b, zeta,p, q,logProposalCorrection);
        //fprintf(stderr,"There p: %0.05f,q:%0.05f,scale:%0.05f,v:%0.05f\n",p,q,scale,param.value);
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

    std::unique_ptr<CModel> CModel::proposeJump(const ProposalScaleMap & scaleMap, std::mt19937 & gen) const {
        //Copy the current Parameter Set
        ParamMap newParams = this->parameters;
        //Update the values of the Copy
        double logHastingsRatio = CModel::ProposeParameterSet(newParams,scaleMap,gen);
        //Return the appropriate Object
        return this->constructAdjacentModel(newParams,logHastingsRatio);
    }

    void CModel::evaluate(const Tree & tree, const StateMap & obs, ctpl::thread_pool & threadPool,std::mt19937 & gen, size_t nSim){
        using namespace std::placeholders;
        size_t nProt = obs.begin()->second.size();
        size_t nTips = tree.getNLeaves();
        size_t nThreads = threadPool.size();
        //Determine the maximum root to tip divergence between root values and observed tips
        this->maxR2TDivA = 0;
        this->maxR2TDivL = 0;
        for(int prot = 0; prot < vInitStates.size(); prot++){
            int rootL = this->vInitStates[prot].length.getMode();
            int rootA = this->vInitStates[prot].abundance.getMode();
            for(const auto & pair : obs){
                double divL = std::abs(pair.second[prot].length.getMode() - rootL)/double(rootL);
                double divA = std::abs(pair.second[prot].abundance.getMode() - rootA)/double(rootA);
                //double div = std::max(divL,divA);
               //fprintf(stderr,"div: %0.04f [%d,%d]\n",div,pair.second[prot].length.getMode(),pair.second[prot].abundance.getMode());
                if(divL > this->maxR2TDivL){
                    this->maxR2TDivL = divL;
                }
                if(divA > this->maxR2TDivA){
                    this->maxR2TDivL = divA;
                }
            }
        }
        logger::Log("Evaluating model with %ld threads ...",logger::DEBUG,nThreads);
        double nLogP = 0;
        if(nThreads > 1){
            std::vector<std::vector<int>> vThreadBlocks = this->partitionEvalBlocks(nThreads);
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

    double CModel::getMaxVectorDist(const GradientMap & gradient) const {
        double maxMultiple = CModel::GoldenSearchBoundMultiple;

        //Determine which parameter is the limiter
        double intervalLength = std::sqrt(double(this->parameters.size() * pow(maxMultiple,2.0)));
        for(const auto & pair : this->parameters){
            //Maximum distance is dist to the bound in the opposite gradient direction
            //Default assume positive slope so head towards lower bound
            double dist = pair.second.value - pair.second.lowerBound;
            if(gradient.at(pair.first) < 0){ // If a negative slope go towards upper bound
                dist = pair.second.upperBound - pair.second.value;
            }
            //If the bound is infinite, and the value is non-zero, set distance to be
            // a multiple of the current value;
            if(std::isinf(dist) && std::abs(pair.second.value) > 0){
                dist = std::abs(pair.second.value) * maxMultiple;
            } else if(pair.second.value == 0){
                dist = maxMultiple; //there is no info on what a big step is, just go a few
            }
            //fprintf(stderr,"%s) maxDist: %0.04f\n", pair.first.c_str(),dist);
            //The maximum interval
            double intervalDist = dist / std::abs(gradient.at(pair.first));
            if(intervalDist < intervalLength){
                intervalLength = intervalDist;
            }
        }
        return intervalLength;
    }

    double CModel::getMinEval() const{
        if(this->nlogJointProbability < 0){
            throw std::logic_error("Attempt to get minimum log Probability of unevaluated model");
        }
        return this->minEvalValue;
    }

    

    std::unique_ptr<CModel> CModel::goldenSearch(const Tree & tree, const StateMap & obs, ctpl::thread_pool & threadPool,std::mt19937 & gen, size_t nSim) const {
        double phi = 1.618033988749895;
        //Establish the gradient (Note: A Regularized Gradient with a magnitude of 1 is
        //used)
        GradientMap gradient = this->estimateGradient(tree,obs,threadPool,gen,nSim);
        if(gradient.empty()){ //The overall slope is zero and there is no changing to be done
            return std::unique_ptr<CModel>(nullptr);
        }
        //for(const auto & pair : gradient){
        //fprintf(stderr,"%s:%0.4f\t",pair.first.c_str(),pair.second);
        //}
        //fprintf(stderr,"\n");
        //Calculate the magnitude of the gradient vector
         
        double intervalLength = this->getMaxVectorDist(gradient); 
        //fprintf(stderr,"intervalLength: %0.04f\n", intervalLength);
        //intervalLength is now the maximum distance to travel in the gradient direction
        
        //f1 is this object
        //f3 is the end point
        //The ratio of f1->f2 / f2->f3 is the golden ratio b/a = φ; a+ b = (f1->f3)
        //  a + b = c; b/a = φ  -> b = aφ; a + aφ = c; a = c / (1 + φ)
        //-1.0 is used to indicate an undefined value as both the x value and the y value
        //are strictly positive
        std::pair<double,double> aProbePoint[4] = {std::make_pair(0.0,this->getNLogP()), //f1
            std::make_pair(intervalLength / (1.0 + phi),-1.0), //f2
            std::make_pair(intervalLength, -1.0) /*f3*/, std::make_pair(-1.0, -1.0)}; //f4

        bool bDone = false;
        double tolerance = intervalLength * CModel::GoldenSearchToleranceProportion;
        while(!bDone) {
            //Define the 4th probe point's x value and mark that it needs to be evaluated
            aProbePoint[3] = std::make_pair(aProbePoint[0].first + (aProbePoint[2].first - aProbePoint[1].first),-1);
            //Just to make things cleaner later we'll make sure that the 4th probe point is
            //always has a higher x than the 2nd probe point;
            if(aProbePoint[3].first < aProbePoint[1].first){
                std::swap(aProbePoint[3],aProbePoint[1]);
            }
            if(aProbePoint[3].first - aProbePoint[0].first < tolerance){
                bDone = true;
            }
            //fprintf(stderr,"Δx = %0.04f\n",aProbePoint[3].first - aProbePoint[0].first);
            //Evaluate Models at probe points as necessary
            for(int i = 0; i < 4; i++){
                if(aProbePoint[i].second >= 0){
                    continue;
                }
                ParamMap newParam = this->parameters;
                    //fprintf(stderr,"[x:%0.04f] ",aProbePoint[i].first);
                for(auto & pair : newParam){
                    pair.second.value -= aProbePoint[i].first * gradient.at(pair.first);
                    if(pair.second.value < pair.second.lowerBound){
                        pair.second.value = pair.second.lowerBound;
                    } else if(pair.second.value > pair.second.upperBound){
                        pair.second.value = pair.second.upperBound;
                    }
                    //fprintf(stderr,"%s) %0.04f -> %0.04f\t",pair.first.c_str(),gradient.at(pair.first),pair.second.value);
                }
                //fprintf(stderr,"\n");
                std::unique_ptr<CModel> pointModel = this->constructAdjacentModel(newParam,0);
                pointModel->evaluate(tree,obs,threadPool,gen,nSim);
                aProbePoint[i].second = pointModel->getNLogP();
            }
            //If the 4th probe point is less than the second, search the right interval
            if(aProbePoint[3].second < aProbePoint[1].second){
                aProbePoint[0] = aProbePoint[1];
                aProbePoint[1] = aProbePoint[3];
            } else { //Search the left interval
                aProbePoint[2] = aProbePoint[3];
            }
        } 

        //Construct the parameters at the best Step
        ParamMap newParam = this->parameters;
        for(auto & pair : newParam){
            pair.second.value -= aProbePoint[1].first * gradient.at(pair.first);
            if(pair.second.value < pair.second.lowerBound){
                pair.second.value = pair.second.lowerBound;
            } else if(pair.second.value > pair.second.upperBound){
                pair.second.value = pair.second.upperBound;
            }
        }
        //We are going to say that the probability of proposing this new model is the same
        //as the parent model, as that choice was random, but from there it is
        //deterministic
        return this->constructAdjacentModel(newParam,this->logHastingsRatio);
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
    }


    //CModel -- concrete, protected
    
    void CModel::determineEvaluationBlocks(){
        using namespace std::placeholders;
        std::vector<int> vIdxs(this->vInitStates.size());
        std::iota(vIdxs.begin(),vIdxs.end(),0);
        //Put the indexes in ascending order of their length
        std::sort(vIdxs.begin(),vIdxs.end(),std::bind(&CModel::evalIsLess,this,_1,_2));
//                [&](int a, int b) -> bool {
//                    return this->vInitStates[a].abundance < this->vInitStates[b].abundance;
//                });
        EvaluationBlock curBlock;
        curBlock.push_back(std::ref(this->vInitStates[vIdxs[0]]));
        for(int i = 0; i < vIdxs.size()-1; i++){
            auto it1 = vIdxs.begin()+i;
            auto it2 = vIdxs.begin()+(i+1);
            //Check if prot[idx[i]]'s length is equal to the subsequent state's length
            if(this->evalIsEqual(*it1,*it2)){
            //if(this->vInitStates[*it1].abundance == this->vInitStates[*it2].abundance){
                //Add the next protein to the current block
                curBlock.push_back(std::ref(this->vInitStates[*it2]));
                vIdxs.erase(it2); //Stop considering this protein
                i--; //Take  step back so that the next iteration compares the current state again
            } else { //The current Block is complete, start a new one
                this->vEvalBlocks.push_back(curBlock);
                curBlock.clear();
                //Start the next block with the next protein
                curBlock.push_back(std::ref(this->vInitStates[*it2]));
            }
        }
        this->vEvalBlocks.push_back(curBlock);
        //Sort the evaluation blocks in descending order by size
        std::sort(this->vEvalBlocks.begin(),this->vEvalBlocks.end(),
            [&](const EvaluationBlock & a, const EvaluationBlock & b) -> bool {
                return b.size() < a.size();
            });
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

    std::vector<std::vector<int>> CModel::partitionEvalBlocks(size_t nThreads) const {
        logger::Log("Paritioning Evaluation Blocks ...",logger::DEBUG+1);
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
        logger::Log("Partitioned into %ld threadBlocks with %ld to %ld proteins each",logger::DEBUG+1,nThreads,min, max);
        return vvEvalBlockIdxs;
    }

    GradientMap CModel::estimateGradient(const Tree & tree, const StateMap & obs, ctpl::thread_pool & threadPool, std::mt19937 & gen, size_t nSim) const {
        double hProp = 0.001; //The step size for gradient estimation; interpret as a proportion of parameter value 
        GradientMap gradient;
        double gradientMagnitude = 0;
        for(const auto & pair : this->parameters){ //For each parameter take small step
            ParamMap newParams = this->parameters;
            //With fixed Priors, ww cannot change the value and the slope is zero
            if(pair.second.lowerBound == pair.second.upperBound){
                gradient[pair.first] = 0.0;
                continue;
            }
            //To make sure the new value always remains within domain, the step size
            //the non-zero minimum of 1.0 (adding hProp), the distance to the nearest
            //domain boundary, the value itself, and the range of the domain
            //be away from the nearest finite bound, or a proportion of the value if both
            //bounds are infinite
            double lowerDist = pair.second.value - pair.second.lowerBound;
            double upperDist = pair.second.upperBound - pair.second.value;
            int sign = 1.0;
            double dist = lowerDist;
            if(lowerDist > upperDist){
                sign = -1.0;
                dist = upperDist;
            }
            std::vector<double> vDist = {dist, pair.second.value, pair.second.upperBound - pair.second.lowerBound};
            double step = 1.0; 
            for(double d : vDist){
                if(d > 0.0 & d < step){
                    step = d;
                }
            }
            //Set the new value for the parameter
            newParams[pair.first].value = pair.second.value + sign * hProp * step;
            //fprintf(stderr,"%s) (v)%0.04f + (s)%d * (h)%0.04f *(t)%0.04f = %0.04f\n",pair.first.c_str(), pair.second.value, sign, hProp, step,newParams[pair.first].value);
            //Generate the model
            std::unique_ptr<CModel> stepModel = this->constructAdjacentModel(newParams,0);
            //Evaluate the model
            stepModel->evaluate(tree,obs,threadPool,gen,nSim);
            //Calculate the slope
            double slope = (stepModel->getNLogP() - this->getNLogP()) / (newParams[pair.first].value - pair.second.value);
            //fprintf(stderr,"%s) Ls:%0.04f Lm:%0.04f, xs:%0.04f, xm:%0.04f\n",pair.first.c_str(),stepModel->getNLogP(),this->getNLogP(),newParams[pair.first].value,pair.second.value);
            //Store the slope
            gradient[pair.first] = slope;
            gradientMagnitude += std::pow(slope,2.0);
        }
        gradientMagnitude = std::sqrt(gradientMagnitude);
        if(gradientMagnitude == 0){ //If all slopes are zero, then the gradient is zero
            //Return an empty map to demonstrate this
            return GradientMap();
        }
        //Rescale the gradient vector so the direction remains the same, but each value now
        //corresponds to the distance traveled in each basis vector for a single step in
        //the gradient direction
        for(auto & pair : gradient){
            pair.second /= gradientMagnitude;
        }
        return gradient;
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
                //fprintf(stderr,"Node %d\n",nodeIdx);
                SVModelStateNode & node = vNodes.atIndex(nodeIdx);
                //Reset the node for each simulation
                node.value.vLength.clear();
                node.value.vAbundance.clear();
                int parentIdx = node.parent;
                SVModelStateNode & parent = vNodes.atIndex(parentIdx);
                double time = node.height;
                this->sampleSimulationNode(evalBlock,node,parent,rootNode,time,gen);
            }
            //Iterate over tips to update counts
            int tipIdx = 0;
            for(const auto & tipPair : obs){
                //First is the label of a tip
                //Second is a vector of all proteins
                auto & node = vNodes.atLabel(tipPair.first);
                int unitIdx = 0;
                for(const InitialModelState & iState: evalBlock){
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

    //CModel -- Virtual, public

    //CStepwiseOUModel -- unique, private

    //CStepwiseOUModel -- overridden, public
    
    //CStepwiseOUModel -- overridden, private

    std::unique_ptr<CModel> CStepwiseOUModel::constructAdjacentModel(ParamMap & newParams, double logHastingsRatio) const {
        return std::unique_ptr<CModel>(new CStepwiseOUModel(this->vInitStates,newParams,logHastingsRatio));
    }

    size_t CStepwiseOUModel::initializeSimulationRootNode(const EvaluationBlock & evalBlock, SVModelStateNode & rootNode) const {
        //Object to keep track of the counts across tips, proteins, and data types
        size_t nProt = 0;
        int length = evalBlock[0].get().length.getMode();
        for(const auto & state : evalBlock){
            int abundance = state.get().abundance.getMode();
            //Every protein has the same length at the root
            //  To have consistent polymorphism, we'll waste the memory and repeat it for
            //  each evalBlock
            rootNode.value.vLength.push_back(length);
            rootNode.value.vAbundance.push_back(abundance);
            nProt += state.get().vProtIdxs.size();
        }
        return nProt;
    }

    void CStepwiseOUModel::sampleSimulationNode(const EvaluationBlock & evalBlock, SVModelStateNode & node, const SVModelStateNode & parent, const SVModelStateNode & root, double time, std::mt19937 & gen) const{
        //Sample the Length Based on the branch length and the parent value
        double lambdaTerm = (parent.value.vLength[0]+1)*this->parameters.at("lambda").value*time;
        double kappaTerm = (parent.value.vLength[0])*this->parameters.at("kappa").value*time;
        int nIns = stats::PoissonQuantile(stats::generate_open_canonical(gen),lambdaTerm);
        int nDel = stats::PoissonQuantile(stats::generate_open_canonical(gen),kappaTerm);
        int length = parent.value.vLength[0] + nIns - nDel;
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
        //Sample each abundance
        for(int protIdx = 0; protIdx < evalBlock.size(); protIdx++){
            //Every protein has the same length at the root
            //  To have consistent polymorphism, we'll waste the memory and repeat it for
            //  each evalBlock
            node.value.vLength.push_back(length);
            double selCoef = this->parameters.at("delta").value*time;
            double pTerm = parent.value.vAbundance[protIdx]*selCoef;
            int lDiff = (node.value.vLength[protIdx] - root.value.vLength[protIdx]);
            double lTerm = this->parameters.at("tau").value * lDiff;
            double rTerm = root.value.vAbundance[protIdx]*(1-selCoef);
            double meanAb = pTerm + rTerm * std::exp(lTerm);
            double driftCoef = this->parameters.at("sigma").value*time;
            int abundance = stats::DiscreteTruncatedNormalQuantile(
                        stats::generate_open_canonical(gen),
                        meanAb, driftCoef,0.0,std::numeric_limits<double>::infinity());
            node.value.vAbundance.push_back(abundance);
        }
    }

    
    //COUStepwiseModel -- overridden, private
    
    std::unique_ptr<CModel> COUStepwiseModel::constructAdjacentModel(ParamMap & newParams, double logHastingsRatio) const{
        return std::unique_ptr<CModel>(new COUStepwiseModel(this->vInitStates,newParams,logHastingsRatio));
    }

    size_t COUStepwiseModel::initializeSimulationRootNode(const EvaluationBlock & evalBlock, SVModelStateNode& rootNode) const{
        //Object to keep track of the counts across tips, proteins, and data types
        size_t nProt = 0;
        int abundance = evalBlock[0].get().abundance.getMode();
        for(const auto & state : evalBlock){
            int length = state.get().length.getMode();
            rootNode.value.vLength.push_back(length);
            //Every protein has the same abundance at the root
            //  To have consistent polymorphism, we'll waste the memory and repeat it for
            //  each evalBlock
            rootNode.value.vAbundance.push_back(abundance);
            nProt += state.get().vProtIdxs.size();
        }
        return nProt;
    }

    void COUStepwiseModel::sampleSimulationNode(const EvaluationBlock & evalBlock, SVModelStateNode & node, const SVModelStateNode & parent, const SVModelStateNode & root, double time, std::mt19937 & gen) const{
        //Sample the abundance
        double selCoef = this->parameters.at("delta").value*time;
        double pTerm = parent.value.vAbundance[0]*selCoef;
        double rTerm = root.value.vAbundance[0]*(1-selCoef);
        double meanAb = pTerm + rTerm;
        double driftCoef = this->parameters.at("sigma").value*time;
        int abundance = stats::DiscreteTruncatedNormalQuantile(
                    stats::generate_open_canonical(gen),
                    meanAb, driftCoef,0.0,std::numeric_limits<double>::infinity());
        for(int protIdx = 0; protIdx < evalBlock.size(); protIdx++){
            node.value.vAbundance.push_back(abundance);
            //Sample the Length Based on the branch length and the parent value
            double fcRoot = (1 + node.value.vAbundance[protIdx]) / (1+ root.value.vAbundance[protIdx]);
            double upsilonTerm = std::exp(fcRoot * this->parameters.at("upsilon").value);
            double lambdaTerm = (parent.value.vLength[protIdx]+1)*this->parameters.at("lambda").value*time;
            lambdaTerm *= upsilonTerm;
            double kappaTerm = (parent.value.vLength[protIdx])*this->parameters.at("kappa").value*time;
            kappaTerm *= kappaTerm;
            int nIns = stats::PoissonQuantile(stats::generate_open_canonical(gen),lambdaTerm);
            int nDel = stats::PoissonQuantile(stats::generate_open_canonical(gen),kappaTerm);
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
            node.value.vLength.push_back(length);
        }
        
    }

}
