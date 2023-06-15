/*### DESCRIPTION ############################################################*/
/*
 Given a newick formatted tree a description of a model, and a set of priors
 Runs an ABC which aims to determine the relationship between LCR length and 
 Transcript abundance
*/

//#include <cstdlib>
#include "AdaptiveParameter.hpp"
#include <algorithm>
#include "Chain.hpp"
#include "ctpl_stl.h"
#include <ctime>
#include <cmath>
#include <fstream>
#include <future>
#include <iostream>
#include "Logging.h"
#include <map>
#include <memory>
#include "optionparser.h"
#include "PriorIO.hpp"
#include <stdexcept>
#include <string>
#include <thread>
#include "Tree.hpp"


/*### DECLARATIONS and GLOBAL VARIABLES ######################################*/

std::string Version = "V0.0.0";

int logger::Verbosity = 2;


//Chain Variables

//Threading
const int MaxThreads = (std::thread::hardware_concurrency()) ? std::thread::hardware_concurrency() : 1;

struct Opts {
    int burnin, simSize, sampleSize, nChains, nThreads, tempHorizon;
    long unsigned int seed;
    double tgtSwpRate,initTempInc;
    std::string obs, prior, resume, tree;
};

enum InputMisMatchCase {VALID, VAR_OBS, PRIOR_OBS, TREE_OBS, MISSING_TIP, NO_OBS};

/*### HELP AND USAGE #########################################################*/

enum optionIndex{UNKNOWN,HELP,BURNIN,CGDISABLE,CPSHORIZON,CPSINIT,CPSOOM,CPSRATE,CSVALPHA,CSVN,MEERROR,MGEPROP,MGSMULT,MGSTOL,NCHAINS,NTHREADS,OBSFILE,PRIORFILE,QUIET,RESUME,SAMPLESIZE,SEED,SIMSIZE,THORIZON,TINIT,TRATE,TREEFILE,VERBOSITY};

struct Arg: public option::Arg{
    static void printError(const char* msg1, const option::Option& opt, const char* msg2){
        fprintf(stderr, "[logger::ERROR] %s", msg1);
        fwrite(opt.name, opt.namelen, 1, stderr);
        fprintf(stderr, "[logger::ERROR] %s", msg2);
    }
    
    static option::ArgStatus Integer(const option::Option& option, bool msg){
        char* endptr = 0;
        if (option.arg != 0 && strtol(option.arg, &endptr, 10)){};
        if (endptr != option.arg && *endptr == 0)
          return option::ARG_OK;
        
        if (msg) printError("Option '", option, "' requires an integer argument\n");
        return option::ARG_ILLEGAL;
    }

    static option::ArgStatus Natural(const option::Option& option, bool msg){
        char* endptr = 0;
        if (option.arg != 0 && strtol(option.arg, &endptr, 10) > 0){};
        if (endptr != option.arg && *endptr == 0)
          return option::ARG_OK;
        
        if (msg) printError("Option '", option, "' requires an natural number argument\n");
        return option::ARG_ILLEGAL;
    }

    static option::ArgStatus String(const option::Option& option, bool msg){
        if (option.arg != 0){return option::ARG_OK;};
        if (msg) printError("Option '", option, "' requires an string argument\n");
        return option::ARG_ILLEGAL;
    }

    static option::ArgStatus Proportion(const option::Option& option, bool msg){
        char* endptr = 0;
        double prop = 0;
        if (option.arg != 0){
           prop = strtod(option.arg, &endptr);
        }
        if (endptr != option.arg && *endptr == 0 && prop >= 0 && prop <= 1){
            return option::ARG_OK;
        }
        if (msg) printError("Option '", option, "' requires a numeric argument between 0 and 1\n");
        return option::ARG_ILLEGAL;
    }

    static option::ArgStatus NonZeroProportion(const option::Option& option, bool msg){
        char* endptr = 0;
        double prop = 0;
        if (option.arg != 0){
           prop = strtod(option.arg, &endptr);
        }
        if (endptr != option.arg && *endptr == 0 && prop > 0 && prop <= 1){
            return option::ARG_OK;
        }
        if (msg) printError("Option '", option, "' requires a real argument between 0 and 1\n");
        return option::ARG_ILLEGAL;
    }

    static option::ArgStatus NonNegativeReal(const option::Option& option, bool msg){
        char* endptr = 0;
        double prop = 0;
        if (option.arg != 0){
           prop = strtod(option.arg, &endptr);
        }
        if (endptr != option.arg && *endptr == 0 && prop >= 0){
            return option::ARG_OK;
        }
        if (msg) printError("Option '", option, "' requires a real argument of at least zero\n");
        return option::ARG_ILLEGAL;
    }

    static option::ArgStatus PositiveReal(const option::Option& option, bool msg){
        char* endptr = 0;
        double prop = 0;
        if (option.arg != 0){
           prop = strtod(option.arg, &endptr);
        }
        if (endptr != option.arg && *endptr == 0 && prop > 0){
            return option::ARG_OK;
        }
        if (msg) printError("Option '", option, "' requires a real argument greater than zero\n");
        return option::ARG_ILLEGAL;
    }
};

const std::map<optionIndex,std::string> usageMessages = {
    {CPSHORIZON,"[" + std::to_string(chain::CChain::GetPSHorizon()) + "] --proposal-scaling-horizon  [1,∞)εZ \tThe number of iterations between updates to scaling factor for parameter proposals. Also affects the horizon for sim-alpha updates, the latter is always greater and coprime to the former."},
    {CPSINIT,"[" + std::to_string(chain::CChain::GetPSInit()) + "] --initial-proposal-scale (0,∞)εR \tThe inital scale factor for parameter proposals."},
    {CPSOOM,"[" + std::to_string(chain::CChain::GetPSOoM()) + "] --proposal-scaling-extrema-magnitude ,ε \tThe bounds on the scale factor for parameter proposals in base 10 orders of magnitude."},
    {CPSRATE,"[" + std::to_string(chain::CChain::GetPSRate()) + "] --target-acceptance-rate (0,1]εR \t.The target rate of acceptance of parameter proposals"},
    {CSVALPHA,"[" + std::to_string(chain::CChain::GetSimVarAlpha()) + "] --initial-sim-variance-alpha (0,1]εR \tThe initial 'confidence level' for evaluated likelihoods; Value adapts to adjust acceptance rates."},
    {CSVN,"[" + std::to_string(chain::CChain::GetSimVarN()) + "] --sim-variance-n [2,∞)εZ \tThe number of full evaluations to perform to estimate the variability between simulations with the same parameter values."},
    {MEERROR,"[" + std::to_string(model::CModel::GetEvalRelError()) +"] --eval-rel-error [0,∞)εR \tThe maximum relative error between observed and simulated data to be considered a match during evaluation."},
    {MGEPROP,"[" + std::to_string(model::CModel::GetGEStepProp()) +"] --gradient-step-proportion (0,∞)εR \tThe size of a step to take when estimating the gradient relative to the current parameter value."},
    {MGSMULT,"[" + std::to_string(model::CModel::GetGSMult()) +"] --golden-search-boundary-mult (0,∞)εR \tWhen minimizing in the gradient direction, the maximum distance to search relative to parameter value."},
    {MGSTOL,"[" + std::to_string(model::CModel::GetGSTolProp()) +"] --golden-search-boundary-mult (0,∞)εR \tWhen minimizing in the gradient direction, searching is terminated once the search interval has shrunk to this proportion of its initial size."}
};

const option::Descriptor usage [] = {
    {UNKNOWN, 0, "","",option::Arg::None, "===Description\n Runs an approximate bayesian calculation with a given prior and its associated model."},
    {UNKNOWN,0, "","",option::Arg::None, "===USAGE\n run_abc [options] --prior file --observations file --tree file"},
    {UNKNOWN,0, "","",option::Arg::None, "===REQUIRED INPUTS"},
    {OBSFILE,0, "o","observations", Arg::String,"--observations, -o path\tA tab delimited file specifying the observed values at the tree tips. The file will have no header and columns of nodeID, abundance, and length."},
    {PRIORFILE,0, "p","prior", Arg::String,"--prior, -p path \tA file specifying the prior(s) to use. Each prior's first line is \">PriorName\". Each subsequent line is in the format key\\tvalue. values may be comma separated lists. Use --priorlist and --default-prior for more information."},
    {TREEFILE,0, "t","tree", Arg::String,"--tree, -t path \tA newick formatted file containing the tree on which to run simulations"},
    {UNKNOWN,0, "","",option::Arg::None, "\tNote: All nodeIDs in the observations must be node labes in the given tree"},
    {UNKNOWN,0, "","",option::Arg::None, "===MCMC OPTIONS"},
    {BURNIN,0,"b","burn-in",Arg::Natural,"[100] --burn-in, -b [1,∞)εZ \tThe number of accepted parameter sets to initially ignore."},
    {NCHAINS,0,"c","chain-count",Arg::Natural,"[1] --chain-count, -c [1,∞)εZ \tThe number of independent chains to run."},
    {SAMPLESIZE,0,"n","samples",Arg::Natural, "[1000] --samples, -n [1,∞)εZ \tThe number of accepted paramter sets at which to stop the calculation and estimate posteriors."},
    {SIMSIZE,0,"s","simulation-count",Arg::Natural, "[200] --simulation-count, -s [1,∞)εZ \tThe number of simulations to run for a generated set of model parameters. Each additional simulation improves the estimation of the distribution of values at tree nodes at the cost of run time."},
    {UNKNOWN,0, "","",option::Arg::None, "===TUNING OPTIONS"},
    {CGDISABLE,0,"","disable-gradient-descent",option::Arg::None," --disable-gradient-descent \tBy default, after a failed proposal, a local gradient is estimated and then a minimum opposite the gradient is found by golden search, disabling will explore less rapidly, but proposals will occur faster."},
    {CPSHORIZON,0,"","proposal-scaling-horizon",Arg::Natural,usageMessages.at(CPSHORIZON).c_str()},
    {CPSINIT,0,"","initial-proposal-scale",Arg::PositiveReal,usageMessages.at(CPSINIT).c_str()},
    {CPSOOM,0,"","proposal-scaling-extrema-magnitude",Arg::PositiveReal,usageMessages.at(CPSOOM).c_str()},
    {CPSRATE,0,"","target-acceptance-rate",Arg::NonZeroProportion,usageMessages.at(CPSRATE).c_str()},
    {CSVALPHA,0,"","initial-sim-variance-alpha",Arg::NonZeroProportion,usageMessages.at(CSVALPHA).c_str()},
    {CSVN,0,"","sim-variance-n",Arg::Natural,usageMessages.at(CSVN).c_str()},
    {MEERROR,0, "","eval-rel-error",Arg::PositiveReal, usageMessages.at(MEERROR).c_str()},
    {MGEPROP,0,"","gradient-step-proportion",Arg::NonNegativeReal, usageMessages.at(MGEPROP).c_str()},
    {MGSMULT,0,"","golden-search-boundary-mult",Arg::NonNegativeReal, usageMessages.at(MGSMULT).c_str()},
    {MGSTOL,0,"","golden-search-tolerance",Arg::NonNegativeReal, usageMessages.at(MGSTOL).c_str()},
    {THORIZON,0,"","swap-horizon",Arg::Natural,"[50] --swap-horizon [1,∞)εZ \tThe number of iterations between attempts to make swaps between chains; also the number of attempted swaps between updates to the temperature increment."},
    {TINIT,0,"","initial-temperature-increment",Arg::PositiveReal,"[0.1] --initial-temperature-increment (0,∞)εR \tThe initial increment in chain temperature between subsequent chains."},
    {TRATE,0,"","target-swap-rate",Arg::NonZeroProportion,"[0.5] --target-swap-rate (0,1]εR \tThe target proportion of proposed swaps which are accepted."},
    {UNKNOWN,0, "","",option::Arg::None, "===GENERAL OPTIONS"},
    {NTHREADS,0, "j","threads",Arg::Natural,"[1] --threads, -j [1,∞)εZ \t The maximum number of threads to use"},
    {RESUME,0,"r","resume-from",Arg::String,"--resume-from, -r path \tA results file from a previous run, will resume with the number of accepted samples from the file; All chains will start from the same point. Adaptive scaling parameters will be reset."},
    {SEED, 0, "", "seed", Arg::Natural, "[none] --seed [1,∞)εZ \tA seed for the random number generator"},
    {HELP, 0, "h", "help", option::Arg::None, " --help, -h \tPrint usage and exit"},
    {QUIET, 0, "q", "quiet", option::Arg::None, " --quiet, -q \tDecrement the level of verbosity."},
    {VERBOSITY, 0, "v", "verbose", option::Arg::None, " --verbose, -v \tIncrement the level of verbosity."},
    {UNKNOWN,0, "","",option::Arg::None, "\tNote: Default level 2: errors and warnings"},
    {0,0,0,0,0,0}
};



Opts ParseOptions(int argc, char ** argv){
    argc -=(argc>0); argv+=(argc>0); //Skips program name
    option::Stats stats(usage,argc,argv);
    option::Option options[stats.options_max], buffer[stats.buffer_max];
    option::Parser parse(usage, argc, argv, options, buffer);

    Opts opts;
    opts.burnin = 100;
    opts.sampleSize = 1000;
    opts.simSize = 200;
    opts.nChains = 1;
    opts.nThreads = 1;
    opts.seed = (long unsigned int)std::time(0);
    opts.tgtSwpRate = 0.5;
    opts.tempHorizon = 50;
    opts.initTempInc = 0.1;

    if(parse.error())
        exit(EXIT_FAILURE);

    if(options[HELP] || !options[OBSFILE] || !options[PRIORFILE] || !options[TREEFILE]){
        option::printUsage(std::cerr,usage);
        exit(EXIT_SUCCESS);
    }

    if(options[UNKNOWN]){
        std::cerr << "Unknown option: " << options[UNKNOWN].name << "\n";
        exit(EXIT_FAILURE);
    }

    if(options[VERBOSITY])
        logger::Verbosity += options[VERBOSITY].count();
    if(options[QUIET])
        logger::Verbosity -= options[QUIET].count();

    size_t cpsHorizon = chain::CChain::GetPSHorizon();
    double cpsInit = chain::CChain::GetPSInit();
    double cpsOoM = chain::CChain::GetPSOoM();
    double cpsRate = chain::CChain::GetPSRate();
    double csvAlpha = chain::CChain::GetSimVarAlpha();
    double csvN = chain::CChain::GetSimVarN();
    double meErr = model::CModel::GetEvalRelError();
    double mgeProp = model::CModel::GetGEStepProp();
    double mgsMult = model::CModel::GetGSMult();
    double mgsTol = model::CModel::GetGSTolProp();

    for(int i = 0; i < parse.optionsCount(); ++i){
        option::Option& opt = buffer[i];
        switch(opt.index()){
            case BURNIN:
                opts.burnin = atoi(opt.arg);
                break;
            case CGDISABLE:
                chain::CChain::DisableGradientDescent();
                break;
            case CPSHORIZON:
                cpsHorizon = atoi(opt.arg);
                break;
            case CPSINIT:
                cpsInit = atof(opt.arg);
                break;
            case CPSOOM:
                cpsOoM = atof(opt.arg);
                break;
            case CPSRATE:
                cpsRate = atof(opt.arg);
                break;
            case CSVALPHA:
                csvAlpha = atof(opt.arg);
                break;
            case CSVN:
                csvN = atoi(opt.arg);
                if(csvN < 2){
                    logger::Log("The simulation variablity handling N value has been reset to its minmum value of 2",logger::WARNING);
                    csvN = 2;
                }
               break;
            case MEERROR:
               meErr = atof(opt.arg);
                break;
            case MGEPROP:
                mgeProp = atof(opt.arg);
                break;
            case MGSMULT:
                mgsMult = atof(opt.arg);
                break;
            case MGSTOL:
                mgsTol = atof(opt.arg);
               break;
            case NCHAINS:
                opts.nChains = atoi(opt.arg);
                break;
            case NTHREADS:
                opts.nThreads = atoi(opt.arg);
                if(opts.nThreads > MaxThreads){
                    logger::Log("Proceeding with maximum sys_threads %d",logger::WARNING,MaxThreads);
                    opts.nThreads = MaxThreads;
                }
                break;
            case OBSFILE:
                opts.obs = opt.arg;
                break;
            case PRIORFILE:
                opts.prior = opt.arg;
                break;
            case RESUME:
                opts.resume = opt.arg;
                break;
            case SAMPLESIZE:
                opts.sampleSize = atoi(opt.arg);
                break;
            case SEED:
                opts.seed = atoi(opt.arg);
                break;
            case SIMSIZE:
                opts.simSize = atoi(opt.arg);
                if(opts.simSize < 1){
                    logger::Log("The simulation has been set to the minimum value of 1",logger::WARNING);
                    opts.simSize = 1;
                }
                break;
            case THORIZON:
                opts.tempHorizon = atoi(opt.arg);
                break;
            case TINIT:
                opts.initTempInc = atof(opt.arg);
                break;
            case TRATE:
                opts.tgtSwpRate = atof(opt.arg);
                break;
            case TREEFILE:
                opts.tree = opt.arg;
                break;
        }
    }

    chain::CChain::TunePS(cpsHorizon,cpsInit,cpsOoM,cpsRate);
    chain::CChain::TuneSimVar(csvAlpha,csvN);
    model::CModel::TuneEvaluation(meErr);
    model::CModel::TuneGradientEstimation(mgeProp);
    model::CModel::TuneGoldenSearch(mgsMult,mgsTol);

    return opts;
}


/*### FUNCTIONS ##############################################################*/

std::string ReadNewickFile(std::string path){
    logger::Log("Reading newick file %s...",logger::INFO,path.c_str());
    std::string newick;
    std::ifstream * file = new std::ifstream;
    file->open(path);
    if(file->fail()){
        throw std::logic_error("File does not exist");
    }
    char symbol;
    while((symbol = file->get()) != ';'){
        if(file->eof()){
            break;
        }
        newick += symbol;
    }
    file->close();
    delete file;
    logger::Log("Newick file read",logger::DEBUG);
    return newick + ';';
}

model::StateMap ReadObservationsFile(std::string path){
    logger::Log("Reading observations file %s...",logger::INFO,path.c_str());
    model::StateMap obs;
    std::ifstream * file = new std::ifstream;
    file->open(path);
    if(file->fail()){
        throw std::logic_error("File does not exist");
    }
    char symbol;
    std::string id, cur;
    model::ModelState state;
    int col = 0;
    while((symbol = file->get()) && !file->eof()){
        switch(symbol){
            case '\t':
               switch(col){
                   case 0:
                       id=cur;
                       break;
                    case 1:
                       state.abundance = CDiscreteFiniteRandomVariable({atoi(cur.c_str())},{1.0});
                       break;
                    case 2:
                       state.length = CDiscreteFiniteRandomVariable({atoi(cur.c_str())},{1.0});
                       break;
               } 
               col++;
               cur.clear();
               break;
            case '\n':
               if(col == 2)
                   state.length = CDiscreteFiniteRandomVariable({atoi(cur.c_str())},{1.0});
               cur.clear();
               col = 0;
               obs[id].push_back(state);
               break;
            default:
               cur += symbol;
        }
    }
    file->close();
    delete file;
    logger::Log("Observations file read",logger::DEBUG,path.c_str());
    return obs;
}

bool ValidateInput(const model::StateMap & obs, const prior::CPrior & prior, const Tree & tree){
    logger::Log("Validating input files...",logger::INFO);
    bool isValid = true;
    size_t nProt = obs.begin()->second.size();
    //Check for variable Sizing;
    for(const auto & x : obs){
        if(!x.second.size()){
            logger::Log("No observations parsed for node %s",logger::ERROR,x.first.c_str());
            isValid = false;
        }
        if(nProt != x.second.size()){
            logger::Log("Number of observed proteins for %s differs from other nodes",logger::ERROR,x.first.c_str());
            isValid = false;
        }
    }
    if(nProt != prior.getNProt()){
        logger::Log("Number of initial states in prior (%lu) must match the observed proteins (%ld) ",logger::ERROR,prior.getNProt(),nProt);
        isValid = false;
    }
    cvTreeNode vLeaves = tree.get_leaves();
    if(obs.size() != vLeaves.size()){
        logger::Log("Number of tips in tree (%ld) must match the number of observed nodes (%ld)",logger::ERROR,vLeaves.size(),obs.size());
        isValid = false;
    }
    for(const TreeNode* node : vLeaves){
        if(obs.find(node->id) == obs.end()){
            logger::Log("No observations provided for leaf node %s",logger::ERROR,node->id.c_str());
            isValid = false;
        }
    }
    logger::Log("Validation complete",logger::DEBUG);
    return isValid;
}

template<typename T>
std::string join(const std::vector<T> & v, std::string delim = ","/*, std::string mark = std::string(), std::string markIdx = 0;*/){
    if(v.empty()){
        return std::string();
    }
    std::string out = std::to_string(v[0]);
    //if(markIdx == 0){
    //    out = mark + out;
    //}
    for(int i = 1;i < v.size(); i++){
        out += delim;
        //if(markIdx == i){
        //    out += mark;
        //}
        out += std::to_string(v[i]);
    }
    return out;
}

std::string join(const std::vector<std::string> & v, std::string delim = ","/*, std::string mark = std::string(), std::string markIdx = 0;*/){
    if(v.empty()){
        return std::string();
    }
    std::string out = v[0];
    //if(markIdx == 0){
    //    out = mark + out;
    //}
    for(int i = 1;i < v.size(); i++){
        out += delim;
        //if(markIdx == i){
        //    out += mark;
        //}
        out += v[i];
    }
    return out;
}    
     
void InitialOutput(size_t nProt,const model::CModel & model){
        for(int i = 0; i < nProt; i++){
            std::cout << "Prot" << i << "\t";
    }
    std::cout << model.getName() << "\t\tnLogP\t";
    std::cout << join(model.getParamNames(),"\t") << "\n" << model << "\n";
}


std::unique_ptr<chain::CChain> ConstructChain(int threadID, int chainID, const prior::CPrior & prior, const Tree & tree, const model::StateMap & obs, unsigned long int seed, size_t nThreads,size_t maxSim){
    logger::Log("Calling chain %d constructor on thread %d ...",logger::DEBUG,chainID,threadID);
    chain::CChain * chainPtr = new chain::CChain(chainID,prior,tree,obs,seed,nThreads,maxSim);
    logger::Log("Constructed: chain %d on thread %d ...",logger::DEBUG,chainID,threadID);
    return std::unique_ptr<chain::CChain>(chainPtr);
}

std::unique_ptr<chain::CChain> ConstructChainFromModelStr(int threadID, int chainID, const prior::CPrior & prior, const Tree & tree, const model::StateMap & obs, unsigned long int seed, size_t nThreads,size_t maxSim, std::string modelStr){
    logger::Log("Calling chain %d constructor on thread %d ...",logger::DEBUG,chainID,threadID);
    chain::CChain * chainPtr = new chain::CChain(chainID,prior,tree,obs,seed,nThreads,maxSim,modelStr);
    logger::Log("Constructed: chain %d on thread %d ...",logger::DEBUG,chainID,threadID);
    return std::unique_ptr<chain::CChain>(chainPtr);
}

void GenerateChains(std::vector<std::unique_ptr<chain::CChain>> & vChains, const prior::CPrior & prior, const Tree & tree, const model::StateMap & obs, unsigned long int seed, size_t nChains, size_t nThreads, size_t maxSim, ctpl::thread_pool & threadPool){
    logger::Log("Generating %ld intial parameter set(s) ...",logger::INFO,nChains);
    //Each chain has its own thread pool
    std::mt19937 gen;
    gen.seed(seed-1);
    auto ptr = prior.GenerateModel(gen);
    size_t nBlocks = ptr->getNEvalBlocks();
    size_t nRemainThreads = nThreads;
    std::vector<std::future<std::unique_ptr<chain::CChain>>> vFutures;
    for(int chain = 0; chain < nChains; chain++){
        size_t nBlockThreads = std::ceil(double(nRemainThreads)/double(nChains - chain));
        nBlockThreads = std::min(nBlockThreads,nBlocks);
        if(nRemainThreads > nBlockThreads){
            nRemainThreads -= nBlockThreads;
        } else {
            nRemainThreads = 1;
        }
        auto future = threadPool.push(ConstructChain,
                chain,std::cref(prior),std::cref(tree),std::cref(obs),seed+chain,nBlockThreads,maxSim); 
        vFutures.push_back(std::move(future));
    }
    //Get the constructed Chains
    for(auto & future : vFutures){
        auto chain = future.get();
        vChains.push_back(std::move(chain));
    }
    logger::Log("Generated intial parameter set(s) ...",logger::DEBUG,nChains);
}


//Reads strings describing models into memory, also spits the resume file back out to cout so that end
//final results file contains the original as well
std::vector<std::string> ParseResultsFile(std::string path){
    logger::Log("Parsing previous results file %s ...",logger::INFO,path.c_str());
    std::ifstream * file = new std::ifstream;
    file->open(path);
    if(file->fail()){
        throw std::logic_error("Resume File does not exist");
    }
    char symbol;
    std::string cur;
    bool bHeader=true;
    std::vector<std::string> vModelStr;
    while((symbol = file->get()) && !file->eof()){
        switch(symbol){
            case '\n':
                std::cout << cur << "\n";
                if(!bHeader){
                    vModelStr.push_back(cur);
                } else{
                    bHeader=false;
                }
                cur.clear();
                break;
            default:
                cur += symbol;
        }
    }
    file->close();
    delete file;
    //Make sure the last entry gets included
    if(!cur.empty() && !bHeader){
        std::cout << cur << "\n";
        vModelStr.push_back(cur);
    }
    return vModelStr;
}

size_t ResumeChains(std::vector<std::unique_ptr<chain::CChain>> & vChains, const prior::CPrior & prior, const Tree & tree, const model::StateMap & obs, unsigned long int seed, size_t nChains, size_t nThreads, size_t maxSim, ctpl::thread_pool & threadPool, std::string resume){
    logger::Log("Generating %ld chains(s) ...",logger::INFO,nChains);
    std::vector<std::string> vModelStr = ParseResultsFile(resume);
    std::mt19937 gen;
    gen.seed(seed-1);
    auto ptr = prior.GenerateModel(gen);
    size_t nBlocks = ptr->getNEvalBlocks();
    size_t nRemainThreads = nThreads;
    std::vector<std::future<std::unique_ptr<chain::CChain>>> vFutures;
    for(int chain = 0; chain < nChains; chain++){
        size_t nBlockThreads = std::ceil(double(nRemainThreads)/double(nChains - chain));
        nBlockThreads = std::min(nBlockThreads,nBlocks);
        if(nRemainThreads > nBlockThreads){
            nRemainThreads -= nBlockThreads;
        } else {
            nRemainThreads = 1;
        }
        auto future = threadPool.push(ConstructChainFromModelStr,
                chain,std::cref(prior),std::cref(tree),std::cref(obs),seed+chain,nBlockThreads,maxSim,vModelStr.back()); 
        vFutures.push_back(std::move(future));
    }
    //Get the constructed Chains
    for(auto & future : vFutures){
        auto chain = future.get();
        vChains.push_back(std::move(chain));
    }
    return vModelStr.size(); 
}

/*### MAIN ###################################################################*/


int main(int argc, char ** argv){
    //TODO: Track down the cause of infinite loops in evaluation
    //Get command line input
    Opts opts = ParseOptions(argc,argv);
    
    //Load Observations
    model::StateMap obs = ReadObservationsFile(opts.obs);
    //Load Tree string
    std::string newick = ReadNewickFile(opts.tree);
    Tree tree = Tree(newick);
    
    //Load in the Prior
    logger::Log("Reading prior file %s ...",logger::INFO,opts.prior.c_str());
    std::unique_ptr<prior::CPrior> prior = prior::PriorIO(opts.prior).next_prior();
    if(!prior){
        logger::Log("No Prior Read from provided file",logger::ERROR);
        return 1;
    }
    //Make sure the files are compatible
    if(!ValidateInput(obs,*prior,tree)){
        logger::Log("Provided prior, tree, and obs files are incompatible",logger::ERROR);
        return 1;
    }
    //Set up randomNumberEngine for the main thread
    std::mt19937 gen;
    gen.seed(opts.seed);
    logger::Log("Main) seed %ld\n",logger::DEBUG,opts.seed);
    std::uniform_int_distribution<> ChainSelect(0,opts.nChains-1);

    //Set up the General ThreadPool which will only handle launching chain iterations
    size_t nChainThreads = std::min(opts.nChains,opts.nThreads);
    ctpl::thread_pool threadPool(nChainThreads);

    //Construct The Chain Objects
    int acceptCount = -1 * opts.burnin;
    std::vector<std::unique_ptr<chain::CChain>> vChains;
    if(opts.resume.empty()){
        GenerateChains(vChains,*prior,tree,obs,opts.seed+1,opts.nChains,opts.nThreads,opts.simSize,threadPool);
        if(vChains[0]->getModel().isFixed()){
            logger::Log("The Prior is Fixed for all parameters, assuming the desired behaviour is evaluation of the logLikelihood of the model defined by the prior\n",logger::WARNING);
            std::cout << vChains[0]->getLastEval() << "±" << vChains[0]->getEvaluationSD() << "\n";
            return 0;
        }
        InitialOutput(prior->getNProt(),vChains[0]->getModel());
    } else {
        acceptCount += ResumeChains(vChains,*prior,tree,obs,opts.seed+1,opts.nChains,opts.nThreads,opts.simSize,threadPool,opts.resume);
    }

    int iteration = 0;
    const std::vector<std::string> & vParamNames = vChains[0]->getModel().getParamNames();
    const size_t nParams = vParamNames.size();

    //Set up Adaptive Chain Heating;
    aparam::CAdaptiveParameter tempIncrement(opts.tgtSwpRate,opts.tempHorizon,opts.initTempInc,0,pow(10,9));

    int nParamCombos = pow(2.0,nParams)-1;

    while(acceptCount < opts.sampleSize){
        uint64_t propFlags = ((++iteration - 1) % nParamCombos) + 1;
        logger::Log("===Iteration %d (%d/%d)",logger::INFO,iteration,acceptCount,opts.sampleSize);
        logger::Log("Proposing changes on 0x%0x ...",logger::INFO,propFlags);

        //Determine which parameters are being modified
        std::unordered_set<std::string> proposalParamSet;
        for(int bit = 0; bit < nParams; bit++){
            if(!((propFlags >> bit) & 1)){
                continue;
            }
            proposalParamSet.insert(vParamNames[bit]);
        }

        //Launch the chains
        std::vector<std::future<bool>> vFutures;
        for(int chain = 0; chain < opts.nChains; chain++){
            double temperature = 1.0 / (1 + tempIncrement.getValue() * chain);
            auto future = threadPool.push(std::bind(&chain::CChain::iterate,vChains[chain].get(),std::placeholders::_1,proposalParamSet,temperature));
            vFutures.push_back(std::move(future));
        }
        //Synchronize the chains
        for(int chain = 0; chain < opts.nChains; chain++){
            if(vFutures[chain].get() && chain == 0){
                std::cout << vChains[0]->getModel() << "\n";
                acceptCount++;
            }
        }
        
        //Attempt Swapping Between Chains
        if(opts.nChains > 1 && iteration % opts.tempHorizon == 0){
            double tempUpdate = 0.0;
            //Shi and Rabosky 2015 (BAMM)
            int chain_i = ChainSelect(gen);
            int chain_j;
            do {chain_j = ChainSelect(gen);} while(chain_j == chain_i);
            if(chain_i > chain_j){
                std::swap(chain_i,chain_j);
            }
            logger::Log("Proposing Swap between chains %d and %d ...",logger::INFO,chain_i,chain_j);
            double temp_i = 1.0 / (1 + tempIncrement.getValue() * chain_i);
            double temp_j = 1.0 / (1 + tempIncrement.getValue() * chain_j);
            double nLogP_i = vChains[chain_i]->getModel().getNLogP();
            double nLogP_j = vChains[chain_j]->getModel().getNLogP();
            double sd_i = vChains[chain_i]->getEvaluationSD(); 
            double sd_j = vChains[chain_j]->getEvaluationSD(); 
            double sd_p = std::sqrt(std::pow(sd_i,2.0) + std::pow(sd_j,2.0));
            double sd_correction = stats::NormalQuantile(0.025) * sd_p;
            //Order is swapped to account for the negation of the logJP
            double logRatio = temp_i * (nLogP_i - nLogP_j + sd_correction);
            logRatio += temp_j * (nLogP_j - nLogP_i - sd_correction);
            double accept = (logRatio >= 0) ? 1.0 : std::exp(logRatio); 
            std::string result = "Reject";
            if(std::generate_canonical<double,10>(gen) < accept){
                vChains[chain_i]->swapModel(*vChains[chain_j]);
                tempUpdate = 1.0;
                result = "Accept";
            }
            if(tempIncrement.update(tempUpdate,gen)){
                logger::Log("Temperature increment updated to %0.04f",tempIncrement.getValue());
            }
            logger::Log("%s with p=%0.1f%%",logger::INFO,result.c_str(),accept*100);
        }
    }
    logger::Log("Done",logger::INFO);
    return 0;
}
