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
#include "Record.hpp"
#include <stdexcept>
#include <string>
#include <thread>
#include "Tree.hpp"


/*### DECLARATIONS and GLOBAL VARIABLES ######################################*/

std::string Version = "V0.0.0";

int logger::Verbosity = 3;

//Threading
const int MaxThreads = (std::thread::hardware_concurrency()) ? std::thread::hardware_concurrency() : 1;

struct Opts {
    Opts() : maxIter(0), maxSampleSize(0), simSize(400), nChains(1), nThreads(1), seed((long unsigned int)std::time(0)), tgtSwpRate(0.5), tempHorizon(50), initTempInc(0.1) {}
    int simSize, maxIter, maxSampleSize, nChains, nThreads, tempHorizon;
    long unsigned int seed;
    double tgtSwpRate,initTempInc;
    std::string obs, prior, resume, tree;
} defaultOpts;

enum InputMisMatchCase {VALID, VAR_OBS, PRIOR_OBS, TREE_OBS, MISSING_TIP, NO_OBS};
enum TerminationStatus {NOTERM, TERMITER, TERMSAMPLE, TERMESS};

/*### STRING HANDLING FUNCTIONS #############################################*/

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

/*### HELP AND USAGE #########################################################*/

enum optionIndex{UNKNOWN,HELP,CPSHORIZON,CPSINIT,CPSOOM,CPSRATE,CSVALPHA,CSVAHORIZ,CSVREHORIZ,CSVN,MAXESS,MAXITER,MEERROR,NCHAINS,NTHREADS,OBSFILE,PRIORFILE,QUIET,RECORDALPHA,RECORDEPSILON,RESUME,SAMPLESIZE,SEED,SIMSIZE,THORIZON,TINIT,TRATE,TREEFILE,VERBOSITY};

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

    static option::ArgStatus Whole(const option::Option& option, bool msg){
        char* endptr = 0;
        if (option.arg != 0 && strtol(option.arg, &endptr, 10) >= 0){};
        if (endptr != option.arg && *endptr == 0)
          return option::ARG_OK;
        
        if (msg) printError("Option '", option, "' requires a whole number argument\n");
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

const std::map<optionIndex,std::string> longOptionNames = {
    {HELP,"help"},
    {CPSHORIZON,"proposal-scaling-horizon"},
    {CPSINIT,"initial-proposal-scale"},
    {CPSOOM,"proposal-scaling-extrema-magnitude"},
    {CPSRATE,"target-acceptance-rate"},
    {CSVALPHA,"initial-sim-variance-alpha"},
    {CSVAHORIZ,"sim-variance-alpha-horizon"},
    {CSVREHORIZ,"sim-variance-evaluation-horizon"},
    {CSVN,"sim-variance-n"},
    {MAXESS,"max-effective-sample-size"},
    {MAXITER,"max-iterations"},
    {MEERROR,"eval-rel-error"},
    {NCHAINS,"chain-count"},
    {NTHREADS,"threads"},
    {OBSFILE,"observations"},
    {PRIORFILE,"prior"},
    {QUIET,"quiet"},
    {RECORDALPHA,"parameter-confidence"},
    {RECORDEPSILON,"parameter-precision"},
    {RESUME,"resume-from"},
    {SAMPLESIZE,"max-samples"},
    {SEED,"seed"},
    {SIMSIZE,"simulation-count"},
    {THORIZON,"swap-horizon"},
    {TINIT,"initial-temperature-increment"},
    {TRATE,"target-swap-rate"},
    {TREEFILE,"tree"},
    {VERBOSITY,"verbose"}
};

const std::map<optionIndex,std::string> usageMessages = {
    {CPSHORIZON,"[" + std::to_string(chain::CChain::GetPSHorizon()) + "] --" + longOptionNames.at(CPSHORIZON) + " [1,∞)εZ \tThe number of iterations between updates to scaling factor for parameter proposals."},
    {CPSINIT,"[" + join(chain::CChain::GetPSInit()) + "] --" + longOptionNames.at(CPSINIT) + " (0,∞)εR \tThe inital scale factor for parameter proposals. Can be either a single value applied to all parameters or a comma separated list which is recycled as necessary. parameter specific scales are processed in alphabetical order of parameter name"},
    {CPSOOM,"[" + std::to_string(chain::CChain::GetPSOoM()) + "] --" + longOptionNames.at(CPSOOM) + " (0,∞)εR \tThe bounds on the scale factor for parameter proposals in base 10 orders of magnitude."},
    {CPSRATE,"[" + std::to_string(chain::CChain::GetPSRate()) + "] --" + longOptionNames.at(CPSRATE) + " (0,1]εR \t.The target rate of acceptance of parameter proposals"},
    {CSVALPHA,"[" + std::to_string(chain::CChain::GetSimVarAlpha()) + "] --" + longOptionNames.at(CSVALPHA) + " (0,1]εR \tThe initial 'confidence level' for evaluated likelihoods; Value adapts to adjust acceptance rates."},
    {CSVAHORIZ,"[" + std::to_string(chain::CChain::GetSimVarAlphaHorizon()) + "] --" + longOptionNames.at(CSVAHORIZ) + " (0,∞)εZ \tThe number of iterations betwwen updates of the alpha value in variance handling. Lower values decrease the time spent at near minima."},
    {CSVREHORIZ,"[" + std::to_string(chain::CChain::GetSimVarReEvalHorizon()) + "] --" + longOptionNames.at(CSVREHORIZ) + " (0,∞)εZ \t The number of iterations between estimates of the current simulation variability. More frequent updates are more accurate, but will significant slow the program."},
    {CSVN,"[" + std::to_string(chain::CChain::GetSimVarN()) + "] --" + longOptionNames.at(CSVN) + " [2,∞)εZ \tThe number of full evaluations to perform to estimate the variability between simulations with the same parameter values."},
    {HELP, " --" + longOptionNames.at(HELP) + ", -h \tPrint usage and exit"},
    {MAXESS ,"[" + std::to_string(record::CRecord::GetMaxESS()) + "] --" + longOptionNames.at(MAXESS) + ", -n [0,∞)εZ \tThe effective sample size at which to stop iterating, 0 indicates no maximum."},  
    {MAXITER,"[" + std::to_string(defaultOpts.maxIter) + "] --" + longOptionNames.at(MAXITER) + " [0,∞)εZ \tThe maximum number of iterations, 0 indicates no maximum."}, 
    {MEERROR,"[" + std::to_string(model::CModel::GetEvalRelError()) +"] --" + longOptionNames.at(MEERROR) + " [0,∞)εR \tThe maximum relative error between observed and simulated data to be considered a match during evaluation."},
    {NCHAINS,"[" + std::to_string(defaultOpts.nChains) + "] --" + longOptionNames.at(NCHAINS) + ", -c [1,∞)εZ \tThe number of independent chains to run."},
    {NTHREADS,"[" + std::to_string(defaultOpts.nThreads) + "] --" + longOptionNames.at(NTHREADS) + ", -j [1,∞)εZ \t The maximum number of threads to use"},
    {PRIORFILE,"--" + longOptionNames.at(PRIORFILE) + ", -p path \tA file specifying the prior(s) to use. Each prior's first line is \">PriorName\". Each subsequent line is in the format key\\tvalue. values may be comma separated lists. Use --priorlist and --default-prior for more information."},
    {OBSFILE,"--" + longOptionNames.at(OBSFILE) + ", -o path\tA tab delimited file specifying the observed values at the tree tips. The file will have no header and columns of nodeID, abundance, and length."},
    {QUIET, " --" + longOptionNames.at(QUIET) + ", -q \tDecrement the level of verbosity."},
    {RECORDALPHA, "[" + std::to_string(record::CRecord::GetAlpha()) + "] --" + longOptionNames.at(RECORDALPHA) + ", -a (0,1]εR \tConfidence level for Parameter estimates."},
    {RECORDEPSILON, "[" + std::to_string(record::CRecord::GetEpsilon()) + "] --" + longOptionNames.at(RECORDEPSILON) + ", -e (0,1]εR \tRelative size of MCMC error to likelihood variance."},
    {RESUME,"--" + longOptionNames.at(RESUME) + ", -r path \tA results file from a previous run, will resume with the number of accepted samples from the file; All chains will start from the same point. Adaptive scaling parameters will be reset."},
    {SEED, "--" + longOptionNames.at(SEED) + " [1,∞)εZ \tA seed for the random number generator"},
    {SIMSIZE, "[" + std::to_string(defaultOpts.simSize) + "] --" + longOptionNames.at(SIMSIZE) + ", -s [1,∞)εZ \tThe number of simulations to run for a generated set of model parameters. Each additional simulation improves the estimation of the distribution of values at tree nodes at the cost of run time."},
    {SAMPLESIZE, "[" + std::to_string(defaultOpts.maxSampleSize) + "] --" + longOptionNames.at(SAMPLESIZE) + ", -N [0,∞)εZ \tThe maximum number of accepted paramter sets (0 indicates no maximum)."},
    {TRATE,"[" + std::to_string(defaultOpts.tgtSwpRate) + "] --" + longOptionNames.at(TRATE) + " (0,1]εR \tThe target proportion of proposed swaps which are accepted."},
    {TINIT,"[" + std::to_string(defaultOpts.initTempInc) + "] --" + longOptionNames.at(TINIT) + " (0,∞)εR \tThe initial increment in chain temperature between subsequent chains."},
    {THORIZON,"[" + std::to_string(defaultOpts.tempHorizon) + "] --" + longOptionNames.at(THORIZON) + " [1,∞)εZ \tThe number of iterations between attempts to make swaps between chains; also the number of attempted swaps between updates to the temperature increment."},
    {TREEFILE,"--" + longOptionNames.at(TREEFILE) + ", -t path \tA newick formatted file containing the tree on which to run simulations"},
    {VERBOSITY, " --" + longOptionNames.at(VERBOSITY) + ", -v \tIncrement the level of verbosity."}
};

const option::Descriptor usage [] = {
    {UNKNOWN, 0, "","",option::Arg::None, "===Description\n Runs an approximate bayesian calculation with a given prior and its associated model."},
    {UNKNOWN,0, "","",option::Arg::None, "===USAGE\n run_abc [options] --prior file --observations file --tree file"},
    {UNKNOWN,0, "","",option::Arg::None, "===REQUIRED INPUTS"},
    {OBSFILE,0, "o",longOptionNames.at(OBSFILE).c_str(), Arg::String,usageMessages.at(OBSFILE).c_str()},
    {PRIORFILE,0, "p",longOptionNames.at(PRIORFILE).c_str(), Arg::String,usageMessages.at(PRIORFILE).c_str()},
    {TREEFILE,0, "t",longOptionNames.at(TREEFILE).c_str(), Arg::String,usageMessages.at(TREEFILE).c_str()},
    {UNKNOWN,0, "","",option::Arg::None, "\tNote: All nodeIDs in the observations must be node labes in the given tree"},
    {UNKNOWN,0, "","",option::Arg::None, "===MCMC OPTIONS"},
    {RECORDALPHA,0,"a",longOptionNames.at(RECORDALPHA).c_str(),Arg::NonZeroProportion,usageMessages.at(RECORDALPHA).c_str()},
    {NCHAINS,0,"c",longOptionNames.at(NCHAINS).c_str(),Arg::Natural,usageMessages.at(NCHAINS).c_str()},
    {RECORDEPSILON,0,"e",longOptionNames.at(RECORDEPSILON).c_str(),Arg::NonZeroProportion,usageMessages.at(RECORDEPSILON).c_str()},
    {MAXESS,0,"n",longOptionNames.at(MAXESS).c_str(),Arg::Whole,usageMessages.at(MAXESS).c_str()},
    {SAMPLESIZE,0,"N",longOptionNames.at(SAMPLESIZE).c_str(),Arg::Whole,usageMessages.at(SAMPLESIZE).c_str()},
    {MAXITER,0,"",longOptionNames.at(MAXITER).c_str(),Arg::Whole,usageMessages.at(MAXITER).c_str()},
    {SIMSIZE,0,"s",longOptionNames.at(SIMSIZE).c_str(),Arg::Natural,usageMessages.at(SIMSIZE).c_str()},
    {UNKNOWN,0, "","",option::Arg::None, "===TUNING OPTIONS"},
    {CPSHORIZON,0,"",longOptionNames.at(CPSHORIZON).c_str(),Arg::Natural,usageMessages.at(CPSHORIZON).c_str()},
    {CPSINIT,0,"",longOptionNames.at(CPSINIT).c_str(),Arg::String,usageMessages.at(CPSINIT).c_str()},
    {CPSOOM,0,"",longOptionNames.at(CPSOOM).c_str(),Arg::PositiveReal,usageMessages.at(CPSOOM).c_str()},
    {CPSRATE,0,"",longOptionNames.at(CPSRATE).c_str(),Arg::NonZeroProportion,usageMessages.at(CPSRATE).c_str()},
    {CSVALPHA,0,"",longOptionNames.at(CSVALPHA).c_str(),Arg::NonZeroProportion,usageMessages.at(CSVALPHA).c_str()},
    {CSVAHORIZ,0,"",longOptionNames.at(CSVAHORIZ).c_str(),Arg::Natural,usageMessages.at(CSVAHORIZ).c_str()},
    {CSVREHORIZ,0,"",longOptionNames.at(CSVREHORIZ).c_str(),Arg::Natural,usageMessages.at(CSVREHORIZ).c_str()},
    {CSVN,0,"",longOptionNames.at(CSVN).c_str(),Arg::Natural,usageMessages.at(CSVN).c_str()},
    {MEERROR,0, "",longOptionNames.at(MEERROR).c_str(),Arg::PositiveReal, usageMessages.at(MEERROR).c_str()},
    {THORIZON,0,"",longOptionNames.at(THORIZON).c_str(),Arg::Natural,usageMessages.at(THORIZON).c_str()},
    {TINIT,0,"",longOptionNames.at(TINIT).c_str(),Arg::PositiveReal,usageMessages.at(TINIT).c_str()},
    {TRATE,0,"",longOptionNames.at(TRATE).c_str(),Arg::NonZeroProportion,usageMessages.at(TRATE).c_str()},
    {UNKNOWN,0, "","",option::Arg::None, "\tNote: Horizons should ideally be co-prime from eachother to minimize synchronized updates"},
    {UNKNOWN,0, "","",option::Arg::None, "===GENERAL OPTIONS"},
    {NTHREADS,0, "j",longOptionNames.at(NTHREADS).c_str(),Arg::Natural,usageMessages.at(NTHREADS).c_str()},
    {RESUME,0,"r",longOptionNames.at(RESUME).c_str(),Arg::String,usageMessages.at(RESUME).c_str()},
    {SEED, 0, "",longOptionNames.at(SEED).c_str(), Arg::Natural,usageMessages.at(SEED).c_str()},
    {HELP, 0, "h",longOptionNames.at(HELP).c_str(), option::Arg::None,usageMessages.at(HELP).c_str()},
    {QUIET, 0, "q",longOptionNames.at(QUIET).c_str(), option::Arg::None,usageMessages.at(QUIET).c_str()},
    {VERBOSITY, 0, "v",longOptionNames.at(VERBOSITY).c_str(), option::Arg::None,usageMessages.at(VERBOSITY).c_str()},
    {UNKNOWN,0, "","",option::Arg::None, "\tNote: Default level 3: errors, warnings, and info"},
    {0,0,0,0,0,0}
};



Opts ParseOptions(int argc, char ** argv){
    argc -=(argc>0); argv+=(argc>0); //Skips program name
    option::Stats stats(usage,argc,argv);
    option::Option options[stats.options_max], buffer[stats.buffer_max];
    option::Parser parse(usage, argc, argv, options, buffer);

    Opts opts;
    

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
    std::vector<double> cpsInit(chain::CChain::GetPSInit());
    double cpsOoM = chain::CChain::GetPSOoM();
    double cpsRate = chain::CChain::GetPSRate();
    double csvAlpha = chain::CChain::GetSimVarAlpha();
    double csvAlphaHorizon = chain::CChain::GetSimVarAlphaHorizon();
    double csvN = chain::CChain::GetSimVarN();
    double csvReEvalHorizon = chain::CChain::GetSimVarReEvalHorizon();
    double meErr = model::CModel::GetEvalRelError();
    double recAlpha = record::CRecord::GetAlpha();
    double recEpsilon = record::CRecord::GetEpsilon();
    double recMaxESS = record::CRecord::GetMaxESS();

    std::string bEnableGradDescentStr = "FALSE";

    for(int i = 0; i < parse.optionsCount(); ++i){
        option::Option& opt = buffer[i];
        size_t pos = 0;
        std::string initStr;
        switch(opt.index()){
            case CPSHORIZON:
                cpsHorizon = atoi(opt.arg);
                break;
            case CPSINIT:
                cpsInit.clear();
                initStr = std::string(opt.arg);
                while((pos = initStr.find(",")) != std::string::npos){
                    cpsInit.push_back(std::stod(initStr.substr(0,pos)));
                    initStr.erase(0, pos+1);
                }
                cpsInit.push_back(std::stod(initStr.substr(0,pos)));
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
            case CSVAHORIZ:
                csvAlphaHorizon = atoi(opt.arg);
                break;
            case CSVN:
                csvN = atoi(opt.arg);
                if(csvN < 2){
                    logger::Log("The simulation variablity handling N value has been reset to its minmum value of 2",logger::WARNING);
                    csvN = 2;
                }
               break;
            case CSVREHORIZ:
               csvReEvalHorizon = atoi(opt.arg);
               break;
            case MAXESS:
               recMaxESS = atoi(opt.arg);
               break;
            case MAXITER:
               opts.maxIter = atoi(opt.arg);
               break;
            case MEERROR:
               meErr = atof(opt.arg);
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
            case RECORDALPHA:
                recAlpha = atof(opt.arg);
                break;
            case RECORDEPSILON:
                recEpsilon = atof(opt.arg);
                break;
            case RESUME:
                opts.resume = opt.arg;
                break;
            case SAMPLESIZE:
                opts.maxSampleSize = atoi(opt.arg);
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
    chain::CChain::TuneSimVar(csvAlpha,csvAlphaHorizon,csvN,csvReEvalHorizon);
    model::CModel::TuneEvaluation(meErr);
    record::CRecord::TuneThreshold(recAlpha,recEpsilon,recMaxESS);

    std::stringstream message;
    message << "\t" << longOptionNames.at(OBSFILE) << ":\t" << opts.obs << "\n";
    message << "\t" << longOptionNames.at(PRIORFILE) << ":\t" << opts.prior << "\n";
    message << "\t" << longOptionNames.at(TREEFILE) << ":\t" << opts.tree << "\n";
    message << "\t" << longOptionNames.at(RECORDALPHA) << ":\t" << recAlpha << "\n";
    message << "\t" << longOptionNames.at(NCHAINS) << ":\t" << opts.nChains << "\n";
    message << "\t" << longOptionNames.at(RECORDEPSILON) << ":\t" << recEpsilon << "\n";
    message << "\t" << longOptionNames.at(MAXESS) << ":\t" << recMaxESS << "\n";
    message << "\t" << longOptionNames.at(SAMPLESIZE) << ":\t" << opts.maxSampleSize << "\n";
    message << "\t" << longOptionNames.at(MAXITER) << ":\t" << opts.maxIter << "\n";
    message << "\t" << longOptionNames.at(SIMSIZE) << ":\t" << opts.simSize << "\n";
    message << "\t" << longOptionNames.at(CPSHORIZON) << ":\t" << cpsHorizon << "\n";
    message << "\t" << longOptionNames.at(CPSINIT) << ":\t" << join(cpsInit) << "\n";
    message << "\t" << longOptionNames.at(CPSOOM) << ":\t" << cpsOoM << "\n";
    message << "\t" << longOptionNames.at(CPSRATE) << ":\t" << cpsRate << "\n";
    message << "\t" << longOptionNames.at(CSVALPHA) << ":\t" << csvAlpha << "\n";
    message << "\t" << longOptionNames.at(CSVAHORIZ) << ":\t" << csvAlphaHorizon << "\n";
    message << "\t" << longOptionNames.at(CSVREHORIZ) << ":\t" << csvReEvalHorizon << "\n";
    message << "\t" << longOptionNames.at(CSVN) << ":\t" << csvN << "\n";
    message << "\t" << longOptionNames.at(MEERROR) << ":\t" << meErr << "\n";
    message << "\t" << longOptionNames.at(THORIZON) << ":\t" << opts.tempHorizon << "\n";
    message << "\t" << longOptionNames.at(TINIT) << ":\t" << opts.initTempInc << "\n";
    message << "\t" << longOptionNames.at(TRATE) << ":\t" << opts.tgtSwpRate << "\n";
    message << "\t" << longOptionNames.at(NTHREADS) << ":\t" << opts.nThreads << "\n";
    message << "\t" << longOptionNames.at(SEED) << ":\t" << opts.seed << "\n";
    if(!opts.resume.empty()){
        message << "\t" << longOptionNames.at(RESUME) << ":\t" << opts.resume << "\n";
    }
    logger::Log("===Parameters\n%s",logger::INFO,message.str().c_str());
    return opts;
}


/*### FUNCTIONS ##############################################################*/

std::string ReadNewickFile(std::string path){
    logger::Log("Reading newick file %s...",logger::INFO,path.c_str());
    std::string newick;
    std::ifstream * file = new std::ifstream;
    file->open(path);
    if(file->fail()){
        throw std::logic_error("Newick file does not exist");
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
        throw std::logic_error("Newick file does not exist");
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

    
     
void InitialOutput(size_t nProt,const model::CModel & model){
        for(int i = 0; i < nProt; i++){
            std::cout << "Prot" << i << "\t";
    }
    std::cout << model.getName() << "\t\tnLogP\t";
    std::cout << join(model.getParamNames(),"\t") << "\n" << model << std::endl;
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
    std::cout << std::flush;
    return vModelStr;
}

void ResumeChains(std::vector<std::unique_ptr<chain::CChain>> & vChains, const prior::CPrior & prior, const Tree & tree, const model::StateMap & obs, unsigned long int seed, size_t nChains, size_t nThreads, size_t maxSim, ctpl::thread_pool & threadPool, std::string resume, record::CRecord & sampleRecord){
    logger::Log("Generating %ld chains(s) ...",logger::INFO,nChains);
    std::vector<std::string> vModelStr = ParseResultsFile(resume);
    std::mt19937 gen;
    gen.seed(seed-1);
    auto ptr = prior.GenerateModel(gen);
    size_t nBlocks = ptr->getNEvalBlocks();
    size_t nRemainThreads = nThreads;
    //Iterate over the model strings and add them to the record
    //  Sips the first model, as it isn't actually an accepted sample
    for(int i = 1; i < vModelStr.size(); i++){
        ptr->setToStr(vModelStr[i]);
        sampleRecord.addSample(*ptr);
    }
    //Construct the chains from the last model
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
}

std::string toBinary(uint64_t flag, size_t bits){
    std::stringstream stream;
    for(int bit = bits; bit >= 0; bit--){
        stream << ((flag >> bit) & 1);
    }
    return stream.str();
}

std::vector<uint64_t> initializePropFlags(const model::CModel & model, std::mt19937 & gen){
    logger::Log("Initializing Proposal Order...",logger::DEBUG);
    const std::vector<std::string> vParamNames = model.getParamNames();
    const model::ParamMap & params = model.getParamMap();
    size_t nParams = params.size();
    int nParamCombos = pow(2.0,nParams)-1;
    std::vector<uint64_t> vPropFlags;
    //Determine which parameters if any are fixed
    uint64_t fixedMask = 0;
    for(int bit = 0; bit < nParams; bit++){
        if(!params.at(vParamNames[bit]).isFixed()) {
            continue;
        }
        fixedMask |= 1 << bit;
    }
    //Get all valid proposal Flags (i.e only proposing on non-fixed parameters)
    for(int i = 0; i < nParamCombos; i++){
        uint64_t flag = i+1;
        if(flag & fixedMask){ //Skip the combo if it includes a fixed parameter
            continue;
        }
        vPropFlags.push_back(flag);
    }
    //Randomize the order of the proposalFlags
    std::shuffle(std::begin(vPropFlags),std::end(vPropFlags),gen);
    //Log the Proposal Flag Vector 
    std::stringstream stream;
    int counter = 1;
    int lineWidth = 80;
    int perLine = lineWidth / (nParams + 6);
    for(auto flag : vPropFlags){
        stream << "0b" << toBinary(flag,nParams-1) << "    ";
        if(counter++ % perLine == 0 && counter <= vPropFlags.size()){
            stream << "\n";
        }
    }
    logger::Log("Proposals will occur in the order:\n%s",logger::INFO,stream.str().c_str());
    return vPropFlags;
}

TerminationStatus getTerminationStatus(size_t iteration, const Opts & opts, record::CRecord & sampleRecord){
    //Check max iterations
    if(opts.maxIter > 0 && iteration >= opts.maxIter){
        return TERMITER;
    }
    //Check max sample size
    if(opts.maxSampleSize > 0 && sampleRecord.size() >= opts.maxSampleSize){
        return TERMSAMPLE;
    }
    if(sampleRecord.isComplete()){
        return TERMESS;
    }
    return NOTERM;
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


    //The record of accepted samples with the list of non-fixed parameters
    record::CRecord *sampleRecord;

    //Construct The Chain Objects
    std::vector<std::unique_ptr<chain::CChain>> vChains;
    if(opts.resume.empty()){
        GenerateChains(vChains,*prior,tree,obs,opts.seed+1,opts.nChains,opts.nThreads,opts.simSize,threadPool);
        std::stringstream stream;
        stream << vChains[0]->getModel() << "\n";
        logger::Log("Initial Model\n%s",logger::DEBUG,stream.str().c_str());
        if(vChains[0]->getModel().isFixed()){
            logger::Log("The Prior is Fixed for all parameters, assuming the desired behaviour is evaluation of the logLikelihood of the model defined by the prior",logger::WARNING);
            vChains[0]->estimateSimulationVariance();
            std::cout << vChains[0]->getLastEval() << "±" << vChains[0]->getEvaluationSD() << std::endl;
            return 0;
        }
        sampleRecord = new record::CRecord(prior->getModelParameterNames(false));
        InitialOutput(prior->getNProt(),vChains[0]->getModel());
    } else {
        sampleRecord = new record::CRecord(prior->getModelParameterNames(false));
        ResumeChains(vChains,*prior,tree,obs,opts.seed+1,opts.nChains,opts.nThreads,opts.simSize,threadPool,opts.resume,*sampleRecord);
    }
    logger::Log("Minimum sample size is %d",logger::INFO,sampleRecord->nextMilestone());

    int iteration = 0;
    const std::vector<std::string> & vParamNames = vChains[0]->getModel().getParamNames();
    const size_t nParams = vParamNames.size();

    //Set up Adaptive Chain Heating;
    aparam::CAdaptiveParameter tempIncrement(opts.tgtSwpRate,opts.tempHorizon,opts.initTempInc,0,pow(10,9));


    //Initialize Vector of bitVectors which specify the proposals to make
    std::vector<uint64_t> vPropFlags = initializePropFlags(vChains[0]->getModel(),gen);

    TerminationStatus termStatus = NOTERM;
    while((termStatus = getTerminationStatus(iteration,opts,*sampleRecord)) == NOTERM){
        int milestone = std::min(opts.maxSampleSize,int(sampleRecord->nextMilestone()));
        if(milestone == 0){
            milestone = int(sampleRecord->nextMilestone());
        }
        uint64_t propFlags = vPropFlags[iteration++ % vPropFlags.size()];
        logger::Log("===Iteration %d (%d/%d+)",logger::INFO,iteration,sampleRecord->size(),milestone);
        logger::Log("Proposing changes on 0b%s ...",logger::INFO,toBinary(propFlags,nParams-1).c_str());

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
                std::cout << vChains[0]->getModel() << std::endl;
                sampleRecord->addSample(vChains[0]->getModel());
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
            //Order is swapped to account for the negation of the logJP
            double logRatio = temp_i * (nLogP_i - nLogP_j);
            logRatio += temp_j * (nLogP_j - nLogP_i);
            double accept = (logRatio >= 0) ? 1.0 : std::exp(logRatio); 
            std::string result = "Reject";
            if(std::generate_canonical<double,10>(gen) < accept){
                vChains[chain_i]->swapModel(*vChains[chain_j]);
                tempUpdate = 1.0;
                result = "Accept";
            }
            if(tempIncrement.update(tempUpdate,gen)){
                logger::Log("Temperature increment updated to %0.04f",logger::INFO,tempIncrement.getValue());
            }
            logger::Log("%s with p=%0.1f%%",logger::INFO,result.c_str(),accept*100);
        }
    }
    std::string exitMessage = "Done";
    logger::LogLevel exitLogLevel = logger::INFO;
    switch(termStatus){
        case TERMITER:
            exitMessage = "Reached User Specified Maximum Iterations";
            exitLogLevel = logger::WARNING;
            break;
        case TERMSAMPLE:
            exitMessage = "Reached User Specified Maximum Samples";
            exitLogLevel = logger::WARNING;
            break;
    }
    logger::Log("%s",exitLogLevel,exitMessage.c_str());
    delete sampleRecord;
    return 0;
}
