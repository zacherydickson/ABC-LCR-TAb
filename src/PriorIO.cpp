#include "Distributions.hpp"
#include "PriorIO.hpp"
#include "Logging.h"
#include <cstdlib>
#include <utility>

namespace prior {

    //#### CON/DESTRUCTION ########################################################
    
    
    //Opens a prior file for reading, automatically parses to the start of the first Prior
    PriorIO::PriorIO(const std::string & file) :
        bEOP (true)
    {
        this->fileHandle = new std::fstream;
        this->fileHandle->open(file);
        if(this->fileHandle->fail()){
            throw std::logic_error("File does not exists");
        }
        //Scan for the start of a prior
        while(!this->fileHandle->eof() && this->fileHandle->get() != '>'){}
    }
    
    PriorIO::~PriorIO(){
        this->fileHandle->close();
        delete this->fileHandle;
        this->fileHandle = NULL;
    }

    //#### STATIC METHODS #########################################################

    bool PriorIO::parseHyperparameterValue(const std::string & curStr, const std::string & name, SParameterPriorLine & line){
        if(curStr.empty() && name.empty()){ //Extra tab between params or at end of line
            return false;
        }
        if(name.empty() || curStr.empty()){
            throw std::invalid_argument("ParameterPrior line has an unnamed hyperparameter value or name without a value");
        }
        char* str_end = NULL;
        const char* cStr = curStr.c_str();
        double value = std::strtod(cStr,&str_end);
        if(*str_end == 0){ //If a valid conversion occurs str_end points to a nul
            line.paramPrior.hyperparameters[name]=value;
        } else if(line.paramPrior.type == stats::Fixed && name == "mu"){
            //Special case that allows parameters to be bound to the value
            //of another parameter
            line.paramPrior.dependency=curStr;
        } else {
            throw std::invalid_argument("ParameterPrior line has a non-numeric hyperparameter value");
        }
        return true;
    }
    
    
    //#### METHODS ################################################################
    
    std::unique_ptr<CPrior> PriorIO::next_prior(){
        if(this->fileHandle->eof()){
            return std::unique_ptr<CPrior>(nullptr);
        }
        std::string name;
        std::fstream *& fh = this->fileHandle;
        char symbol;
        while(!fh->eof() && (symbol = fh->get()) != '\n'){
            name += symbol;
        }
        this->bEOP = false;
        model::ModelType type;
        try{
            type = model::str2ModelType(name);
        } catch (std::invalid_argument & e){
            throw std::invalid_argument("Prior read from file requests an unknown model");
        }
        ParameterPriorMap mParamPriors;
        while(!this->bEOP){
            SParameterPriorLine line = this->parse_line();
            if(!this->bEOP){
                mParamPriors[line.name] = line.paramPrior;
            }
        }
        //Iterate over and check for dependencies
        for(auto & pair : mParamPriors){
            if(!pair.second.dependency.empty()){
                if(mParamPriors.count(pair.second.dependency) == 0){
                    logger::Log("Encountered unknown parameter dependency (%s) for %s",
                            logger::ERROR,pair.second.dependency.c_str(),pair.first.c_str());
                    throw std::invalid_argument("Unknown parameter dependency");
                }
                const SParameterPriorSpecification & dependency = mParamPriors.at(pair.second.dependency);
                pair.second.type = dependency.type;
                pair.second.hyperparameters = dependency.hyperparameters;
            }
        }
        return std::unique_ptr<CPrior>(new CPrior(type,mParamPriors));
    }


    
    //Line Format
    //paramName\tDisttype\tHyperparamName0:HyperparamValue0\t...HyperparamNameN:HyperparamValueN\n
    SParameterPriorLine PriorIO::parse_line(){
        std::string name, curStr;
        bool bType = false;
        bool bLineComplete = false;
        char symbol;
        std::fstream *& fh = this->fileHandle;
        SParameterPriorLine line;
        while(!bLineComplete && !fh->eof() && (symbol = fh->get()) != '>'){
            switch(symbol){
                case '\t': //Separator between paramName, DistType, and each hyperparameter
                    if(line.name.empty()){ //First column is Name
                        if(curStr.empty()){
                            throw std::invalid_argument("ParameterPrior Line does not have a parameter name");
                        }
                        line.name = curStr;
                    } else if(!bType){ //The second column is the distribution Type
                        try {
                            line.paramPrior.type = stats::str2DistributionType(curStr);
                        } catch (std::invalid_argument & e){
                            throw std::invalid_argument("ParameterPrior line has an unrecognized distribution type");
                        }
                        bType = true;
                    } else { //curStr now contains the value of a hyperparameter
                        if(!PriorIO::parseHyperparameterValue(curStr,name,line)){
                            break;
                        }
                        name.clear();
                    }
                    curStr.clear();
                    break;
                case ':': //curStr contains the name of a hyperparameter
                    name = curStr;
                    curStr.clear();
                    break;
                case '\n'://Make sure the last hyperparameter value gets stored
                    if(curStr.empty() && line.name.empty()){ //Empty Line;Skip
                        break;
                    }
                    bLineComplete = true;
                    PriorIO::parseHyperparameterValue(curStr,name,line);
                    break;
                default:
                    curStr += symbol;
            }
        }
        if(!bLineComplete){ //If the finish wasn't on a new line it was at a > or eof
            this->bEOP = true;
        }
        return line;
    }

}
