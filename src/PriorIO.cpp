#include "Distributions.hpp"
#include "PriorIO.hpp"
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
                        if(curStr.empty() && name.empty()){ //Extra tab between params
                            break;
                        }
                        if(name.empty() || curStr.empty()){
                            throw std::invalid_argument("ParameterPrior line has an unnamed hyperparameter value or name without a value");
                        }
                        line.paramPrior.hyperparameters[name]=std::strtod(curStr.c_str(),NULL);
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
                    if(curStr.empty() && name.empty()){ //Extra tab at end of line
                        break;
                    }
                    if(name.empty() || curStr.empty()){
                        throw std::invalid_argument("ParameterPrior line ends with an unnamed hyperparameter value or a name without a value");
                    }
                    line.paramPrior.hyperparameters[name]=std::strtod(curStr.c_str(),NULL);
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
