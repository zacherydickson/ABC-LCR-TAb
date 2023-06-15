#ifndef RANDOMVARIABLE_HPP_
#define RANDOMVARIABLE_HPP_

#include <unordered_map>
#include <limits>
#include <vector>

typedef std::unordered_map<int,double> ProbabilityMap;

class CDiscreteFiniteRandomVariable {
    //Con-/Desctruction
    public:
        CDiscreteFiniteRandomVariable(){}
        CDiscreteFiniteRandomVariable(std::vector<int> vVal, std::vector<double> vProb);
    //Members
    protected:
        ProbabilityMap probMassMap;
        int mode;
    //Methods
    public:
        double PMF(int x) const;
        int getMode() const {return mode;}
        const ProbabilityMap & getProbMap() const {return probMassMap;}
    //Operators
    public:
        CDiscreteFiniteRandomVariable & operator=(const CDiscreteFiniteRandomVariable & rhs){
            this->probMassMap = rhs.probMassMap;
            this->mode = rhs.mode;
            return *this;
        }
        bool operator==(const CDiscreteFiniteRandomVariable & rhs) const{
            if(this->probMassMap.size() != rhs.probMassMap.size()){
                return false;
            }
            for(const auto & pair : this->probMassMap){
                if(rhs.probMassMap.count(pair.first) == 0){
                    return false;
                }
                if(rhs.probMassMap.at(pair.first) != pair.second){
                    return false;
                }
            }
            return true;
        }
        bool operator<(const CDiscreteFiniteRandomVariable & rhs) const{
            if(this->mode < rhs.mode){
                return true;
            }
            return false;
        }
};

struct SProbabilityPair{
    int x;
    double p;
};

class CProbabilityEstimate{
    //Con-/Desctruction
    public:
        CProbabilityEstimate(double tolerance) :
            CProbabilityEstimate(tolerance,
                    std::numeric_limits<int>::min(),
                    std::numeric_limits<int>::max(),
                    std::numeric_limits<double>::quiet_NaN()) {}
        CProbabilityEstimate(double tolerance,double mode) :
            CProbabilityEstimate(tolerance,
                    std::numeric_limits<int>::min(),
                    std::numeric_limits<int>::max(),
                    mode) {}
        CProbabilityEstimate(double tolerance, int domainMin) :
            CProbabilityEstimate(tolerance,domainMin,
                    std::numeric_limits<int>::max(),
                    std::numeric_limits<double>::quiet_NaN()) {}
        CProbabilityEstimate(double tolerance, int domainMin, double mode) :
            CProbabilityEstimate(tolerance,domainMin,
                    std::numeric_limits<int>::max(),
                    mode) {}
        CProbabilityEstimate(double tolerance,int domainMin,int domainMax) :
            CProbabilityEstimate(tolerance,domainMin,domainMax,
                    std::numeric_limits<double>::quiet_NaN()) {};
        CProbabilityEstimate(double tolerance,int domainMin,int domainMax, double mode);
    //Members
    public:
        const double mode;
        const double tolerance;
        const int domainMin;
        const int domainMax;
    protected:
        bool bComplete;
        bool bModeCovered;
        double estimate;
        SProbabilityPair lower;
        SProbabilityPair max;
        SProbabilityPair upper;
    //Methods
    public:
        bool addTerm(int x, double p, bool bLogP = false);
        void forceComplete(){bComplete = true;};
        double getEstimate() const;
        bool isComplete()const {return bComplete;}
        int nextX() const{return this->nextX(max.x);}
        int nextX(int hint) const;
};



#endif //RANDOMVARIABLE_HPP_


