#include <cmath>
#include "RandomVariable.hpp"
#include <stdexcept>

CDiscreteFiniteRandomVariable::CDiscreteFiniteRandomVariable(std::vector<int> vVal, std::vector<double> vProb){
    if(vVal.size() != vProb.size()){
        throw std::invalid_argument("Attempt to Construct CDiscreteRandomVariable with incongruent value and probability vectors\n");
    }
    int maxP = -1;
    for(int i = 0; i < vVal.size(); i++){
        this->probMassMap[vVal[i]] = vProb[i];
        if(vProb[i] > maxP){
            this->mode = vVal[i];
            maxP = vProb[i];
        }
    }
}

double CDiscreteFiniteRandomVariable::PMF(int x) const{
    if(this->probMassMap.count(x)){
        return this->probMassMap.at(x);
    }
    return 0.0;
}

//CProbabilityEstimate
CProbabilityEstimate::CProbabilityEstimate(double tolerance, int domainMin, int domainMax,double mode) :
    tolerance(tolerance), estimate(0), bComplete(false), bModeCovered(false),
    domainMin(domainMin), domainMax(domainMax), mode(mode)
{
    lower.x = domainMax;
    upper.x = domainMin;
    lower.p = upper.p = max.p = -1;
}

bool CProbabilityEstimate::addTerm(int x, double p, bool bLogP){
    double tmp = (bLogP) ? p/std::log(10) : std::log10(p);
   //fprintf(stderr,"addTerm(%d,%0.4f)\n",x,tmp);
    double estError = 0;
    if(bLogP){
        double tmp = std::log(this->estimate);
        tmp += std::log(1+std::exp(p-tmp));
        if(this->estimate != 0 || std::isfinite(p)){
            estError = std::exp(p-tmp);
        }
        this->estimate = std::exp(tmp);
    } else {
        this->estimate += p;
        if(p != 0 || this->estimate != 0){
            estError = p/this->estimate;
        }
    }
    SProbabilityPair & l = this->lower;
    SProbabilityPair & m = this->max;
    SProbabilityPair & u = this->upper;
    if(x > u.x || u.x < l.x){
        u.x = x;
        u.p = p;
    }
    if(x < l.x || u.x < l.x){
        l.x = x;
        l.p = p;
    }
    if(p > m.p){
        m.x = x;
        m.p = p;
    }
    if(m.p > l.p && m.p > u.p){
        this->bModeCovered = true;
    }
    if(m.x == this->domainMin && u.x > m.x || m.x == this->domainMax && l.x < m.x){
        this->bModeCovered = true;
    }
    if(!std::isnan(this->mode) && l.x <= this->mode && u.x >= this->mode){
        this->bModeCovered = true;
    }
    if(this->bModeCovered && estError < this->tolerance){
        this->bComplete = true;
    } if(l.x <= this->domainMin && u.x >= this->domainMax){
        this->bComplete = true;
    }
    //fprintf(stderr,"E=%0.5f(ε:%0.04f),l[%d,%0.4f] m[%d,%0.4f], u[%d,%0.4f] (%d,%d)\n",std::log10(this->estimate),std::log10(estError),l.x,std::log10(l.p),m.x,std::log10(m.p),u.x,std::log10(u.p),bModeCovered,bComplete);
    return bComplete;
}

double CProbabilityEstimate::getEstimate() const{
    if(!bComplete){
        throw std::logic_error("Attempt to get an incomplete probability estimate");
    }
    return this->estimate;
}

//Hint is a guess as to where the mode might be to reduce search time in 0 probability
//areas
int CProbabilityEstimate::nextX(int hint) const{
    const SProbabilityPair & l = this->lower;
    const SProbabilityPair & m = this->max;
    const SProbabilityPair & u = this->upper;
    std::vector<int> vNexts;
    if(l.x > this->domainMin){
        vNexts.push_back(l.x - 1);
    }
    if(u.x < this->domainMax){
        vNexts.push_back(u.x + 1);
    }
    if(vNexts.size() == 0){ //The domain is covered and the estimate can't be improved
        throw std::logic_error("Attempt to get the next x for an exact probability estimate");
    }
   //fprintf(stderr,"l.x:%d, u.x:%d, m.x:%d\n",l.x,u.x,m.x);
   //fprintf(stderr,"l.p:%d, u.p:%d, m.p:%d\n",l.p,u.p,m.p);
    if(vNexts.size() == 1){ //Only one next is possible, go with it
       //fprintf(stderr,"nextX[0]: %d\n",vNexts[0]);
        return vNexts[0];
    }
    //Strategy: (Assumes a smooth unimodal curve)
    //  if mode is covered
    //      go to closer Bound
    //  if mode is uncovered
    //      go to the bound with the greater probability
    int idx = 0; //Default to going lower
    int modeX = (this->bModeCovered) ? m.x : hint;
    if(l.p < u.p){
        idx = 1;
    } else if(std::abs(modeX - l.x) > std::abs(modeX - u.x)){
        idx = 1;
    }
   //fprintf(stderr,"nextX: %d\n",vNexts[idx]);
    return vNexts[idx];
}


