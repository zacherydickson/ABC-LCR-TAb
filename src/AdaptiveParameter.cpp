#include "AdaptiveParameter.hpp"
#include "Logging.h"
#include <cmath>

namespace aparam {

    CAdaptiveParameter::CAdaptiveParameter(double target, size_t horizon, double initValue, double minValue, double maxValue, double scaleFactor) :
        target(target), horizon(horizon), factor(scaleFactor),
        counter(0), rate(0), rate0(target),
        value(initValue), maxValue(maxValue), minValue(minValue)
    {
        if(minValue > maxValue){
            //fprintf(stderr,"%0.4f\t%0.4f\n",minValue,maxValue);
            throw std::invalid_argument("Attempt construct CAdaptiveParameter with invalid bounds");
        }
    }
    
    
    //For a simple acceptance Rate the increment should be either 0 or 1 (a fail or a pass)
    bool CAdaptiveParameter::update(double increment, std::mt19937 & gen){
        rate = (counter * rate + increment) / double(counter + 1.0);
        if(++counter < horizon){
            return false;
        }
        logger::Log("Updating Adaptive Parameter: target: %0.4f, rate: %0.4f, rate0: %0.4f, value: %0.4f",logger::DEBUG+3,target,rate,rate0,value);
        //All that is important is the sign 
        double offset = rate - target; //neg if below target pos otherwise
        double delta = (rate == rate0) ? offset : rate - rate0; //neg if neg change, pos if pos, offset otherwise 
        //neg * neg = pos -> decrease value
        //pos * pos = pos -> increase value
        //neg * pos = neg -> hold value
        counter = 0;
        rate0 = rate;
        rate = 0;
        if(offset * delta >= 0){ //The offset and the rate are in the same direction
            double f = (offset > 0) ? factor : (1 / factor); 
            double p = std::generate_canonical<double,10>(gen);
            //Weighted average of itself and the scaled amount
            value *= 1 + p*(f - 1);
            if(value > maxValue){
                value = maxValue;
            } else if(value < minValue){
                value = minValue;
            }
            return true;
        }
        return false;
    }

}
