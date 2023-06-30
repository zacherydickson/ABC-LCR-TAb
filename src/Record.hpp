#ifndef __RECORD__HPP
#define __RECORD__HPP

#include <map>
#include "Model.hpp"
#include <string>

namespace record {

    typedef std::map<std::string,double> ValueMap;
    typedef std::map<std::String,double> ValueVectorMap;


    //TODO: Change implementation to judge burnin just off of the nlogP
    class CRecord{
        //Cons/Destruction
        CSample() {}
        //Members
        protected:
            ValueVectorMap vSamples;
        //Static Methods
        protected:
            size_t estimateBurnin(const std::vector<double> & values) const;
            double estimateEffectiveSampleSize(const std::vector<double> & values) const;
            size_t kneedle(const std::vector<double> & values);
            std::vector<double> lowess(const std::vector<double> & values, double localityFactor);
        //Methods
        public:
            void addSample(const model::CModel & model) {addSample(model.getParamMap());}
            void addSample(const model::ParamMap & params);
            ValueMap estimateEffectiveSampleSizes();
        protected:
    };

}


#endif //__RECORD__HPP


