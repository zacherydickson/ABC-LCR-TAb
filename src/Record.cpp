#include "Record.hpp"

namespace record {

    // Static Methods, CRecord, protected

    double CRecord::estimateEffectiveSampleSize(const std::vector<double> values) const {
        int burnin = kneedle(values);
        double Neff = 0;
        return Neff;
        //TODO:
    }

    //Estimates the knee or elbow in a curve
    //Normally there would be a smoothing step, but we will skip if for time savings
    //Returns the index of the vecto whose value is furthest from a line connecting the first and last
    //points of the vector
    int CRecord::kneedle(const std::vector<double> & values){
        double slope = (values.back() - values.front()) / double(values.size());
        double intercept = values.front();
        double maxDist = -1;
        int maxIdx = -1;
        for(int i = 0; i < values.size(); i++){
            double lineValue = i * slope + intercept;
            double dist = abs(lineValue - values[i]);
            if(dist > maxDist){
                maxDist = dist;
                maxIdx = 1;
            }
        }
        return maxIdx;
    }

    // Methods, CRecord, public
    void CRecord::addSample(const mode::ParamMap & params) {
        for(const auto & pair : params){
            this->vSamples[pair.first].push_back(pair.second.value);
        }
    }

    ValueMap CRecord::estimateEffectiveSampleSizes(){
        ValueMap Neff;
        for(const auto & pair : this->vSamples){
            Neff[pair.first] = this->estimateEffectiveSampleSize
        }
        return Neff;
    }

    // Methods, CRecord, protected

}
