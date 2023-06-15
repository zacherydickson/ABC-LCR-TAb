#ifndef ADAPTIVE_PARAMETER_HPP
#define ADAPTIVE_PARAMETER_HPP

#include <limits>
#include <random>
#include <vector>

namespace aparam {

    class CAdaptiveParameter{
        //An adaptive paramter that attempts to keep the rate at some target by altering
        // its value
        // A scaleFactor defines the behaviour if the rate is above target
        //  i.e 0.5 means cut in half if above target, and double if below
        //  i.e 3 means triple if above target, and cut to a third if below
        //All scaling has a random component where the updated value will be between
        //  itself and the scaled amount
        //Con-Destruction
        public:
            CAdaptiveParameter(double targetRate, size_t horizon, double initValue,
                    double minValue = -std::numeric_limits<double>::infinity(),
                    double maxValue = std::numeric_limits<double>::infinity(),
                    double scaleFactor = 2);
        //Static Members
        public:
        //Members
        public:
            const size_t horizon;
            const double target;
            const double factor;
        protected:
            int counter;
            double rate;
            double rate0;
            double value;
            double minValue;
            double maxValue;
        //Methods
        public:
            const double & getValue() const {return this->value;}
            const double & getMax() const {return this->maxValue;}
            const double & getMin() const {return this->minValue;}
            bool update(double increment, std::mt19937 & gen);
    };
}



#endif //ADAPTIVE_PARAMETER_HPP



