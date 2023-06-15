#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <atomic>

namespace sim {

    struct SPolicy {
        SPolicy(int nsim, int maxsim) {nSim=nsim;maxSim=maxsim,nZero=0;nTotal=0;}
        int nSim;
        int maxSim;
        std::atomic<int> nZero;
        std::atomic<int> nTotal;
    };

} //namespace sim


#endif //SIMULATION_HPP


