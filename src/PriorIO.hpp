#ifndef PRIORIO_HPP
#define PRIORIO_HPP

#include <fstream>
#include <map>
#include <memory>
#include "Prior.hpp"
#include <string>
#include <vector>

namespace prior {

    struct SParameterPriorLine{
        std::string name;
        SParameterPriorSpecification paramPrior;
    };

    class PriorIO {
        public:
            PriorIO(const std::string & file);
            ~PriorIO();
        //Members
        private:
            std::fstream* fileHandle;
            bool bEOP;
            std::string lastName;
        //Methods
        private:
            size_t parse_pairs();
            SParameterPriorLine parse_line();
        public:
            std::unique_ptr<CPrior> next_prior();
    };

}

#endif


