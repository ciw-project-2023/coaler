//
// Created by niklas on 12/18/23.
//

#ifndef COALER_MATCHERMCS_H
#define COALER_MATCHERMCS_H

#include "Matcher.hpp"
#include "GraphMol/RDKitBase.h"

namespace coaler::core {
    class MatcherMCS : public Matcher {
    public:
        MatcherMCS(RDKit::MOL_SPTR_VECT mols);
    private:
        std::optional<RDKit::MCSResult> findMCS(RDKit::MOL_SPTR_VECT mols);
    };
}

#endif //COALER_MATCHERMCS_H
