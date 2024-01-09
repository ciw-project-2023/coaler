//
// Created by niklas on 1/9/24.
//

#include "UniquePoseID.hpp"
#include <unordered_set>
#include <boost/functional/hash.hpp>

#ifndef COALER_UNITQUEPOSESET_HPP
#define COALER_UNITQUEPOSESET_HPP
namespace coaler::multialign {
    using UniquePoseSet = std::unordered_set<UniquePoseID, UniquePoseIdentifierHash>;
}

#endif  // COALER_UNITQUEPOSESET_HPP
