//
// Created by niklas on 1/9/24.
//

#include <boost/functional/hash.hpp>
#include <unordered_set>

#include "UniquePoseID.hpp"

#ifndef COALER_UNIQUEPOSESET_HPP
#define COALER_UNIQUEPOSESET_HPP
namespace coaler::multialign {
    using UniquePoseSet = std::unordered_set<UniquePoseID, UniquePoseIdentifierHash>;
}

#endif  // COALER_UNIQUEPOSESET_HPP
