#pragma once

#include <boost/functional/hash.hpp>
#include <unordered_set>

#include "UniquePoseID.hpp"

namespace coaler::multialign {
    using UniquePoseSet = std::unordered_set<UniquePoseID, UniquePoseIdentifierHash>;
}
