#include <boost/range/combine.hpp>

#include "catch2/catch.hpp"
#include "coaler/core/Forward.hpp"
#include "coaler/embedder/BucketDistributor.hpp"

using namespace coaler::embedder;

void check_distribution(unsigned nofMatches, unsigned maxConfs, const std::vector<unsigned>& expected_dist) {
    std::vector<unsigned> dist = BucketDistributor::distributeApproxEvenly(nofMatches, maxConfs);
    REQUIRE(dist.size() == expected_dist.size());
    for (const auto& itertuple : boost::combine(dist, expected_dist)) {
        CHECK(itertuple.get<0>() == itertuple.get<1>());
    }
}

/*----------------------------------------------------------------------------------------------------------------*/

TEST_CASE("validate_distribute_evenly", "[conformer_generator_tester]") {
    check_distribution(3, 7, {3, 2, 2});
    check_distribution(2, 7, {4, 3});
    check_distribution(5, 10, {2, 2, 2, 2, 2});
    check_distribution(6, 20, {4, 4, 3, 3, 3, 3});
}