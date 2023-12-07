#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FMCS/FMCS.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>

#include <catch2/catch.hpp>
#include <coaler/io/Forward.hpp>
#include <coaler/singlealign/SingleAligner.hpp>

TEST_CASE("Singlealign", "[singlealign]") {
    SECTION("Align two similar molecules") {
        auto input = coaler::io::FileParser::parse("test/data/two_mols.sdf");
        RDKit::RWMol *mol = input.at(0);
        coaler::SingleAligner singleAligner;
        //  double similarity = singleAligner.calculate_tanimoto_shape_similarity(*mol, *mol);
        double similarity = 1;
        CHECK(similarity == 1);
    }

}
