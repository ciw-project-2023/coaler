//
// Created by malte on 12/4/23.
//

#pragma once

#include "Forward.hpp"

namespace coaler::core {

    class Core {
      public:

        /**
         * constructor for object of class Core
         * @param molecules the molecules a core will be calculated for
         * @param coreType type of core calculated: MCS or Murcko Scaffold
         */
        Core(const std::vector<RDKit::RWMol>& molecules, const coreType coreType);

        /**
         * getter function for the type of scaffold the core is calculated on
         * @return coreType of core
         */
        [[nodiscard]] coreType getScaffoldType() const;

        /**
         * getter function for the core
         * @return core as ROMol
         */
        [[nodiscard]] CoreAsMol getCore() const;


      private:
        /**
         * calculates the MCS of the molecules
         * @return MCS as ROMol
         */
        CoreAsMol calculateCoreMcs() const;

        /**
         * calculates the Murcko scaffold of the molecules
         * @return Murcko Scaffold as ROMol
         */
        CoreAsMol calculateCoreMurcko() const;

        CoreAsMol m_core;
        std::vector<RDKit::RWMol> m_molecules;
        coreType m_coreType;
    };
}