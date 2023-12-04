//
// Created by malte on 12/4/23.
//

#pragma once

#include "Forward.hpp"

namespace coaler::core {

    class Core {
      public:
        Core(const std::vector<RDKit::RWMol>& molecules, const coreType coreType);

        [[nodiscard]] coreType getScaffoldType() const;

        [[nodiscard]] CoreAsMol getCore() const;



      private:
        bool calculateCore() const;

        CoreAsMol m_core;
        std::vector<RDKit::RWMol> m_molecules{};
        coreType m_coreType;

    };


}



