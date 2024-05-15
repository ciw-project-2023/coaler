#pragma once
#include "../LigandAlignmentAssembly.hpp"
#include "Forward.hpp"

namespace coaler::multialign {

    /**
     * The AssemblyIDManager class provides functionality for the management of assembly IDs.
     */
    class AssemblyIDManager {
      public:
        AssemblyIDManager() = default;

        /**
         * Checks if the assembly is not previously considered.
         * @param assembly
         * @return
         */
        bool isAssemblyNew(const LigandAlignmentAssembly &assembly);

      private:
        std::vector<size_t> m_existing_assembly_hashes{};
    };
}  // namespace coaler::multialign
