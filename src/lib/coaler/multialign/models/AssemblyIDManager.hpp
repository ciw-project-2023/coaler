#pragma once
#include "../LigandAlignmentAssembly.hpp"
#include "Forward.hpp"

namespace coaler::multialign {

    /**
     * @brief This class is responsible for managing Assembly IDs.
     */
    class AssemblyIDManager {
      public:
        AssemblyIDManager() = default;

        /**
         * Checks if the assembly is not previously considered.
         * @param assembly
         * @return
         */
        bool is_assembly_new(const LigandAlignmentAssembly &assembly);

      private:
        std::vector<size_t> m_existing_assembly_hashes{};
    };
}  // namespace coaler::multialign