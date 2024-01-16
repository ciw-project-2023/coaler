#include "AssemblyIDManager.hpp"

namespace coaler::multialign {

    namespace {
        std::size_t calculateAssemblyHash(const coaler::multialign::LigandAlignmentAssembly &assembly) {
            std::vector<size_t> uniquePosIdHashs;
            for (auto idPair : assembly.getAssemblyMapping()) {
                LigandID ligandId = idPair.first;
                PoseID posId = idPair.second;
                UniquePoseID uPoseId(ligandId, posId);
                uniquePosIdHashs.emplace_back(UniquePoseIdentifierHash()(uPoseId));
            }

            std::hash<size_t> hasher;
            size_t assemblyHash = 0;
            for (size_t value : uniquePosIdHashs) {
                assemblyHash ^= hasher(value) + 0x9e3779b9 + (assemblyHash << 6) + (assemblyHash >> 2);
            }
            return assemblyHash;
        }
    }  // namespace

    bool AssemblyIDManager::isAssemblyNew(const coaler::multialign::LigandAlignmentAssembly &assembly) {
        size_t assemblyHash = calculateAssemblyHash(assembly);
        auto it = std::find(m_existing_assembly_hashes.begin(), m_existing_assembly_hashes.end(), assemblyHash);

        // Check if assembly combination is unique
        if (it != m_existing_assembly_hashes.end()) {
            return false;
        }

        // add new unique pose id
        m_existing_assembly_hashes.emplace_back(assemblyHash);
        return true;
    }

}  // namespace coaler::multialign
