#include "AssemblyIDManager.hpp"

namespace coaler::multialign {

    namespace {
        std::size_t calculateAssemblyHash(const coaler::multialign::LigandAlignmentAssembly &assembly) {
            std::vector<size_t> uniquePosIdHashs;
            for (auto id_pair : assembly.getAssemblyMapping()) {
                LigandID ligandId = id_pair.first;
                PoseID posId = id_pair.second;
                UniquePoseID uPoseId(ligandId, posId);
                uniquePosIdHashs.emplace_back(UniquePoseIdentifierHash()(uPoseId));
            }

            std::hash<size_t> hasher;
            size_t assembly_hash = 0;
            for (size_t value : uniquePosIdHashs) {
                assembly_hash ^= hasher(value) + 0x9e3779b9 + (assembly_hash << 6) + (assembly_hash >> 2);
            }
            return assembly_hash;
        }
    }  // namespace

    bool AssemblyIDManager::isAssemblyNew(const coaler::multialign::LigandAlignmentAssembly &assembly) {
        size_t assembly_hash = calculateAssemblyHash(assembly);
        auto it = std::find(m_existing_assembly_hashes.begin(), m_existing_assembly_hashes.end(), assembly_hash);

        // Check if assembly combination is unique
        if (it != m_existing_assembly_hashes.end()) {
            return false;
        }

        // add new unique pose id
        m_existing_assembly_hashes.emplace_back(assembly_hash);
        return true;
    }

}  // namespace coaler::multialign
