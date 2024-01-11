#include "AssemblyIDManager.hpp"

namespace coaler::multialign {

    namespace {
        std::size_t calculate_assembly_hash(const coaler::multialign::LigandAlignmentAssembly &assembly) {
            std::vector<size_t> unique_pos_id_hashs;
            for (auto id_pair : assembly.getAssemblyMapping()) {
                LigandID ligand_id = id_pair.first;
                PoseID pos_id = id_pair.second;
                UniquePoseID u_pose_id(ligand_id, pos_id);
                unique_pos_id_hashs.emplace_back(UniquePoseIdentifierHash()(u_pose_id));
            }

            std::hash<size_t> hasher;
            size_t assembly_hash = 0;
            for (size_t value : unique_pos_id_hashs) {
                assembly_hash ^= hasher(value) + 0x9e3779b9 + (assembly_hash << 6) + (assembly_hash >> 2);
            }
            return assembly_hash;
        }
    }  // namespace

    bool AssemblyIDManager::is_assembly_new(const coaler::multialign::LigandAlignmentAssembly &assembly) {
        size_t assembly_hash = calculate_assembly_hash(assembly);
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