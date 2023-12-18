//
// Created by malte on 11/27/23.
//

#include "UniquePoseID.hpp"

namespace coaler::multialign {

    UniquePoseID::UniquePoseID(LigandID ligandId, PoseID ligandInternalPoseId)
        : m_ligandId(ligandId), m_ligandInternalPoseId(ligandInternalPoseId) {}

    /*----------------------------------------------------------------------------------------------------------------*/

    bool UniquePoseID::operator==(const UniquePoseID &other) const {
        return this->m_ligandId == other.m_ligandId && this->m_ligandInternalPoseId == other.m_ligandInternalPoseId;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    bool UniquePoseID::operator!=(const UniquePoseID &other) const { return !(*this == other); }

    /*----------------------------------------------------------------------------------------------------------------*/

    bool UniquePoseID::operator<(const UniquePoseID &other) const {
        if (this->m_ligandId == other.m_ligandId) {
            return this->m_ligandInternalPoseId < other.m_ligandInternalPoseId;
        }
        return this->m_ligandId < other.m_ligandId;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    bool UniquePoseID::operator>(const UniquePoseID &other) const {
        if (this->m_ligandId == other.m_ligandId) {
            return this->m_ligandInternalPoseId > other.m_ligandInternalPoseId;
        }
        return this->m_ligandId > other.m_ligandId;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    [[maybe_unused]] std::string UniquePoseID::toString() const noexcept {
        return std::to_string(m_ligandId) + "-" + std::to_string(m_ligandInternalPoseId);
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    LigandID UniquePoseID::getLigandId() const { return m_ligandId; }

    /*----------------------------------------------------------------------------------------------------------------*/

    PoseID UniquePoseID::getLigandInternalPoseId() const { return m_ligandInternalPoseId; }

    /*----------------------------------------------------------------------------------------------------------------*/

    std::size_t UniquePoseIdentifierHash::operator()(const UniquePoseID &uniquePoseId) const {
        const std::string key = uniquePoseId.toString();
        return std::hash<std::string>{}(key);
    }

}  // namespace coaler::multialign
