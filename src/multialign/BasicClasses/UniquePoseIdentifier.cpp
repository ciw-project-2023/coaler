//
// Created by malte on 11/27/23.
//

#include "UniquePoseIdentifier.hpp"

namespace MultiAlign
{

    UniquePoseIdentifier::UniquePoseIdentifier(LigandID ligandId,
                         PoseID ligandInternalPoseId)
        : m_ligandId(ligandId), m_ligandInternalPoseId(ligandInternalPoseId){}

    /*----------------------------------------------------------------------------------------------------------------*/

    bool UniquePoseIdentifier::operator==(const UniquePoseIdentifier &other) const
    {
        return this->m_ligandId == other.m_ligandId && this->m_ligandInternalPoseId == other.m_ligandInternalPoseId;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    bool UniquePoseIdentifier::operator!=(const UniquePoseIdentifier &other) const
    {
        return !(*this == other);
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    bool UniquePoseIdentifier::operator<(const UniquePoseIdentifier &other) const
    {
        if(this->m_ligandId == other.m_ligandId)
        {
            return this->m_ligandInternalPoseId <
                   other.m_ligandInternalPoseId;
        }
        return this->m_ligandId < other.m_ligandId;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    bool UniquePoseIdentifier::operator>(const UniquePoseIdentifier &other) const
    {
        if(this->m_ligandId == other.m_ligandId)
        {
            return this->m_ligandInternalPoseId >
                   other.m_ligandInternalPoseId;
        }
        return this->m_ligandId > other.m_ligandId;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    [[maybe_unused]] std::string UniquePoseIdentifier::toString() const noexcept
    {
        return std::to_string(m_ligandId) + "-" + std::to_string(m_ligandInternalPoseId);
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    LigandID UniquePoseIdentifier::getLigandId() const
    {
            return m_ligandId;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    PoseID UniquePoseIdentifier::getLigandInternalPoseId() const
    {
            return m_ligandInternalPoseId;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    std::size_t UniquePoseIdentifierHash::operator()(const UniquePoseIdentifier& uniquePoseId) const
    {
        const std::size_t seed = 0;
        const std::string key = uniquePoseId.toString();
        return std::hash<std::string>{}(key);
    }

}
