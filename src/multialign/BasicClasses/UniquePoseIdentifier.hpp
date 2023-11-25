//
// Created by chris on 11/16/23.
//

#pragma once
#include "Forward.hpp"

namespace MultiAlign
{
/**
 * Identifier for Conformers using the id of their ligand and the
 * internal pose id.
 */
struct UniquePoseIdentifier //TODO move implementation to cpp
{
    UniquePoseIdentifier()
    {};

    /*----------------------------------------------------------------------------------------------------------------*/

    UniquePoseIdentifier(LigandID ligandId,
                         PoseID ligandInternalPoseId)
                         : m_ligandId(ligandId)
                         , m_ligandInternalPoseId(ligandInternalPoseId)
    {
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    [[nodiscard]] std::string toString() const noexcept{
        return std::to_string(m_ligandId) + "-" + std::to_string(m_ligandInternalPoseId);
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    //TODO wtf
    bool operator!=(const UniquePoseIdentifier &other) const{
        return this->m_ligandId != other.m_ligandId
        || this->m_ligandInternalPoseId != other.m_ligandInternalPoseId;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    bool operator==(const UniquePoseIdentifier &other) const{
        return !(*this != other);
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    bool operator>(const UniquePoseIdentifier &other) const
    {
        if(this->m_ligandId == other.m_ligandId)
        {
            return this->m_ligandInternalPoseId >
                   other.m_ligandInternalPoseId;
        }
        return this->m_ligandId > other.m_ligandId;
    }

    /*----------------------------------------------------------------------------------------------------------------*/

    LigandID m_ligandId;
    PoseID m_ligandInternalPoseId;
};

struct UniquePoseIdentifierHash
{
    std::size_t operator()(const UniquePoseIdentifier& id) const
    {
        std::size_t seed = 0;
        std::string key = std::to_string(id.m_ligandId) + "-" +  std::to_string(id.m_ligandInternalPoseId);
        return std::hash<std::string>{}(key);
    }
};
}