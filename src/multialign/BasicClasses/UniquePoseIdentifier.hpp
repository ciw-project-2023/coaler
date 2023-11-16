//
// Created by chris on 11/16/23.
//

#pragma once
#include "Forward.hpp"

struct UniquePoseIdentifier
{
public:
    UniquePoseIdentifier(LigandID ligandId,
                         PoseID ligandInternalPoseId)
                         : m_ligandId(ligandId)
                         , m_ligandInternalPoseId(ligandInternalPoseId)
    {
    }

private:
    LigandID m_ligandId;
    PoseID m_ligandInternalPoseId;
};

struct UniquePoseIdentifierHash
{
    std::size_t operator()(const UniquePoseIdentifier& id) const
    {
        std::size_t seed = 0;
        std::string key = std::to_string(id.m_ligandId) + "-" +  std::to_string(id.m_ligandInternalPoseId);
        return std::hash<std::string>(key);
    }};
};
