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
    UniquePoseIdentifier() = default;

    UniquePoseIdentifier(LigandID ligandId,
                         PoseID ligandInternalPoseId);

    /*----------------------------------------------------------------------------------------------------------------*/

    bool operator==(const UniquePoseIdentifier &other) const;

    bool operator!=(const UniquePoseIdentifier &other) const;

    bool operator<(const UniquePoseIdentifier &other) const;

    bool operator>(const UniquePoseIdentifier &other) const;

    [[maybe_unused]] [[nodiscard]] std::string toString() const noexcept;

    [[nodiscard]] LigandID getLigandId() const;

    [[nodiscard]] PoseID getLigandInternalPoseId() const;

    /*----------------------------------------------------------------------------------------------------------------*/

  private:

    LigandID m_ligandId;
    PoseID m_ligandInternalPoseId;
};

struct UniquePoseIdentifierHash
{
    std::size_t operator()(const UniquePoseIdentifier& uniquePoseId) const;

    friend class UniquePoseIdentifier;
};
}