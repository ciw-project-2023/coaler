#pragma once

#include <string>

#include "Alias.hpp"

namespace coaler::multialign {
    /**
     * Identifier for Conformers using the id of their ligand and the
     * internal pose id.
     */
    struct UniquePoseID  // TODO move implementation to cpp
    {
        UniquePoseID() = default;

        UniquePoseID(LigandID ligandId, PoseID ligandInternalPoseId);

        /*----------------------------------------------------------------------------------------------------------------*/

        bool operator==(const UniquePoseID &other) const;

        bool operator!=(const UniquePoseID &other) const;

        bool operator<(const UniquePoseID &other) const;

        bool operator>(const UniquePoseID &other) const;

        [[nodiscard]] std::string toString() const noexcept;

        [[nodiscard]] LigandID getLigandId() const;

        [[nodiscard]] PoseID getLigandInternalPoseId() const;

        /*----------------------------------------------------------------------------------------------------------------*/

      private:
        LigandID m_ligandId;
        PoseID m_ligandInternalPoseId;
    };

    struct UniquePoseIdentifierHash {
        std::size_t operator()(const UniquePoseID &uniquePoseId) const;

        friend class UniquePoseID;
    };
}  // namespace coaler::multialign
