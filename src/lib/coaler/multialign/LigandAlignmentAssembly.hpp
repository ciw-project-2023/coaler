#pragma once
#include <unordered_map>

#include "models/Forward.hpp"

namespace coaler::multialign {

    /**
     * An alignment of a set of ligands. Contains one pose for each ligand.
     */
    class LigandAlignmentAssembly {
      public:
        explicit LigandAlignmentAssembly(const std::unordered_map<LigandID, PoseID>& initialAssembly);

        /**
         * Exchange the associated pose for a given Ligand
         * @param ligandId The ligand whose pose is to be swapped.
         * @param newPoseId The new pose.
         *
         * @note if the ligand isnt present in the assembly its added with the pose assigned to it
         */
        void swapPoseForLigand(LigandID ligandId, PoseID newPoseId);

        /**
         * @param ligandId Ligand to get the associated pose for.
         * @return The Pose that is associated with the @p ligandId.
         */
        [[nodiscard]] PoseID getPoseOfLigand(LigandID ligandId) const;

        /**
         * increase missing ligands count by one
         */
        void incrementMissingLigandsCount();

        /**
         * decrease missing ligands count by one
         */
        void decrementMissingLigandsCount();

        /**
         * @return The missing ligand count for the assembly.
         */
        [[nodiscard]] unsigned getMissingLigandsCount() const noexcept;

        // NOLINTBEGIN(cppcoreguidelines-non-private-member-variables-in-classes)
        // NOLINTBEGIN(misc-non-private-member-variables-in-classes)
        std::unordered_map<LigandID, PoseID> getAssemblyMapping() const noexcept;
        // NOLINTEND(cppcoreguidelines-non-private-member-variables-in-classes)
        // NOLINTEND(misc-non-private-member-variables-in-classes)

      private:
        void setMissingLigandsCount(unsigned count);

        LigandAlignmentAssembly() = default;

        bool insertLigandPose(LigandID ligand, PoseID pose);

        std::unordered_map<LigandID, PoseID> m_assembly;
        unsigned m_missingLigandsCount{0};

        friend class StartingAssemblyGenerator;
        friend class MultiAligner;
    };

}  // namespace coaler::multialign
