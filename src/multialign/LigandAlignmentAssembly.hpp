//
// Created by chris on 11/5/23.
//

#pragma once
#include "Forward.hpp"
#include <unordered_map>

namespace MultiAlign {

    /**
     * An alignment of a set of ligands. Contains one pose for each ligand.
     */
    class LigandAlignmentAssembly {

    public:
        explicit LigandAlignmentAssembly(
            const std::unordered_map<LigandID, PoseID>& initialAssembly);

        /**
         * Exchange the associated pose for a given Ligand
         * @param ligandId The ligand whose pose is to be swapped.
         * @param newPoseId The new pose.
         */
        void swapPoseForLigand(LigandID ligandId,
                               PoseID newPoseId);

        /**
         * @param ligandId Ligand to get the associated pose for.
         * @return The Pose that is associated with the @p ligandId.
         */
        PoseID getPoseOfLigand(LigandID ligandId);

        /**
         * increase missing ligands count by one
         */
        void incrementMissingLigandsCount();

    private:
        LigandAlignmentAssembly()=default;

        bool insertLigandPose(LigandID ligand, PoseID pose);

        std::unordered_map<LigandID, PoseID> m_assembly;
        unsigned m_missingLigandsCount{0};

        friend class StartingAssemblyGenerator;
    };

}