//
// Created by chris on 11/5/23.
//

#pragma once
#include "Forward.hpp"
#include <unordered_map>

namespace MultiAlign {

    class LigandAlignmentAssembly {

    public:
        explicit LigandAlignmentAssembly(
            const std::unordered_map<LigandID, PoseID>& initialAssembly);

        void swapPoseForLigand(LigandID ligandId,
                               PoseID newPoseId);

        PoseID getPoseOfLigand(LigandID ligandId);
    private:
        std::unordered_map<LigandID, PoseID> m_assembly;
    };

}