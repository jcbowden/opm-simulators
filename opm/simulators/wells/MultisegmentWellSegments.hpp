/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef OPM_MULTISEGMENTWELL_SEGMENTS_HEADER_INCLUDED
#define OPM_MULTISEGMENTWELL_SEGMENTS_HEADER_INCLUDED

#include <opm/simulators/wells/MultisegmentWellPrimaryVariables.hpp>

#include <vector>

namespace Opm
{

class UnitSystem;
class WellInterfaceGeneric;

template<typename FluidSystem, typename Indices, typename Scalar>
class MultisegmentWellSegments
{
    using PrimaryVariables = MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>;
    using EvalWell = typename PrimaryVariables::EvalWell;

public:
    MultisegmentWellSegments(const int numSegments,
                             WellInterfaceGeneric& well);

    void computeFluidProperties(const EvalWell& temperature,
                                const EvalWell& saltConcentration,
                                const PrimaryVariables& primary_variables,
                                int pvt_region_index,
                                DeferredLogger& deferred_logger);

    //! \brief Update upwinding segments.
    void updateUpwindingSegments(const PrimaryVariables& primary_variables);

    EvalWell getHydroPressureLoss(const int seg) const;

    EvalWell getSurfaceVolume(const EvalWell& temperature,
                              const EvalWell& saltConcentration,
                              const PrimaryVariables& primary_variables,
                              const int pvt_region_index,
                              const int seg_idx) const;

    EvalWell getFrictionPressureLoss(const int seg) const;

    // pressure drop for Spiral ICD segment (WSEGSICD)
    EvalWell pressureDropSpiralICD(const int seg) const;

    // pressure drop for Autonomous ICD segment (WSEGAICD)
    EvalWell pressureDropAutoICD(const int seg,
                                 const UnitSystem& unit_system) const;

    // pressure drop for sub-critical valve (WSEGVALV)
    EvalWell pressureDropValve(const int seg) const;

    // pressure loss due to acceleration
    EvalWell accelerationPressureLoss(const int seg) const;

    // TODO: trying to use the information from the Well opm-parser as much
    // as possible, it will possibly be re-implemented later for efficiency reason.

    // the completions that is related to each segment
    // the completions's ids are their index in the vector well_index_, well_cell_
    // This is also assuming the order of the completions in Well is the same with
    // the order of the completions in wells.
    // it is for convinience reason. we can just calcuate the inforation for segment once then using it for all the perofrations
    // belonging to this segment
    std::vector<std::vector<int>> perforations_;

    // depth difference between the segment and the perforation
    // or in another way, the depth difference between the perforation and
    // the segment the perforation belongs to
    std::vector<double> perforation_depth_diffs_;

    // the inlet segments for each segment. It is for convenience and efficiency reason
    std::vector<std::vector<int>> inlets_;

    std::vector<double> depth_diffs_;

    // the densities of segment fluids
    // we should not have this member variable
    std::vector<EvalWell> densities_;

    // the mass rate of the segments
    std::vector<EvalWell> mass_rates_;

    // the viscosity of the segments
    std::vector<EvalWell> viscosities_;

    // the upwinding segment for each segment based on the flow direction
    std::vector<int> upwinding_segments_;

    std::vector<std::vector<EvalWell>> phase_densities_;
    std::vector<std::vector<EvalWell>> phase_fractions_;
    std::vector<std::vector<EvalWell>> phase_viscosities_;

private:
    const WellInterfaceGeneric& well_;
};

}

#endif // OPM_MULTISEGMENTWELL_SEGMENTS_HEADER_INCLUDED
