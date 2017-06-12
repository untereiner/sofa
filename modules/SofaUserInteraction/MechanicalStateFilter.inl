/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, development version     *
*                (c) 2006-2015 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
//
// C++ Models: MechanicalStateFilter
//
// Description:
//
//
// Author:
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef SOFA_COMPONENT_CONTROLLER_MECHANICALSTATEFILTER_INL
#define SOFA_COMPONENT_CONTROLLER_MECHANICALSTATEFILTER_INL

#include <SofaUserInteraction/MechanicalStateFilter.h>

namespace sofa
{
namespace component
{
namespace controller
{
template <class DataTypes>
MechanicalStateFilter<DataTypes>::MechanicalStateFilter()
    : filteredIndices(initData(&filteredIndices, "filteredIndices", "indices of the constrained joints"))
    , range(initData(&range, "range", "range of freedom beetwen a max e min value for each articulation in index"))
    , controlledIndices(initData(&controlledIndices, "controlledIndices", "indices of the controlled joints"))
    , rotationAngle1(initData(&rotationAngle1, "rotationAngle1", "angles of the controlled joints"))
    , rotationAngle2(initData(&rotationAngle2, "rotationAngle2", "angles of the controlled joints"))
    , baseRotation(initData(&baseRotation, "baseRotation", "rotation angle of the base (first joint)"))
    , baseTranslation(initData(&baseTranslation, "baseTranslation", "translation of the base (first joint)"))
{

}

// -------------------------- Init() --------------------------------
template <class DataTypes>
void MechanicalStateFilter<DataTypes>::init()
{
    using core::behavior::MechanicalState;
    mState = dynamic_cast<MechanicalState<DataTypes> *> (this->getContext()->getMechanicalState());

    f_listening.setValue(true); // force event handling to avoid adding listening='true' in XML scene description
    //rotationAngle1 = 0.0;
    //rotationAngle2 = 0.0;
    if (!mState)
        serr << "MechanicalStateFilter has no binding MechanicalState" << sendl;
}

// -------------------------- onBeginAnimationStep() --------------------------------
template <class DataTypes>
void MechanicalStateFilter<DataTypes>::onBeginAnimationStep(const double /*dt*/)
{
    using sofa::defaulttype::Vec;
    double scaleFactor = 0.1;
    if (mState)
    {
        // helper::WriteAccessor<Data<VecCoord> > x0_write = *this->mState->write(core::VecCoordId::restPosition());
        helper::WriteAccessor<Data<VecCoord> > x_write = *this->mState->write(core::VecCoordId::position());
        helper::WriteAccessor<Data<VecDeriv> > v_write = *this->mState->write(core::VecDerivId::velocity());

        // set base rotation
        x_write[0].x() = baseTranslation.getValue() / 180.0 * M_PI;
        v_write[0].x() = 0;

        // set base translation
        x_write[1].x() = baseRotation.getValue() / 180.0 * M_PI;
        v_write[1].x() = 0;

        // set rotation for specific joints
        for (unsigned int i = 2*controlledIndices.getValue()[0]+2; i < x_write.size(); i+=2)
        {
            std::cout << "Angle of joint: " << i << " set to: " << (rotationAngle1.getValue() / 180.0 * M_PI) << std::endl;
            x_write[i].x() = scaleFactor*(rotationAngle1.getValue() / 180.0 * M_PI);
            v_write[i].x() = 0;
         }

        for (unsigned int i = 2*controlledIndices.getValue()[1]+3; i < x_write.size(); i+=2)
        {
            std::cout << "Angle of joint: " << i << " set to: " << (rotationAngle2.getValue() / 180.0 * M_PI) << std::endl;
            x_write[i].x() = scaleFactor*(rotationAngle2.getValue() / 180.0 * M_PI);
            v_write[i].x() = 0;
         }

        // filter angles out of range
        for (unsigned int i = 0; i < filteredIndices.getValue().size(); i++)
        {
            //std::cout << "Joint: " << filteredIndices.getValue()[i] << " = " <<  x_write[filteredIndices.getValue()[i]].x() * 180.0 / M_PI << " [Min: " << range.getValue()[2 * i] << " Max: " << range.getValue()[2 * i + 1] << "]" << std::endl;
            unsigned int constrainedIndex = filteredIndices.getValue()[i];
            if (constrainedIndex > 0) {
                //Min
                if (((x_write[constrainedIndex].x() - x_write[constrainedIndex-1].x()) * 180 / M_PI) < range.getValue()[2*i])
                {
                    std::cout << "    >>>> joint " << constrainedIndex << " blocked at " << (range.getValue()[2*i]) << std::endl;
                    x_write[constrainedIndex].x() = (range.getValue()[2*i]) * M_PI / 180.0;
                    v_write[constrainedIndex].x() = 0;
                }

                //Max
                else if (((x_write[constrainedIndex].x() - x_write[constrainedIndex-1].x()) * 180 / M_PI) > range.getValue()[2*i+1])
                {
                    std::cout << "    >>>> joint " << constrainedIndex << " blocked at " << (range.getValue()[2*i+1]) << std::endl;
                    x_write[constrainedIndex].x() = range.getValue()[2*i+1] * M_PI / 180.0;
                    v_write[constrainedIndex].x() = 0;
                }
            }
        }
    }

    sofa::simulation::Node *node = static_cast<sofa::simulation::Node*> (this->getContext());
    sofa::simulation::MechanicalPropagatePositionAndVelocityVisitor mechaVisitor(core::MechanicalParams::defaultInstance()); mechaVisitor.execute(node);
    // sofa::simulation::UpdateMappingVisitor updateVisitor(core::ExecParams::defaultInstance()); updateVisitor.execute(node);
}


} // namespace controller

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_CONTROLLER_MECHANICALSTATEFILTER_INL
