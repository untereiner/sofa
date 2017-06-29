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
// MechanicalRigidStateFilter
//
// Description: limits the angular motion of a rigid link
// It can be used in an articulated object or beam structure for instance
//
//
// Author: Stephane Cotin (2016)
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef SOFA_COMPONENT_CONTROLLER_MECHANICALRIGIDSTATEFILTER_INL
#define SOFA_COMPONENT_CONTROLLER_MECHANICALRIGIDSTATEFILTER_INL

#include <SofaUserInteraction/MechanicalRigidStateFilter.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/helper/Quater.h>
#include <sofa/simulation/Node.h>
#include <sofa/simulation/MechanicalVisitor.h>
#include <sofa/simulation/UpdateMappingVisitor.h>

namespace sofa
{

	namespace component
	{

		namespace controller
		{

			template <class DataTypes>
            MechanicalRigidStateFilter<DataTypes>::MechanicalRigidStateFilter()
                : indices(initData(&indices, "indices", "indices of the constrained nodes"))
                , range(initData(&range, "range", "range of freedom (max; min) values for each articulation"))
			{
                rotationAngle = 0.0;
			}

			template <class DataTypes>
            void MechanicalRigidStateFilter<DataTypes>::init()
			{
				using core::behavior::MechanicalState;
				mState = dynamic_cast<MechanicalState<DataTypes> *> (this->getContext()->getMechanicalState());
				f_listening.setValue(true);
				if (!mState)
                    serr << "MechanicalRigidStateFilter has no binding MechanicalState" << sendl;
			}

			template <class DataTypes>
            void MechanicalRigidStateFilter<DataTypes>::onBeginAnimationStep(const double /*dt*/)
			{
				using sofa::defaulttype::Vec;
                defaulttype::Vec<3,double> axis;
                double phi, alpha;
                sofa::defaulttype::Quat q_lim;
				if (mState)
				{
                    helper::WriteAccessor<Data<VecCoord> > x_write = *this->mState->write(core::VecCoordId::position());
                    helper::ReadAccessor<Data<VecCoord> > x0_read = *this->mState->read(core::VecCoordId::restPosition());
                    helper::WriteAccessor<Data<VecDeriv> > v_write = *this->mState->write(core::VecDerivId::velocity());

                    rotationAngle += 0.5 / 180.0 * M_PI;
                    axis.x()=0.0;
                    axis.y()=0.0;
                    axis.z()=1.0;
                    q_lim.axisToQuat(axis, rotationAngle);

                    // set rotation for specific joints
                    x_write[5].getOrientation() = q_lim;

                    for (unsigned int i = 0; i < indices.getValue().size(); i++)
					{
                        // compute deviation from initial orientation
                        unsigned int nodeIndex = indices.getValue()[i];
                        sofa::defaulttype::Quat q = x_write[nodeIndex].getOrientation();
                        sofa::defaulttype::Quat q0 = x0_read[nodeIndex].getOrientation();

                        q.quatToAxis(axis, phi);
                        phi = phi * 180.0 / M_PI;
                        std::cout << "main axis of rotation" << axis << "   - angle: " << phi << std::endl;

                        //Min
                        if (phi > range.getValue()[2 * i])
                        {
                            alpha = range.getValue()[2 * i] * M_PI / 180.0;
                            q_lim.axisToQuat(axis, alpha);
                            x_write[nodeIndex].getOrientation() = q_lim;
                        }
                        //Max
                        if (phi < range.getValue()[2 * i + 1])
						{
                            alpha = range.getValue()[2 * i + 1] * M_PI / 180.0;
                            q_lim.axisToQuat(axis, alpha);
                            x_write[nodeIndex].getOrientation() = q_lim;
                        }
                    }
				}

				sofa::simulation::Node *node = static_cast<sofa::simulation::Node*> (this->getContext());
				sofa::simulation::MechanicalPropagatePositionAndVelocityVisitor mechaVisitor(core::MechanicalParams::defaultInstance()); mechaVisitor.execute(node);
				sofa::simulation::UpdateMappingVisitor updateVisitor(core::ExecParams::defaultInstance()); updateVisitor.execute(node);
			}

			template <class DataTypes>
            core::behavior::MechanicalState<DataTypes> *MechanicalRigidStateFilter<DataTypes>::getMechanicalState() const
			{
				return mState;
			}

			template <class DataTypes>
            void MechanicalRigidStateFilter<DataTypes>::setMechanicalState(core::behavior::MechanicalState<DataTypes> *_mState)
			{
				mState = _mState;
			}


			template <class DataTypes>
            unsigned int MechanicalRigidStateFilter<DataTypes>::getIndex() const
			{
				return index.getValue();
			}


			template <class DataTypes>
            void MechanicalRigidStateFilter<DataTypes>::setIndex(const unsigned int _index)
			{
				index.setValue(_index);
			}

		} // namespace controller

	} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_CONTROLLER_MECHANICALRIGIDSTATEFILTER_INL
