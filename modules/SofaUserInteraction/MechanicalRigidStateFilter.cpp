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
// C++ Implementation : MechanicalRigidStateFilter
//
// Description:
//
//
// Author: ???????????????????????????????????????????????????????
//
// Copyright: See COPYING file that comes with this distribution
//
//
#define SOFA_COMPONENT_CONTROLLER_MECHANICALRIGIDSTATEFILTER_CPP
#include <SofaUserInteraction/MechanicalRigidStateFilter.inl>
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>


namespace sofa
{

namespace component
{

namespace controller
{

using namespace sofa::defaulttype;

SOFA_DECL_CLASS(MechanicalRigidStateFilter)

// Register in the Factory
int MechanicalRigidStateFilterClass = core::RegisterObject("Filter articulation position in a mechanical object")
#ifndef SOFA_FLOAT
//.add< MechanicalRigidStateFilter<Vec3dTypes> >()
//.add< MechanicalRigidStateFilter<Vec2dTypes> >()
        .add< MechanicalRigidStateFilter<Rigid3dTypes> >()
        
#endif
#ifndef SOFA_DOUBLE
//.add< MechanicalRigidStateFilter<Vec3fTypes> >()
//.add< MechanicalRigidStateFilter<Vec2fTypes> >()
        .add< MechanicalRigidStateFilter<Rigid3fTypes> >()
#endif
        ;

#ifndef SOFA_FLOAT
//template class SOFA_USER_INTERACTION_API MechanicalRigidStateFilter<Vec3dTypes>;
//template class SOFA_USER_INTERACTION_API MechanicalRigidStateFilter<Vec2dTypes>;
template class SOFA_USER_INTERACTION_API MechanicalRigidStateFilter<Rigid3dTypes>;

#endif
#ifndef SOFA_DOUBLE
//template class SOFA_USER_INTERACTION_API MechanicalRigidStateFilter<Vec3fTypes>;
//template class SOFA_USER_INTERACTION_API MechanicalRigidStateFilter<Vec2fTypes>;
template class SOFA_USER_INTERACTION_API MechanicalRigidStateFilter<Rigid3fTypes>;

#endif


} // namespace controller

} // namespace component

} // namespace sofa
