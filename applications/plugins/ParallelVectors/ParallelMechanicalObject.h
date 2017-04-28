/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, development version     *
*                (c) 2006-2017 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_PARALLELMECHANICALOBJECT_H
#define SOFA_PARALLELMECHANICALOBJECT_H

#include "ParallelTypes.h"
#include <SofaBaseMechanics/MechanicalObject.h>

namespace sofa
{

namespace component
{

namespace container
{

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_CONTAINER_MECHANICALOBJECT_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_BASE_MECHANICS_API MechanicalObject<defaulttype::ParallelVec3dTypes>;
extern template class SOFA_BASE_MECHANICS_API MechanicalObject<defaulttype::ParallelVec2dTypes>;
extern template class SOFA_BASE_MECHANICS_API MechanicalObject<defaulttype::ParallelVec1dTypes>;
extern template class SOFA_BASE_MECHANICS_API MechanicalObject<defaulttype::ParallelVec6dTypes>;
//extern template class SOFA_BASE_MECHANICS_API MechanicalObject<defaulttype::ParallelRigid3dTypes>;
//extern template class SOFA_BASE_MECHANICS_API MechanicalObject<defaulttype::ParallelRigid2dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_BASE_MECHANICS_API MechanicalObject<defaulttype::ParallelVec3fTypes>;
extern template class SOFA_BASE_MECHANICS_API MechanicalObject<defaulttype::ParallelVec2fTypes>;
extern template class SOFA_BASE_MECHANICS_API MechanicalObject<defaulttype::ParallelVec1fTypes>;
extern template class SOFA_BASE_MECHANICS_API MechanicalObject<defaulttype::ParallelVec6fTypes>;
//extern template class SOFA_BASE_MECHANICS_API MechanicalObject<defaulttype::ParallelRigid3fTypes>;
//extern template class SOFA_BASE_MECHANICS_API MechanicalObject<defaulttype::ParallelRigid2fTypes>;
#endif
//extern template class SOFA_BASE_MECHANICS_API MechanicalObject<defaulttype::ParallelLaparoscopicRigid3Types>;
#endif



}

} // namespace component

} // namespace sofa



#endif //SOFA_PARALLELMECHANICALOBJECT_H
