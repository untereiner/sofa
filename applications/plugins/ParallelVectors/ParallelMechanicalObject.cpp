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
#ifndef SOFA_PARALLELMECHANICALOBJECT_CPP
#define SOFA_PARALLELMECHANICALOBJECT_CPP

//#include <SofaBaseMechanics/MechanicalObject.inl>
#include "ParallelTypes.h"
#include "ParallelMechanicalObject.h"
//#include <sofa/helper/Quater.h>
#include <sofa/core/ObjectFactory.h>


namespace sofa
{

namespace component
{

namespace container
{

using namespace core::behavior;
using namespace defaulttype;


SOFA_DECL_CLASS(MechanicalObject)

int ParallelMechanicalObjectClass = core::RegisterObject("parallel mechanical state vectors")
#ifdef SOFA_FLOAT
        .add< MechanicalObject<ParallelVec3fTypes> >(true) // default template
#else
        .add< MechanicalObject<ParallelVec3dTypes> >(true) // default template
#ifndef SOFA_DOUBLE
        .add< MechanicalObject<ParallelVec3fTypes> >()
#endif
#endif
#ifndef SOFA_FLOAT
        .add< MechanicalObject<ParallelVec2dTypes> >()
        .add< MechanicalObject<ParallelVec1dTypes> >()
        .add< MechanicalObject<ParallelVec6dTypes> >()
//        .add< MechanicalObject<ParallelRigid3dTypes> >()
//        .add< MechanicalObject<ParallelRigid2dTypes> >()
#endif
#ifndef SOFA_DOUBLE
        .add< MechanicalObject<ParallelVec2fTypes> >()
        .add< MechanicalObject<ParallelVec1fTypes> >()
        .add< MechanicalObject<ParallelVec6fTypes> >()
//        .add< MechanicalObject<ParallelRigid3fTypes> >()
//        .add< MechanicalObject<ParallelRigid2fTypes> >()
#endif
//        .add< MechanicalObject<ParallelLaparoscopicRigid3Types> >()
        ;


} // namespace container

} // namespace component

} // namespace sofa

#endif // SOFA_PARALLELMECHANICALOBJECT_CPP
