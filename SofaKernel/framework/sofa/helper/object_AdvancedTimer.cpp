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

#define SOFA_HELPER_OBJECT_ADVANCEDTIMER_CPP

#include <sofa/helper/object_AdvancedTimer.inl>
#include <sofa/core/ObjectFactory.h>


namespace sofa
{

namespace helper
{


SOFA_DECL_CLASS(Object_AdvancedTimer)

// Register in factory

int Object_AdvancedTimerClass = core::RegisterObject("Object used for advencedTimer python usage")
        .add<Object_AdvancedTimer>();

} // Namespace helper

} // Namespace sofa
