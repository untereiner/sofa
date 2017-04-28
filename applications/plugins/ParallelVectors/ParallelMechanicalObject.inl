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
#ifndef SOFA_PARALLELMECHANICALOBJECT_INL
#define SOFA_PARALLELMECHANICALOBJECT_INL

#include <SofaBaseMechanics/MechanicalObject.inl>
#include "ParallelMechanicalObject.h"


namespace sofa
{

namespace component
{

namespace container
{

//template <class DataTypes>
//bool MechanicalObject<DataTypes>::addBBox(SReal* minBBox, SReal* maxBBox)
//{
//    // participating to bbox only if it is drawn
//    if( !showObject.getValue() ) return false;
//
//    static const unsigned spatial_dimensions = std::min( (unsigned)DataTypes::spatial_dimensions, 3u );
//
//    const VecCoord& x = read(core::ConstVecCoordId::position())->getValue();
//    for( std::size_t i=0; i<x.size(); i++ )
//    {
//        defaulttype::Vec<3,Real> p;
//        DataTypes::get( p[0], p[1], p[2], x[i] );
//
//        for( unsigned int j=0 ; j<spatial_dimensions; ++j )
//        {
//            if(p[j]<minBBox[j]) minBBox[j]=p[j];
//            if(p[j]>maxBBox[j]) maxBBox[j]=p[j];
//        }
//    }
//    return true;
//}
//
//
//template <class DataTypes>
//void MechanicalObject<DataTypes>::computeBBox(const core::ExecParams* params, bool onlyVisible)
//{
//    // participating to bbox only if it is drawn
//    if( onlyVisible && !showObject.getValue() ) return;
//    Inherited::computeBBox( params );
//}
//
//template <class DataTypes>
//bool MechanicalObject<DataTypes>::isIndependent() const
//{
//    return static_cast<const simulation::Node*>(this->getContext())->mechanicalMapping.empty();
//}


} // namespace container

} // namespace component

} // namespace sofa

#endif  //SOFA_PARALLELMECHANICALOBJECT_INL

