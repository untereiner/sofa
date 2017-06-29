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
// C++ Interface: MechanicalStateFilter
//
// Description:
//
//
// Author: S. Cotin (Inria)
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef SOFA_COMPONENT_CONTROLLER_MECHANICALSTATEFILTER_H
#define SOFA_COMPONENT_CONTROLLER_MECHANICALSTATEFILTER_H
#include "config.h"

#include <SofaUserInteraction/Controller.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <SofaBaseTopology/TopologySubsetData.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/simulation/Node.h>
#include <sofa/core/objectmodel/KeypressedEvent.h>
#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/simulation/AnimateEndEvent.h>
#include <sofa/simulation/MechanicalVisitor.h>
#include <sofa/simulation/UpdateMappingVisitor.h>

namespace sofa
{

namespace component
{

namespace controller
{


/**
 * @brief MechanicalStateFilter Class
 *
 * Provides functionalities for controlling and limiting the angular motion
 * of an articulated rigid body system
 */
template<class DataTypes>
class MechanicalStateFilter : public Controller
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(MechanicalStateFilter,DataTypes),Controller);
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord    Coord   ;
    typedef typename DataTypes::Deriv    Deriv   ;
    typedef typename Coord::value_type   Real    ;
    typedef sofa::helper::vector<int>    VecInt  ;
    typedef sofa::helper::vector<float>  VecFloat;

    Data<VecInt> filteredIndices;
    Data<VecFloat> range;
    Data<VecInt> controlledIndices;
    Data <Real> rotationAngle1;
    Data <Real> rotationAngle2;
    Data <Real> baseRotation;
    Data <Real> baseTranslation;

protected:
    /**
     * @brief Default Constructor.
     */
    MechanicalStateFilter();

    /**
     * @brief Default Destructor.
     */
    virtual ~MechanicalStateFilter() {}

    // Data< unsigned int > index; ///< Controlled DOF index.
    core::behavior::MechanicalState<DataTypes> *mState; ///< Controlled MechanicalState.
    sofa::defaulttype::Vector3 position;

public:
    /**
     * @brief SceneGraph callback initialization method.
     */
    void init();

  
    /**
     * @brief Begin Animation event callback.
     */
    void onBeginAnimationStep(const double dt);

    //@}

    /**
     * @brief End Animation event callback.
     */
    //void onEndAnimationStep(const double dt);

    //@}

    /**
     * @brief Event (keyborad, mouse) handler callback.
     */
//    void handleEvent(core::objectmodel::Event *event);

    //@}

    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const MechanicalStateFilter<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_CONTROLLER_MECHANICALSTATEFILTER_CPP)
#ifndef SOFA_FLOAT
//extern template class SOFA_USER_INTERACTION_API MechanicalStateFilter<defaulttype::Vec3dTypes>;
//extern template class SOFA_USER_INTERACTION_API MechanicalStateFilter<defaulttype::Vec2dTypes>;
extern template class SOFA_USER_INTERACTION_API MechanicalStateFilter<defaulttype::Vec1dTypes>;


#endif
#ifndef SOFA_DOUBLE
//extern template class SOFA_USER_INTERACTION_API MechanicalStateFilter<defaulttype::Vec3fTypes>;
//extern template class SOFA_USER_INTERACTION_API MechanicalStateFilter<defaulttype::Vec2fTypes>;
extern template class SOFA_USER_INTERACTION_API MechanicalStateFilter<defaulttype::Vec1fTypes>;

#endif
#endif

} // namespace controller

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_CONTROLLER_MECHANICALSTATEFILTER_H
