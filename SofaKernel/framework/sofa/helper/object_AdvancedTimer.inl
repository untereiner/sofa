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

#ifndef SOFA_HELPER_OBJECT_ADVANCEDTIMER_INL
#define SOFA_HELPER_OBJECT_ADVANCEDTIMER_INL

#include <sofa/helper/object_AdvancedTimer.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/core/objectmodel/BaseObject.h>

namespace sofa
{

namespace helper
{

    Object_AdvancedTimer::Object_AdvancedTimer()
    {
        timerId = "Animate";
        timerTimeStamp = 1000;
        AdvancedTimer::setEnabled(timerId, true);
        AdvancedTimer::setInterval(timerId, timerTimeStamp);
    }


    Object_AdvancedTimer::Object_AdvancedTimer(std::string id, bool enabled, int timestamp)
    {
        timerId = id;
        timerTimeStamp = timestamp;
        AdvancedTimer::setEnabled(id, enabled);
        AdvancedTimer::setInterval(id, timestamp);
    }


    std::string Object_AdvancedTimer::getID()
    {
        return timerId;
    }


    std::string Object_AdvancedTimer::setID(std::string newID, bool keepOldTimerActive)
    {
        if(keepOldTimerActive)
        {
            std::string oldTimerID = timerId;
            timerId = newID;
            AdvancedTimer::setEnabled(timerId, true);
            AdvancedTimer::setInterval(timerId, timerTimeStamp);
            return oldTimerID;
        }
        else
        {
            AdvancedTimer::setEnabled(timerId, false);
            timerId = newID;
            AdvancedTimer::setEnabled(timerId, true);
            AdvancedTimer::setInterval(timerId, timerTimeStamp);
            return NULL;
        }
    }



    void Object_AdvancedTimer::startTimer()
    {
        AdvancedTimer::setEnabled(timerId, true);
        AdvancedTimer::setInterval(timerId, timerTimeStamp);
    }


    void Object_AdvancedTimer::stopTimer()
    {
        AdvancedTimer::setEnabled(timerId, false);
    }



} // helper

} // sofa

#endif
