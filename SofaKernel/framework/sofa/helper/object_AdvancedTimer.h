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

/**
 * Class : Object-AdvancedTimer
 * Brief : This class is used to create an object AdvacedTimer in the scene graph.
 *         It uses the Binding_BaseContext::createObject to be available in sofaPython.
 * Method : startTimer - return void. Used to start the timer after creating it.
 **/

#ifndef SOFA_HELPER_OBJECT_ADVANCEDTIMER_H
#define SOFA_HELPER_OBJECT_ADVANCEDTIMER_H

#include <sofa/core/objectmodel/BaseObject.h>


namespace sofa
{

namespace helper
{

class Object_AdvancedTimer : public virtual core::objectmodel::BaseContext
{
// -----------------------
// Methods
public:
    SOFA_CLASS(Object_AdvancedTimer, sofa::helper::Object_AdvancedTimer);
    std::string timerId = "";
    int timerTimeStamp = 0;

// -----------------------
// Methods
public:
    // Default constructor with id = "Animate", enabled = true and timestamp = 1000
    Object_AdvancedTimer();


    /**
     * @brief Object_AdvancedTimer : constructor with paramters
     * @param id, std::string - name of the timer
     * @param enabled, bool - boolean passed to setEnabled method
     * @param timestamp, int - time stamp for the timer
     **/
    Object_AdvancedTimer(std::string id, bool enabled, int timestamp);


    /**
     * @brief getID : return the id of the actual timer
     * @return std::string, the id of the actual timer
     **/
    std::string getID();


    /**
     * @brief setID Set the new id of the timer. Stop the old one unless keepOldTimerActive is set to true.
     * @param keepOldTimerActive bool - it equal to true, then the old timer is kept active and its id is return.
     * @return id of the old timer if keepOldTimerActive is egual to true else null
     */
    std::string setID(std::string newID, bool keepOldTimerActive);


    /**
     * @brief startTimer Start the timer with timerId id and with timerTimeStamp time stamp
     **/
    void startTimer();


    /**
     * @brief stopTimer Stop the timer with timeID id
     **/
    void stopTimer();

    // TO-DO :
    // Une méthode permettant de récupérer uniquement le champ voulu dans un fichier pour en faire un graphe
    // Une méthode permettant de stopper le timer, de le relancer, etc...

};

} // helper

} // sofa

#endif
