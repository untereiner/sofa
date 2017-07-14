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

#include "Binding_AdvancedTimer.h"

using namespace sofa::helper;


/**
 * Method : AdvancedTimer_clear
 * Desc   : Wrapper for python usage. Clear the timer.
 * Param  : PyObject*, self - Object of the python script
 * Return : Py_RETURN_NONE
 */
extern "C" PyObject * AdvancedTimer_clear(PyObject *self, PyObject * /*args*/)
{
    AdvancedTimer::clear();  // Method call
    Py_RETURN_NONE;
}


/**
 * Method : AdvancedTimer_isEnabled
 * Desc   : Wrapper for python usage. Return if the timer is enable or not.
 * Param  : PyObject*, self - Object of the python script
 * Param  : PyObject*, args - given arguments to apply to the method
 * Return : Py_RETURN_NONE
 */
extern "C" PyObject * AdvancedTimer_isEnabled(PyObject *self, PyObject *args)
{
    unsigned int id = 0;
    bool answer = false;

    if(!PyArg_ParseTuple(args, "I", id))
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }

    answer = AdvancedTimer::isEnabled(id);  // Method call

    if(answer)
    {
        Py_RETURN_TRUE;
    }
    Py_RETURN_FALSE;
}


/**
 * Method : AdvancedTimer_setEnabled
 * Desc   : Wrapper for python usage. /!\ Need to pass an int in arguments insteed of a bool in the python script.
 * Param  : PyObject*, self - Object of the python script
 * Param  : PyObject*, args - given arguments to apply to the method
 * Return : Py_RETURN_NONE
 */
extern "C" PyObject * AdvancedTimer_setEnabled(PyObject *self, PyObject *args)
{
    std::string id = "";
    int tempBool = 0;
    bool val = false;

    if(!PyArg_ParseTuple(args, "si", &id, &tempBool))
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }

    if(tempBool == 1)
    {
        val = true;
    }

    AdvancedTimer::setEnabled(id, val);  // Method call
    Py_RETURN_NONE;
}


/**
 * Method : AdvancedTimer_getInterval
 * Desc   : Wrapper for python usage.
 * Param  : PyObject*, self - Object of the python script
 * Param  : PyObject*, args - given arguments to apply to the method
 * Return : Py_RETURN_NONE
 */
extern "C" PyObject * AdvancedTimer_getInterval(PyObject *self, PyObject *args)
{
    std::string id = "";
    int answer = 0;

    if(!PyArg_ParseTuple(args, "s", &id))
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }

    answer = AdvancedTimer::getInterval(id);  // Method call

    return PyInt_FromLong(static_cast<long int>(answer));
}


/**
 * Method : AdvancedTimer_setInterval
 * Desc   : Wrapper for python usage.
 * Param  : PyObject*, self - Object of the python script
 * Param  : PyObject*, args - given arguments to apply to the method
 * Return : Py_RETURN_NONE
 */
extern "C" PyObject * AdvancedTimer_setInterval(PyObject *self, PyObject *args)
{
    std::string id = "";
    int newValue = 0;

    if(!PyArg_ParseTuple(args, "si", &id, &newValue))
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }

    AdvancedTimer::setInterval(id, newValue);  // Method call

    Py_RETURN_NONE;
}


/**
 * Method : AdvancedTimer_begin
 * Desc   : Wrapper for python usage.
 * Param  : PyObject*, self - Object of the python script
 * Param  : PyObject*, args - given arguments to apply to the method
 * Return : Py_RETURN_NONE
 */
extern "C" PyObject * AdvancedTimer_begin(PyObject *self, PyObject *args)
{
    std::string id = "";

    if(!PyArg_ParseTuple(args, "s", &id))
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }

    AdvancedTimer::begin(id);  // Method call

    Py_RETURN_NONE;
}


/**
 * Method : AdvancedTimer_end
 * Desc   : Wrapper for python usage.
 * Param  : PyObject*, self - Object of the python script
 * Param  : PyObject*, args - given arguments to apply to the method
 * Return : Py_RETURN_NONE
 */
extern "C" PyObject * AdvancedTimer_end(PyObject *self, PyObject *args)
{
    std::string id = "";

    if(!PyArg_ParseTuple(args, "s", &id))
    {
        PyErr_BadArgument();
        Py_RETURN_NONE;
    }

    AdvancedTimer::end(id);  // Method call

    Py_RETURN_NONE;
}



SP_CLASS_METHODS_BEGIN(AdvancedTimer)
SP_CLASS_METHOD(AdvancedTimer, clear)
SP_CLASS_METHOD(AdvancedTimer, isEnabled)
SP_CLASS_METHOD(AdvancedTimer, setEnabled)
SP_CLASS_METHOD(AdvancedTimer, getInterval)
SP_CLASS_METHOD(AdvancedTimer, setInterval)
SP_CLASS_METHOD(AdvancedTimer, begin)
SP_CLASS_METHOD(AdvancedTimer, end)
SP_CLASS_METHODS_END

//SP_CLASS_TYPE_SPTR(AdvancedTimer, AdvancedTimer, AdvancedTimer)

