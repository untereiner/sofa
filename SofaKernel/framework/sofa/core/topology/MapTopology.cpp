/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, development version     *
*                (c) 2006-2016 INRIA, USTL, UJF, CNRS, MGH                    *
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
#include <sofa/core/topology/MapTopology.h>

namespace sofa
{

namespace core
{

namespace topology
{

void MapTopology::cleanup()
{
	Inherit1::cleanup();
}

MapTopology::MapTopology() :
	Inherit1(),
	d_initPoints (initData(&d_initPoints, "position", "Initial position of points")),
	d_triangle(initData(&d_triangle, "triangles", "List of triangle indices")),
	d_quad(initData(&d_quad, "quads", "List of quad indices")),
	d_tetra(initData(&d_tetra, "tetrahedra", "List of tetrahedron indices")),
	d_hexa(initData(&d_hexa, "hexahedra", "List of hexahedron indices")),
	mech_state_(initLink("mstate", "mechanical state linked to the topology"))
{

}

void MapTopology::init()
{
	Inherit1::init();
}

void MapTopology::reinit()
{
	Inherit1::reinit();
}

void MapTopology::reset()
{
	Inherit1::reset();
}

void MapTopology::bwdInit()
{
	Inherit1::bwdInit();
}

} // namespace topology

} // namespace core

} // namespace sofa
