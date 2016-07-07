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
#include <SofaBaseTopology/VolumeTopologyContainer.h>

#include <cgogn/io/volume_import.h>
#include <cgogn/geometry/types/eigen.h>

namespace sofa
{

namespace component
{

namespace topology
{

void VolumeTopologyContainer::initFromMeshLoader()
{
	helper::ReadAccessor< Data< VecCoord > > m_position = d_initPoints;
	helper::ReadAccessor< Data< helper::vector< Tetra > > > m_tetra = d_tetra;
	helper::ReadAccessor< Data< helper::vector< Hexa > > > m_hexa = d_hexa;

	cgogn::io::VolumeImport<Topo_Traits::MapTraits> volume_import;
	volume_import.set_nb_vertices(m_position.size());

	auto* pos_att = volume_import.template position_attribute<Eigen::Vector3d>();
	for(std::size_t i = 0ul, end = m_position.size(); i < end ; ++i)
	{
		const auto id = volume_import.insert_line_vertex_container();
		const auto& src = m_position[i];
		auto& dest = pos_att->operator [](id);
		dest[0] = src[0];
		dest[1] = src[1];
		dest[2] = src[2];
	}

	for(const Tetra& t : m_tetra.ref())
		volume_import.add_tetra(*pos_att, t[0], t[1], t[2], t[3], true);
	for(const Hexa& t : m_hexa.ref())
		volume_import.add_hexa(*pos_att, t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], true);

	volume_import.create_map(topology_);
}

void VolumeTopologyContainer::init()
{
	topology_.clear_and_remove_attributes();
	Inherit1::init();
	initFromMeshLoader();
}

void VolumeTopologyContainer::bwdInit()
{
	Inherit1::bwdInit();
}

void VolumeTopologyContainer::reinit()
{
	Inherit1::reinit();
}

void VolumeTopologyContainer::reset()
{
	topology_.clear_and_remove_attributes();
	Inherit1::reset();
}

void VolumeTopologyContainer::cleanup()
{
	Inherit1::cleanup();
}


} // namespace topology

} // namespace component

} // namespace sofa
