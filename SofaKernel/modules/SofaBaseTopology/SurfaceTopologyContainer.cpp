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
#include <SofaBaseTopology/SurfaceTopologyContainer.h>
#include <sofa/core/ObjectFactory.h>
#include <cgogn/io/surface_import.h>

namespace sofa
{

namespace component
{

namespace topology
{

SOFA_DECL_CLASS(SurfaceTopologyContainer)
int SurfaceTopologyContainerClass = core::RegisterObject("Surface topology container")
        .add< SurfaceTopologyContainer >()
        ;

SurfaceTopologyContainer::SurfaceTopologyContainer()
{

}

SurfaceTopologyContainer::~SurfaceTopologyContainer()
{

}

void SurfaceTopologyContainer::initFromMeshLoader()
{
	helper::ReadAccessor< Data< VecCoord > > m_position = d_initPoints;
	helper::ReadAccessor< Data< helper::vector< TriangleIds > > > m_tri = d_triangle;
	helper::ReadAccessor< Data< helper::vector< QuadIds > > > m_quad = d_quad;

	cgogn::io::SurfaceImport<Topo_Traits::MapTraits> surface_import;
	surface_import.reserve(m_tri.size() + m_quad.size());

	auto* pos_att = surface_import.template position_attribute<Eigen::Vector3d>();
	for(std::size_t i = 0ul, end = m_position.size(); i < end ; ++i)
	{
		const auto id = surface_import.insert_line_vertex_container();
		const auto& src = m_position[i];
		auto& dest = pos_att->operator [](id);
		dest[0] = src[0];
		dest[1] = src[1];
		dest[2] = src[2];
	}

	for(const TriangleIds& t : m_tri.ref())
		surface_import.add_triangle(t[0], t[1], t[2]);
	for(const QuadIds& q : m_quad.ref())
		surface_import.add_quad(q[0], q[1], q[2], q[3]);

	surface_import.create_map(topology_);
}

void SurfaceTopologyContainer::init()
{
	topology_.clear_and_remove_attributes();
	Inherit1::init();
	initFromMeshLoader();

	if (!this->edge_dofs_.is_valid())
	{
		this->edge_dofs_ = this->template add_attribute<helper::fixed_array<unsigned int, 2>, Edge::ORBIT>("edge_dofs_");
		this->parallel_foreach_cell([&](Edge e, cgogn::uint32)
		{
			this->edge_dofs_[e.dart] = helper::fixed_array<unsigned int, 2>(get_dof(Vertex(e.dart)), get_dof(Vertex(phi2(e.dart))));
		});
	}


	if (!this->face_dofs_.is_valid())
	{
		this->face_dofs_ = this->template add_attribute<helper::fixed_array<unsigned int, 4>, Face::ORBIT>("face_dofs_");
		this->parallel_foreach_cell([&](Face f, cgogn::uint32)
		{
			auto & dofs = this->face_dofs_[f.dart];
			unsigned int i = 0u;
			foreach_incident_vertex(f, [&](Vertex v)
			{
				dofs[i++] = get_dof(v);
			});
		});
	}
}

void SurfaceTopologyContainer::bwdInit()
{
	Inherit1::bwdInit();
}

void SurfaceTopologyContainer::reinit()
{
	Inherit1::reinit();
}

void SurfaceTopologyContainer::reset()
{
//	topology_.clear_and_remove_attributes(); // reset() seems to be called after init() at the beginning of the scene (?!)
	Inherit1::reset();
}

void SurfaceTopologyContainer::cleanup()
{
	Inherit1::cleanup();
}

} // namespace topology

} // namespace component

} // namespace sofa

