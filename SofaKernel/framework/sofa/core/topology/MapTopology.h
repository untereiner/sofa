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
*                              SOFA :: Framework                              *
*                                                                             *
* Authors: The SOFA Team (see Authors.txt)                                    *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_CORE_TOPOLOGY_MAPTOPOLOGY_H
#define SOFA_CORE_TOPOLOGY_MAPTOPOLOGY_H

#include <functional>
#include <sofa/core/topology/BaseTopology.h>
#include <sofa/core/State.h>
#include <cgogn/core/cmap/cmap3.h>

namespace sofa
{

namespace core
{

namespace topology
{

struct CGOGN_Traits
{
	using index_type = cgogn::Dart;
	using MapTraits = cgogn::DefaultMapTraits;
	using Topology2 = cgogn::CMap2<MapTraits>;
	using Topology3 = cgogn::CMap3<MapTraits>;

	using Vertex2	= Topology2::Vertex;
	using Edge2		= Topology2::Edge;
	using Face2		= Topology2::Face;
	using Volume2	= Topology2::Volume;

	using Vertex3	= Topology3::Vertex;
	using Edge3		= Topology3::Edge;
	using Face3		= Topology3::Face;
	using Volume3	= Topology3::Volume;

	template<typename T>
	using Attribute_T = cgogn::Attribute_T<MapTraits,T>;

	struct Vertex
	{
		inline Vertex() : id_() {}
		inline Vertex(index_type id) : id_(id) {}
		index_type id_;
		inline operator index_type() const
		{
			return id_;
		}
	};

	struct Edge
	{
		inline Edge() : id_() {}
		inline Edge(index_type id) : id_(id) {}
		index_type id_;
		inline operator index_type() const
		{
			return id_;
		}
	};

	struct Face
	{
		inline Face() : id_() {}
		inline Face(index_type id) : id_(id) {}
		index_type id_;
		inline operator index_type() const
		{
			return id_;
		}
	};

	struct Volume
	{
		inline Volume() : id_() {}
		inline Volume(index_type id) : id_(id) {}
		index_type id_;
		inline operator index_type() const
		{
			return id_;
		}
	};
};

class SOFA_CORE_API MapTopology : public TopologyContainer
{
public:
	SOFA_CLASS(MapTopology,core::topology::TopologyContainer);

	using Topo_Traits = CGOGN_Traits;
	using Dart = cgogn::Dart;
	using Orbit = cgogn::Orbit;
	template<typename T>
	using Attribute_T = CGOGN_Traits::Attribute_T<T>;

	using Vertex = Topo_Traits::Vertex;
	using Edge = Topo_Traits::Edge;
	using Face = Topo_Traits::Face;
	using Volume = Topo_Traits::Volume;

	using Vec3Types = sofa::defaulttype::Vec3Types;
	using VecCoord = Vec3Types::VecCoord;

	MapTopology();

	virtual void foreach_vertex(std::function<void(Vertex)> const &) = 0;

	virtual void foreach_edge(std::function<void(Edge)> const &) = 0;

	virtual void foreach_face(std::function<void(Face)> const &) = 0;

	virtual void foreach_volume(std::function<void(Volume)> const &) = 0;

	virtual void foreach_incident_vertex_of_edge(Edge /*edge_id*/, std::function<void(Vertex)> const & /*func*/) = 0;

	virtual void foreach_incident_vertex_of_face(Face /*face_id*/, std::function<void(Vertex)> const & /*func*/) = 0;

	virtual void foreach_incident_vertex_of_volume(Volume /*w*/, std::function<void(Vertex)> const & /*func*/) = 0;

	virtual void foreach_incident_edge_of_face(Face /*f_id*/, std::function<void(Edge)> const & /*func*/) = 0;

	virtual void foreach_incident_edge_of_volume(Volume /*vol_id*/, std::function<void(Edge)> const & /*func*/) = 0;

	virtual void foreach_incident_face_of_volume(Volume /*vol_id*/, std::function<void(Face)> const & /*func*/) = 0;

	virtual void init() override;
	virtual void bwdInit() override;
	virtual void reinit() override;
	virtual void reset() override;
	virtual void cleanup() override;

protected:
	virtual void initFromMeshLoader() = 0;

	Data< VecCoord > d_initPoints;
	Data< helper::vector< Triangle > > d_triangle;
	Data< helper::vector< Quad > > d_quad;
	Data< helper::vector< Tetra > > d_tetra;
	Data< helper::vector< Hexa > > d_hexa;

	SingleLink< MapTopology, core::State< Vec3Types >, BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK > mech_state_;
private:
	Attribute_T<Dart> first_vertex_of_edge;
	Attribute_T<Dart> first_vertex_of_face;
	Attribute_T<Dart> first_vertex_of_volume;
};


} // namespace topology

} // namespace core

} // namespace sofa

#endif // SOFA_CORE_TOPOLOGY_MAPTOPOLOGY_H
