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
#include <cgogn/core/basic/dart_marker.h>
#include <cgogn/core/basic/cell_marker.h>
namespace sofa
{

namespace core
{

namespace cm_topology
{

/// The enumeration used to give unique identifiers to Topological objects.
enum TopologyObjectType
{
	VERTEX,
	EDGE,
	FACE,
	VOLUME
};

} // namespace cm_topology

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

class SOFA_CORE_API MapTopology : virtual public TopologyContainer
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
	~MapTopology() override;

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
	Data< helper::vector< core::topology::Topology::Edge > > d_edge;
	Data< helper::vector< Triangle > > d_triangle;
	Data< helper::vector< Quad > > d_quad;
	Data< helper::vector< Tetra > > d_tetra;
	Data< helper::vector< Hexa > > d_hexa;

	SingleLink< MapTopology, core::State< Vec3Types >, BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK > mech_state_;
private:
	Attribute_T<Dart> first_vertex_of_edge_;
	Attribute_T<Dart> first_vertex_of_face_;
	Attribute_T<Dart> first_vertex_of_volume_;

	// compatibility
protected:

	Attribute_T<EdgesAroundVertex>			m_edgesAroundVertex;
	Attribute_T<EdgesInQuad>				m_edgesInQuad;
	Attribute_T<QuadsAroundEdge>			m_quadsAroundEdge;
	Attribute_T<QuadsAroundVertex>			m_quadsAroundVertex;
	Attribute_T<QuadsAroundVertex>			m_orientedQuadsAroundVertex;
	Attribute_T<EdgesAroundVertex>			m_orientedEdgesAroundVertex;
	Attribute_T<TrianglesInTetrahedron>		m_trianglesInTetrahedron;
	Attribute_T<EdgesInHexahedron>			m_edgesInHexahedron;
	Attribute_T<EdgesInTetrahedron>			m_edgesInTetrahedron;
	Attribute_T<QuadsInHexahedron>			m_quadsInHexahedron;
	Attribute_T<TetrahedraAroundVertex>		m_tetrahedraAroundVertex;
	Attribute_T<TetrahedraAroundEdge>		m_tetrahedraAroundEdge;
	Attribute_T<TetrahedraAroundTriangle>	m_tetrahedraAroundTriangle;
	Attribute_T<HexahedraAroundVertex>		m_hexahedraAroundVertex;
	Attribute_T<HexahedraAroundEdge>		m_hexahedraAroundEdge;
	Attribute_T<HexahedraAroundQuad>		m_hexahedraAroundQuad;
//	Attribute_T<core::topology::Topology::Edge> edges_;

	// compatibility
public:
	virtual bool hasPos() const override;
	virtual SReal getPX(int) const override;
	virtual SReal getPY(int) const override;
	virtual SReal getPZ(int) const override;

	virtual const SeqEdges&getEdges() override;
	virtual const SeqTriangles&getTriangles() override;
	virtual const SeqQuads&getQuads() override;
	virtual const SeqTetrahedra& getTetrahedra() override;
	virtual const SeqHexahedra& getHexahedra() override;

	virtual const EdgesAroundVertex&getEdgesAroundVertex(PointID i) override;
	virtual const EdgesInTriangle&getEdgesInTriangle(TriangleID i) override;
	virtual const EdgesInQuad&getEdgesInQuad(QuadID i) override;
	virtual const TrianglesAroundVertex&getTrianglesAroundVertex(PointID i) override;
	virtual const TrianglesAroundEdge&getTrianglesAroundEdge(EdgeID i) override;
	virtual const QuadsAroundVertex&getQuadsAroundVertex(PointID i) override;
	virtual const QuadsAroundEdge&getQuadsAroundEdge(EdgeID i) override;
	virtual const VerticesAroundVertex getVerticesAroundVertex(PointID i) override;
	virtual const sofa::helper::vector<index_type> getElementAroundElement(index_type elem) override;
	virtual const sofa::helper::vector<index_type> getElementAroundElements(sofa::helper::vector<index_type> elems) override;
	virtual int getEdgeIndex(PointID v1, PointID v2) override;
	virtual int getTriangleIndex(PointID v1, PointID v2, PointID v3) override;
	virtual int getQuadIndex(PointID v1, PointID v2, PointID v3, PointID v4) override;
	virtual int getVertexIndexInTriangle(const Triangle& t, PointID vertexIndex) const override;
	virtual int getEdgeIndexInTriangle(const EdgesInTriangle& t, EdgeID edgeIndex) const override;
	virtual int getVertexIndexInQuad(const Quad& t, PointID vertexIndex) const override;
	virtual int getEdgeIndexInQuad(const EdgesInQuad& t, EdgeID edgeIndex) const override;
	virtual void clear() override;
	virtual void addPoint(SReal px, SReal py, SReal pz) override;
	virtual void addEdge(int a, int b) override;
	virtual void addTriangle(int a, int b, int c) override;
	virtual void addQuad(int a, int b, int c, int d) override;
	virtual bool checkConnexity() override;
	virtual unsigned int getNumberOfConnectedComponent() override;
	virtual const sofa::helper::vector<index_type> getConnectedElement(index_type elem) override;
	virtual void reOrientateTriangle(TriangleID id) override;
};


} // namespace topology

} // namespace core

} // namespace sofa

#endif // SOFA_CORE_TOPOLOGY_MAPTOPOLOGY_H
