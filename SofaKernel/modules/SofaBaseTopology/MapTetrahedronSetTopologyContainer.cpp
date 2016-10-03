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

#include <SofaBaseTopology/MapTetrahedronSetTopologyContainer.h>
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace topology
{

SOFA_DECL_CLASS(MapTetrahedronSetTopologyContainer)
int MapTetrahedronSetTopologyContainerClass = core::RegisterObject("Tetrahedron set topology container backward compatibility")
		.add< MapTetrahedronSetTopologyContainer >()
		;

MapTetrahedronSetTopologyContainer::MapTetrahedronSetTopologyContainer():
	Inherit1(),
	map_(nullptr)
{}

MapTetrahedronSetTopologyContainer::~MapTetrahedronSetTopologyContainer()
{

}

int MapTetrahedronSetTopologyContainer::getNbPoints() const
{
	return map_->nb_cells<Vertex::ORBIT>();
}

void MapTetrahedronSetTopologyContainer::init()
{
	this->getContext()->get(map_);
	if (!map_)
		return;

	map_->add_attribute(edges_around_vertex_, "edges_around_vertex");
	map_->add_attribute(edges_in_triangle_, "edges_in_triangle");
	map_->add_attribute(triangles_around_vertex_, "triangles_around_vertex");
	map_->add_attribute(triangles_around_edge_, "triangles_around_edge");
	map_->add_attribute(oriented_triangles_around_vertex_, "oriented_triangles_around_vertex");
	map_->add_attribute(tetrahedra_around_vertex_, "tetrahedra_around_vertex");
	map_->add_attribute(tetrahedra_around_edge_, "tetrahedra_around_edge");
	map_->add_attribute(tetrahedra_around_triangle_, "tetrahedra_around_triangle");
	map_->add_attribute(edges_in_tetrahedron_, "edges_in_tetrahedron");
	map_->add_attribute(triangles_in_tetrahedron_, "triangles_in_tetrahedron");
	Inherit1::init();
}

void MapTetrahedronSetTopologyContainer::bwdInit()
{
	Inherit1::bwdInit();
}

void MapTetrahedronSetTopologyContainer::reinit()
{
	Inherit1::reinit();
}

void MapTetrahedronSetTopologyContainer::reset()
{
	Inherit1::reset();
}

void MapTetrahedronSetTopologyContainer::cleanup()
{
	if (edges_in_triangle_.is_valid())
		map_->remove_attribute(edges_in_triangle_);
	if (triangles_around_vertex_.is_valid())
		map_->remove_attribute(triangles_around_vertex_);
	if (triangles_around_edge_.is_valid())
		map_->remove_attribute(triangles_around_edge_);
	if (oriented_triangles_around_vertex_.is_valid())
		map_->remove_attribute(oriented_triangles_around_vertex_);
	if (tetrahedra_around_vertex_.is_valid())
		map_->remove_attribute(tetrahedra_around_vertex_);
	if (tetrahedra_around_edge_.is_valid())
		map_->remove_attribute(tetrahedra_around_edge_);
	if (tetrahedra_around_triangle_.is_valid())
		map_->remove_attribute(tetrahedra_around_triangle_);
	if (edges_in_tetrahedron_.is_valid())
		map_->remove_attribute(edges_in_tetrahedron_);
	if (triangles_in_tetrahedron_.is_valid())
		map_->remove_attribute(triangles_in_tetrahedron_);

	Inherit1::cleanup();
}

void MapTetrahedronSetTopologyContainer::draw(const sofa::core::visual::VisualParams*)
{
}

bool MapTetrahedronSetTopologyContainer::load(const char* filename)
{
}

const EdgeSetTopologyContainer::SeqEdges&MapTetrahedronSetTopologyContainer::getEdges()
{
	return map_->getEdges();
}

const TriangleSetTopologyContainer::SeqTriangles&MapTetrahedronSetTopologyContainer::getTriangles()
{
	return map_->getTriangles();
}

const TetrahedronSetTopologyContainer::SeqTetrahedra&MapTetrahedronSetTopologyContainer::getTetrahedra()
{
	return map_->getTetrahedra();
}

int MapTetrahedronSetTopologyContainer::getNbEdges()
{
	return map_->nb_cells<Edge::ORBIT>();
}

int MapTetrahedronSetTopologyContainer::getNbTriangles()
{
	return map_->nb_cells<Face::ORBIT>();
}

int MapTetrahedronSetTopologyContainer::getNbTetrahedra()
{
	return map_->template nb_cells<Volume::ORBIT>();
}

int MapTetrahedronSetTopologyContainer::getNbTetras()
{
	return getNbTetrahedra();
}

const EdgeSetTopologyContainer::EdgesAroundVertex&MapTetrahedronSetTopologyContainer::getEdgesAroundVertex(TetrahedronSetTopologyContainer::PointID i)
{
	if(!hasEdgesAroundVertex())
		createEdgesAroundVertexArray();

	return edges_around_vertex_[i];
}

const TriangleSetTopologyContainer::EdgesInTriangle&MapTetrahedronSetTopologyContainer::getEdgesInTriangle(TetrahedronSetTopologyContainer::TriangleID i)
{
}

const TetrahedronSetTopologyContainer::EdgesInTetrahedron&MapTetrahedronSetTopologyContainer::getEdgesInTetrahedron(TetrahedronSetTopologyContainer::TetraID i)
{
}

const TriangleSetTopologyContainer::TrianglesAroundVertex&MapTetrahedronSetTopologyContainer::getTrianglesAroundVertex(TetrahedronSetTopologyContainer::PointID i)
{
}

const TriangleSetTopologyContainer::TrianglesAroundEdge&MapTetrahedronSetTopologyContainer::getTrianglesAroundEdge(TetrahedronSetTopologyContainer::EdgeID i)
{
}

const TetrahedronSetTopologyContainer::TrianglesInTetrahedron&MapTetrahedronSetTopologyContainer::getTrianglesInTetrahedron(TetrahedronSetTopologyContainer::TetraID i)
{
}

const TetrahedronSetTopologyContainer::TetrahedraAroundVertex&MapTetrahedronSetTopologyContainer::getTetrahedraAroundVertex(TetrahedronSetTopologyContainer::PointID i)
{
}

const TetrahedronSetTopologyContainer::TetrahedraAroundEdge&MapTetrahedronSetTopologyContainer::getTetrahedraAroundEdge(TetrahedronSetTopologyContainer::EdgeID i)
{
}

const TetrahedronSetTopologyContainer::TetrahedraAroundTriangle&MapTetrahedronSetTopologyContainer::getTetrahedraAroundTriangle(TetrahedronSetTopologyContainer::TriangleID i)
{
}

const sofa::core::topology::BaseMeshTopology::VerticesAroundVertex MapTetrahedronSetTopologyContainer::getVerticesAroundVertex(TetrahedronSetTopologyContainer::PointID i)
{
}

const sofa::helper::vector<sofa::core::topology::Topology::index_type> MapTetrahedronSetTopologyContainer::getElementAroundElement(sofa::core::topology::Topology::index_type elem)
{
}

const sofa::helper::vector<sofa::core::topology::Topology::index_type> MapTetrahedronSetTopologyContainer::getElementAroundElements(sofa::helper::vector<sofa::core::topology::Topology::index_type> elems)
{
}

int MapTetrahedronSetTopologyContainer::getEdgeIndex(TetrahedronSetTopologyContainer::PointID v1, TetrahedronSetTopologyContainer::PointID v2)
{
}

int MapTetrahedronSetTopologyContainer::getTriangleIndex(TetrahedronSetTopologyContainer::PointID v1, TetrahedronSetTopologyContainer::PointID v2, TetrahedronSetTopologyContainer::PointID v3)
{
}

int MapTetrahedronSetTopologyContainer::getTetrahedronIndex(TetrahedronSetTopologyContainer::PointID v1, TetrahedronSetTopologyContainer::PointID v2, TetrahedronSetTopologyContainer::PointID v3, TetrahedronSetTopologyContainer::PointID v4)
{
}

int MapTetrahedronSetTopologyContainer::getVertexIndexInTriangle(const TetrahedronSetTopologyContainer::Triangle& t, TetrahedronSetTopologyContainer::PointID vertexIndex) const
{
}

int MapTetrahedronSetTopologyContainer::getEdgeIndexInTriangle(const TriangleSetTopologyContainer::EdgesInTriangle& t, TetrahedronSetTopologyContainer::EdgeID edgeIndex) const
{
}

int MapTetrahedronSetTopologyContainer::getVertexIndexInTetrahedron(const TetrahedronSetTopologyContainer::Tetra& t, TetrahedronSetTopologyContainer::PointID vertexIndex) const
{
}

int MapTetrahedronSetTopologyContainer::getEdgeIndexInTetrahedron(const TetrahedronSetTopologyContainer::EdgesInTetrahedron& t, TetrahedronSetTopologyContainer::EdgeID edgeIndex) const
{
}

int MapTetrahedronSetTopologyContainer::getTriangleIndexInTetrahedron(const TetrahedronSetTopologyContainer::TrianglesInTetrahedron& t, TetrahedronSetTopologyContainer::TriangleID triangleIndex) const
{
}

 core::topology::Topology::Edge MapTetrahedronSetTopologyContainer::getLocalEdgesInTetrahedron(const TetrahedronSetTopologyContainer::PointID i) const
{
}

TetrahedronSetTopologyContainer::Triangle MapTetrahedronSetTopologyContainer::getLocalTrianglesInTetrahedron(const TetrahedronSetTopologyContainer::PointID i) const
{
}

void MapTetrahedronSetTopologyContainer::clear()
{
}

void MapTetrahedronSetTopologyContainer::addPoint(SReal px, SReal py, SReal pz)
{
}

void MapTetrahedronSetTopologyContainer::addEdge(int a, int b)
{
}

void MapTetrahedronSetTopologyContainer::addTriangle(int a, int b, int c)
{
}

void MapTetrahedronSetTopologyContainer::addQuad(int a, int b, int c, int d)
{
}

void MapTetrahedronSetTopologyContainer::addTetra(int a, int b, int c, int d)
{
}

void MapTetrahedronSetTopologyContainer::addHexa(int a, int b, int c, int d, int e, int f, int g, int h)
{
}

bool MapTetrahedronSetTopologyContainer::checkConnexity()
{
	return true;
}

unsigned int MapTetrahedronSetTopologyContainer::getNumberOfConnectedComponent()
{
	return 0u;
}

const sofa::helper::vector<sofa::core::topology::Topology::index_type> MapTetrahedronSetTopologyContainer::getConnectedElement(sofa::core::topology::Topology::index_type elem)
{
	return sofa::helper::vector<sofa::core::topology::Topology::index_type>();
}

void MapTetrahedronSetTopologyContainer::reOrientateTriangle(TetrahedronSetTopologyContainer::TriangleID id)
{
}

const sofa::helper::vector<TetrahedronSetTopologyContainer::TriangleID>&MapTetrahedronSetTopologyContainer::getTrianglesOnBorder()
{
	static const sofa::helper::vector<TetrahedronSetTopologyContainer::TriangleID> empty;
	return empty;
}

const sofa::helper::vector<TetrahedronSetTopologyContainer::EdgeID>&MapTetrahedronSetTopologyContainer::getEdgesOnBorder()
{
	static const sofa::helper::vector<TetrahedronSetTopologyContainer::EdgeID> empty;
	return empty;
}

const sofa::helper::vector<TetrahedronSetTopologyContainer::PointID>&MapTetrahedronSetTopologyContainer::getPointsOnBorder()
{
	static const sofa::helper::vector<TetrahedronSetTopologyContainer::PointID> empty;
	return empty;
}

unsigned int MapTetrahedronSetTopologyContainer::getNumberOfElements() const
{
	return map_->template nb_cells<Volume::ORBIT>();
}

bool MapTetrahedronSetTopologyContainer::checkTopology() const
{
	return true;
}

bool MapTetrahedronSetTopologyContainer::hasEdgesAroundVertex() const
{
	return edges_around_vertex_.is_valid();
}

void MapTetrahedronSetTopologyContainer::createEdgesAroundVertexArray()
{
	if (!edges_around_vertex_.is_valid())
		map_->add_attribute(edges_around_vertex_, "edges_around_vertex_");

	map_->parallel_foreach_cell([&](Vertex v, uint32_t)
	{
		auto& edges = edges_around_vertex_[v];
		edges.clear();
		map_->foreach_incident_edge(v, [&](Edge e)
		{
			edges.push_back(map_->embedding(e));
		});
	});
}

void MapTetrahedronSetTopologyContainer::createEdgesInTetrahedronArray()
{
	if (!edges_in_tetrahedron_.is_valid())
		map_->add_attribute(edges_in_tetrahedron_,"edges_in_tetrahedron_");
}

void MapTetrahedronSetTopologyContainer::createTrianglesInTetrahedronArray()
{
	if (!triangles_in_tetrahedron_.is_valid())
		map_->add_attribute(triangles_in_tetrahedron_,"triangles_in_tetrahedron_");
}

void MapTetrahedronSetTopologyContainer::createTetrahedraAroundVertexArray()
{
	if (!tetrahedra_around_vertex_.is_valid())
		map_->add_attribute(tetrahedra_around_vertex_,"tetrahedra_around_vertex_");
}

void MapTetrahedronSetTopologyContainer::createTetrahedraAroundEdgeArray()
{
	if (!tetrahedra_around_edge_.is_valid())
		map_->add_attribute(tetrahedra_around_edge_,"tetrahedra_around_edge_");
}

void MapTetrahedronSetTopologyContainer::createTetrahedraAroundTriangleArray()
{
	if (!tetrahedra_around_triangle_.is_valid())
		map_->add_attribute(tetrahedra_around_triangle_,"tetrahedra_around_triangle_");
}

} // namespace topology
} // namespace component
} // namespace sofa

