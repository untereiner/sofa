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

#include <cgogn/io/surface_import.h>


namespace sofa
{

namespace component
{

namespace topology
{

void SurfaceTopologyContainer::initFromMeshLoader()
{
	helper::ReadAccessor< Data< VecCoord > > m_position = d_initPoints;
	helper::ReadAccessor< Data< helper::vector< Triangle > > > m_tri = d_triangle;
	helper::ReadAccessor< Data< helper::vector< Quad> > > m_quad = d_quad;

	cgogn::io::SurfaceImport<Topo_Traits::MapTraits> surface_import;
	surface_import.set_nb_vertices(m_position.size());

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

	for(const Triangle& t : m_tri.ref())
		surface_import.add_triangle(t[0], t[1], t[2]);
	for(const Quad& q : m_quad.ref())
		surface_import.add_quad(q[0], q[1], q[2], q[3]);

	surface_import.create_map(topology_);
}

void SurfaceTopologyContainer::init()
{
	topology_.clear_and_remove_attributes();
	Inherit1::init();
	initFromMeshLoader();
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
	topology_.clear_and_remove_attributes();
	Inherit1::reset();
}

void SurfaceTopologyContainer::cleanup()
{
	Inherit1::cleanup();
}

bool SurfaceTopologyContainer::hasPos() const
{
}

SReal SurfaceTopologyContainer::getPX(int) const
{
}

SReal SurfaceTopologyContainer::getPY(int) const
{
}

SReal SurfaceTopologyContainer::getPZ(int) const
{
}

const sofa::core::topology::BaseMeshTopology::SeqEdges&SurfaceTopologyContainer::getEdges()
{
}

const sofa::core::topology::BaseMeshTopology::SeqTriangles&SurfaceTopologyContainer::getTriangles()
{
}

const sofa::core::topology::BaseMeshTopology::SeqQuads&SurfaceTopologyContainer::getQuads()
{
}

const sofa::core::topology::BaseMeshTopology::SeqTetrahedra&SurfaceTopologyContainer::getTetrahedra()
{
}

const sofa::core::topology::BaseMeshTopology::SeqHexahedra&SurfaceTopologyContainer::getHexahedra()
{
}

int SurfaceTopologyContainer::getNbTriangles()
{
}

int SurfaceTopologyContainer::getNbQuads()
{
}

int SurfaceTopologyContainer::getNbTetrahedra()
{
}

int SurfaceTopologyContainer::getNbHexahedra()
{
}

const SurfaceTopologyContainer::Edge SurfaceTopologyContainer::getEdge(core::topology::Topology::EdgeID i)
{
}

const core::topology::Topology::Triangle SurfaceTopologyContainer::getTriangle(core::topology::Topology::TriangleID i)
{
}

const core::topology::Topology::Quad SurfaceTopologyContainer::getQuad(core::topology::Topology::QuadID i)
{
}

const core::topology::Topology::Tetra SurfaceTopologyContainer::getTetrahedron(core::topology::Topology::TetraID i)
{
}

const core::topology::Topology::Hexa SurfaceTopologyContainer::getHexahedron(core::topology::Topology::HexaID i)
{
}

int SurfaceTopologyContainer::getNbTetras()
{
}

int SurfaceTopologyContainer::getNbHexas()
{
}

core::topology::Topology::Tetra SurfaceTopologyContainer::getTetra(core::topology::Topology::TetraID i)
{
}

core::topology::Topology::Hexa SurfaceTopologyContainer::getHexa(core::topology::Topology::HexaID i)
{
}

const sofa::core::topology::BaseMeshTopology::SeqTetrahedra&SurfaceTopologyContainer::getTetras()
{
}

const sofa::core::topology::BaseMeshTopology::SeqHexahedra&SurfaceTopologyContainer::getHexas()
{
}

const sofa::core::topology::BaseMeshTopology::EdgesAroundVertex&SurfaceTopologyContainer::getEdgesAroundVertex(core::topology::Topology::PointID i)
{
}

const sofa::core::topology::BaseMeshTopology::EdgesInTriangle&SurfaceTopologyContainer::getEdgesInTriangle(core::topology::Topology::TriangleID i)
{
}

const sofa::core::topology::BaseMeshTopology::EdgesInQuad&SurfaceTopologyContainer::getEdgesInQuad(core::topology::Topology::QuadID i)
{
}

const sofa::core::topology::BaseMeshTopology::EdgesInTetrahedron&SurfaceTopologyContainer::getEdgesInTetrahedron(core::topology::Topology::TetraID i)
{
}

const sofa::core::topology::BaseMeshTopology::EdgesInHexahedron&SurfaceTopologyContainer::getEdgesInHexahedron(core::topology::Topology::HexaID i)
{
}

const sofa::core::topology::BaseMeshTopology::TrianglesAroundVertex&SurfaceTopologyContainer::getTrianglesAroundVertex(core::topology::Topology::PointID i)
{
}

const sofa::core::topology::BaseMeshTopology::TrianglesAroundEdge&SurfaceTopologyContainer::getTrianglesAroundEdge(core::topology::Topology::EdgeID i)
{
}

const sofa::core::topology::BaseMeshTopology::TrianglesInTetrahedron&SurfaceTopologyContainer::getTrianglesInTetrahedron(core::topology::Topology::TetraID i)
{
}

const sofa::core::topology::BaseMeshTopology::QuadsAroundVertex&SurfaceTopologyContainer::getQuadsAroundVertex(core::topology::Topology::PointID i)
{
}

const sofa::core::topology::BaseMeshTopology::QuadsAroundEdge&SurfaceTopologyContainer::getQuadsAroundEdge(core::topology::Topology::EdgeID i)
{
}

const sofa::core::topology::BaseMeshTopology::QuadsInHexahedron&SurfaceTopologyContainer::getQuadsInHexahedron(core::topology::Topology::HexaID i)
{
}

const sofa::core::topology::BaseMeshTopology::TetrahedraAroundVertex&SurfaceTopologyContainer::getTetrahedraAroundVertex(core::topology::Topology::PointID i)
{
}

const sofa::core::topology::BaseMeshTopology::TetrahedraAroundEdge&SurfaceTopologyContainer::getTetrahedraAroundEdge(core::topology::Topology::EdgeID i)
{
}

const sofa::core::topology::BaseMeshTopology::TetrahedraAroundTriangle&SurfaceTopologyContainer::getTetrahedraAroundTriangle(core::topology::Topology::TriangleID i)
{
}

const sofa::core::topology::BaseMeshTopology::HexahedraAroundVertex&SurfaceTopologyContainer::getHexahedraAroundVertex(core::topology::Topology::PointID i)
{
}

const sofa::core::topology::BaseMeshTopology::HexahedraAroundEdge&SurfaceTopologyContainer::getHexahedraAroundEdge(core::topology::Topology::EdgeID i)
{
}

const sofa::core::topology::BaseMeshTopology::HexahedraAroundQuad&SurfaceTopologyContainer::getHexahedraAroundQuad(core::topology::Topology::QuadID i)
{
}

const sofa::core::topology::BaseMeshTopology::VerticesAroundVertex SurfaceTopologyContainer::getVerticesAroundVertex(core::topology::Topology::PointID i)
{
}

const sofa::helper::vector<core::topology::Topology::index_type> SurfaceTopologyContainer::getElementAroundElement(core::topology::Topology::index_type elem)
{
}

const sofa::helper::vector<core::topology::Topology::index_type> SurfaceTopologyContainer::getElementAroundElements(sofa::helper::vector<core::topology::Topology::index_type> elems)
{
}

int SurfaceTopologyContainer::getEdgeIndex(core::topology::Topology::PointID v1, core::topology::Topology::PointID v2)
{
}

int SurfaceTopologyContainer::getTriangleIndex(core::topology::Topology::PointID v1, core::topology::Topology::PointID v2, core::topology::Topology::PointID v3)
{
}

int SurfaceTopologyContainer::getQuadIndex(core::topology::Topology::PointID v1, core::topology::Topology::PointID v2, core::topology::Topology::PointID v3, core::topology::Topology::PointID v4)
{
}

int SurfaceTopologyContainer::getTetrahedronIndex(core::topology::Topology::PointID v1, core::topology::Topology::PointID v2, core::topology::Topology::PointID v3, core::topology::Topology::PointID v4)
{
}

int SurfaceTopologyContainer::getHexahedronIndex(core::topology::Topology::PointID v1, core::topology::Topology::PointID v2, core::topology::Topology::PointID v3, core::topology::Topology::PointID v4, core::topology::Topology::PointID v5, core::topology::Topology::PointID v6, core::topology::Topology::PointID v7, core::topology::Topology::PointID v8)
{
}

int SurfaceTopologyContainer::getVertexIndexInTriangle(const core::topology::Topology::Triangle& t, core::topology::Topology::PointID vertexIndex) const
{
}

int SurfaceTopologyContainer::getEdgeIndexInTriangle(const sofa::core::topology::BaseMeshTopology::EdgesInTriangle& t, core::topology::Topology::EdgeID edgeIndex) const
{
}

int SurfaceTopologyContainer::getVertexIndexInQuad(const core::topology::Topology::Quad& t, core::topology::Topology::PointID vertexIndex) const
{
}

int SurfaceTopologyContainer::getEdgeIndexInQuad(const sofa::core::topology::BaseMeshTopology::EdgesInQuad& t, core::topology::Topology::EdgeID edgeIndex) const
{
}

int SurfaceTopologyContainer::getVertexIndexInTetrahedron(const core::topology::Topology::Tetra& t, core::topology::Topology::PointID vertexIndex) const
{
}

int SurfaceTopologyContainer::getEdgeIndexInTetrahedron(const sofa::core::topology::BaseMeshTopology::EdgesInTetrahedron& t, core::topology::Topology::EdgeID edgeIndex) const
{
}

int SurfaceTopologyContainer::getTriangleIndexInTetrahedron(const sofa::core::topology::BaseMeshTopology::TrianglesInTetrahedron& t, core::topology::Topology::TriangleID triangleIndex) const
{
}

int SurfaceTopologyContainer::getVertexIndexInHexahedron(const core::topology::Topology::Hexa& t, core::topology::Topology::PointID vertexIndex) const
{
}

int SurfaceTopologyContainer::getEdgeIndexInHexahedron(const sofa::core::topology::BaseMeshTopology::EdgesInHexahedron& t, core::topology::Topology::EdgeID edgeIndex) const
{
}

int SurfaceTopologyContainer::getQuadIndexInHexahedron(const sofa::core::topology::BaseMeshTopology::QuadsInHexahedron& t, core::topology::Topology::QuadID quadIndex) const
{
}

SurfaceTopologyContainer::Edge SurfaceTopologyContainer::getLocalEdgesInTetrahedron(const core::topology::Topology::PointID i) const
{
}

core::topology::Topology::Triangle SurfaceTopologyContainer::getLocalTrianglesInTetrahedron(const core::topology::Topology::PointID i) const
{
}

SurfaceTopologyContainer::Edge SurfaceTopologyContainer::getLocalEdgesInHexahedron(const core::topology::Topology::PointID i) const
{
}

core::topology::Topology::Quad SurfaceTopologyContainer::getLocalQuadsInHexahedron(const core::topology::Topology::PointID i) const
{
}

void SurfaceTopologyContainer::clear()
{
}

void SurfaceTopologyContainer::addPoint(SReal px, SReal py, SReal pz)
{
}

void SurfaceTopologyContainer::addEdge(int a, int b)
{
}

void SurfaceTopologyContainer::addTriangle(int a, int b, int c)
{
}

void SurfaceTopologyContainer::addQuad(int a, int b, int c, int d)
{
}

void SurfaceTopologyContainer::addTetra(int a, int b, int c, int d)
{
}

void SurfaceTopologyContainer::addHexa(int a, int b, int c, int d, int e, int f, int g, int h)
{
}

bool SurfaceTopologyContainer::checkConnexity()
{
}

unsigned int SurfaceTopologyContainer::getNumberOfConnectedComponent()
{
}

const sofa::helper::vector<core::topology::Topology::index_type> SurfaceTopologyContainer::getConnectedElement(core::topology::Topology::index_type elem)
{
}

void SurfaceTopologyContainer::reOrientateTriangle(core::topology::Topology::TriangleID id)
{
}

} // namespace topology

} // namespace component

} // namespace sofa

