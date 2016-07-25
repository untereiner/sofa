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
	return true;
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
	return d_edge.getValue();
}

const sofa::core::topology::BaseMeshTopology::SeqTriangles&SurfaceTopologyContainer::getTriangles()
{
	return d_triangle.getValue();
}

const sofa::core::topology::BaseMeshTopology::SeqQuads&SurfaceTopologyContainer::getQuads()
{
	return d_quad.getValue();
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

const sofa::core::topology::BaseMeshTopology::TrianglesAroundVertex&SurfaceTopologyContainer::getTrianglesAroundVertex(core::topology::Topology::PointID i)
{
}

const sofa::core::topology::BaseMeshTopology::TrianglesAroundEdge&SurfaceTopologyContainer::getTrianglesAroundEdge(core::topology::Topology::EdgeID i)
{
}

const sofa::core::topology::BaseMeshTopology::QuadsAroundVertex&SurfaceTopologyContainer::getQuadsAroundVertex(core::topology::Topology::PointID i)
{
}

const sofa::core::topology::BaseMeshTopology::QuadsAroundEdge&SurfaceTopologyContainer::getQuadsAroundEdge(core::topology::Topology::EdgeID i)
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

