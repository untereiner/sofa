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

#include <SofaBaseTopology/MapTriangleSetTopologyContainer.h>
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace topology
{

SOFA_DECL_CLASS(MapTriangleSetTopologyContainer)
int MapTriangleSetTopologyContainerClass = core::RegisterObject("Triangle set topology container backward compatibility")
        .add< MapTriangleSetTopologyContainer >()
        ;


MapTriangleSetTopologyContainer::MapTriangleSetTopologyContainer() :
	Inherit1(),
	map_(nullptr)
{
}

MapTriangleSetTopologyContainer::~MapTriangleSetTopologyContainer()
{
}

void MapTriangleSetTopologyContainer::init()
{
	Inherit1::init();
	this->getContext()->get(map_);
	if (!map_)
		return;
	map_->add_attribute(m_edgesInTriangle, "edgesInTriangle");
	map_->add_attribute(m_trianglesAroundVertex, "trianglesAroundVertex");
	map_->add_attribute(m_trianglesAroundEdge, "trianglesAroundEdge");
	map_->add_attribute(m_orientedTrianglesAroundVertex, "orientedTrianglesAroundVertex");
	this->nbPoints.setValue(map_->nb_cells<Vertex::ORBIT>());
}

void MapTriangleSetTopologyContainer::bwdInit()
{
    Inherit1::bwdInit();
}

void MapTriangleSetTopologyContainer::reinit()
{
    Inherit1::reinit();
}

void MapTriangleSetTopologyContainer::reset()
{
    Inherit1::reset();
}

void MapTriangleSetTopologyContainer::cleanup()
{
    Inherit1::cleanup();
}

const EdgeSetTopologyContainer::SeqEdges&MapTriangleSetTopologyContainer::getEdges()
{
    return map_->getEdges();
}

const TriangleSetTopologyContainer::SeqTriangles&MapTriangleSetTopologyContainer::getTriangles()
{
        return map_->getTriangles();
}

const core::topology::BaseMeshTopology::SeqQuads&MapTriangleSetTopologyContainer::getQuads()
{
    return Inherit1::getQuads();
}

const core::topology::BaseMeshTopology::SeqTetrahedra&MapTriangleSetTopologyContainer::getTetrahedra()
{
    return Inherit1::getTetrahedra();
}

const core::topology::BaseMeshTopology::SeqHexahedra&MapTriangleSetTopologyContainer::getHexahedra()
{
    return Inherit1::getHexahedra();
}

int MapTriangleSetTopologyContainer::getNbEdges()
{
	return map_->nb_cells<Edge::ORBIT>();
}

int MapTriangleSetTopologyContainer::getNbTriangles()
{
	return map_->nb_cells<Face::ORBIT>();
}

int MapTriangleSetTopologyContainer::getNbQuads()
{
    return 0;
}

int MapTriangleSetTopologyContainer::getNbTetrahedra()
{
    return 0;
}

int MapTriangleSetTopologyContainer::getNbHexahedra()
{
    return 0;
}

const TriangleSetTopologyContainer::Edge MapTriangleSetTopologyContainer::getEdge(TriangleSetTopologyContainer::EdgeID i)
{
    return map_->getEdges()[i];
}

const TriangleSetTopologyContainer::Triangle MapTriangleSetTopologyContainer::getTriangle(TriangleSetTopologyContainer::TriangleID i)
{
    return map_->getTriangles()[i];
}

const core::topology::Topology::Quad MapTriangleSetTopologyContainer::getQuad(core::topology::Topology::QuadID i)
{
    return Inherit1::getQuad(i);
}

const core::topology::Topology::Tetra MapTriangleSetTopologyContainer::getTetrahedron(core::topology::Topology::TetraID i)
{
    return Inherit1::getTetrahedron(i);
}

const core::topology::Topology::Hexa MapTriangleSetTopologyContainer::getHexahedron(core::topology::Topology::HexaID i)
{
    return Inherit1::getHexahedron(i);
}

int MapTriangleSetTopologyContainer::getNbTetras()
{
    return 0;
}

int MapTriangleSetTopologyContainer::getNbHexas()
{
    return 0;
}

core::topology::Topology::Tetra MapTriangleSetTopologyContainer::getTetra(core::topology::Topology::TetraID i)
{
    return Inherit1::getTetra(i);
}

core::topology::Topology::Hexa MapTriangleSetTopologyContainer::getHexa(core::topology::Topology::HexaID i)
{
	return Inherit1::getHexa(i);
}

const core::topology::BaseMeshTopology::SeqTetrahedra&MapTriangleSetTopologyContainer::getTetras()
{
	return Inherit1::getTetras();
}

const core::topology::BaseMeshTopology::SeqHexahedra&MapTriangleSetTopologyContainer::getHexas()
{
	return Inherit1::getHexas();
}

const EdgeSetTopologyContainer::EdgesAroundVertex&MapTriangleSetTopologyContainer::getEdgesAroundVertex(TriangleSetTopologyContainer::PointID i)
{
	return m_edgesAroundVertex[i];
}

const TriangleSetTopologyContainer::EdgesInTriangle&MapTriangleSetTopologyContainer::getEdgesInTriangle(TriangleSetTopologyContainer::TriangleID i)
{
	return m_edgesInTriangle[i];
}

const core::topology::BaseMeshTopology::EdgesInQuad&MapTriangleSetTopologyContainer::getEdgesInQuad(core::topology::Topology::QuadID i)
{
	return Inherit1::getEdgesInQuad(i);
}

const core::topology::BaseMeshTopology::EdgesInTetrahedron&MapTriangleSetTopologyContainer::getEdgesInTetrahedron(core::topology::Topology::TetraID i)
{
	return Inherit1::getEdgesInTetrahedron(i);
}

const core::topology::BaseMeshTopology::EdgesInHexahedron&MapTriangleSetTopologyContainer::getEdgesInHexahedron(core::topology::Topology::HexaID i)
{
	return Inherit1::getEdgesInHexahedron(i);
}

const TriangleSetTopologyContainer::TrianglesAroundVertex&MapTriangleSetTopologyContainer::getTrianglesAroundVertex(TriangleSetTopologyContainer::PointID i)
{
	return m_trianglesAroundVertex[i];
}

const TriangleSetTopologyContainer::TrianglesAroundEdge&MapTriangleSetTopologyContainer::getTrianglesAroundEdge(TriangleSetTopologyContainer::EdgeID i)
{
	return m_trianglesAroundEdge[i];
}

const core::topology::BaseMeshTopology::TrianglesInTetrahedron&MapTriangleSetTopologyContainer::getTrianglesInTetrahedron(core::topology::Topology::TetraID i)
{
	return Inherit1::getTrianglesInTetrahedron(i);
}

const core::topology::BaseMeshTopology::QuadsAroundVertex&MapTriangleSetTopologyContainer::getQuadsAroundVertex(TriangleSetTopologyContainer::PointID i)
{
	return Inherit1::getQuadsAroundVertex(i);
}

const core::topology::BaseMeshTopology::QuadsAroundEdge&MapTriangleSetTopologyContainer::getQuadsAroundEdge(TriangleSetTopologyContainer::EdgeID i)
{
	return Inherit1::getQuadsAroundEdge(i);
}

const core::topology::BaseMeshTopology::QuadsInHexahedron&MapTriangleSetTopologyContainer::getQuadsInHexahedron(core::topology::Topology::HexaID i)
{
	return Inherit1::getQuadsInHexahedron(i);
}

const core::topology::BaseMeshTopology::TetrahedraAroundVertex&MapTriangleSetTopologyContainer::getTetrahedraAroundVertex(TriangleSetTopologyContainer::PointID i)
{
	return Inherit1::getTetrahedraAroundVertex(i);
}

const core::topology::BaseMeshTopology::TetrahedraAroundEdge&MapTriangleSetTopologyContainer::getTetrahedraAroundEdge(TriangleSetTopologyContainer::EdgeID i)
{
	return Inherit1::getTetrahedraAroundEdge(i);
}

const core::topology::BaseMeshTopology::TetrahedraAroundTriangle&MapTriangleSetTopologyContainer::getTetrahedraAroundTriangle(TriangleSetTopologyContainer::TriangleID i)
{
	return Inherit1::getTetrahedraAroundTriangle(i);
}

const core::topology::BaseMeshTopology::HexahedraAroundVertex&MapTriangleSetTopologyContainer::getHexahedraAroundVertex(TriangleSetTopologyContainer::PointID i)
{
	return Inherit1::getHexahedraAroundVertex(i);
}

const core::topology::BaseMeshTopology::HexahedraAroundEdge&MapTriangleSetTopologyContainer::getHexahedraAroundEdge(TriangleSetTopologyContainer::EdgeID i)
{
	return Inherit1::getHexahedraAroundEdge(i);
}

const core::topology::BaseMeshTopology::HexahedraAroundQuad&MapTriangleSetTopologyContainer::getHexahedraAroundQuad(core::topology::Topology::QuadID i)
{
	return Inherit1::getHexahedraAroundQuad(i);
}

const core::topology::BaseMeshTopology::VerticesAroundVertex MapTriangleSetTopologyContainer::getVerticesAroundVertex(TriangleSetTopologyContainer::PointID i)
{
	return Inherit1::getVerticesAroundVertex(i);
}

int MapTriangleSetTopologyContainer::getEdgeIndex(TriangleSetTopologyContainer::PointID v1, TriangleSetTopologyContainer::PointID v2)
{
	return map_->getEdgeIndex(v1, v2);
}

int MapTriangleSetTopologyContainer::getTriangleIndex(TriangleSetTopologyContainer::PointID v1, TriangleSetTopologyContainer::PointID v2, TriangleSetTopologyContainer::PointID v3)
{
	return map_->getTriangleIndex(v1, v2, v3);
}

int MapTriangleSetTopologyContainer::getQuadIndex(TriangleSetTopologyContainer::PointID, TriangleSetTopologyContainer::PointID, TriangleSetTopologyContainer::PointID,
												  TriangleSetTopologyContainer::PointID)
{
	return -1;
}

int MapTriangleSetTopologyContainer::getTetrahedronIndex(TriangleSetTopologyContainer::PointID, TriangleSetTopologyContainer::PointID,
														 TriangleSetTopologyContainer::PointID, TriangleSetTopologyContainer::PointID)
{
	return -1;
}

int MapTriangleSetTopologyContainer::getHexahedronIndex(TriangleSetTopologyContainer::PointID , TriangleSetTopologyContainer::PointID , TriangleSetTopologyContainer::PointID ,
														TriangleSetTopologyContainer::PointID , TriangleSetTopologyContainer::PointID , TriangleSetTopologyContainer::PointID ,
														TriangleSetTopologyContainer::PointID , TriangleSetTopologyContainer::PointID)
{
	return -1;
}

int MapTriangleSetTopologyContainer::getVertexIndexInTriangle(const TriangleSetTopologyContainer::Triangle& t, TriangleSetTopologyContainer::PointID vertexIndex) const
{
	return map_->getVertexIndexInTriangle(t, vertexIndex);
}

int MapTriangleSetTopologyContainer::getEdgeIndexInTriangle(const TriangleSetTopologyContainer::EdgesInTriangle& t, TriangleSetTopologyContainer::EdgeID edgeIndex) const
{
	return map_->getEdgeIndexInTriangle(t,edgeIndex);
}

int MapTriangleSetTopologyContainer::getVertexIndexInQuad(const core::topology::Topology::Quad& , TriangleSetTopologyContainer::PointID) const
{
	return -1;
}

int MapTriangleSetTopologyContainer::getEdgeIndexInQuad(const core::topology::BaseMeshTopology::EdgesInQuad& , TriangleSetTopologyContainer::EdgeID ) const
{
	return -1;
}

int MapTriangleSetTopologyContainer::getVertexIndexInTetrahedron(const core::topology::Topology::Tetra& , TriangleSetTopologyContainer::PointID ) const
{
	return -1;
}

int MapTriangleSetTopologyContainer::getEdgeIndexInTetrahedron(const core::topology::BaseMeshTopology::EdgesInTetrahedron& , TriangleSetTopologyContainer::EdgeID ) const
{
	return -1;
}

int MapTriangleSetTopologyContainer::getTriangleIndexInTetrahedron(const core::topology::BaseMeshTopology::TrianglesInTetrahedron& , TriangleSetTopologyContainer::TriangleID ) const
{
	return -1;
}

int MapTriangleSetTopologyContainer::getVertexIndexInHexahedron(const core::topology::Topology::Hexa& , TriangleSetTopologyContainer::PointID ) const
{
	return -1;
}

int MapTriangleSetTopologyContainer::getEdgeIndexInHexahedron(const core::topology::BaseMeshTopology::EdgesInHexahedron& , TriangleSetTopologyContainer::EdgeID ) const
{
	return -1;
}

int MapTriangleSetTopologyContainer::getQuadIndexInHexahedron(const core::topology::BaseMeshTopology::QuadsInHexahedron& , core::topology::Topology::QuadID ) const
{
	return -1;
}

TriangleSetTopologyContainer::Edge MapTriangleSetTopologyContainer::getLocalEdgesInTetrahedron(const TriangleSetTopologyContainer::PointID i) const
{
	return Inherit1::getLocalEdgesInTetrahedron(i);
}

TriangleSetTopologyContainer::Triangle MapTriangleSetTopologyContainer::getLocalTrianglesInTetrahedron(const TriangleSetTopologyContainer::PointID i) const
{
	return Inherit1::getLocalTrianglesInTetrahedron(i);
}

TriangleSetTopologyContainer::Edge MapTriangleSetTopologyContainer::getLocalEdgesInHexahedron(const TriangleSetTopologyContainer::PointID i) const
{
	return Inherit1::getLocalEdgesInHexahedron(i);
}

core::topology::Topology::Quad MapTriangleSetTopologyContainer::getLocalQuadsInHexahedron(const TriangleSetTopologyContainer::PointID i) const
{
	return Inherit1::getLocalQuadsInHexahedron(i);
}

void MapTriangleSetTopologyContainer::clear()
{
	Inherit1::clear();
}

void MapTriangleSetTopologyContainer::addEdge(int, int)
{
//	map_->addEdge(a,b);
}

void MapTriangleSetTopologyContainer::addTriangle(int, int, int)
{
//	map_->addTriangle(a,b,c);
}

void MapTriangleSetTopologyContainer::addQuad(int, int, int, int)
{
//	map_->addQuad(a,b,c,d);
}

void MapTriangleSetTopologyContainer::addTetra(int, int, int, int)
{
//	map_->addTetra(a,b,c,d);
}

void MapTriangleSetTopologyContainer::addHexa(int, int, int, int, int, int, int, int)
{
//	map_->addHexa(a,b,c,d,e,f,g,h);
}

bool MapTriangleSetTopologyContainer::checkConnexity()
{
//	return map_->checkConnexity();
}

unsigned int MapTriangleSetTopologyContainer::getNumberOfConnectedComponent()
{
	return map_->getNumberOfConnectedComponent();
}

int MapTriangleSetTopologyContainer::getRevision() const
{
//	return map_->getRevision();
	return -1;
}

void MapTriangleSetTopologyContainer::reOrientateTriangle(TriangleSetTopologyContainer::TriangleID id)
{
	return map_->reOrientateTriangle(id);
}

void MapTriangleSetTopologyContainer::updateTopologyEngineGraph()
{
	Inherit1::updateTopologyEngineGraph();
//	map_->updateTopologyEngineGraph();
}

unsigned int MapTriangleSetTopologyContainer::getNumberOfElements() const
{
	return map_->getNumberOfConnectedComponent();
}

bool MapTriangleSetTopologyContainer::checkTopology() const
{
	// TODO : use map_->checkTopology();
	return Inherit1::checkTopology();
//	return map_->checkTopology();
}

void MapTriangleSetTopologyContainer::createEdgeSetArray()
{
//	map_->
}

const TriangleSetTopologyContainer::VecTriangleID MapTriangleSetTopologyContainer::getConnectedElement(TriangleSetTopologyContainer::TriangleID elem)
{
	return map_->getConnectedElement(elem);
}

const TriangleSetTopologyContainer::VecTriangleID MapTriangleSetTopologyContainer::getElementAroundElement(TriangleSetTopologyContainer::TriangleID elem)
{
	return map_->getElementAroundElement(elem);
}

const TriangleSetTopologyContainer::VecTriangleID MapTriangleSetTopologyContainer::getElementAroundElements(TriangleSetTopologyContainer::VecTriangleID elems)
{
	return map_->getElementAroundElements(elems);
}

void MapTriangleSetTopologyContainer::createTriangleSetArray()
{
	 // TODO
}

void MapTriangleSetTopologyContainer::createEdgesInTriangleArray()
{
	// TODO
}

void MapTriangleSetTopologyContainer::createTrianglesAroundVertexArray()
{
		// NOTHING TODO
}

void MapTriangleSetTopologyContainer::createTrianglesAroundEdgeArray()
{
		// NOTHING TODO
}

TriangleSetTopologyContainer::TrianglesAroundVertex&MapTriangleSetTopologyContainer::getTrianglesAroundVertexForModification(const TriangleSetTopologyContainer::PointID vertexIndex)
{
	return m_trianglesAroundVertex[vertexIndex];
}

TriangleSetTopologyContainer::TrianglesAroundEdge&MapTriangleSetTopologyContainer::getTrianglesAroundEdgeForModification(const TriangleSetTopologyContainer::EdgeID edgeIndex)
{
	return m_trianglesAroundEdge[edgeIndex];
}


} // namespace topology
} // namespace component
} // namespace sofa


