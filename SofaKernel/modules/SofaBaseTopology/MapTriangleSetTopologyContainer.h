#ifndef MAPTRIANGLESETTOPOLOGYCONTAINER_H
#define MAPTRIANGLESETTOPOLOGYCONTAINER_H

#include <SofaBaseTopology/TriangleSetTopologyContainer.h>
#include <SofaBaseTopology/SurfaceTopologyContainer.h>

namespace sofa
{

namespace component
{

namespace topology
{

class SOFA_BASE_TOPOLOGY_API MapTriangleSetTopologyContainer : public TriangleSetTopologyContainer
{
    friend class TriangleSetTopologyModifier;
public:
	SOFA_CLASS(MapTriangleSetTopologyContainer,TriangleSetTopologyContainer);
	template<typename T>
	using Attribute_T = core::topology::MapTopology::Attribute_T<T>;
	using Orbit = SurfaceTopologyContainer::Orbit;
	template <typename T, Orbit ORBIT>
	using Attribute = SurfaceTopologyContainer::Attribute<T,ORBIT>;
	using Vertex = SurfaceTopologyContainer::Vertex;
	using Edge = SurfaceTopologyContainer::Edge;
	using Face = SurfaceTopologyContainer::Face;

	MapTriangleSetTopologyContainer();
	virtual ~MapTriangleSetTopologyContainer() override;

    // BaseObject interface
public:
    virtual void init() override;
    virtual void bwdInit() override;
    virtual void reinit() override;
    virtual void reset() override;
    virtual void cleanup() override;

    // BaseMeshTopology interface
public:
    virtual const SeqEdges&getEdges() override;
    virtual const SeqTriangles&getTriangles() override;
    virtual const SeqQuads&getQuads() override;
    virtual const SeqTetrahedra&getTetrahedra() override;
    virtual const SeqHexahedra&getHexahedra() override;
    virtual int getNbEdges() override;
    virtual int getNbTriangles() override;
    virtual int getNbQuads() override;
    virtual int getNbTetrahedra() override;
    virtual int getNbHexahedra() override;
    virtual const core::topology::Topology::Edge getEdge(EdgeID i) override;
    virtual const Triangle getTriangle(TriangleID i) override;
    virtual const Quad getQuad(QuadID i) override;
    virtual const Tetra getTetrahedron(TetraID i) override;
    virtual const Hexa getHexahedron(HexaID i) override;
    virtual int getNbTetras() override;
    virtual int getNbHexas() override;
    virtual Tetra getTetra(TetraID i) override;
    virtual Hexa getHexa(HexaID i) override;
    virtual const SeqTetrahedra&getTetras() override;
    virtual const SeqHexahedra&getHexas() override;
    virtual const EdgesAroundVertex&getEdgesAroundVertex(PointID i) override;
    virtual const EdgesInTriangle&getEdgesInTriangle(TriangleID i) override;
    virtual const EdgesInQuad&getEdgesInQuad(QuadID i) override;
    virtual const EdgesInTetrahedron&getEdgesInTetrahedron(TetraID i) override;
    virtual const EdgesInHexahedron&getEdgesInHexahedron(HexaID i) override;
    virtual const TrianglesAroundVertex&getTrianglesAroundVertex(PointID i) override;
    virtual const TrianglesAroundEdge&getTrianglesAroundEdge(EdgeID i) override;
    virtual const TrianglesInTetrahedron&getTrianglesInTetrahedron(TetraID i) override;
    virtual const QuadsAroundVertex&getQuadsAroundVertex(PointID i) override;
    virtual const QuadsAroundEdge&getQuadsAroundEdge(EdgeID i) override;
    virtual const QuadsInHexahedron&getQuadsInHexahedron(HexaID i) override;
    virtual const TetrahedraAroundVertex&getTetrahedraAroundVertex(PointID i) override;
    virtual const TetrahedraAroundEdge&getTetrahedraAroundEdge(EdgeID i) override;
    virtual const TetrahedraAroundTriangle&getTetrahedraAroundTriangle(TriangleID i) override;
    virtual const HexahedraAroundVertex&getHexahedraAroundVertex(PointID i) override;
    virtual const HexahedraAroundEdge&getHexahedraAroundEdge(EdgeID i) override;
    virtual const HexahedraAroundQuad&getHexahedraAroundQuad(QuadID i) override;
    virtual const VerticesAroundVertex getVerticesAroundVertex(PointID i) override;
    virtual int getEdgeIndex(PointID v1, PointID v2) override;
    virtual int getTriangleIndex(PointID v1, PointID v2, PointID v3) override;
    virtual int getQuadIndex(PointID v1, PointID v2, PointID v3, PointID v4) override;
    virtual int getTetrahedronIndex(PointID v1, PointID v2, PointID v3, PointID v4) override;
    virtual int getHexahedronIndex(PointID v1, PointID v2, PointID v3, PointID v4, PointID v5, PointID v6, PointID v7, PointID v8) override;
    virtual int getVertexIndexInTriangle(const Triangle& t, PointID vertexIndex) const override;
    virtual int getEdgeIndexInTriangle(const EdgesInTriangle& t, EdgeID edgeIndex) const override;
    virtual int getVertexIndexInQuad(const Quad& t, PointID vertexIndex) const override;
    virtual int getEdgeIndexInQuad(const EdgesInQuad& t, EdgeID edgeIndex) const override;
    virtual int getVertexIndexInTetrahedron(const Tetra& t, PointID vertexIndex) const override;
    virtual int getEdgeIndexInTetrahedron(const EdgesInTetrahedron& t, EdgeID edgeIndex) const override;
    virtual int getTriangleIndexInTetrahedron(const TrianglesInTetrahedron& t, TriangleID triangleIndex) const override;
    virtual int getVertexIndexInHexahedron(const Hexa& t, PointID vertexIndex) const override;
    virtual int getEdgeIndexInHexahedron(const EdgesInHexahedron& t, EdgeID edgeIndex) const override;
    virtual int getQuadIndexInHexahedron(const QuadsInHexahedron& t, QuadID quadIndex) const override;
    virtual core::topology::Topology::Edge getLocalEdgesInTetrahedron(const PointID i) const override;
    virtual Triangle getLocalTrianglesInTetrahedron(const PointID i) const override;
    virtual core::topology::Topology::Edge getLocalEdgesInHexahedron(const PointID i) const override;
    virtual Quad getLocalQuadsInHexahedron(const PointID i) const override;
    virtual void clear() override;
    virtual void addEdge(int a, int b) override;
    virtual void addTriangle(int a, int b, int c) override;
    virtual void addQuad(int a, int b, int c, int d) override;
    virtual void addTetra(int a, int b, int c, int d) override;
    virtual void addHexa(int a, int b, int c, int d, int e, int f, int g, int h) override;
    virtual bool checkConnexity() override;
    virtual unsigned int getNumberOfConnectedComponent() override;
    virtual int getRevision() const override;
    virtual void reOrientateTriangle(TriangleID id) override;
    virtual const sofa::helper::vector<TriangleID>&getTrianglesOnBorder() override;
    virtual const sofa::helper::vector<EdgeID>&getEdgesOnBorder() override;
    virtual const sofa::helper::vector<PointID>&getPointsOnBorder() override;

    // TopologyContainer interface
protected:
    virtual void updateTopologyEngineGraph() override;

    // PointSetTopologyContainer interface
public:
    virtual unsigned int getNumberOfElements() const override;
    virtual bool checkTopology() const override;

    // EdgeSetTopologyContainer interface
protected:
    virtual void createEdgeSetArray() override;

    // TriangleSetTopologyContainer interface
public:
    virtual const VecTriangleID getConnectedElement(TriangleID elem) override;
    virtual const VecTriangleID getElementAroundElement(TriangleID elem) override;
    virtual const VecTriangleID getElementAroundElements(VecTriangleID elems) override;

protected:
    virtual void createTriangleSetArray() override;
    virtual void createEdgesInTriangleArray() override;
    virtual void createTrianglesAroundVertexArray() override;
    virtual void createTrianglesAroundEdgeArray() override;
    virtual TrianglesAroundVertex&getTrianglesAroundVertexForModification(const PointID vertexIndex) override;
    virtual TrianglesAroundEdge&getTrianglesAroundEdgeForModification(const EdgeID edgeIndex) override;

private:
	SurfaceTopologyContainer*						map_;
	Attribute<EdgesInTriangle, Face::ORBIT>			m_edgesInTriangle;
	Attribute<TrianglesAroundVertex, Vertex::ORBIT>	m_trianglesAroundVertex;
	Attribute<TrianglesAroundEdge, Edge::ORBIT>		m_trianglesAroundEdge;
	Attribute<TrianglesAroundVertex, Vertex::ORBIT>	m_orientedTrianglesAroundVertex;
};

} // namespace topology

} // namespace component

} // namespace sofa

#endif // MAPTRIANGLESETTOPOLOGYCONTAINER_H
