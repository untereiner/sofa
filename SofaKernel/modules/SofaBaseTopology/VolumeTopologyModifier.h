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
#ifndef SOFA_COMPONENT_TOPOLOGY_VOLUMETOPOLOGYMODIFIER_H
#define SOFA_COMPONENT_TOPOLOGY_VOLUMETOPOLOGYMODIFIER_H
#include "config.h"

#include <sofa/core/topology/CMBaseTopology.h>
#include <SofaBaseTopology/VolumeTopologyContainer.h>
#include <cgogn/modeling/algos/tetrahedralization.h>
#include <sofa/simulation/StateChangeVisitor.h>

namespace sofa
{

namespace component
{

namespace topology
{

/**
* A class that modifies the topology by adding and removing volumes
*/
class SOFA_BASE_TOPOLOGY_API VolumeTopologyModifier : public sofa::core::cm_topology::TopologyModifier
{
public:
	using Map = VolumeTopologyContainer::Map;
	using BaseVertex = core::topology::MapTopology::Vertex;
	using BaseEdge = core::topology::MapTopology::Edge;
	using BaseFace = core::topology::MapTopology::Face;
	using BaseVolume = core::topology::MapTopology::Volume;
	using Vertex = VolumeTopologyContainer::Vertex;
	using Edge = VolumeTopologyContainer::Edge;
	using Face = VolumeTopologyContainer::Face;
	using Volume = VolumeTopologyContainer::Volume;
	using MechanicalState = core::behavior::MechanicalState<defaulttype::Vec3Types>;
    SOFA_CLASS(VolumeTopologyModifier,sofa::core::cm_topology::TopologyModifier);


//    typedef core::topology::BaseMeshTopology::TetraID TetraID;
//    typedef core::topology::BaseMeshTopology::Tetra Tetra;
//    typedef core::topology::BaseMeshTopology::SeqTetrahedra SeqTetrahedra;
//    typedef core::topology::BaseMeshTopology::TetrahedraAroundVertex TetrahedraAroundVertex;
//    typedef core::topology::BaseMeshTopology::TetrahedraAroundEdge TetrahedraAroundEdge;
//    typedef core::topology::BaseMeshTopology::TetrahedraAroundTriangle TetrahedraAroundTriangle;
//    typedef core::topology::BaseMeshTopology::EdgesInTetrahedron EdgesInTetrahedron;
//    typedef core::topology::BaseMeshTopology::TrianglesInTetrahedron TrianglesInTetrahedron;
//    typedef Tetra Tetrahedron;


    Data< bool > removeIsolated; ///< Controlled DOF index.
protected:
    VolumeTopologyModifier()
        : Inherit1()
        , removeIsolated( initData(&removeIsolated,true, "removeIsolated", "remove Isolated dof") )
    {}

    virtual ~VolumeTopologyModifier() override {}
public:
    virtual void init() override;

    virtual void reinit() override;

	void addDOF(Vertex v, const helper::vector<Vertex>& ancestors, const helper::vector<double>& coeffs);

    /// \brief function to propagate topological change events by parsing the list of topologyEngines linked to this topology.
	void propagateTopologicalEngineChanges();

    /** \brief add a set of tetrahedra
    @param tetrahedra an array of vertex indices describing the tetrahedra to be created
    */
	virtual void addTetrahedra(const sofa::helper::vector< BaseVolume > &vols);

    /** \brief add a set of tetrahedra
    @param quads an array of vertex indices describing the tetrahedra to be created
    @param ancestors for each tetrahedron to be created provides an array of tetrahedron ancestors (optional)
    @param baryCoefs for each tetrahedron provides the barycentric coordinates (sum to 1) associated with each ancestor (optional)
    *
    */
	virtual void addTetrahedra(const sofa::helper::vector< BaseVolume > &tetrahedra,
			const sofa::helper::vector< sofa::helper::vector< BaseVolume > > & ancestors,
            const sofa::helper::vector< sofa::helper::vector< double > >& baryCoefs) ;


    /** \brief Sends a message to warn that some tetrahedra were added in this topology.
    *
    * \sa addTetrahedraProcess
    */
    void addTetrahedraWarning(
			const sofa::helper::vector< BaseVolume >& tetrahedraList);

    /** \brief Sends a message to warn that some tetrahedra were added in this topology.
    *
    * \sa addTetrahedraProcess
    */
    void addTetrahedraWarning(
			const sofa::helper::vector< BaseVolume >& tetrahedraList,
			const sofa::helper::vector< sofa::helper::vector< BaseVolume> > & ancestors,
            const sofa::helper::vector< sofa::helper::vector< double > >& baryCoefs);

    /** \brief Add a tetrahedron.
    *
    */
	void addTetrahedronProcess(BaseVolume e);

    /** \brief Actually Add some tetrahedra to this topology.
    *
    * \sa addTetrahedraWarning
    */
	virtual void addTetrahedraProcess(const sofa::helper::vector< BaseVolume > &tetrahedra);

    /** \brief Sends a message to warn that some tetrahedra are about to be deleted.
    *
    * \sa removeTetrahedraProcess
    *
    * Important : parameter indices is not const because it is actually sorted from the highest index to the lowest one.
    */
	void removeTetrahedraWarning( const sofa::helper::vector<BaseVolume> &tetrahedra);

    /** \brief Remove a subset of tetrahedra
    *
    * Elements corresponding to these points are removed form the mechanical object's state vectors.
    *
    * Important : some structures might need to be warned BEFORE the points are actually deleted, so always use method removeEdgesWarning before calling removeEdgesProcess.
    * \sa removeTetrahedraWarning
    * @param removeIsolatedItems if true remove isolated triangles, edges and vertices
    */
	virtual void removeTetrahedraProcess( const sofa::helper::vector<BaseVolume> &indices,
            const bool removeIsolatedItems=false);

    /** \brief Actually Add some triangles to this topology.
    *
    * \sa addTrianglesWarning
    */
//    virtual void addTrianglesProcess(const sofa::helper::vector< Triangle > &triangles);

//    /** \brief Remove a subset of triangles
//    *
//    * Important : some structures might need to be warned BEFORE the points are actually deleted, so always use method removeEdgesWarning before calling removeEdgesProcess.
//    * @param removeIsolatedEdges if true isolated edges are also removed
//    * @param removeIsolatedPoints if true isolated vertices are also removed
//    */
//    virtual void removeTrianglesProcess(const sofa::helper::vector<unsigned int> &indices,
//            const bool removeIsolatedEdges=false,
//            const bool removeIsolatedPoints=false);

//    /** \brief Add some edges to this topology.
//    *
//    * \sa addEdgesWarning
//    */
//    virtual void addEdgesProcess(const sofa::helper::vector< Edge > &edges);

    /** \brief Remove a subset of edges
    *
    * Important : some structures might need to be warned BEFORE the points are actually deleted, so always use method removeEdgesWarning before calling removeEdgesProcess.
    * \sa removeEdgesWarning
    *
    * Important : parameter indices is not const because it is actually sorted from the highest index to the lowest one.
    * @param removeIsolatedItems if true remove isolated vertices
    */
//    virtual void removeEdgesProcess( const sofa::helper::vector<unsigned int> &indices,
//            const bool removeIsolatedItems=false);

    /** \brief Add some points to this topology.
    *
    * \sa addPointsWarning
    */
//    virtual void addPointsProcess(const unsigned int nPoints);

    /** \brief Remove a subset of points
    *
    * Elements corresponding to these points are removed form the mechanical object's state vectors.
    *
    * Important : some structures might need to be warned BEFORE the points are actually deleted, so always use method removePointsWarning before calling removePointsProcess.
    * \sa removePointsWarning
    * Important : the points are actually deleted from the mechanical object's state vectors iff (removeDOF == true)
    */
//    virtual void removePointsProcess(const sofa::helper::vector<unsigned int> &indices, const bool removeDOF = true);

    /** \brief Reorder this topology.
    *
    * Important : the points are actually renumbered in the mechanical object's state vectors iff (renumberDOF == true)
    * \see MechanicalObject::renumberValues
    */
//    virtual void renumberPointsProcess( const sofa::helper::vector<unsigned int> &index,
//            const sofa::helper::vector<unsigned int> &/*inv_index*/,
//            const bool renumberDOF = true);

    /** \brief Remove a set  of tetrahedra
    @param tetrahedra an array of tetrahedron indices to be removed (note that the array is not const since it needs to be sorted)
    *
    */
//    virtual void removeTetrahedra(const sofa::helper::vector< unsigned int >& tetrahedraIds);

    /** \brief Generic method to remove a list of items.
    */
//    virtual void removeItems(const sofa::helper::vector<Volume> &items);

    /** \brief  Removes all tetrahedra in the ball of center "ind_ta" and of radius dist(ind_ta, ind_tb)
    */
//    void RemoveTetraBall(unsigned int ind_ta, unsigned int ind_tb);

    /** \brief Generic method for points renumbering
    */
//    virtual void renumberPoints( const sofa::helper::vector<unsigned int> &/*index*/,
//            const sofa::helper::vector<unsigned int> &/*inv_index*/);


	Vertex split1to4(Volume w, helper::vector< double > coeff = helper::vector< double >(), helper::vector< Vertex > ancestorsVertices = helper::vector< Vertex >());

	Vertex split1to3(Face f, helper::vector< double > coeff = helper::vector< double >(), helper::vector< Vertex > ancestorsVertices = helper::vector< Vertex >());

	Edge swap23(Face f);

	Vertex trianguleFace(Face f, helper::vector< Vertex > ancestorPoints = helper::vector< Vertex >(), sofa::helper::vector< double > coeffs = sofa::helper::vector< double >(), bool sendAddWarning = true, bool sendRemovalWarning = true );

	Face splitVolume(const std::vector<cgogn::Dart>& edges);

	Edge splitFace(cgogn::Dart e, cgogn::Dart f);

	helper::vector<Vertex> deleteVolume(Volume w);

	Vertex edgeBissection(Edge e, helper::vector< Vertex > ancestorPoints = helper::vector< Vertex >(), helper::vector< double > coeff = helper::vector< double >(), bool sendAddWarning = true, bool sendRemovalWarning = true);

	void updateTetrahedraAroundVertexAttributeInFF(Vertex );

	std::vector<Volume> swap32genOptimized(Edge e);

	std::vector<std::pair<Vertex, Vertex>> unsewVolumes(Face f);

	inline Map& getMap() { return m_container->topology_; }
private:
	VolumeTopologyContainer* 	m_container;
	MechanicalState* m_state;
};

} // namespace topology

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_TOPOLOGY_VOLUMETOPOLOGYMODIFIER_H