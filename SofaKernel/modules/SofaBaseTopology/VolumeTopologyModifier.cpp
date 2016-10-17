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
#include <SofaBaseTopology/VolumeTopologyModifier.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/topology/CMTopologyChange.h>
#include <SofaBaseTopology/VolumeTopologyContainer.h>
#include <algorithm>
#include <functional>
#include <iostream>
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace topology
{
SOFA_DECL_CLASS(VolumeTopologyModifier)
int VolumeTopologyModifierClass = core::RegisterObject("Volume set topology modifier")
        .add< VolumeTopologyModifier >();

using namespace std;
using namespace sofa::defaulttype;
using namespace sofa::core::topology;

//const unsigned int edgesInTetrahedronArray[6][2] = {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}};


void VolumeTopologyModifier::init()
{
	Inherit1::init();

	this->getContext()->get(m_container);
}

void VolumeTopologyModifier::reinit()
{
	Inherit1::reinit();
}


void VolumeTopologyModifier::addTetrahedra(const sofa::helper::vector<Volume> &vols)
{

    /// effectively add triangles in the topology container
    addTetrahedraProcess(vols);

    // add topology event in the stack of topological events
    addTetrahedraWarning(vols);

    // inform other objects that the edges are already added
    propagateTopologicalChanges();
}


void VolumeTopologyModifier::addTetrahedra(const sofa::helper::vector<Volume>& tetrahedra,
        const sofa::helper::vector<sofa::helper::vector<Volume> >& ancestors,
        const sofa::helper::vector<sofa::helper::vector<double> > &baryCoefs)
{
    /// effectively add triangles in the topology container
    addTetrahedraProcess(tetrahedra);

    // add topology event in the stack of topological events
    addTetrahedraWarning(tetrahedra, ancestors, baryCoefs);

    // inform other objects that the edges are already added
    propagateTopologicalChanges();
}

void VolumeTopologyModifier::addTetrahedronProcess(Volume t)
{
 // TODO
}


void VolumeTopologyModifier::addTetrahedraProcess(const sofa::helper::vector< Volume > &tetrahedra)
{
    for (unsigned int i = 0; i < tetrahedra.size(); ++i)
    {
        addTetrahedronProcess(tetrahedra[i]);
    }
}


void VolumeTopologyModifier::addTetrahedraWarning(
		const sofa::helper::vector< Volume >& tetrahedraList)
{
	//    m_container->setTetrahedronTopologyToDirty(); // TODO
	// Warning that tetrahedra just got created
	core::cm_topology::VolumesAdded *e = new core::cm_topology::VolumesAdded(tetrahedraList);
	addTopologyChange(e);
}


void VolumeTopologyModifier::addTetrahedraWarning(
		const sofa::helper::vector< Volume >& tetrahedraList,
		const sofa::helper::vector< sofa::helper::vector< Volume> > & ancestors,
		const sofa::helper::vector< sofa::helper::vector< double > >& baryCoefs)
{
	//    m_container->setTetrahedronTopologyToDirty(); // TODO
	// Warning that tetrahedra just got created
	core::cm_topology::VolumesAdded *e = new core::cm_topology::VolumesAdded(tetrahedraList, ancestors, baryCoefs);
	addTopologyChange(e);
}


void VolumeTopologyModifier::removeTetrahedraWarning( sofa::helper::vector<Volume> &tetrahedra )
{
//    m_container->setTetrahedronTopologyToDirty(); // TODO
    /// sort vertices to remove in a descendent order
//    std::sort( tetrahedra.begin(), tetrahedra.end(), std::greater<unsigned int>() );

    // Warning that these edges will be deleted
	core::cm_topology::VolumesRemoved *e=new core::cm_topology::VolumesRemoved(tetrahedra);
	addTopologyChange(e);
}

void VolumeTopologyModifier::removeTetrahedraProcess( const sofa::helper::vector<Volume> &indices,
        const bool removeIsolatedItems)
{
 // TODO
}

//void VolumeTopologyModifier::addPointsProcess(const unsigned int nPoints)
//{
//    // start by calling the parent's method.
//    TriangleSetTopologyModifier::addPointsProcess( nPoints );

//    if(m_container->hasTetrahedraAroundVertex())
//        m_container->m_tetrahedraAroundVertex.resize( m_container->getNbPoints() );
//}

//void VolumeTopologyModifier::addEdgesProcess(const sofa::helper::vector< Edge > &edges)
//{
//    // start by calling the parent's method.
//    TriangleSetTopologyModifier::addEdgesProcess( edges );

//    if(m_container->hasTetrahedraAroundEdge())
//        m_container->m_tetrahedraAroundEdge.resize( m_container->getNumberOfEdges() );
//}

//void VolumeTopologyModifier::addTrianglesProcess(const sofa::helper::vector< Triangle > &triangles)
//{
//    // start by calling the parent's method.
//    TriangleSetTopologyModifier::addTrianglesProcess( triangles );
//    if(m_container->hasTetrahedraAroundTriangle())
//        m_container->m_tetrahedraAroundTriangle.resize( m_container->getNumberOfTriangles() );
//}

//void VolumeTopologyModifier::removePointsProcess(const sofa::helper::vector<unsigned int> &indices,
//        const bool removeDOF)
//{
//    if(m_container->hasTetrahedra())
//    {
//        if(!m_container->hasTetrahedraAroundVertex())
//        {
//            m_container->createTetrahedraAroundVertexArray();
//        }

//        helper::WriteAccessor< Data< sofa::helper::vector<Tetrahedron> > > m_tetrahedron = m_container->d_tetrahedron;
//        unsigned int lastPoint = m_container->getNbPoints() - 1;
//        for (unsigned int i=0; i<indices.size(); ++i, --lastPoint)
//        {
//            // updating the edges connected to the point replacing the removed one:
//            // for all edges connected to the last point
//            for (sofa::helper::vector<unsigned int>::iterator itt=m_container->m_tetrahedraAroundVertex[lastPoint].begin();
//                    itt!=m_container->m_tetrahedraAroundVertex[lastPoint].end(); ++itt)
//            {
//                unsigned int vertexIndex = m_container->getVertexIndexInTetrahedron(m_tetrahedron[(*itt)],lastPoint);
//                m_tetrahedron[(*itt)][vertexIndex]=indices[i];
//            }

//            // updating the edge shell itself (change the old index for the new one)
//            m_container->m_tetrahedraAroundVertex[ indices[i] ] = m_container->m_tetrahedraAroundVertex[ lastPoint ];
//        }

//        m_container->m_tetrahedraAroundVertex.resize( m_container->m_tetrahedraAroundVertex.size() - indices.size() );
//    }

//    // Important : the points are actually deleted from the mechanical object's state vectors iff (removeDOF == true)
//    // call the parent's method.
//    TriangleSetTopologyModifier::removePointsProcess(  indices, removeDOF );
//}

//void VolumeTopologyModifier::removeEdgesProcess( const sofa::helper::vector<unsigned int> &indices,
//        const bool removeIsolatedItems)
//{
//    if(!m_container->hasEdges()) // this method should only be called when edges exist
//        return;

//    if (m_container->hasEdgesInTetrahedron())
//    {
//        if(!m_container->hasTetrahedraAroundEdge())
//            m_container->createTetrahedraAroundEdgeArray();

//        unsigned int lastEdge = m_container->getNumberOfEdges() - 1;
//        for (unsigned int i=0; i<indices.size(); ++i, --lastEdge)
//        {
//            for (sofa::helper::vector<unsigned int>::iterator itt=m_container->m_tetrahedraAroundEdge[lastEdge].begin();
//                    itt!=m_container->m_tetrahedraAroundEdge[lastEdge].end(); ++itt)
//            {
//                unsigned int edgeIndex=m_container->getEdgeIndexInTetrahedron(m_container->m_edgesInTetrahedron[(*itt)],lastEdge);
//                m_container->m_edgesInTetrahedron[(*itt)][edgeIndex]=indices[i];
//            }

//            // updating the edge shell itself (change the old index for the new one)
//            m_container->m_tetrahedraAroundEdge[ indices[i] ] = m_container->m_tetrahedraAroundEdge[ lastEdge ];
//        }

//        m_container->m_tetrahedraAroundEdge.resize( m_container->m_tetrahedraAroundEdge.size() - indices.size() );
//    }

//    // call the parent's method.
//    TriangleSetTopologyModifier::removeEdgesProcess( indices, removeIsolatedItems );
//}

//void VolumeTopologyModifier::removeTrianglesProcess( const sofa::helper::vector<unsigned int> &indices,
//        const bool removeIsolatedEdges,
//        const bool removeIsolatedPoints)
//{
//    if(!m_container->hasTriangles()) // this method should only be called when triangles exist
//        return;

//    if (m_container->hasTrianglesInTetrahedron())
//    {
//        if(!m_container->hasTetrahedraAroundTriangle())
//            m_container->createTetrahedraAroundTriangleArray();

//        size_t lastTriangle = m_container->m_tetrahedraAroundTriangle.size() - 1;
//        for (unsigned int i = 0; i < indices.size(); ++i, --lastTriangle)
//        {
//            for (sofa::helper::vector<unsigned int>::iterator itt=m_container->m_tetrahedraAroundTriangle[lastTriangle].begin();
//                    itt!=m_container->m_tetrahedraAroundTriangle[lastTriangle].end(); ++itt)
//            {
//                unsigned int triangleIndex=m_container->getTriangleIndexInTetrahedron(m_container->m_trianglesInTetrahedron[(*itt)],lastTriangle);
//                m_container->m_trianglesInTetrahedron[(*itt)][triangleIndex] = indices[i];
//            }

//            // updating the triangle shell itself (change the old index for the new one)
//            m_container->m_tetrahedraAroundTriangle[ indices[i] ] = m_container->m_tetrahedraAroundTriangle[ lastTriangle ];
//        }
//        m_container->m_tetrahedraAroundTriangle.resize( m_container->m_tetrahedraAroundTriangle.size() - indices.size() );
//    }

//    // call the parent's method.
//    TriangleSetTopologyModifier::removeTrianglesProcess( indices, removeIsolatedEdges, removeIsolatedPoints );
//}

//void VolumeTopologyModifier::renumberPointsProcess( const sofa::helper::vector<unsigned int> &index,
//        const sofa::helper::vector<unsigned int> &inv_index,
//        const bool renumberDOF)
//{
//    if(m_container->hasTetrahedra())
//    {
//        helper::WriteAccessor< Data< sofa::helper::vector<Tetrahedron> > > m_tetrahedron = m_container->d_tetrahedron;
//        if(m_container->hasTetrahedraAroundVertex())
//        {
//            sofa::helper::vector< sofa::helper::vector< unsigned int > > tetrahedraAroundVertex_cp = m_container->m_tetrahedraAroundVertex;
//            for (unsigned int i = 0; i < index.size(); ++i)
//            {
//                m_container->m_tetrahedraAroundVertex[i] = tetrahedraAroundVertex_cp[ index[i] ];
//            }
//        }

//        for (unsigned int i=0; i<m_tetrahedron.size(); ++i)
//        {
//            m_tetrahedron[i][0]  = inv_index[ m_tetrahedron[i][0]  ];
//            m_tetrahedron[i][1]  = inv_index[ m_tetrahedron[i][1]  ];
//            m_tetrahedron[i][2]  = inv_index[ m_tetrahedron[i][2]  ];
//            m_tetrahedron[i][3]  = inv_index[ m_tetrahedron[i][3]  ];
//        }
//    }

//    // call the parent's method.
//    TriangleSetTopologyModifier::renumberPointsProcess( index, inv_index, renumberDOF );
//}

//void VolumeTopologyModifier::removeTetrahedra(const sofa::helper::vector<unsigned int> &tetrahedraIds)
//{
//    sofa::helper::vector<unsigned int> tetrahedraIds_filtered;
//    for (unsigned int i = 0; i < tetrahedraIds.size(); i++)
//    {
//        if( tetrahedraIds[i] >= m_container->getNumberOfTetrahedra())
//            std::cout << "Error: VolumeTopologyModifier::removeTetrahedra: tetrahedra: "<< tetrahedraIds[i] <<" is out of bound and won't be removed." << std::endl;
//        else
//            tetrahedraIds_filtered.push_back(tetrahedraIds[i]);
//    }

//    removeTetrahedraWarning(tetrahedraIds_filtered);

//    // inform other objects that the triangles are going to be removed
//    propagateTopologicalChanges();

//    // now destroy the old tetrahedra.
//    removeTetrahedraProcess(tetrahedraIds_filtered ,true);

//    m_container->checkTopology();

//    m_container->addRemovedTetraIndex(tetrahedraIds_filtered);
//}

//void VolumeTopologyModifier::removeItems(const sofa::helper::vector<Volume>& items)
//{
//    removeTetrahedra(items);
//}

//void VolumeTopologyModifier::renumberPoints( const sofa::helper::vector<unsigned int> &index,
//        const sofa::helper::vector<unsigned int> &inv_index)
//{
//    /// add the topological changes in the queue
//    renumberPointsWarning(index, inv_index);
//    // inform other objects that the triangles are going to be removed
//    propagateTopologicalChanges();
//    // now renumber the points
//    renumberPointsProcess(index, inv_index);

//    m_container->checkTopology();
//}


//void VolumeTopologyModifier::propagateTopologicalEngineChanges()
//{
//    if (m_container->beginChange() == m_container->endChange()) return; // nothing to do if no event is stored

//    if (!m_container->isTetrahedronTopologyDirty()) // tetrahedron Data has not been touched
//        return TriangleSetTopologyModifier::propagateTopologicalEngineChanges();

//    sofa::helper::list <sofa::core::topology::TopologyEngine *>::iterator it;

//    for ( it = m_container->m_enginesList.begin(); it!=m_container->m_enginesList.end(); ++it)
//    {
//        sofa::core::topology::TopologyEngine* topoEngine = (*it);
//        if (topoEngine->isDirty())
//        {
//#ifndef NDEBUG
//            std::cout << "VolumeTopologyModifier::performing: " << topoEngine->getName() << std::endl;
//#endif
//            topoEngine->update();
//        }
//    }

//    m_container->cleanTetrahedronTopologyFromDirty();
//    TriangleSetTopologyModifier::propagateTopologicalEngineChanges();
//}

} // namespace topology

} // namespace component

} // namespace sofa

