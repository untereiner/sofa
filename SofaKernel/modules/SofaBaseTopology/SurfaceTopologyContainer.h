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
#ifndef SOFA_COMPONENT_TOPOLOGY_SURFACETOPOLOGYCONTAINER_H
#define SOFA_COMPONENT_TOPOLOGY_SURFACETOPOLOGYCONTAINER_H

#include "config.h"
#include <sofa/core/topology/MapTopology.h>
#include <cgogn/io/surface_import.h>
namespace sofa
{

namespace component
{

namespace topology
{

class SOFA_BASE_TOPOLOGY_API SurfaceImport : public cgogn::io::SurfaceImport<cgogn::DefaultMapTraits>
{
protected:
	virtual bool import_file_impl(const std::string& filename) override;
};

class SOFA_BASE_TOPOLOGY_API SurfaceTopologyContainer : public core::topology::MapTopology
{
    SOFA_CLASS(SurfaceTopologyContainer, core::topology::MapTopology);
    using Topology = Topo_Traits::Topology2;
    using Vertex = Topo_Traits::Vertex2;
    using Edge = Topo_Traits::Edge2;
    using Face = Topo_Traits::Face2;
    using Volume = Topo_Traits::Volume2;
    using BaseVertex = Inherit1::Vertex;
    using BaseEdge   = Inherit1::Edge;
    using BaseFace   = Inherit1::Face;
    using BaseVolume = Inherit1::Volume;

    using CellCache = cgogn::CellCache<Topology>;


    // MapTopology interface
public:

    virtual void foreach_vertex(const std::function<void (BaseVertex)>& func) override
    {
	topology_.foreach_cell([&](Vertex v) { func((v.dart));});
    }
    virtual void foreach_edge(const std::function<void (BaseEdge)>& func) override
    {
	topology_.foreach_cell([&](Edge e) { func((e.dart));});
    }
    virtual void foreach_face(const std::function<void (BaseFace)>& func) override
    {
	topology_.foreach_cell([&](Face f) { func((f.dart));});
    }
    virtual void foreach_volume(const std::function<void (BaseVolume)>& func) override
    {
	topology_.foreach_cell([&](Volume w) { func((w.dart));});
    }

    template<typename FUNC>
    inline void foreach_cell(const FUNC& f)
    {
        topology_.foreach_cell(f);
    }


    virtual void foreach_incident_vertex_of_edge(BaseEdge e, const std::function<void (BaseVertex)>& func) override
    {
	topology_.foreach_incident_vertex(Edge(e.id_),[&func](Vertex v) { func(v.dart); });
    }
    template<typename FUNC>
    inline void foreach_incident_vertex(Edge e,const FUNC& f)
    {
	topology_.foreach_incident_vertex(e,f);
    }


    virtual void foreach_incident_vertex_of_face(BaseFace f, const std::function<void (BaseVertex)>& func) override
    {
	topology_.foreach_incident_vertex(Face(f.id_),[&func](Vertex v) { func(v.dart); });
    }
    template<typename FUNC>
    inline void foreach_incident_vertex(Face f,const FUNC& func)
    {
	topology_.foreach_incident_vertex(f,func);
    }


    virtual void foreach_incident_vertex_of_volume(BaseVolume w, const std::function<void (BaseVertex)>& func) override
    {
	topology_.foreach_incident_vertex(Volume(w.id_),[&func](Vertex v) { func(v.dart); });
    }
    template<typename FUNC>
    inline void foreach_incident_vertex(Volume w,const FUNC& func)
    {
	topology_.foreach_incident_vertex(w,func);
    }

    virtual void foreach_incident_edge_of_face(BaseFace f, const std::function<void (BaseEdge)>& func) override
    {
	topology_.foreach_incident_edge(Face(f.id_),[&func](Edge e) { func(e.dart); });
    }
    template<typename FUNC>
    inline void foreach_incident_edge(Face f,const FUNC& func)
    {
        topology_.foreach_incident_edge(f,func);
    }


    virtual void foreach_incident_edge_of_volume(BaseVolume w, const std::function<void (BaseEdge)>& func) override
    {
	topology_.foreach_incident_edge(Volume(w.id_),[&func](Edge e) { func(e.dart); });
    }
    template<typename FUNC>
    inline void foreach_incident_edge(Volume vol,const FUNC& func)
    {
        topology_.foreach_incident_edge(vol,func);
    }


    virtual void foreach_incident_face_of_volume(BaseVolume w, const std::function<void (BaseFace)>& func) override
    {
	topology_.foreach_incident_face(Volume(w.id_),[&func](Face f) { func(f.dart); });
    }
    template<typename FUNC>
    inline void foreach_incident_face(Volume vol,const FUNC& func)
    {
        topology_.foreach_incident_face(vol,func);
    }

public:
	virtual void init() override;
	virtual void bwdInit() override;
	virtual void reinit() override;
	virtual void reset() override;
	virtual void cleanup() override;

protected:
	virtual void initFromMeshLoader() override;
private:
    Topology topology_;
    CellCache cache_;
};

} // namespace topology

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_TOPOLOGY_SURFACETOPOLOGYCONTAINER_H
