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
#ifndef SOFA_COMPONENT_TOPOLOGY_VOLUMETOPOLOGYCONTAINER_H
#define SOFA_COMPONENT_TOPOLOGY_VOLUMETOPOLOGYCONTAINER_H

#include "config.h"
#include <sofa/core/topology/MapTopology.h>
#include <cgogn/io/map_export.h>

namespace sofa
{

namespace component
{

namespace topology
{

//fw decl of VolumeTopologyModifier
class VolumeTopologyModifier;

class SOFA_BASE_TOPOLOGY_API VolumeTopologyContainer : public core::topology::MapTopology
{
	friend class VolumeTopologyModifier;
	public:

	SOFA_CLASS(VolumeTopologyContainer, core::topology::MapTopology);
	using Topology = Topo_Traits::Topology3;
	template<Orbit ORBIT>
	using Cell = cgogn::Cell<ORBIT>;
	using Vertex = Topo_Traits::Vertex3;
	using Edge = Topo_Traits::Edge3;
	using Face = Topo_Traits::Face3;
	using Volume = Topo_Traits::Volume3;
	using BaseVertex = Inherit1::Vertex;
	using BaseEdge   = Inherit1::Edge;
	using BaseFace   = Inherit1::Face;
	using BaseVolume = Inherit1::Volume;
	using Map = Topology;

	using CellCache = cgogn::CellCache<Map>;
	using DartMarker = cgogn::DartMarker<Map>;
	template<Orbit ORBIT>
	using CellMarker = cgogn::CellMarker<Map, ORBIT>;

	VolumeTopologyContainer();
	~VolumeTopologyContainer() override;

	virtual unsigned int getNbPoints() const override
	{
		return this->nb_cells<Vertex::ORBIT>();
	}

	template<Orbit ORBIT>
	inline uint32_t embedding(cgogn::Cell<ORBIT> c) const
	{
		return topology_.embedding(c);
	}

	// attributes
	template<typename T, typename CellType>
	inline Attribute<T,CellType::ORBIT> add_attribute(const std::string& attribute_name)
	{
		return topology_.add_attribute<T,CellType>(attribute_name);
	}

	template<typename T, typename CellType>
	inline Attribute<T,CellType::ORBIT> get_attribute(const std::string& attribute_name)
	{
		return topology_.get_attribute<T,CellType>(attribute_name);
	}

	template<typename T, typename CellType>
	inline Attribute<T,CellType::ORBIT> check_attribute(const std::string& attribute_name)
	{
		if (topology_.has_attribute(CellType::ORBIT,attribute_name))
			return topology_.get_attribute<T,CellType>(attribute_name);
		else
			return topology_.add_attribute<T,CellType>(attribute_name);
	}

	template<typename T, Orbit ORB>
	void remove_attribute(const Attribute<T,ORB>& attribute)
	{
		topology_.remove_attribute(attribute);
	}

	// MapTopology interface

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

	template<typename FUNC>
	inline void parallel_foreach_cell(const FUNC& f)
	{
		topology_.parallel_foreach_cell(f);
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


	template<typename FUNC>
	inline void foreach_incident_edge(Vertex v,const FUNC& func)
	{
		topology_.foreach_incident_edge(v,func);
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

	template<typename FUNC>
	inline void foreach_incident_volume(Vertex v,const FUNC& func)
	{
		topology_.foreach_incident_volume(v,func);
	}

	template<typename FUNC>
	inline void foreach_incident_volume(Edge e,const FUNC& func)
	{
		topology_.foreach_incident_volume(e,func);
	}

	template<typename FUNC>
	inline void foreach_incident_volume(Face f,const FUNC& func)
	{
		topology_.foreach_incident_volume(f,func);
	}

	template<typename FUNC>
	inline void foreach_incident_face(Vertex v,const FUNC& func)
	{
		topology_.foreach_incident_face(v,func);
	}

	template<typename FUNC>
	inline void foreach_incident_face(Edge e,const FUNC& func)
	{
		topology_.foreach_incident_face(e,func);
	}

	template<Orbit ORBIT>
	inline unsigned int nb_cells() const
	{
		return topology_.nb_cells<ORBIT>();
	}

	inline unsigned int get_dof(Vertex v)
	{
		return topology_.embedding(v);
	}

	inline const helper::fixed_array<unsigned int, 2>& get_dofs(Edge e) const
	{
		return this->edge_dofs_[e.dart];
	}

	inline const helper::vector<unsigned int>& get_dofs(Face f) const
	{
		return this->face_dofs_[f.dart];
	}

	inline const helper::vector<unsigned int>& get_dofs(Volume w) const
	{
		return this->volume_dofs_[w.dart];
	}

	template<Orbit ORBIT>
	inline bool same_cell(Cell<ORBIT> c1, Cell<ORBIT> c2)
	{
		return topology_.same_cell(c1, c2);
	}

protected:
	virtual void initFromMeshLoader() override;

	virtual void createEdgesInTriangleArray() override
	{

	}
	virtual void createTrianglesAroundVertexArray() override
	{

	}
	virtual void createTrianglesAroundEdgeArray() override;

	virtual void createEdgesInQuadArray() override;

	virtual void createEdgesAroundVertexArray() override;

	virtual void createEdgesInTetrahedronArray() override;

	virtual void createTrianglesInTetrahedronArray() override;

	virtual void createTetrahedraAroundTriangleArray() override;

	virtual void createTriangleSetArray() override;

	// BaseObject interface
public:
	virtual void init() override;
	virtual void bwdInit() override;
	virtual void reinit() override;
	virtual void reset() override;
	virtual void cleanup() override;
	virtual void draw(const core::visual::VisualParams*) override;

	virtual void exportMesh(const std::string& filename) override
	{
		cgogn::io::export_volume(topology_, cgogn::io::ExportOptions(filename, std::make_pair(cgogn::Orbit(Vertex::ORBIT), std::string("position"))));
	}

public:
	inline Dart phi1(Dart d) { return topology_.phi1(d); }
	inline Dart phi_1(Dart d) { return topology_.phi_1(d); }
	inline Dart phi2(Dart d) { return topology_.phi2(d); }
	inline Dart phi3(Dart d) { return topology_.phi3(d); }

	inline bool isBoundaryVertex(Vertex v) { return topology_.is_incident_to_boundary(v); }
	inline bool isBoundaryEdge(Edge e) { return topology_.is_incident_to_boundary(e); }
	inline bool isBoundaryFace(Face f) { return topology_.is_incident_to_boundary(f); }


	inline Face findBoundaryFaceOfEdge(Edge e) { return topology_.boundary_face_of_edge(e); }
	inline Face findBoundaryFaceOfVertex(Vertex v) { return topology_.boundary_face_of_Vertex(v); }

private:
	Topology topology_;
	std::unique_ptr<CellCache> cache_;
};

} // namespace topology

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_TOPOLOGY_VOLUMETOPOLOGYCONTAINER_H
