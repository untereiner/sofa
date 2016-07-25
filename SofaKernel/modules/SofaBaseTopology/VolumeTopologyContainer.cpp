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
#include <SofaBaseTopology/VolumeTopologyContainer.h>
#include <sofa/core/ObjectFactory.h>
#include <cgogn/io/volume_import.h>
#include <cgogn/geometry/types/eigen.h>

namespace sofa
{

namespace component
{

namespace topology
{

SOFA_DECL_CLASS(VolumeTopologyContainer)
int VolumeTopologyContainerClass = core::RegisterObject("Volume topology container")
        .add< VolumeTopologyContainer >()
        ;

VolumeTopologyContainer::VolumeTopologyContainer()
{

}

VolumeTopologyContainer::~VolumeTopologyContainer()
{

}

void VolumeTopologyContainer::initFromMeshLoader()
{
	helper::ReadAccessor< Data< VecCoord > > m_position = d_initPoints;
	helper::ReadAccessor< Data< helper::vector< Tetra > > > m_tetra = d_tetra;
	helper::ReadAccessor< Data< helper::vector< Hexa > > > m_hexa = d_hexa;

	cgogn::io::VolumeImport<Topo_Traits::MapTraits> volume_import;
	volume_import.set_nb_vertices(m_position.size());

	auto* pos_att = volume_import.template position_attribute<Eigen::Vector3d>();
	for(std::size_t i = 0ul, end = m_position.size(); i < end ; ++i)
	{
		const auto id = volume_import.insert_line_vertex_container();
		const auto& src = m_position[i];
		auto& dest = pos_att->operator [](id);
		dest[0] = src[0];
		dest[1] = src[1];
		dest[2] = src[2];
	}

	for(const Tetra& t : m_tetra.ref())
		volume_import.add_tetra(*pos_att, t[0], t[1], t[2], t[3], true);
	for(const Hexa& h : m_hexa.ref())
		volume_import.add_hexa(*pos_att, h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7], true);

	volume_import.create_map(topology_);
}

void VolumeTopologyContainer::init()
{
	topology_.clear_and_remove_attributes();
	Inherit1::init();
	initFromMeshLoader();
}

void VolumeTopologyContainer::bwdInit()
{
	Inherit1::bwdInit();
}

void VolumeTopologyContainer::reinit()
{
	Inherit1::reinit();
}

void VolumeTopologyContainer::reset()
{
	topology_.clear_and_remove_attributes();
	Inherit1::reset();
}

void VolumeTopologyContainer::cleanup()
{
	Inherit1::cleanup();
}

void VolumeTopologyContainer::draw(const core::visual::VisualParams* /*vparams*/)
{
//	TriangleSetGeometryAlgorithms<DataTypes>::draw(vparams);

//    const VecCoord& coords =mech_state_->read(core::ConstVecCoordId::position())->getValue();

//    // Draw Tetra
////    if (d_drawTetrahedra.getValue())
////    {
//	if (vparams->displayFlags().getShowWireFrame())
//		vparams->drawTool()->setPolygonMode(0, true);
//	const sofa::defaulttype::Vec4f& color_tmp = d_drawColorTetrahedra.getValue();
//	defaulttype::Vec4f color4(color_tmp[0] - 0.2f, color_tmp[1] - 0.2f, color_tmp[2] - 0.2f, 1.0);

//	const sofa::helper::vector<Tetrahedron> &tetraArray = this->m_topology->getTetrahedra();
//	std::vector<defaulttype::Vector3>   pos;
//	pos.reserve(tetraArray.size()*4u);

//	for (unsigned int i = 0; i<tetraArray.size(); ++i)
//	{
//		const Tetrahedron& tet = tetraArray[i];
//		for (unsigned int j = 0u; j<4u; ++j)
//		{
//			pos.push_back(defaulttype::Vector3(DataTypes::getCPos(coords[tet[j]])));
//		}
//	}

//	const float& scale = d_drawScaleTetrahedra.getValue();

//	if (scale >= 1.0 && scale < 0.001)
//		vparams->drawTool()->drawTetrahedra(pos, color4);
//	else
//		vparams->drawTool()->drawScaledTetrahedra(pos, color4, scale);

//	if (vparams->displayFlags().getShowWireFrame())
//		vparams->drawTool()->setPolygonMode(0, false);
////    }
}

} // namespace topology

} // namespace component

} // namespace sofa
