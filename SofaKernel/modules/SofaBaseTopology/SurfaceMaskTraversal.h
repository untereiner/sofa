#ifndef SOFA_COMPONENT_TOPOLOGY_SURFACEMASKTRAVERSAL_H
#define SOFA_COMPONENT_TOPOLOGY_SURFACEMASKTRAVERSAL_H

#include <sofa/core/objectmodel/BaseObject.h>
#include <SofaBaseTopology/SurfaceTopologyContainer.h>

namespace sofa
{

namespace component
{

namespace topology
{


class SurfaceMaskTraversal : public virtual core::objectmodel::BaseObject
{

public:
	SOFA_CLASS(SurfaceMaskTraversal, core::objectmodel::BaseObject);
	//SOFA_BASE_CAST_IMPLEMENTATION(SurfaceMaskTraversal)

protected:
	SurfaceMaskTraversal():BaseObject() {}
	virtual ~SurfaceMaskTraversal()
	{}

public:
	//virtual void init();
	//bool select(cgogn::Dart d);
};


class FixedConstraintMask : public SurfaceMaskTraversal
{

public:
	SOFA_CLASS(FixedConstraintMask, SurfaceMaskTraversal);
	//SOFA_BASE_CAST_IMPLEMENTATION(SurfaceMaskTraversal)

protected:
	FixedConstraintMask():BaseObject() {}
	virtual ~FixedConstraintMask()
	{}

//public:
//	bool select(cgogn::Dart d) override;
};


} //end namespace topology

} //end namespace component

} //end namespace sofa

#endif // SOFA_COMPONENT_TOPOLOGY_SURFACEMASKTRAVERSAL_H
