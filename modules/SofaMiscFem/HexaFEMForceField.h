/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 MGH, INRIA, USTL, UJF, CNRS                    *
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
#ifndef SOFA_COMPONENT_FORCEFIELD_HEXA_FEM_FF_H
#define SOFA_COMPONENT_FORCEFIELD_HEXA_FEM_FF_H

#include <sofa/core/behavior/ForceField.h>
#include "ElementFEMForceField.h"
#include <SofaBaseTopology/TopologyData.h>
#include <sofa/helper/vector.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/helper/map.h>


#include <SofaOpenglVisual/OglColorMap.h>
#include <numeric>

namespace sofa
{

namespace component
{

namespace forcefield
{

struct HexaFF_Order
{
	static const unsigned int ind[8];
};

template < class DataTypes>
class HexaFEMForceField : public virtual ElementFEMForceField< DataTypes, 8 >
{

public:
	SOFA_CLASS(SOFA_TEMPLATE(HexaFEMForceField, DataTypes),
			   SOFA_TEMPLATE2(ElementFEMForceField, DataTypes, 8));

	enum
	{
		LARGE = 0, ///< Symbol of corotational large displacements tetrahedron solver based on a QR decomposition    -> Nesme et al
		/// 2005 "Efficient, Physically Plausible Finite Elements"
		POLAR = 1, ///< Symbol of corotational large displacements tetrahedron solver based on a polar decomposition -> Muller et al
		/// 2004 "Interactive Virtual Materials"
	};

	typedef typename Inherit1::Coord Coord;
	typedef typename Inherit1::Deriv Deriv;
	typedef typename Inherit1::Real Real;

	typedef typename Inherit1::VecCoord VecCoord;
	typedef VecCoord Vector;
	typedef typename Inherit1::VecDeriv VecDeriv;
	typedef typename Inherit1::VecReal VecReal;

	typedef typename Inherit1::DataVecDeriv DataVecDeriv;
	typedef typename Inherit1::DataVecCoord DataVecCoord;
	//  typedef typename cgogn_plugin::VolumeIHMTopology Topology;
	typedef typename Inherit1::Topology Topology;
	using Vertex = typename Inherit1::Vertex;
	using Volume = typename Inherit1::Volume;
	using Dart = cgogn::Dart;
	template<typename T>
	using VolumeAttribute = typename Inherit1::template VolumeAttribute<T>;
	template<typename T>
	using VertexAttribute = typename Inherit1::template VertexAttribute<T>;
	using MarkerVolume = typename Inherit1::template CellMarker<Volume::ORBIT>;

	/// @name Per element (tetrahedron) data
	/// @{

	/// Displacement vector (deformation of the 4 corners of a tetrahedron)
	typedef typename Inherit1::Displacement Displacement;

	/// Material stiffness matrix (E)
	typedef typename Inherit1::MaterialStiffnessMatrix MaterialStiffnessMatrix;

	/// Strain-displacement matrix (B, 12 / 24 lines, 6 columns)
	typedef typename Inherit1::StrainDisplacementMatrix StrainDisplacementMatrix;

	/// Rotation matrix (R)
	typedef typename Inherit1::TransformationMatrix TransformationMatrix;

	/// Stiffness matrix (Ke)
	typedef typename Inherit1::StiffnessMatrix StiffnessMatrix;

	/// @}

	/*
  /// Information stored for each element (tetrahedron)
  class TetrahedronInformation {
  public:
	/// Material stiffness matrix of a tetrahedron (E)
	MaterialStiffnessMatrix stressStrainMatrix;

	/// Strain-displacement matrix transposed (B)
	StrainDisplacementMatrix strainDisplacementMatrix;

	/// (Rigid) Rotation matrix (R)
	TransformationMatrix rotation;

	/// For large, polar and SVD method
	helper::fixed_array<Coord, 4> rotatedInitialElements;

	/// For SVD method
	TransformationMatrix initialTransformation;

	TetrahedronInformation() {}

	/// Output stream
	inline friend std::ostream &operator<<(
		std::ostream &os, const TetrahedronInformation & ) {
	  return os;
	}

	/// Input stream
	inline friend std::istream &operator>>(std::istream &in,
										   TetrahedronInformation & ) {
	  return in;
	}
  };
  */

	/// @name Full system matrix assembly support
	/// @{

	typedef std::pair< int, SReal > Col_Value;
	typedef helper::vector< Col_Value > CompressedValue;
	typedef helper::vector< CompressedValue > CompressedMatrix;
	typedef unsigned int Index;

	CompressedMatrix _stiffnesses;

	/// @}

protected:
	/// This class is made to take care transparently of all topology changes that might
	/// happen (non exhaustive list: Edges added, removed, fused, renumbered).
	class HexahedronHandler : public core::cm_topology::TopologyElementHandler< core::topology::MapTopology::Volume >
	{
		using Inherit = core::cm_topology::TopologyElementHandler< core::topology::MapTopology::Volume >;
	public:
		HexahedronHandler(HexaFEMForceField< DataTypes >* ff) :Inherit(), hexaCoroFF(ff)
		{
		}

		virtual void applyCreateFunction(unsigned int tetraDartIndex);
		virtual void applyCreateFunction(const sofa::helper::vector< unsigned int >& tetraIndices);

	protected:
		HexaFEMForceField< DataTypes >* hexaCoroFF;
	};
	VolumeAttribute< StiffnessMatrix > m_stiffnessElementMatrix;
	VolumeAttribute< SReal > m_elemLambda;
	VolumeAttribute< SReal > m_elemMu;
	VolumeAttribute< defaulttype::Vec< 6, Real > > m_vStrainPerElement;
	VertexAttribute< defaulttype::Vec< 6, Real > > m_averagedStrainPerNode;
	VertexAttribute< defaulttype::Vec< 6, Real > > m_stdrdDevStrainPerNode;
	VolumeAttribute< SReal > m_vonMisesPerElement;
	VertexAttribute< SReal > m_averagedVonMisesPerNode;
	VertexAttribute< SReal > m_stdrdDevVonMisesPerNode;
	VertexAttribute< defaulttype::Vec4f > m_stressPerNodeColors;

	VolumeAttribute< helper::fixed_array< StrainDisplacementMatrix, 8 > > m_strainDisplMatrices;
	VolumeAttribute< Eigen::MatrixXd > m_shapeFunction;

	VolumeAttribute< SReal > m_sxx;


	Data< bool > d_computePhuocVonMises;

	Data< int > d_computeVonMisesStress;
	Data< float > d_vonMisesStressAlpha;
	Data< helper::vector< defaulttype::Vec3f > > d_vonMisesStressColorPerVertex;
	//  Data< helper::vector< SReal > > d_vonMisesStress; // per node
	//  Data< bool > d_showVonMisesStressPerNode;
	SingleLink< HexaFEMForceField, component::visualmodel::OglColorMap, BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK > l_vonMisesStressColorMap;


	helper::vector<helper::vector<unsigned int> > neighborIDVector;

	//  Data<helper::vector<Vec3f> > d_vonMisesStressColors;

	Real prevMaxStress;

	//  Data< std::string > d_showStressColorMap;
	//  component::visualmodel::ColorMap::SPtr m_showStressColorMapReal;



	int m_method;
	defaulttype::Mat< 8, 3, int > _coef;

	HexahedronHandler* m_hexaedronHandler;
	core::cm_topology::TopologyEngine* m_topoEngine;

public:
	HexaFEMForceField();

	virtual ~HexaFEMForceField();

	virtual void init();
	virtual void reinit();
	void draw(const core::visual::VisualParams* vparams);

	void computeVonMisesStress();

	void setMethod(int val)
	{
		this->m_method = val;
	}
	void setMstate(core::behavior::BaseMechanicalState* state)
	{
		this->mstate.set(dynamic_cast< core::behavior::MechanicalState< DataTypes >* >(state));
	}
	void setTopology(Topology* topology)
	{
		Inherit1::setTopology(topology);
	}

protected:
	virtual void addForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_f, const DataVecCoord& d_x,
						  const DataVecDeriv& d_v);

	virtual void addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_df, const DataVecDeriv& d_dx);

	virtual void addKToMatrix(sofa::defaulttype::BaseMatrix* m, SReal kFactor, unsigned int& offset);

	// specific to the type of element
	virtual void computeElementStrainDisplacementMatrix(StrainDisplacementMatrix& B, Coord a, Coord b, Coord c, Coord d, Coord e,
														Coord f, Coord g, Coord h);

	// Divide by volume of the element
	void computeMaterialStiffnessMatrix(Volume w, const Index& a, const Index& b, const Index& c, const Index& d,
										const Index& e, const Index& f, const Index& g, const Index& h);

	void computeMaterialStiffnessMatrix(MaterialStiffnessMatrix& E, const Index& a, const Index& b, const Index& c, const Index& d,
										const Index& e, const Index& f, const Index& g, const Index& h,
										double localStiffnessFactor = 1);

	//  void computeElementForce( Displacement &F, const Displacement &Udepl, const MaterialStiffnessMatrix &M, const Coord &a,
	//  const Coord &b, const Coord &c, const Coord &d, const Coord &e, const Coord &f, const Coord &g, const Coord &h, double
	//  stiffnessFactor = 1.0);

	void computeElementStiffnessMatrix(Volume w, const MaterialStiffnessMatrix& E, const Coord& a, const Coord& b,
									   const Coord& c, const Coord& d, const Coord& e, const Coord& f, const Coord& g,
									   const Coord& h, const double stiffnessFactor = 1.0);

	/*
  void computeElementStiffnessMatrix(StiffnessMatrix &S,
									 StiffnessMatrix &SR,
									 const MaterialStiffnessMatrix &E,
									 const StrainDisplacementMatrix &B,
									 const TransformationMatrix &Rot);
  */

	/*void computeElementForce(Displacement &F,
					const Displacement &U,
					const MaterialStiffnessMatrix &E,
					const StrainDisplacementMatrix &B);

  void computeElementForce(Displacement &F,
					const Displacement &U,
					const MaterialStiffnessMatrix &E,
					const StrainDisplacementMatrix &B,
					double fact);*/

	defaulttype::Mat< 3, 3, SReal > integrateStiffness(int signx0, int signy0, int signz0, int signx1, int signy1, int signz1, const SReal u,
													   const SReal v, const SReal w, const defaulttype::Mat< 3, 3, SReal >& J_1);

public:
	////////////// large displacements method
	void initLarge(Volume w, const Index& a, const Index& b, const Index& c, const Index& d, const Index& e, const Index& f,
				   const Index& g, const Index& h);
	void computeRotationLarge(TransformationMatrix& r, Coord& edgex, Coord& edgey);
	void accumulateForceLarge(Vector& force, const Vector& p, Volume w, const Index& a, const Index& b, const Index& c,
							  const Index& d, const Index& e, const Index& f, const Index& g, const Index& h);
	void applyStiffnessCorotational(Vector& force, const Vector& x, Volume w, const Index& a, const Index& b, const Index& c,
									const Index& d, const Index& e, const Index& f, const Index& g, const Index& h,
									double fact = 1.0);

	////////////// polar decomposition method
	void initPolar(Volume w, const Index& a, const Index& b, const Index& c, const Index& d, const Index& e, const Index& f,
				   const Index& g, const Index& h);
	void computeRotationPolar(TransformationMatrix& r, Coord& edgex, Coord& edgey, Coord& edgez);
	void accumulateForcePolar(Vector& force, const Vector& p, Volume w, const Index& a, const Index& b, const Index& c,
							  const Index& d, const Index& e, const Index& f, const Index& g, const Index& h);
};


#ifdef SOFA_EXTERN_TEMPLATE

#if !defined(SOFA_FLOAT) && !defined(SOFA_COMPONENT_FORCEFIELD_HEXA_FEM_FF_CPP)
extern template class SOFA_MISC_FEM_API HexaFEMForceField< defaulttype::Vec3dTypes >;
#endif

#if !defined(SOFA_DOUBLE) && !defined(SOFA_COMPONENT_FORCEFIELD_HEXA_FEM_FF_CPP)
extern template class SOFA_MISC_FEM_API HexaFEMForceField< defaulttype::Vec3fTypes >;
#endif

#endif

} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_HEXA_FEM_FF_H
