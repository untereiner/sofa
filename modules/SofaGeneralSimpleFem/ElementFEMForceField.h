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
#ifndef SOFA_COMPONENT_FORCEFIELD_ELEMENT_FEM_FF_H
#define SOFA_COMPONENT_FORCEFIELD_ELEMENT_FEM_FF_H

#include <sofa/core/behavior/ForceField.h>
#include <SofaBaseTopology/TopologyData.h>
#include <sofa/helper/vector.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/helper/map.h>
#include <SofaBaseTopology/VolumeTopologyContainer.h>
#include <SofaBaseTopology/CMTopologyEngine.h>
#include <SofaBaseTopology/CMTopologyDataHandler.h>
#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>


namespace sofa
{

namespace component
{

namespace forcefield
{

template < class DataTypes, int NBVert > class ElementFEMForceField : public core::behavior::ForceField< DataTypes >
{

public:
	using Topology = component::topology::VolumeTopologyContainer;
	using Vertex = Topology::Vertex;
	using Volume = Topology::Volume;
//  typedef typename boost::mpl::if_< boost::is_same< DIM, VolumeMap >, VolumeTopology, VolumeIHMTopology >::type Topology;
  SOFA_ABSTRACT_CLASS(SOFA_TEMPLATE2(ElementFEMForceField, DataTypes, NBVert),
                      SOFA_TEMPLATE(sofa::core::behavior::ForceField, DataTypes));
  typedef typename DataTypes::Coord Coord;
  typedef typename DataTypes::Deriv Deriv;
  typedef typename Coord::value_type Real;

  typedef typename DataTypes::VecCoord VecCoord;
  typedef VecCoord Vector;
  typedef typename DataTypes::VecDeriv VecDeriv;
  typedef typename DataTypes::VecReal VecReal;

  typedef core::objectmodel::Data< VecDeriv > DataVecDeriv;
  typedef core::objectmodel::Data< VecCoord > DataVecCoord;
//  typedef typename cgogn_plugin::VolumeIHMTopology Topology;
  using MarkerVolume = Topology::CellMarker<Volume::ORBIT>;
  template<cgogn::Orbit ORBIT>
  using CellMarker = Topology::CellMarker<ORBIT>;
  typedef typename Topology::Map Map;
  template<typename T>
  using VertexAttribute = Topology::Attribute<T, Vertex::ORBIT>;
  template<typename T>
  using VolumeAttribute = Topology::Attribute<T, Volume::ORBIT>;

protected:
  enum
  {
    SMALL = 0, ///< Symbol of small displacements tetrahedron solver
    LARGE = 1, ///< Symbol of corotational large displacements tetrahedron solver based on a QR decomposition    -> Nesme et al
    /// 2005 "Efficient, Physically Plausible Finite Elements"
    POLAR = 2, ///< Symbol of corotational large displacements tetrahedron solver based on a polar decomposition -> Muller et al
    /// 2004 "Interactive Virtual Materials"
    SVD = 3 ///< Symbol of corotational large displacements tetrahedron solver based on a SVD decomposition   -> inspired from
    /// Irving et al 2004 "Invertible Finite Element for Robust Simulation of Large Deformation"
  };

  /// @name Per element (tetrahedron) data
  /// @{

  /// Displacement vector (deformation of the 4 corners of a tetrahedron)
  typedef defaulttype::VecNoInit< 3 * NBVert, Real > Displacement;

  // Material stiffness matrix
  typedef defaulttype::Mat< 6, 6, Real > MaterialStiffnessMatrix;

  /// Strain-displacement matrix (12 lines, 6 columns)
  typedef defaulttype::Mat< 3 * NBVert, 6, Real > StrainDisplacementMatrix;

  // Rigid transformation (rotation) matrix
  typedef defaulttype::MatNoInit< 3, 3, Real > TransformationMatrix;

  /// Stiffness matrix
  typedef defaulttype::Mat< 3 * NBVert, 3 * NBVert, Real > StiffnessMatrix;

  /// @}

  /// Information stored for each element (tetrahedron)
  class TetrahedronInformation
  {
  public:
    /// Material stiffness matrix of a tetrahedron (D)
    MaterialStiffnessMatrix stressStrainMatrix;

    /// Strain-displacement matrix transposed (Bt)
    StrainDisplacementMatrix strainDisplacementMatrix;

    /// (Rigid) Rotation matrix (R)
    TransformationMatrix rotation;

    /// For large, polar and SVD method
    helper::fixed_array< Coord, 4 > rotatedInitialElements;

    /// For SVD method
    TransformationMatrix initialTransformation;

    TetrahedronInformation()
    {
    }

    /// Output stream
    inline friend std::ostream& operator<<(std::ostream& os, const TetrahedronInformation&)
    {
      return os;
    }

    /// Input stream
    inline friend std::istream& operator>>(std::istream& in, TetrahedronInformation&)
    {
      return in;
    }
  };

protected:
  /// Topology of the object
  Topology* m_topology;
  // CMTopologyEngine  *m_topoEngine;

public:
	VolumeAttribute< MaterialStiffnessMatrix > m_materielStiffness;
	VolumeAttribute< StrainDisplacementMatrix > m_strainDisplacement;
	VolumeAttribute< helper::fixed_array< Coord, 8 > > m_rotatedInitialElements;
	VolumeAttribute< TransformationMatrix > m_rotation;
	VolumeAttribute< TransformationMatrix > m_initialTransformation;
	VolumeAttribute< SReal > m_localStiffnessFactor;

protected:
  Data< std::string > f_method;   /// Computation method of the displacements
  Data< Real > f_poissonRatio;    /// Incompressibility
  Data< VecReal > f_youngModulus; /// Stiffness
  Data< bool > f_updateStiffnessMatrix;
  Data< bool > f_assembling;
  Data< bool > f_drawing;
  Data< bool > f_displayWholeVolume;
  Data< VecReal > f_localStiffnessFactor;
  Data< bool > f_nonlinearStiffness;

  //  VecReal youngModulusArray;

public:
  ElementFEMForceField();

  virtual ~ElementFEMForceField()
  {
  }

  virtual void init() = 0;
  virtual void reinit() = 0;
  virtual void draw(const core::visual::VisualParams* vparams) = 0;

  Topology* getTopology()
  {
    return this->m_topology;
  }

  void setTopology(Topology* topology)
  {
    this->m_topology = topology;
  }

  void setPoissonRatio(const Real& val)
  {
    this->f_poissonRatio.setValue(val);
  }

  void setYoungModulus(const VecReal& val)
  {
    this->f_youngModulus.setValue(val);
  }

  void setLocalStiffnessFactor(const VecReal& val)
  {
    this->f_localStiffnessFactor.setValue(val);
  }

  void setUpdateStiffnessMatrix(const bool& val)
  {
    this->f_updateStiffnessMatrix.setValue(val);
  }

  void setNonlinearStiffness(const bool& val)
  {
    this->f_nonlinearStiffness.setValue(val);
  }

  void setComputeGlobalMatrix(const bool& val)
  {
    this->f_assembling.setValue(val);
  }

protected:
  virtual void addForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_f, const DataVecCoord& d_x,
                        const DataVecDeriv& d_v) = 0;

  virtual void addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_df,
                         const DataVecDeriv& d_dx) = 0;

  virtual void addKToMatrix(sofa::defaulttype::BaseMatrix* m, SReal kFactor, unsigned int& offset) = 0;

  virtual void computeElementStrainDisplacementMatrix(StrainDisplacementMatrix& B, sofa::helper::vector< Coord > vecCoords)
  {
  }

  virtual void computeMaterialStiffnessMatrix(MaterialStiffnessMatrix& E, double localStiffnessFactor = 1);

public:
  void computeElementStiffnessMatrix(StiffnessMatrix& S, StiffnessMatrix& SR, const MaterialStiffnessMatrix& E,
                                     const StrainDisplacementMatrix& B, const TransformationMatrix& Rot);

  /* virtual void computeElementForce(Displacement &F,
                                    const Displacement &U,
                                    const MaterialStiffnessMatrix &E,
                                    const StrainDisplacementMatrix &B) = 0;*/

  virtual void computeElementForce(Displacement& F, const Displacement& U, const MaterialStiffnessMatrix& E,
                                   const StrainDisplacementMatrix& B, double fact = 1);

  // BaseForceField interface
public:
  virtual SReal getPotentialEnergy(const core::MechanicalParams *mparams) const
  {
//      NOT_IMPLEMENTED_METHOD;
      return SReal(0);
  }
  virtual SReal getPotentialEnergy(const core::MechanicalParams* /*mparams*/, const DataVecCoord& ) const
  {
//      NOT_IMPLEMENTED_METHOD;
      return SReal(0);
  }

  virtual std::string getTemplateName() const
  {
      return templateName(this);
  }

  static std::string templateName(const ElementFEMForceField< DataTypes, NBVert>* = NULL);


};



} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_ELEMENT_FEM_FF_H
