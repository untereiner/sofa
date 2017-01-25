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
#ifndef SOFA_COMPONENT_FORCEFIELD_ELEMENT_FEM_FF_INL
#define SOFA_COMPONENT_FORCEFIELD_ELEMENT_FEM_FF_INL

#include <SofaMiscFem/config.h>
#include <sofa/core/behavior/ForceField.inl>
#include "ElementFEMForceField.h"
#include <SofaBaseTopology/VolumeTopologyContainer.h>
#include <sofa/core/visual/VisualParams.h>
#include <SofaBaseTopology/GridTopology.h>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/simulation/Simulation.h>
#include <sofa/helper/decompose.h>
#include <sofa/helper/gl/template.h>
#include <SofaBaseTopology/TopologyData.inl>
#include <iostream>

namespace sofa
{

namespace component
{

namespace forcefield
{

template < class DataTypes, int NBVert >
std::string ElementFEMForceField< DataTypes, NBVert >::templateName(const ElementFEMForceField<DataTypes, NBVert>* )
{
    return DataTypes::Name();
}

template < class DataTypes, int NBVert >
ElementFEMForceField< DataTypes, NBVert >::ElementFEMForceField()
    : f_method(initData(&f_method, std::string("large"), "method", "\"small\", \"large\" (by QR) or \"polar\" displacements")),
      f_poissonRatio(core::objectmodel::BaseObject::initData(&f_poissonRatio, (Real)0.45, "poissonRatio", "FEM Poisson Ratio")),
      f_youngModulus(initData(&f_youngModulus, "youngModulus", "FEM Young Modulus")),
      f_localStiffnessFactor(core::objectmodel::BaseObject::initData(
          &f_localStiffnessFactor, "localStiffnessFactor", "Allow specification of different stiffness per element. If there are "
                                                           "N element and M values are specified, the youngModulus factor for "
                                                           "element i would be localStiffnessFactor[i*M/N]")),
      f_nonlinearStiffness(initData(&f_nonlinearStiffness, bool(false), "nonlinearStiffness", "non linear bla bla"))

      ,
      f_updateStiffnessMatrix(
          core::objectmodel::BaseObject::initData(&f_updateStiffnessMatrix, false, "updateStiffnessMatrix", "")),
      f_assembling(core::objectmodel::BaseObject::initData(&f_assembling, false, "computeGlobalMatrix", "")),
      f_drawing(initData(&f_drawing, true, "drawing", " draw the forcefield if true")), m_topology(NULL)
{
  this->addAlias(&f_assembling, "assembling");
  f_poissonRatio.setWidget("poissonRatio");
}

template < class DataTypes, int NBVert > void ElementFEMForceField< DataTypes, NBVert >::init()
{
  if (!m_topology)
  {
    m_topology = this->getContext()->template get< Topology>();
  }

  if (!this->mstate.get())
  {
    Inherit1::init();
  }
}

template < class DataTypes, int NBVert > void ElementFEMForceField< DataTypes, NBVert >::reinit()
{
}

template < class DataTypes, int NBVert >
void ElementFEMForceField< DataTypes, NBVert >::computeMaterialStiffnessMatrix(MaterialStiffnessMatrix& E,
                                                                               double localStiffnessFactor)
{
  const Real youngModulus = (Real)localStiffnessFactor;
  const Real poissonRatio = f_poissonRatio.getValue();

  /* Material Stiffness Matrix E[6 x 6]
  -----------------------------------------------------------
  |E[0][0]   E[0][1]   E[0][2]                              |
  |E[1][0]   E[1][1]   E[1][2]                              |
  |E[2][0]   E[2][1]   E[2][2]                              |
  |                              E[3][3]                    |
  |                                        E[4][4]          |
  |                                                  E[5][5]|
  -----------------------------------------------------------*/

  // ZER0s
  E[0][3] = E[0][4] = E[0][5] = 0;
  E[1][3] = E[1][4] = E[1][5] = 0;
  E[2][3] = E[2][4] = E[2][5] = 0;
  E[3][0] = E[3][1] = E[3][2] = E[3][4] = E[3][5] = 0;
  E[4][0] = E[4][1] = E[4][2] = E[4][3] = E[4][5] = 0;
  E[5][0] = E[5][1] = E[5][2] = E[5][3] = E[5][4] = 0;

  E[0][0] = E[1][1] = E[2][2] = 1;
  E[0][1] = E[0][2] = E[1][0] = E[1][2] = E[2][0] = E[2][1] = poissonRatio / (1 - poissonRatio);

  E[3][3] = E[4][4] = E[5][5] = (1 - 2 * poissonRatio) / (2 * (1 - poissonRatio));

  E *= (youngModulus * (1 - poissonRatio)) / ((1 + poissonRatio) * (1 - 2 * poissonRatio));
}

template < class DataTypes, int NBVert >
void ElementFEMForceField< DataTypes, NBVert >::computeElementStiffnessMatrix(StiffnessMatrix& S, StiffnessMatrix& SR,
                                                                              const MaterialStiffnessMatrix& E,
                                                                              const StrainDisplacementMatrix& B,
                                                                              const TransformationMatrix& Rot)
{
  /// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX///
  /// Pourquoi R*B*E*Bt*Rt et non pas Rt*Bt*E*B*R ??
  defaulttype::MatNoInit< 6, 3 * NBVert, Real > Bt;
  Bt.transpose(B);

  defaulttype::MatNoInit< 3 * NBVert, 3 * NBVert, Real > BEBt;
  BEBt = B * E * Bt;

  defaulttype::MatNoInit< 3 * NBVert, 3 * NBVert, Real > R, Rt;
  R.clear();
  Rt.clear();

  /* Rot[3 x 3] Rotation Matrix of the barycentric coordinate of the Element (extracted from stretching by polar decomposition)
   * R[3*NBVert x 3*NBVert] Rotation Matrix of the Element
  -------------------------------------------------------------------------------------------------------------------------|
  |Rot[0][0]  Rot[0][1]  Rot[0][2]                                                                                         |
  |Rot[1][0]  Rot[1][1]  Rot[1][2]                                                                                         |
  |Rot[2][0]  Rot[2][1]  Rot[2][2]                                                                                         |
  |                                 Rot[0][0]  Rot[0][1]  Rot[0][2]                                                        |
  |                                 Rot[1][0]  Rot[1][1]  Rot[1][2]                                                        |
  |                                 Rot[2][0]  Rot[2][1]  Rot[2][2]                                                        |
  |                                                                          ....                                          |
  |                                                                          ....                                          |
  |                                                                          ....                                          |
  |                                                                                         Rot[0][0]  Rot[0][1]  Rot[0][2]|
  |                                                                                         Rot[1][0]  Rot[1][1]  Rot[1][2]|
  |                                                                                         Rot[2][0]  Rot[2][1]  Rot[2][2]|
  |------------------------------------------------------------------------------------------------------------------------*/
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      for (int k = 0; k < NBVert; ++k)
      {
        R[i + 3 * k][j + 3 * k] = Rot[i][j];
        Rt[i + 3 * k][j + 3 * k] = Rot[j][i];
      }
    }
  }

  S = R * BEBt;
  SR = S * Rt;
}

template < class DataTypes, int NBVert >
void ElementFEMForceField< DataTypes, NBVert >::computeElementForce(Displacement& F, const Displacement& U,
                                                                    const MaterialStiffnessMatrix& E,
                                                                    const StrainDisplacementMatrix& B, double fact)
{
  /// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX///
  /// Pourquoi B*E*Bt*U et non pas Bt*E*B*U ??
  // F = B*E*Bt*U

  /* We have these zeros:
   E[6 x 6]
  -----------------------------------------------------------
  |                              E[0][3]   E[0][4]   E[0][5]|
  |                              E[1][3]   E[1][4]   E[1][5]|
  |                              E[2][3]   E[2][4]   E[2][5]|
  |E[3][0]   E[3][1]   E[3][2]             E[3][4]   E[3][5]|
  |E[4][0]   E[4][1]   E[4][2]   E[4][3]             E[4][5]|
  |E[5][0]   E[5][1]   E[5][2]   E[5][3]   E[5][4]          |
  -----------------------------------------------------------

   B[6 x 3*NBVert]
  ------------------------------------------------------------
  |          B[0][1]   B[0][2]             B[0][4]            |
  |B[1][0]             B[1][2]                        B[1][5] |
  |B[2][0]   B[2][1]             B[2][3]                      |
  |          B[3][1]   B[3][2]             B[3][4]            |
  |B[4][0]             B[4][2]                        B[4][5] |
  |B[5][0]   B[5][1]             B[5][3]                      |
  |            ...       ...                 ...              |
  |  ...                 ...                            ...   |
  |  ...       ...                 ...                        |
  |          B[N-2][1] B[N-2][2]           B[N-2][4]          |
  |B[N-1][0]           B[N-1][2]                     B[N-1][5]|
  |B[N][0]   B[N][1]             B[N][3]                      |
  ------------------------------------------------------------
  */

  defaulttype::Vec< 6, Real > BtU;
  for (unsigned int i = 0; i < NBVert; i++)
  {
    BtU[0] += B[3 * i][0] * U[3 * i] /*+  B[3*i+1][0]*U[3*i+1]  +  B[3*i+2][0]*U[3*i+2]*/;
    BtU[1] += /*B[ 3*i ][1]*U[ 3*i ]  +*/ B[3 * i + 1][1] * U[3 * i + 1] /*+  B[3*i+2][1]*U[3*i+2]*/;
    BtU[2] += /*B[ 3*i ][2]*U[ 3*i ]  +  B[3*i+1][2]*U[3*i+1]  +*/ B[3 * i + 2][2] * U[3 * i + 2];
    BtU[3] += B[3 * i][3] * U[3 * i] + B[3 * i + 1][3] * U[3 * i + 1] /*+  B[3*i+2][3]*U[3*i+2]*/;
    BtU[4] += /*B[ 3*i ][4]*U[ 3*i ]  +*/ B[3 * i + 1][4] * U[3 * i + 1] + B[3 * i + 2][4] * U[3 * i + 2];
    BtU[5] += B[3 * i][5] * U[3 * i] + /*B[3*i+1][5]*U[3*i+1]  +*/ B[3 * i + 2][5] * U[3 * i + 2];
  }

  defaulttype::VecNoInit< 6, Real > EBtU;
  EBtU[0] = E[0][0] * BtU[0] + E[0][1] * BtU[1] + E[0][2] * BtU[2];
  /*E[0][3]*BtU[3] +  E[0][4]*BtU[4] +  E[0][5]*BtU[5]*/
  EBtU[1] = E[1][0] * BtU[0] + E[1][1] * BtU[1] + E[1][2] * BtU[2];
  /*E[1][3]*BtU[3] +  E[1][4]*BtU[4] +  E[1][5]*BtU[5]*/
  EBtU[2] = E[2][0] * BtU[0] + E[2][1] * BtU[1] + E[2][2] * BtU[2];
  /*E[2][3]*BtU[3] +  E[2][4]*BtU[4] +  E[2][5]*BtU[5]*/
  EBtU[3] =                                  /*E[3][0]*BtU[0] +  E[3][1]*BtU[1] +  E[3][2]*BtU[2]+*/
      E[3][3] * BtU[3];                      /*E[3][4]*BtU[4] +  E[3][5]*BtU[5]*/
  EBtU[4] =                                  /*E[4][0]*BtU[0] +  E[4][1]*BtU[1] +  E[4][2]*BtU[2]+*/
      /*E[4][3]*BtU[3] +*/ E[4][4] * BtU[4]; /*E[4][5]*BtU[5]*/
  EBtU[5] =                                  /*E[5][0]*BtU[0] +  E[5][1]*BtU[1] +  E[5][2]*BtU[2]+*/
      /*E[5][3]*BtU[3] +  E[5][4]*BtU[4] +*/ E[5][5] * BtU[5];

  EBtU *= fact;

  for (unsigned int i = 0; i < NBVert; i++)
  {
    F[3 * i] = B[3 * i][0] * EBtU[0] + /*B[ 3*i ][1]*EBtU[1]  +  B[ 3*i ][2]*EBtU[2]+*/
               B[3 * i][3] * EBtU[3] + /*B[ 3*i ][4]*EBtU[4]  +*/ B[3 * i][5] * EBtU[5];
    F[3 * i + 1] = /*B[3*i+1][0]*EBtU[0]  +*/ B[3 * i + 1][1] * EBtU[1] + /*B[3*i+1][2]*EBtU[2]+*/
                   B[3 * i + 1][3] * EBtU[3] + B[3 * i + 1][4] * EBtU[4]; /*B[3*i+1][5]*EBtU[5]*/
    F[3 * i + 2] = /*B[3*i+2][0]*EBtU[0]  +  B[3*i+2][1]*EBtU[1]  +*/ B[3 * i + 2][2] * EBtU[2] +
                   /*B[3*i+2][3]*EBtU[3]  +*/ B[3 * i + 2][4] * EBtU[4] + B[3 * i + 2][5] * EBtU[5];
  }
}

} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_ELEMENT_FEM_FF_INL
