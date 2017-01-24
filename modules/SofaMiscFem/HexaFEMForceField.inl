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
#ifndef SOFA_COMPONENT_FORCEFIELD_HEXA_FEM_FF_INL
#define SOFA_COMPONENT_FORCEFIELD_HEXA_FEM_FF_INL

#include "ElementFEMForceField.inl"
#include "HexaFEMForceField.h"
#include <sofa/core/visual/VisualParams.h>
#include <SofaBaseTopology/GridTopology.h>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/simulation/Simulation.h>
#include <sofa/helper/decompose.h>
#include <sofa/helper/gl/template.h>
#include <SofaBaseTopology/TopologyData.inl>
#include <SofaBaseTopology/CMTopologyEngine.h>

namespace sofa
{

namespace component
{

namespace forcefield
{

template < class DataTypes>
HexaFEMForceField< DataTypes>::HexaFEMForceField() :    Inherit1()
  , m_hexaedronHandler(NULL)
  , d_computePhuocVonMises(initData(&d_computePhuocVonMises,false,"phuocPowaa","tu l'as ou tu l'as pas"))
  , d_computeVonMisesStress(initData(&d_computeVonMisesStress,int(2),"computeVonMisesStress","von Mises bla bla"))
  , d_vonMisesStressColorPerVertex(initData(&d_vonMisesStressColorPerVertex, "vonMisesStressColorPerVertex", "Vector of colors describing the VonMises stress"))
  , d_vonMisesStressAlpha(initData(&d_vonMisesStressAlpha, 1.0f, "vonMisesStressAlpha", "Alpha color component for vonMises visualisation (transparency > 0.0 to 1.0 < opacity)"))
  , l_vonMisesStressColorMap(initLink("vonMisesStressColorMap", "Path to the ColorMap component on scene to "))

//  , d_showStressColorMap(initData(&d_showStressColorMap,"showStressColorMap", "Color map used to show stress values"))
//  , m_showStressColorMapReal(sofa::core::objectmodel::New< sofa::component::visualmodel::ColorMap >())

{
  this->addAlias(&this->f_assembling, "assembling");
  this->f_poissonRatio.setWidget("poissonRatio");
  m_hexaedronHandler = new HexahedronHandler(this);

  _coef[0][0] = -1;
  _coef[1][0] = 1;
  _coef[2][0] = 1;
  _coef[3][0] = -1;
  _coef[4][0] = -1;
  _coef[5][0] = 1;
  _coef[6][0] = 1;
  _coef[7][0] = -1;
  _coef[0][1] = -1;
  _coef[1][1] = -1;
  _coef[2][1] = 1;
  _coef[3][1] = 1;
  _coef[4][1] = -1;
  _coef[5][1] = -1;
  _coef[6][1] = 1;
  _coef[7][1] = 1;
  _coef[0][2] = -1;
  _coef[1][2] = -1;
  _coef[2][2] = -1;
  _coef[3][2] = -1;
  _coef[4][2] = 1;
  _coef[5][2] = 1;
  _coef[6][2] = 1;
  _coef[7][2] = 1;
}

template < class DataTypes> HexaFEMForceField< DataTypes >::~HexaFEMForceField()
{
  delete m_hexaedronHandler;
}

template < class DataTypes>
void HexaFEMForceField< DataTypes >::HexahedronHandler::applyCreateFunction(unsigned int hexadronDartIndex)
{
  if (hexaCoroFF)
  {
    const typename Topology::Volume w{Dart(hexadronDartIndex)};

    //        hexaCoroFF->m_initialTransformation[w] = TransformationMatrix();
    //        if (hexaCoroFF->m_localStiffnessFactor.isValid())
    //        {
    //            hexaCoroFF->m_localStiffnessFactor[w] = SReal();
    //        }
    //        hexaCoroFF->m_materielStiffness[w] = MaterialStiffnessMatrix();
    //        hexaCoroFF->m_rotatedInitialElements[w] = helper::fixed_array<Coord, 8>();
    //        hexaCoroFF->m_rotation[w] = TransformationMatrix();
    //        hexaCoroFF->m_strainDisplacement[w] =  StrainDisplacementMatrix();

    const auto& t = hexaCoroFF->getTopology()->get_dofs(w);

    switch (hexaCoroFF->m_method)
    {
    case LARGE:
      hexaCoroFF->computeMaterialStiffnessMatrix(w, t[HexaFF_Order::ind[0]], t[HexaFF_Order::ind[1]], t[HexaFF_Order::ind[2]],
                                                    t[HexaFF_Order::ind[3]], t[HexaFF_Order::ind[4]], t[HexaFF_Order::ind[5]],
                                                    t[HexaFF_Order::ind[6]], t[HexaFF_Order::ind[7]]);
      hexaCoroFF->initLarge(w, t[HexaFF_Order::ind[0]], t[HexaFF_Order::ind[1]], t[HexaFF_Order::ind[2]], t[HexaFF_Order::ind[3]],
                               t[HexaFF_Order::ind[4]], t[HexaFF_Order::ind[5]], t[HexaFF_Order::ind[6]], t[HexaFF_Order::ind[7]]);
      break;
    case POLAR:
    default:
      hexaCoroFF->computeMaterialStiffnessMatrix(w, t[HexaFF_Order::ind[0]], t[HexaFF_Order::ind[1]], t[HexaFF_Order::ind[2]],
                                                    t[HexaFF_Order::ind[3]], t[HexaFF_Order::ind[4]], t[HexaFF_Order::ind[5]],
                                                    t[HexaFF_Order::ind[6]], t[HexaFF_Order::ind[7]]);
      hexaCoroFF->initPolar(w, t[HexaFF_Order::ind[0]], t[HexaFF_Order::ind[1]], t[HexaFF_Order::ind[2]], t[HexaFF_Order::ind[3]],
                               t[HexaFF_Order::ind[4]], t[HexaFF_Order::ind[5]], t[HexaFF_Order::ind[6]], t[HexaFF_Order::ind[7]]);
      break;
    }
    //    tetraCoroFF->reinit();
  }
}

template < class DataTypes>
void HexaFEMForceField< DataTypes >::HexahedronHandler::applyCreateFunction(
    const sofa::helper::vector< unsigned int >& hexadronDartIndices)
{
  std::cout << "in applyCreateFunction< multiple >: " << hexadronDartIndices.size() << " new hexahedra !" << std::endl;
  if (hexaCoroFF)
  {
    const unsigned int size = hexadronDartIndices.size();
    for (unsigned int i = 0u; i < size; ++i)
    {
      Volume w(Dart(hexadronDartIndices[i]));
      hexaCoroFF->m_initialTransformation[w] = TransformationMatrix();
      if (hexaCoroFF->m_localStiffnessFactor.is_valid())
      {
          hexaCoroFF->m_localStiffnessFactor[w] = SReal();
      }
      hexaCoroFF->m_materielStiffness[w] = MaterialStiffnessMatrix();
      hexaCoroFF->m_rotatedInitialElements[w] = helper::fixed_array<Coord, 8>();
      hexaCoroFF->m_strainDisplMatrices[w] = helper::fixed_array< sofa::defaulttype::Mat< 24, 6, Real >, 8 >();
      hexaCoroFF->m_rotation[w] = TransformationMatrix();
      hexaCoroFF->m_strainDisplacement[w] =  StrainDisplacementMatrix();
      hexaCoroFF->m_stiffnessElementMatrix[w] = StiffnessMatrix();

      const auto& t = hexaCoroFF->getTopology()->get_dofs(w);
      switch (hexaCoroFF->m_method)
      {
      case LARGE:
        hexaCoroFF->computeMaterialStiffnessMatrix(w, t[HexaFF_Order::ind[0]], t[HexaFF_Order::ind[1]], t[HexaFF_Order::ind[2]],
                                                      t[HexaFF_Order::ind[3]], t[HexaFF_Order::ind[4]], t[HexaFF_Order::ind[5]],
                                                      t[HexaFF_Order::ind[6]], t[HexaFF_Order::ind[7]]);
        hexaCoroFF->initLarge(w, t[HexaFF_Order::ind[0]], t[HexaFF_Order::ind[1]], t[HexaFF_Order::ind[2]], t[HexaFF_Order::ind[3]],
                                 t[HexaFF_Order::ind[4]], t[HexaFF_Order::ind[5]], t[HexaFF_Order::ind[6]], t[HexaFF_Order::ind[7]]);
//        hexaCoroFF->computeMaterialStiffnessMatrix(w, t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7]);
//        hexaCoroFF->initLarge(w, t[0], t[1], t[2], t[3], t[4], t[5], t[7], t[6]);
        break;
      case POLAR:
      default:
        hexaCoroFF->computeMaterialStiffnessMatrix(w, t[HexaFF_Order::ind[0]], t[HexaFF_Order::ind[1]], t[HexaFF_Order::ind[2]],
                                                      t[HexaFF_Order::ind[3]], t[HexaFF_Order::ind[4]], t[HexaFF_Order::ind[5]],
                                                      t[HexaFF_Order::ind[6]], t[HexaFF_Order::ind[7]]);
        hexaCoroFF->initPolar(w, t[HexaFF_Order::ind[0]], t[HexaFF_Order::ind[1]], t[HexaFF_Order::ind[2]], t[HexaFF_Order::ind[3]],
                           t[HexaFF_Order::ind[4]], t[HexaFF_Order::ind[5]], t[HexaFF_Order::ind[6]], t[HexaFF_Order::ind[7]]);
//        hexaCoroFF->computeMaterialStiffnessMatrix(w, t[0], t[1], t[3], t[2], t[4], t[5], t[6], t[7]);
//        hexaCoroFF->initPolar(w, t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7]);
        break;
      }
    }
  }
}

template < class DataTypes>
void HexaFEMForceField< DataTypes >::init()
{
  if (!this->mstate.get() || !this->getTopology())
  {
    Inherit1::init();
  }

  if (this->f_method.getValue() == "large")
  {
    this->setMethod(LARGE);
  }
  else
  {
    this->setMethod(POLAR);
  }

  //TODOOO
  //  this->m_topoEngine = new component::cm_topology::TopologyEngineImpl(this->getTopology(), m_hexaedronHandler);

  // volume FF attributes (links between Sofa atributes and CGOGN topoloogy, useful to retrieve attributes of specific Darts
  // during the simulation)
  this->m_materielStiffness =       this->getTopology()->template add_attribute< MaterialStiffnessMatrix, Volume >           ("hexaFemFF__MaterielStiffness");
  this->m_strainDisplacement =      this->getTopology()->template add_attribute< StrainDisplacementMatrix, Volume >          ("hexaFemFF__StrainDisplacement");
  this->m_strainDisplMatrices =     this->getTopology()->template add_attribute< helper::fixed_array< StrainDisplacementMatrix, 8>, Volume >  ("hexaFemFF__StrainDispMatrices");
  this->m_rotatedInitialElements =  this->getTopology()->template add_attribute< helper::fixed_array< Coord, 8 >, Volume >   ("hexaFemFF__RotatedInitialElements");
  this->m_rotation =                this->getTopology()->template add_attribute< TransformationMatrix, Volume >              ("hexaFemFF__Rotation");
  this->m_initialTransformation =   this->getTopology()->template add_attribute< TransformationMatrix, Volume >              ("hexaFemFF__InitialTransformation");
  this->m_localStiffnessFactor =    this->getTopology()->template add_attribute< SReal, Volume >                             ("hexaFemFF__LocalStiffnessFactor");
  this->m_stiffnessElementMatrix =  this->getTopology()->template add_attribute< StiffnessMatrix, Volume >                   ("hexaFemFF__StiffnessMatrix");

  this->m_elemLambda =              this->getTopology()->template add_attribute<SReal, Volume>                       ("hexaFemFF__ElemLambda");
  this->m_elemMu =                  this->getTopology()->template add_attribute<SReal, Volume>                       ("hexaFemFF__ElemMu");
  this->m_vStrainPerElement =       this->getTopology()->template add_attribute<defaulttype::Vec<6,Real>, Volume>    ("hexaFemFF__StrainPerElement");
  this->m_averagedStrainPerNode =   this->getTopology()->template add_attribute<defaulttype::Vec<6,Real>, Vertex>    ("hexaFemFF__AverageStrainPerNode");
  this->m_stdrdDevStrainPerNode =   this->getTopology()->template add_attribute<defaulttype::Vec<6,Real>, Vertex>    ("hexaFemFF__StdDevStrainPerNode");
  this->m_vonMisesPerElement =      this->getTopology()->template add_attribute<SReal, Volume>                       ("hexaFemFF__VonMisesPerElement");
  this->m_averagedVonMisesPerNode = this->getTopology()->template add_attribute<SReal, Vertex>                       ("hexaFemFF__VonMisesPerNode");
  this->m_stdrdDevVonMisesPerNode = this->getTopology()->template add_attribute<SReal, Vertex>                       ("hexaFemFF__StdDevVonMisesPerNode");
  this->m_stressPerNodeColors =     this->getTopology()->template add_attribute<defaulttype::Vec4f, Vertex>          ("hexaFemFF__StressPerNodeColors");

  this->m_shapeFunction =           this->getTopology()->template add_attribute<Eigen::MatrixXd, Volume>             ("hexaFemFF__ShapeFunction");


  this->m_sxx =                     this->getTopology()->template add_attribute< SReal, Volume >                             ("hexaFemFF__Sxx");


  if (this->getTopology()->getNbHexahedra() == 0)
  {
    serr << "ERROR(HexaFEMForceField): no hexahedron found in the CMTopology." << sendl;
    return;
  }

  reinit(); // compute per-element stiffness matrices and other precomputed values

}

template < class DataTypes>
void HexaFEMForceField< DataTypes >::reinit()
{
  const helper::ReadAccessor<Data<VecCoord> >& X0 =  this->mstate->read(core::VecCoordId::restPosition());

  neighborIDVector.clear();
  neighborIDVector.resize(X0.size());

  for(unsigned int i=0; i<X0.size(); i++)
  {
      // Same, but set to -1, to detect how many points are really surrounding
      neighborIDVector[i].resize(57);
      neighborIDVector[i][0] = i;
      for(unsigned int j=1; j<57; j++)
      {
          neighborIDVector[i][j] = -1;
      }
  }


//  typedef typename  Topology::TraversorVolumes TraversorVolumes;
//  TraversorVolumes traVol(this->getTopology()->newTraversorVolumes());
  this->m_topology->foreach_cell([this,&X0](Volume w)
  {
	  if (this->getTopology()->is_hexa(w))
	  {
		  const auto& h = this->getTopology()->get_dofs(w);
		  switch (m_method)
		  {
			  case LARGE:
				  this->computeMaterialStiffnessMatrix(w, h[HexaFF_Order::ind[0]], h[HexaFF_Order::ind[1]], h[HexaFF_Order::ind[2]],
						  h[HexaFF_Order::ind[3]], h[HexaFF_Order::ind[4]], h[HexaFF_Order::ind[5]],
						  h[HexaFF_Order::ind[6]], h[HexaFF_Order::ind[7]]);
				  this->initLarge(w, h[HexaFF_Order::ind[0]], h[HexaFF_Order::ind[1]], h[HexaFF_Order::ind[2]], h[HexaFF_Order::ind[3]],
						  h[HexaFF_Order::ind[4]], h[HexaFF_Order::ind[5]], h[HexaFF_Order::ind[6]], h[HexaFF_Order::ind[7]]);
				  break;
			  default:
			  case POLAR:
				  this->computeMaterialStiffnessMatrix(w, h[HexaFF_Order::ind[0]], h[HexaFF_Order::ind[1]], h[HexaFF_Order::ind[2]],
						  h[HexaFF_Order::ind[3]], h[HexaFF_Order::ind[4]], h[HexaFF_Order::ind[5]],
						  h[HexaFF_Order::ind[6]], h[HexaFF_Order::ind[7]]);
				  this->initPolar(w, h[HexaFF_Order::ind[0]], h[HexaFF_Order::ind[1]], h[HexaFF_Order::ind[2]], h[HexaFF_Order::ind[3]],
						  h[HexaFF_Order::ind[4]], h[HexaFF_Order::ind[5]], h[HexaFF_Order::ind[6]], h[HexaFF_Order::ind[7]]);
				  break;
		  }

		  Real dX = X0[h[1]][0] - X0[h[0]][0];
		  Real dY = X0[h[3]][1] - X0[h[0]][1];
		  Real dZ = X0[h[4]][2] - X0[h[0]][2];

		  Eigen::MatrixXd matVert(8,3);

		  matVert(0,0) = -dY * dZ / 4.0 / (dX*dY*dZ);
		  matVert(0,1) = -dX * dZ / 4.0 / (dX*dY*dZ);
		  matVert(0,2) = -dY * dX / 4.0 / (dX*dY*dZ);

		  matVert(1,0) =  dY * dZ / 4.0 / (dX*dY*dZ);
		  matVert(1,1) = -dX * dZ / 4.0 / (dX*dY*dZ);
		  matVert(1,2) = -dY * dX / 4.0 / (dX*dY*dZ);

		  matVert(2,0) =  dY * dZ / 4.0 / (dX*dY*dZ);
		  matVert(2,1) =  dX * dZ / 4.0 / (dX*dY*dZ);
		  matVert(2,2) = -dY * dX / 4.0 / (dX*dY*dZ);

		  matVert(3,0) = -dY * dZ / 4.0 / (dX*dY*dZ);
		  matVert(3,1) =  dX * dZ / 4.0 / (dX*dY*dZ);
		  matVert(3,2) = -dY * dX / 4.0 / (dX*dY*dZ);

		  matVert(4,0) = -dY * dZ / 4.0 / (dX*dY*dZ);
		  matVert(4,1) = -dX * dZ / 4.0 / (dX*dY*dZ);
		  matVert(4,2) =  dY * dX / 4.0 / (dX*dY*dZ);

		  matVert(5,0) =  dY * dZ / 4.0 / (dX*dY*dZ);
		  matVert(5,1) = -dX * dZ / 4.0 / (dX*dY*dZ);
		  matVert(5,2) =  dY * dX / 4.0 / (dX*dY*dZ);

		  matVert(6,0) =  dY * dZ / 4.0 / (dX*dY*dZ);
		  matVert(6,1) =  dX * dZ / 4.0 / (dX*dY*dZ);
		  matVert(6,2) =  dY * dX / 4.0 / (dX*dY*dZ);

		  matVert(7,0) = -dY * dZ / 4.0 / (dX*dY*dZ);
		  matVert(7,1) =  dX * dZ / 4.0 / (dX*dY*dZ);
		  matVert(7,2) =  dY * dX / 4.0 / (dX*dY*dZ);


		  this->m_shapeFunction[w] = matVert;
		  //      std::cout << m_shapeFunction[w] << std::endl;

		  //        std::cout<<  "-------------"  << std::endl;
		  //        std::cout<<  "vol = " << dX*dY*dZ  << std::endl;
		  //        std::cout<<  dX  << std::endl;
		  //        std::cout<<  dY<< std::endl;
		  //        std::cout<<  dZ  << std::endl;
		  //        std::cout<<  matVert  << std::endl;
		  //        std::cout<<  "-------------"  << std::endl;
	  }
  });

  if (d_computeVonMisesStress.getValue() > 0)
  {
      /// initialization of data structures for vonMises stress computations
      prevMaxStress = -1.0;

      helper::WriteAccessor<Data<helper::vector<defaulttype::Vec3f> > > vonMisesStressColors(d_vonMisesStressColorPerVertex);
      vonMisesStressColors.clear();
      vonMisesStressColors.resize(this->getTopology()->getNbPoints());

      if (l_vonMisesStressColorMap == NULL)
      {
          l_vonMisesStressColorMap.set( sofa::core::objectmodel::New< component::visualmodel::OglColorMap>() );
      }
  }

}

template< class DataTypes>
void HexaFEMForceField< DataTypes >::computeVonMisesStress()
{
    if (d_computeVonMisesStress.getValue() < 1)
    {
        return;
    }

    helper::ReadAccessor<Data<VecCoord> > X =  this->mstate->readPositions();
    helper::ReadAccessor<Data<VecCoord> > X0 =  this->mstate->readRestPositions();


///////////////////begin hugo///////
    VecCoord U2;
    U2.resize(X.size());
    for (size_t i = 0; i < X0.size(); i++)
        U2[i] = X[i] - X0[i];
//////////////////end hugo//////////

    defaulttype::Mat< 24, 1, Real > U;
    helper::fixed_array<double,8> vonMisesStressAtIP;
    helper::fixed_array<double,8> sXXAtIP;
    defaulttype::Mat< 6, 24, Real > Bt;
	this->m_topology->foreach_cell([&,this](Volume w)
	{
		const auto& hexa = this->getTopology()->get_dofs(w) ;
		if(!d_computePhuocVonMises.getValue())
		{

			Eigen::Matrix3d gradU;
			const Eigen::MatrixXd& shf = this->m_shapeFunction[w];

			/// compute gradU
			for (size_t k = 0; k < 3; k++)
			{
				for (size_t l = 0; l < 3; l++)
				{
					gradU(k,l) = 0.0;
					for (size_t m = 0; m < 8; m++)
						gradU(k,l) += shf(m,l) * U2[hexa[m]][k];
				}
			}

			Eigen::Matrix3d strainPerElement = ((Real)0.5)*(gradU + gradU.transpose() + gradU.transpose()*gradU);

			for (size_t i = 0; i < 3; i++)
				m_vStrainPerElement[w][i] = strainPerElement(i,i);
			m_vStrainPerElement[w][3] = strainPerElement(1,2);
			m_vStrainPerElement[w][4] = strainPerElement(0,2);
			m_vStrainPerElement[w][5] = strainPerElement(0,1);

			Real lambda = m_elemLambda[w];
			Real mu = m_elemMu[w];

			/// stress
			defaulttype::Vec<6,Real> s;// = m_vStrainPerElement[w];

			Real traceStrain = 0.0;
			for (size_t k = 0; k < 3; k++) {
				traceStrain += m_vStrainPerElement[w][k];
				s[k] = m_vStrainPerElement[w][k]*2*mu;
			}

			for (size_t k = 3; k < 6; k++)
				s[k] = m_vStrainPerElement[w][k]*2*mu;

			for (size_t k = 0; k < 3; k++)
				s[k] += lambda*traceStrain;

			m_sxx[w] = s[2];

			m_vonMisesPerElement[w] = helper::rsqrt(s[0]*s[0] + s[1]*s[1] + s[2]*s[2] - s[0]*s[1] - s[1]*s[2] - s[2]*s[0] + 3*s[3]*s[3] + 3*s[4]*s[4] + 3*s[5]*s[5]);
		}
		else
		{
			for (size_t i = 0; i < 8; i++)
			{
				U[3*i][0]   = (X[hexa[i]] - X0[hexa[i]])[0];
				U[3*i+1][0] = (X[hexa[i]] - X0[hexa[i]])[1];
				U[3*i+2][0] = (X[hexa[i]] - X0[hexa[i]])[2];
			}

			for (int integrationPointId = 0; integrationPointId < 8; ++integrationPointId)
			{
				Bt.transpose(this->m_strainDisplMatrices[w][integrationPointId]);

				defaulttype::Mat< 6, 1, Real > strain = Bt * U;

				for (size_t i = 0; i < 6; i++)
					m_vStrainPerElement[w][i] = strain[i][0];

				Real lambda = m_elemLambda[w];
				Real mu = m_elemMu[w];

				/// stress
				defaulttype::Vec<6,Real> s;// = m_vStrainPerElement[w];

				Real traceStrain = 0.0;
				for (size_t k = 0; k < 3; k++) {
					traceStrain += m_vStrainPerElement[w][k];
					s[k] = m_vStrainPerElement[w][k]*2*mu;
				}

				for (size_t k = 3; k < 6; k++)
					s[k] = m_vStrainPerElement[w][k]*2*mu;

				for (size_t k = 0; k < 3; k++)
					s[k] += lambda*traceStrain;

				vonMisesStressAtIP[integrationPointId] = helper::rsqrt(s[0]*s[0] + s[1]*s[1] + s[2]*s[2] - s[0]*s[1] - s[1]*s[2] - s[2]*s[0] + 3*s[3]*s[3] + 3*s[4]*s[4] + 3*s[5]*s[5]);
				sXXAtIP[integrationPointId] = s[2];
			}
			m_sxx[w] = std::accumulate(sXXAtIP.begin(),sXXAtIP.end(),0.0)/8.0;

			m_vonMisesPerElement[w] = std::accumulate(vonMisesStressAtIP.begin(),vonMisesStressAtIP.end(),0.0)/8.0;
		}

		if (m_vonMisesPerElement[w] < 1e-10)
			m_vonMisesPerElement[w] = 0.0;
	});

    defaulttype::Vec<6,Real> subStrainPerNodePerElem;
    defaulttype::Vec<6,Real> stdDevStrain;

    SReal subStrainPerNodePerElemReal;
    SReal stdDevStrainReal;


    /// compute the values of vonMises stress in nodes
    SReal maxVonMisesPerNode = 0.0f;
    SReal minVMN = (SReal)1e20, maxVMN = (SReal)-1e20;

	this->m_topology->foreach_cell([&,this](Vertex v)
	{
		m_averagedVonMisesPerNode[v] = 0.0;
		m_averagedStrainPerNode[v] = defaulttype::Vec<6,Real>();
		int nbrVol = 0;
		this->m_topology->foreach_incident_volume(v, [&,this](Volume w2)
		{
			m_averagedVonMisesPerNode[v] += m_vonMisesPerElement[w2];
			m_averagedStrainPerNode[v] += m_vStrainPerElement[w2];
			++nbrVol;
		});
		if(nbrVol>0)
		{
			m_averagedVonMisesPerNode[v] /= SReal(nbrVol);
			m_averagedStrainPerNode[v] /= SReal(nbrVol);
		}

		if(maxVonMisesPerNode < m_averagedVonMisesPerNode[v])
		{
			maxVonMisesPerNode = m_averagedVonMisesPerNode[v];
		}

		stdDevStrain = defaulttype::Vec<6,SReal>();
		subStrainPerNodePerElem = defaulttype::Vec<6,SReal>();
		stdDevStrainReal = 0.0f;
		subStrainPerNodePerElemReal = 0.0f;
		this->m_topology->foreach_incident_volume(v, [&,this](Volume w2)
		{
			stdDevStrain = (m_averagedStrainPerNode[v] - m_vStrainPerElement[w2]);
			std::transform( stdDevStrain.begin(), stdDevStrain.end(), stdDevStrain.begin(), stdDevStrain.begin(), std::multiplies<double>());
			subStrainPerNodePerElem += stdDevStrain;

			stdDevStrainReal = ( m_averagedVonMisesPerNode[v] - m_vonMisesPerElement[w2] );
			subStrainPerNodePerElemReal += stdDevStrainReal * stdDevStrainReal;
		});

		if(nbrVol>0)
		{
			subStrainPerNodePerElem /= nbrVol;
			std::transform(subStrainPerNodePerElem.begin(), subStrainPerNodePerElem.end(), subStrainPerNodePerElem.begin(), (SReal(*)(SReal)) sqrt);
			m_stdrdDevStrainPerNode[v] = subStrainPerNodePerElem;

			subStrainPerNodePerElemReal /=nbrVol;
			m_stdrdDevVonMisesPerNode[v] = sqrt(subStrainPerNodePerElemReal);

		}
		else {
			m_stdrdDevStrainPerNode[v] = defaulttype::Vec<6,Real>();
			m_stdrdDevVonMisesPerNode[v] = 0.0f;
		}

		minVMN = ( m_stdrdDevVonMisesPerNode[v] < minVMN) ? m_stdrdDevVonMisesPerNode[v] : minVMN;
		maxVMN = ( m_stdrdDevVonMisesPerNode[v] > maxVMN) ? m_stdrdDevVonMisesPerNode[v] : maxVMN;
	});

    minVMN /=maxVonMisesPerNode;
    maxVMN /=maxVonMisesPerNode;

    /// vonMises stress
    Real minVM = (Real)1e20, maxVM = (Real)-1e20;

	this->m_topology->foreach_cell([&,this](Volume w)
	{
		minVM = (m_vonMisesPerElement[w] < minVM) ? m_vonMisesPerElement[w] : minVM;
		maxVM = (m_vonMisesPerElement[w] > maxVM) ? m_vonMisesPerElement[w] : maxVM;
	});

	if (maxVM < prevMaxStress)
		maxVM = prevMaxStress;


    helper::WriteAccessor<Data<helper::vector<defaulttype::Vec3f> > > vonMisesStressColors(d_vonMisesStressColorPerVertex);
    vonMisesStressColors.clear();
    vonMisesStressColors.resize(this->getTopology()->getNbPoints());



	sofa::helper::ColorMap::evaluator<SReal> evalColor = l_vonMisesStressColorMap->getEvaluator(minVMN, maxVMN);
	this->m_topology->foreach_cell([&,this](Vertex v)
	{
        defaulttype::Vec4f col = evalColor(m_stdrdDevVonMisesPerNode[v]/maxVonMisesPerNode); //*vM[i]);
        defaulttype::Vec3f col3(col[0], col[1], col[2]);
        m_stressPerNodeColors[v] = col;

        vonMisesStressColors[this->getTopology()->get_dof(v)] += col3;
    });

//    std::cout<<"   -- deflexion X[12] = "<<X[12]<<std::endl;
}

template < class DataTypes>
void HexaFEMForceField< DataTypes >::addForce(const core::MechanicalParams* /* mparams */, DataVecDeriv& d_f,
											  const DataVecCoord& d_x, const DataVecDeriv& /* d_v */)
{
	VecDeriv& f = *d_f.beginEdit();
	const VecCoord& p = d_x.getValue();

	switch (m_method)
	{
		case LARGE:
		{
			this->m_topology->foreach_cell([&,this](Volume w)
			{
				assert(this->getTopology()->is_hexa(w));
				{
					const auto& t = this->getTopology()->get_dofs(w);
					accumulateForceLarge(f, p, w, t[HexaFF_Order::ind[0]], t[HexaFF_Order::ind[1]], t[HexaFF_Order::ind[2]],
							t[HexaFF_Order::ind[3]], t[HexaFF_Order::ind[4]], t[HexaFF_Order::ind[5]], t[HexaFF_Order::ind[6]],
							t[HexaFF_Order::ind[7]]);
				}
			});
			break;
		}
		case POLAR:
		{
			this->m_topology->foreach_cell([&,this](Volume w)
			{
				if (this->getTopology()->is_hexa(w))
				{
					const auto& t = this->getTopology()->get_dofs(w);
					accumulateForcePolar(f, p, w, t[HexaFF_Order::ind[0]], t[HexaFF_Order::ind[1]], t[HexaFF_Order::ind[2]],
							t[HexaFF_Order::ind[3]], t[HexaFF_Order::ind[4]], t[HexaFF_Order::ind[5]], t[HexaFF_Order::ind[6]],
							t[HexaFF_Order::ind[7]]);
				}
			});
			break;
		}
	}
	d_f.endEdit();

	this->computeVonMisesStress();
}

template < class DataTypes>
void HexaFEMForceField< DataTypes >::addDForce(const core::MechanicalParams* mparams, DataVecDeriv& d_df,
                                               const DataVecDeriv& d_dx)
{

	VecDeriv& df = *d_df.beginEdit();
	const VecDeriv& dx = d_dx.getValue();

	Real kFactor = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());

	this->m_topology->foreach_cell([&,this](Volume w)
	{
		assert(this->getTopology()->is_hexa(w));
		{
			const auto& t = this->getTopology()->get_dofs(w);
			applyStiffnessCorotational(df, dx, w, t[HexaFF_Order::ind[0]], t[HexaFF_Order::ind[1]], t[HexaFF_Order::ind[2]],
					t[HexaFF_Order::ind[3]], t[HexaFF_Order::ind[4]], t[HexaFF_Order::ind[5]],
					t[HexaFF_Order::ind[6]], t[HexaFF_Order::ind[7]], kFactor);
		}
	});

	d_df.endEdit();
}

template < class DataTypes>
void HexaFEMForceField< DataTypes >::addKToMatrix(sofa::defaulttype::BaseMatrix* mat, SReal kFactor, unsigned int& offset)
{
    unsigned int i, j, n1, n2, row, column, ROW, COLUMN;
    Index node1, node2;

	this->m_topology->foreach_cell([&,this](Volume w)
	{
        assert(this->getTopology()->is_hexa(w));
        {
            // this->Inherit1::computeElementStiffnessMatrix(BEBt,tmp,this->m_materielStiffness[w],
            // this->m_strainDisplacement[w],this->m_rotation[w]);

            const auto& t = this->getTopology()->get_dofs(w);

            TransformationMatrix& Rot = this->m_rotation[w];

            // find index of node 1
            for (n1 = 0u; n1 < 8u; ++n1)
            {
                node1 = t[n1];

                for (i = 0u; i < 3u; ++i)
                {
                    ROW = offset + 3u * node1 + i;
                    row = 3u * n1 + i;

                    // find index of node 2
                    for (n2 = 0u; n2 < 8u; ++n2)
                    {
                        node2 = t[n2];
                        const StiffnessMatrix& sm = this->m_stiffnessElementMatrix[w];
                        const defaulttype::Mat< 3, 3, Real >& tmp =
                                Rot.multTranspose(defaulttype::Mat< 3, 3, Real >(Coord(sm[3 * n1 + 0][3 * n2 + 0],
                                                  sm[3 * n1 + 0][3 * n2 + 1],
                                sm[3u * n1 + 0u][3u * n2 + 2u]),
                                Coord(sm[3 * n1 + 1u][3u * n2 + 0u],
                                sm[3u * n1 + 1u][3u * n2 + 1u],
                                sm[3u * n1 + 1u][3u * n2 + 2u]),
                                Coord(sm[3u * n1 + 2u][3u * n2 + 0u],
                                sm[3u * n1 + 2u][3u * n2 + 1u],
                                sm[3u * n1 + 2u][3u * n2 + 2u]))) *
                                Rot;

                        for (j = 0u; j < 3u; j++)
                        {
                            COLUMN = offset + 3u * node2 + j;
                            column = 3u * n2 + j;
                            mat->add(ROW, COLUMN, -tmp[i][j] * kFactor);
                        }
                    }
                }
            }
        }
    });
}

template < class DataTypes>
void HexaFEMForceField< DataTypes >::computeElementStrainDisplacementMatrix(StrainDisplacementMatrix& B, Coord a, Coord b,
                                                                            Coord c, Coord d, Coord e, Coord f, Coord g, Coord h)
{

  /* B[6 x 24 Strain-Displacement Matrix    !! Mat<L,C> !!
      ------------------------------------------------------------
      |B[0][0]                       B[0][3]             B[0][5] |
      |          B[1][1]             B[1][3]   B[1][4]           |
      |                    B[2][2]             B[2][4]   B[2][5] |
      |B[3][0]                       B[3][3]             B[3][5] |
      |          B[4][1]             B[4][3]   B[4][4]           |
      |                    B[5][2]             B[5][4]   B[5][5] |
      |  ...                           ...                 ...   |
      |            ...                 ...       ...             |
      |                      ...                 ...       ...   |
      |B[21][0]                      B[21][3]            B[21][5]|
      |          B[22][1]            B[22][3]  B[22][4]          |
      |                    B[23][2]            B[23][4]  B[23][5]|
      ---------------------------------------------------------- */

  // compute the elements of B
  // which are Nn_dx, Nn_dy, and Nn_dz for n = 0...7

  // B = [ B0 B1 B2 B3 B4 B5 B6 B7 ]T
  // where Bn (n = 0,1,2,3,4,5,6,7) are
  //      | bn  0  0 |T
  //      |  0 cn  0 |
  // Bn = |  0  0 dn |
  //      | cn bn  0 |
  //      |  0 dn cn |
  //      | dn  0 bn |
  // where
  // | bn |            | Nn_dxi   |
  // | cn | = J_n^{-1} | Nn_deta  |
  // | dn |            | Nn_dzeta |
}

template < class DataTypes>
void HexaFEMForceField< DataTypes >::computeElementStiffnessMatrix(Volume w, const MaterialStiffnessMatrix& M, const Coord& a,
                                                                   const Coord& b, const Coord& c, const Coord& d, const Coord& e,
                                                                   const Coord& f, const Coord& g, const Coord& h,
                                                                   const double stiffnessFactor)
{
  using defaulttype::Mat;
  Mat< 24, 24, Real > K;
  K.fill(0.0);

  Mat< 3, 3, Real > J;   // J[i][j] = dXi/dxj
  Mat< 3, 3, Real > J_1; // J_1[i][j] = dxi/dXj
  Mat< 3, 3, Real > J_1t;
  Real detJ = (Real)1.0;
  // check if the hexaedra is a parallelepiped
  Coord lx = b - a;
  Coord ly = d - a;
  Coord lz = e - a;
  bool isParallel = false;
  if ((d + lx - c).norm() < lx.norm() * 0.001 && (a + lz - e).norm() < lz.norm() * 0.001 &&
      (b + lz - f).norm() < lz.norm() * 0.001 && (c + lz - g).norm() < lz.norm() * 0.001 &&
      (d + lz - h).norm() < lz.norm() * 0.001)
  {
    isParallel = true;
    for (int k = 0; k < 3; ++k)
    {
      J[k][0] = lx[k] / 2;
      J[k][1] = ly[k] / 2;
      J[k][2] = lz[k] / 2;
    }
    detJ = defaulttype::determinant(J);
    J_1.invert(J);
    J_1t.transpose(J_1);
  }

  const Real U = M[0][0];
  const Real V = M[0][1];
  const Real W = M[3][3];
  int integrationPointId = 0;
  const double inv_sqrt3 = 1.0 / sqrt(3.0);
  for (int gx1 = -1; gx1 <= 1; gx1 += 2)
  {
    for (int gx2 = -1; gx2 <= 1; gx2 += 2)
    {
      for (int gx3 = -1; gx3 <= 1; gx3 += 2)
      {
        double x1 = gx1 * inv_sqrt3;
        double x2 = gx2 * inv_sqrt3;
        double x3 = gx3 * inv_sqrt3;
        if (!isParallel)
        {
          for (int k = 0; k < 3; ++k)
          {
            J[k][0] = (Real)((b[k] - a[k]) * (1 - x2) * (1 - x3) / 8 + (c[k] - d[k]) * (1 + x2) * (1 - x3) / 8 +
                             (f[k] - e[k]) * (1 - x2) * (1 + x3) / 8 + (g[k] - h[k]) * (1 + x2) * (1 + x3) / 8);
            J[k][1] = (Real)((d[k] - a[k]) * (1 - x1) * (1 - x3) / 8 + (c[k] - b[k]) * (1 + x1) * (1 - x3) / 8 +
                             (h[k] - e[k]) * (1 - x1) * (1 + x3) / 8 + (g[k] - f[k]) * (1 + x1) * (1 + x3) / 8);
            J[k][2] = (Real)((e[k] - a[k]) * (1 - x1) * (1 - x2) / 8 + (f[k] - b[k]) * (1 + x1) * (1 - x2) / 8 +
                             (g[k] - c[k]) * (1 + x1) * (1 + x2) / 8 + (h[k] - d[k]) * (1 - x1) * (1 + x2) / 8);
          }
          detJ = defaulttype::determinant(J);
          J_1.invert(J);
          J_1t.transpose(J_1);
        }

        Real qx[8];
        Real qy[8];
        Real qz[8];
        StrainDisplacementMatrix strainDisplacementMatrix;
        for (int i = 0; i < 8; ++i)
        {
          // Ni = 1/8 (1+_coef[i][0]x1)(1+_coef[i][1]x2)(1+_coef[i][2]x3)
          // qxi = dNi/dx = dNi/dx1 dx1/dx + dNi/dx2 dx2/dx + dNi/dx3 dx3/dx
          Real dNi_dx1 = (Real)((_coef[i][0]) * (1 + _coef[i][1] * x2) * (1 + _coef[i][2] * x3) / 8.0);
          Real dNi_dx2 = (Real)((1 + _coef[i][0] * x1) * (_coef[i][1]) * (1 + _coef[i][2] * x3) / 8.0);
          Real dNi_dx3 = (Real)((1 + _coef[i][0] * x1) * (1 + _coef[i][1] * x2) * (_coef[i][2]) / 8.0);
          qx[i] = dNi_dx1 * J_1[0][0] + dNi_dx2 * J_1[1][0] + dNi_dx3 * J_1[2][0];
          qy[i] = dNi_dx1 * J_1[0][1] + dNi_dx2 * J_1[1][1] + dNi_dx3 * J_1[2][1];
          qz[i] = dNi_dx1 * J_1[0][2] + dNi_dx2 * J_1[1][2] + dNi_dx3 * J_1[2][2];
          //                    qx[i] = dNi_dx1;
          //                    qy[i] = dNi_dx2;
          //                    qz[i] = dNi_dx3;

          strainDisplacementMatrix[3*i+0][0] = qx[i];
          strainDisplacementMatrix[3*i+1][0] = (Real)0;
          strainDisplacementMatrix[3*i+2][0] = (Real)0;
          strainDisplacementMatrix[3*i+0][1] = (Real)0;
          strainDisplacementMatrix[3*i+1][1] = qy[i];
          strainDisplacementMatrix[3*i+2][1] = (Real)0;
          strainDisplacementMatrix[3*i+0][2] = (Real)0;
          strainDisplacementMatrix[3*i+1][2] = (Real)0;
          strainDisplacementMatrix[3*i+2][2] = qz[i];
          strainDisplacementMatrix[3*i+0][3] = qy[i];
          strainDisplacementMatrix[3*i+1][3] = qx[i];
          strainDisplacementMatrix[3*i+2][3] = (Real)0;
          strainDisplacementMatrix[3*i+0][4] = (Real)0;
          strainDisplacementMatrix[3*i+1][4] = qz[i];
          strainDisplacementMatrix[3*i+2][4] = qy[i];
          strainDisplacementMatrix[3*i+0][5] = qz[i];
          strainDisplacementMatrix[3*i+1][5] = (Real)0;
          strainDisplacementMatrix[3*i+2][5] = qx[i];
        }
        this->m_strainDisplMatrices[w][integrationPointId] = strainDisplacementMatrix;
        for (int i = 0; i < 8; ++i)
        {
          Mat< 6, 3, Real > MBi;
          MBi[0][0] = U * qx[i];
          MBi[0][1] = V * qy[i];
          MBi[0][2] = V * qz[i];
          MBi[1][0] = V * qx[i];
          MBi[1][1] = U * qy[i];
          MBi[1][2] = V * qz[i];
          MBi[2][0] = V * qx[i];
          MBi[2][1] = V * qy[i];
          MBi[2][2] = U * qz[i];
          MBi[3][0] = W * qy[i];
          MBi[3][1] = W * qx[i];
          MBi[3][2] = (Real)0;
          MBi[4][0] = (Real)0;
          MBi[4][1] = W * qz[i];
          MBi[4][2] = W * qy[i];
          MBi[5][0] = W * qz[i];
          MBi[5][1] = (Real)0;
          MBi[5][2] = W * qx[i];
          for (int j = i; j < 8; ++j)
          {
            Mat< 3, 3, Real > k; // k = BjtMBi
            k[0][0] = qx[j] * MBi[0][0] + qy[j] * MBi[3][0] + qz[j] * MBi[5][0];
            k[0][1] = qx[j] * MBi[0][1] + qy[j] * MBi[3][1] /*+ qz[j]*MBi[5][1]*/;
            k[0][2] = qx[j] * MBi[0][2] /*+ qy[j]*MBi[3][2]*/ + qz[j] * MBi[5][2];

            k[1][0] = qy[j] * MBi[1][0] + qx[j] * MBi[3][0] /*+ qz[j]*MBi[4][0]*/;
            k[1][1] = qy[j] * MBi[1][1] + qx[j] * MBi[3][1] + qz[j] * MBi[4][1];
            k[1][2] = qy[j] * MBi[1][2] /*+ qx[j]*MBi[3][2]*/ + qz[j] * MBi[4][2];

            k[2][0] = qz[j] * MBi[2][0] /*+ qy[j]*MBi[4][0]*/ + qx[j] * MBi[5][0];
            k[2][1] = qz[j] * MBi[2][1] + qy[j] * MBi[4][1] /*+ qx[j]*MBi[5][1]*/;
            k[2][2] = qz[j] * MBi[2][2] + qy[j] * MBi[4][2] + qx[j] * MBi[5][2];

            k = J_1t * k * J_1;

            k *= detJ;
            for (int m = 0; m < 3; ++m)
            {
              for (int l = 0; l < 3; ++l)
              {
                K[i * 3 + m][j * 3 + l] += k[l][m];
              }
            }
          }
        }
        ++integrationPointId;
      }
    }
  }
  for (int i = 0; i < 24; ++i)
  {
    for (int j = i + 1; j < 24; ++j)
    {
      K[j][i] = K[i][j];
    }
  }

  // Mat33 J_1; // only accurate for orthogonal regular hexa
  J_1.fill(0.0);
  Coord l = g - a;
  J_1[0][0] = 2.0f / l[0];
  J_1[1][1] = 2.0f / l[1];
  J_1[2][2] = 2.0f / l[2];

  Real vol = ((b - a).norm() * (d - a).norm() * (e - a).norm());
  vol /= 8.0; // ???

  K.clear();

  for (int i = 0; i < 8; ++i)
  {
    Mat< 3, 3, Real > k = integrateStiffness(_coef[i][0], _coef[i][1], _coef[i][2], _coef[i][0], _coef[i][1], _coef[i][2],
                                             M[0][0], M[0][1], M[3][3], J_1) *
                          vol;

    for (int m = 0; m < 3; ++m)
    {
      for (int l = 0; l < 3; ++l)
      {
        K[i * 3 + m][i * 3 + l] += k[m][l];
      }
    }

    for (int j = i + 1; j < 8; ++j)
    {
      Mat< 3, 3, Real > k = integrateStiffness(_coef[i][0], _coef[i][1], _coef[i][2], _coef[j][0], _coef[j][1], _coef[j][2],
                                               M[0][0], M[0][1], M[3][3], J_1) *
                            vol;

      for (int m = 0; m < 3; ++m)
      {
        for (int l = 0; l < 3; ++l)
        {
          K[i * 3 + m][j * 3 + l] += k[m][l];
        }
      }
    }
  }

  for (int i = 0; i < 24; ++i)
  {
    for (int j = i + 1; j < 24; ++j)
    {
      K[j][i] = K[i][j];
    }
  }

  this->m_stiffnessElementMatrix[w] = K * (Real)stiffnessFactor;
}

template < class DataTypes>
defaulttype::Mat< 3, 3, SReal > HexaFEMForceField< DataTypes >::integrateStiffness(int signx0, int signy0, int signz0, int signx1, int signy1,
                                                                      int signz1, const SReal u, const SReal v, const SReal w,
                                                                      const defaulttype::Mat< 3, 3, SReal >& J_1)
{
  using defaulttype::Mat;
  Mat< 3, 3, Real > K;

  Real t1 = J_1[0][0] * J_1[0][0];    // m^-2            (J0J0             )
  Real t2 = t1 * signx0;              // m^-2            (J0J0    sx0      )
  Real t3 = (Real)(signy0 * signz0);  //                 (           sy0sz0)
  Real t4 = t2 * t3;                  // m^-2            (J0J0    sx0sy0sz0)
  Real t5 = w * signx1;               // kg m^-4 s^-2    (W       sx1      )
  Real t6 = (Real)(signy1 * signz1);  //                 (           sy1sz1)
  Real t7 = t5 * t6;                  // kg m^-4 s^-2    (W       sx1sy1sz1)
  Real t10 = t1 * signy0;             // m^-2            (J0J0       sy0   )
  Real t12 = w * signy1;              // kg m^-4 s^-2    (W          sy1   )
  Real t13 = t12 * signz1;            // kg m^-4 s^-2    (W          sy1sz1)
  Real t16 = t2 * signz0;             // m^-2            (J0J0    sx0   sz0)
  Real t17 = u * signx1;              // kg m^-4 s^-2    (U       sx1      )
  Real t18 = t17 * signz1;            // kg m^-4 s^-2    (U       sx1   sz1)
  Real t21 = t17 * t6;                // kg m^-4 s^-2    (U       sx1sy1sz1)
  Real t24 = t2 * signy0;             // m^-2            (J0J0    sx0sy0   )
  Real t25 = t17 * signy1;            // kg m^-4 s^-2    (U       sx1sy1   )
  Real t28 = t5 * signy1;             // kg m^-4 s^-2    (W       sx1sy1   )
  Real t32 = w * signz1;              // kg m^-4 s^-2    (W             sz1)
  Real t37 = t5 * signz1;             // kg m^-4 s^-2    (W       sx1   sz1)
  Real t43 = J_1[0][0] * signx0;      // m^-1            (J0      sx0      )
  Real t45 = v * J_1[1][1];           // kg m^-5 s^-2    (VJ1              )
  Real t49 = J_1[0][0] * signy0;      // m^-1            (J0         sy0   )
  Real t50 = t49 * signz0;            // m^-1            (J0         sy0sz0)
  Real t51 = w * J_1[1][1];           // kg m^-5 s^-2    (WJ1              )
  Real t52 = (Real)(signx1 * signz1); //                 (        sx1   sz1)
  Real t53 = t51 * t52;               // kg m^-5 s^-2    (WJ1     sx1   sz1)
  Real t56 = t45 * signy1;            // kg m^-5 s^-2    (VJ1        sy1   )
  Real t64 = v * J_1[2][2];           // kg m^-5 s^-2    (VJ2              )
  Real t68 = w * J_1[2][2];           // kg m^-5 s^-2    (WJ2              )
  Real t69 = (Real)(signx1 * signy1); //                 (        sx1sy1   )
  Real t70 = t68 * t69;               // kg m^-5 s^-2    (WJ2     sx1sy1   )
  Real t73 = t64 * signz1;            // kg m^-5 s^-2    (VJ2           sz1)
  Real t81 = J_1[1][1] * signy0;      // m^-1            (J1         sy0   )
  Real t83 = v * J_1[0][0];           // kg m^-5 s^-2    (VJ0              )
  Real t87 = J_1[1][1] * signx0;      // m^-1            (J1      sx0      )
  Real t88 = t87 * signz0;            // m^-1            (J1      sx0   sz0)
  Real t89 = w * J_1[0][0];           // kg m^-5 s^-2    (WJ0              )
  Real t90 = t89 * t6;                // kg m^-5 s^-2    (WJ0        sy1sz1)
  Real t93 = t83 * signx1;            // kg m^-5 s^-2    (VJ0     sx1      )
  Real t100 = J_1[1][1] * J_1[1][1];  // m^-2            (J1J1             )
  Real t101 = t100 * signx0;          // m^-2            (J1J1    sx0      )
  Real t102 = t101 * t3;              // m^-2            (J1J1    sx0sy0sz0)
  Real t110 = t100 * signy0;          // m^-2            (J1J1       sy0   )
  Real t111 = t110 * signz0;          // m^-2            (J1J1       sy0sz0)
  Real t112 = u * signy1;             // kg m^-4 s^-2    (U          sy1   )
  Real t113 = t112 * signz1;          // kg m^-4 s^-2    (U          sy1sz1)
  Real t116 = t101 * signy0;          // m^-2            (J1J1    sx0sy0   )
  Real t144 = J_1[2][2] * signy0;     // m^-1            (J2         sy0   )
  Real t149 = J_1[2][2] * signx0;     // m^-1            (J2      sx0      )
  Real t150 = t149 * signy0;          // m^-1            (J2      sx0sy0   )
  Real t153 = J_1[2][2] * signz0;     // m^-1            (J2            sz0)
  Real t172 = J_1[2][2] * J_1[2][2];  // m^-2            (J2J2             )
  Real t173 = t172 * signx0;          // m^-2            (J2J2    sx0      )
  Real t174 = t173 * t3;              // m^-2            (J2J2    sx0sy0sz0)
  Real t177 = t173 * signz0;          // m^-2            (J2J2    sx0   sz0)
  Real t180 = t172 * signy0;          // m^-2            (J2J2       sy0   )
  Real t181 = t180 * signz0;          // m^-2            (J2J2       sy0sz0)
  // kg m^-6 s^-2
  K[0][0] = (float)(t4 * t7 / 36.0 + t10 * signz0 * t13 / 12.0 + t16 * t18 / 24.0 + t4 * t21 / 72.0 + t24 * t25 / 24.0 +
                    t24 * t28 / 24.0 + t1 * signz0 * t32 / 8.0 + t10 * t12 / 8.0 + t16 * t37 / 24.0 + t2 * t17 / 8.0);

  K[0][1] = (float)(t43 * signz0 * t45 * t6 / 24.0 + t50 * t53 / 24.0 + t43 * t56 / 8.0 + t49 * t51 * signx1 / 8.0);

  K[0][2] =
      (float)(t43 * signy0 * t64 * t6 / 24.0 + t50 * t70 / 24.0 + t43 * t73 / 8.0 + J_1[0][0] * signz0 * t68 * signx1 / 8.0);

  K[1][0] = (float)(t81 * signz0 * t83 * t52 / 24.0 + t88 * t90 / 24.0 + t81 * t93 / 8.0 + t87 * t89 * signy1 / 8.0);

  K[1][1] = (float)(t102 * t7 / 36.0 + t102 * t21 / 72.0 + t101 * signz0 * t37 / 12.0 + t111 * t113 / 24.0 + t116 * t28 / 24.0 +
                    t100 * signz0 * t32 / 8.0 + t111 * t13 / 24.0 + t116 * t25 / 24.0 + t110 * t112 / 8.0 + t101 * t5 / 8.0);

  K[1][2] =
      (float)(t87 * signy0 * t64 * t52 / 24.0 + t88 * t70 / 24.0 + t81 * t73 / 8.0 + J_1[1][1] * signz0 * t68 * signy1 / 8.0);

  K[2][0] = (float)(t144 * signz0 * t83 * t69 / 24.0 + t150 * t90 / 24.0 + t153 * t93 / 8.0 + t149 * t89 * signz1 / 8.0);

  K[2][1] = (float)(t149 * signz0 * t45 * t69 / 24.0 + t150 * t53 / 24.0 + t153 * t56 / 8.0 + t144 * t51 * signz1 / 8.0);

  K[2][2] =
      (float)(t174 * t7 / 36.0 + t177 * t37 / 24.0 + t181 * t13 / 24.0 + t174 * t21 / 72.0 + t173 * signy0 * t28 / 12.0 +
              t180 * t12 / 8.0 + t181 * t113 / 24.0 + t177 * t18 / 24.0 + t172 * signz0 * u * signz1 / 8.0 + t173 * t5 / 8.0);

  return K;
}

template < class DataTypes>
void HexaFEMForceField< DataTypes >::computeMaterialStiffnessMatrix(Volume w, const Index& a, const Index& b, const Index& c,
                                                                    const Index& d, const Index& e, const Index& f,
                                                                    const Index& g, const Index& h)
{
  unsigned i = this->getTopology()->embedding(w);
  const VecReal& localStiffnessFactor = this->f_localStiffnessFactor.getValue();

  computeMaterialStiffnessMatrix(
      this->m_materielStiffness[w], a, b, c, d, e, f, g, h,
      this->f_youngModulus.getValue()[0] *
          (localStiffnessFactor.empty()
               ? 1.0f
               : localStiffnessFactor[i * localStiffnessFactor.size() / this->getTopology()->template nb_cells<Volume::ORBIT>()]));

  if (d_computeVonMisesStress.getValue() >0) {
      m_elemLambda[w] = this->m_materielStiffness[w][0][1];
      m_elemMu[w] = this->m_materielStiffness[w][3][3];
  }
}

template < class DataTypes>
void HexaFEMForceField< DataTypes >::computeMaterialStiffnessMatrix(MaterialStiffnessMatrix& E, const Index& a, const Index& b,
                                                                    const Index& c, const Index& d, const Index& e,
                                                                    const Index& f, const Index& g, const Index& h,
                                                                    double localStiffnessFactor)
{
  Inherit1::computeMaterialStiffnessMatrix(E, localStiffnessFactor);
}

//////////////////////////////////////////////////////////////////////
////////////////////  large displacements method  ////////////////////
//////////////////////////////////////////////////////////////////////

template < class DataTypes>
void HexaFEMForceField< DataTypes >::initLarge(Volume w, const Index& a, const Index& b, const Index& c, const Index& d,
                                               const Index& e, const Index& f, const Index& g, const Index& h)
{
  helper::ReadAccessor< Data< VecCoord > > X0 = this->mstate->read(core::VecCoordId::restPosition());
  helper::fixed_array< Coord, 8 >& rie = this->m_rotatedInitialElements[w];

  TransformationMatrix R_0_1;

  Coord horizontal = (X0[b] - X0[a] + X0[c] - X0[d] + X0[f] - X0[e] + X0[g] - X0[h]) * .25;
  Coord vertical = (X0[d] - X0[a] + X0[c] - X0[b] + X0[h] - X0[e] + X0[g] - X0[f]) * .25;
  computeRotationLarge(R_0_1, horizontal, vertical);

  rie[0] = R_0_1 * X0[a];
  rie[1] = R_0_1 * X0[b];
  rie[2] = R_0_1 * X0[c];
  rie[3] = R_0_1 * X0[d];

  rie[4] = R_0_1 * X0[e];
  rie[5] = R_0_1 * X0[f];
  rie[6] = R_0_1 * X0[g];
  rie[7] = R_0_1 * X0[h];

  this->m_initialTransformation[w] = this->m_rotation[w] = R_0_1;

  this->computeElementStiffnessMatrix(w, this->m_materielStiffness[w], rie[0], rie[1], rie[2], rie[3], rie[4], rie[5], rie[6],
                                      rie[7]);
}

template < class DataTypes>
inline void HexaFEMForceField< DataTypes >::computeRotationLarge(TransformationMatrix& r, Coord& edgex, Coord& edgey)
{
  // frame:
  // first axis is first edge of tetra
  // second axis orthogonal to the first on the plane of the two first edges
  // third axis orthogonal to first and second
  edgex.normalize();
  edgey.normalize();

  Coord edgez = cross(edgex, edgey);
  edgez.normalize();

  edgey = cross(edgez, edgex);
  edgey.normalize();

  r[0][0] = edgex[0];
  r[0][1] = edgex[1];
  r[0][2] = edgex[2];
  r[1][0] = edgey[0];
  r[1][1] = edgey[1];
  r[1][2] = edgey[2];
  r[2][0] = edgez[0];
  r[2][1] = edgez[1];
  r[2][2] = edgez[2];
}

template < class DataTypes>
void HexaFEMForceField< DataTypes >::accumulateForceLarge(Vector& force, const Vector& p, Volume w, const Index& a,
                                                          const Index& b, const Index& c, const Index& d, const Index& e,
                                                          const Index& f, const Index& g, const Index& h)
{
    TransformationMatrix& rot = this->m_rotation[w];

  // Rotation matrix (deformed and displaced Tetrahedron/world)
  Coord horizontal;
  horizontal = (p[b] - p[a] + p[c] - p[d] + p[f] - p[e] + p[g] - p[h]) * .25;
  Coord vertical;
  vertical = (p[d] - p[a] + p[c] - p[b] + p[h] - p[e] + p[g] - p[f]) * .25;
  computeRotationLarge(rot, horizontal, vertical);
  // positions of the deformed and displaced Tetrahedron in its frame
  helper::fixed_array< Coord, 8 > deformed;
  deformed[0] = rot * p[a];
  deformed[1] = rot * p[b];
  deformed[2] = rot * p[c];
  deformed[3] = rot * p[d];
  deformed[4] = rot * p[e];
  deformed[5] = rot * p[f];
  deformed[6] = rot * p[g];
  deformed[7] = rot * p[h];

  // displacement
  Displacement U;
  for (int k = 0; k < 8; ++k)
  {
    int indice = k * 3;
    for (int l = 0; l < 3; ++l)
      U[indice + l] = this->m_rotatedInitialElements[w][k][l] - deformed[k][l];
  }

  Displacement F;
  if (this->f_updateStiffnessMatrix.getValue())
  {
    computeElementStiffnessMatrix(w, this->m_materielStiffness[w], deformed[0], deformed[1], deformed[2], deformed[3],
                                  deformed[4], deformed[5], deformed[6], deformed[7]);
  }

  // compute element force
  F = this->m_stiffnessElementMatrix[w] * U;

  force[a] += rot.multTranspose(Deriv(F[0], F[1], F[2]));
  force[b] += rot.multTranspose(Deriv(F[3], F[4], F[5]));
  force[c] += rot.multTranspose(Deriv(F[6], F[7], F[8]));
  force[d] += rot.multTranspose(Deriv(F[9], F[10], F[11]));
  force[e] += rot.multTranspose(Deriv(F[12], F[13], F[14]));
  force[f] += rot.multTranspose(Deriv(F[15], F[16], F[17]));
  force[g] += rot.multTranspose(Deriv(F[18], F[19], F[20]));
  force[h] += rot.multTranspose(Deriv(F[21], F[22], F[23]));
}

template < class DataTypes>
void HexaFEMForceField< DataTypes >::applyStiffnessCorotational(Vector& df, const Vector& dx, Volume w, const Index& a,
                                                                const Index& b, const Index& c, const Index& d, const Index& e,
                                                                const Index& f, const Index& g, const Index& h, double kFact)
{

  Displacement U;
  Coord x_2;
  TransformationMatrix& rot = this->m_rotation[w];
  x_2 = rot * dx[a];
  U[0] = x_2[0];
  U[1] = x_2[1];
  U[2] = x_2[2];

  x_2 = rot * dx[b];
  U[3] = x_2[0];
  U[4] = x_2[1];
  U[5] = x_2[2];

  x_2 = rot * dx[c];
  U[6] = x_2[0];
  U[7] = x_2[1];
  U[8] = x_2[2];

  x_2 = rot * dx[d];
  U[9] = x_2[0];
  U[10] = x_2[1];
  U[11] = x_2[2];

  x_2 = rot * dx[e];
  U[12] = x_2[0];
  U[13] = x_2[1];
  U[14] = x_2[2];

  x_2 = rot * dx[f];
  U[15] = x_2[0];
  U[16] = x_2[1];
  U[17] = x_2[2];

  x_2 = rot * dx[g];
  U[18] = x_2[0];
  U[19] = x_2[1];
  U[20] = x_2[2];

  x_2 = rot * dx[h];
  U[21] = x_2[0];
  U[22] = x_2[1];
  U[23] = x_2[2];

  Displacement F;
  F = this->m_stiffnessElementMatrix[w] * U;

  df[a] -= rot.multTranspose(Deriv(F[0], F[1], F[2])) * kFact;
  df[b] -= rot.multTranspose(Deriv(F[3], F[4], F[5])) * kFact;
  df[c] -= rot.multTranspose(Deriv(F[6], F[7], F[8])) * kFact;
  df[d] -= rot.multTranspose(Deriv(F[9], F[10], F[11])) * kFact;
  df[e] -= rot.multTranspose(Deriv(F[12], F[13], F[14])) * kFact;
  df[f] -= rot.multTranspose(Deriv(F[15], F[16], F[17])) * kFact;
  df[g] -= rot.multTranspose(Deriv(F[18], F[19], F[20])) * kFact;
  df[h] -= rot.multTranspose(Deriv(F[21], F[22], F[23])) * kFact;
}

//////////////////////////////////////////////////////////////////////
////////////////////  polar decomposition method  ////////////////////
//////////////////////////////////////////////////////////////////////

template < class DataTypes>
void HexaFEMForceField< DataTypes >::initPolar(Volume w, const Index& a, const Index& b, const Index& c, const Index& d,
                                               const Index& e, const Index& f, const Index& g, const Index& h)
{
  helper::ReadAccessor< Data< VecCoord > > X0 = this->mstate->read(core::VecCoordId::restPosition());
  TransformationMatrix R_0_1;

  Coord horizontal;
  horizontal = (X0[b] - X0[a] + X0[c] - X0[d] + X0[f] - X0[e] + X0[g] - X0[h]) * .25;
  Coord vertical;
  vertical = (X0[d] - X0[a] + X0[c] - X0[b] + X0[h] - X0[e] + X0[g] - X0[f]) * .25;
  Coord normal2handv;
  normal2handv = (X0[e] - X0[a] + X0[f] - X0[b] + X0[h] - X0[d] + X0[g] - X0[c]) * .25;
  computeRotationPolar(R_0_1, horizontal, vertical, normal2handv);

  this->m_initialTransformation[w] = this->m_rotation[w] = R_0_1;

  helper::fixed_array< Coord, 8 >& rie = this->m_rotatedInitialElements[w];
  rie[0] = R_0_1 * X0[a];
  rie[1] = R_0_1 * X0[b];
  rie[2] = R_0_1 * X0[c];
  rie[3] = R_0_1 * X0[d];

  rie[4] = R_0_1 * X0[e];
  rie[5] = R_0_1 * X0[f];
  rie[6] = R_0_1 * X0[g];
  rie[7] = R_0_1 * X0[h];

  this->computeElementStiffnessMatrix(w, this->m_materielStiffness[w], rie[0], rie[1], rie[2], rie[3], rie[4], rie[5], rie[6],
                                      rie[7]);
}

template < class DataTypes>
inline void HexaFEMForceField< DataTypes >::computeRotationPolar(TransformationMatrix& r, Coord& edgex, Coord& edgey,
                                                                 Coord& edgez)
{
  TransformationMatrix A;

  // Coord Edge =(Xs[1]-Xs[0] + Xs[2]-Xs[3] + Xs[5]-Xs[4] + Xs[6]-Xs[7])*.25;
  A[0][0] = edgex[0];
  A[0][1] = edgex[1];
  A[0][2] = edgex[2];

  // Edge = (Xs[3]-Xs[0] + Xs[2]-Xs[1] + Xs[7]-Xs[4] + Xs[6]-Xs[5])*.25;
  A[1][0] = edgey[0];
  A[1][1] = edgey[1];
  A[1][2] = edgey[2];

  // Edge = (Xs[4]-Xs[0] + Xs[5]-Xs[1] + Xs[7]-Xs[3] + Xs[6]-Xs[2])*.25;
  A[2][0] = edgez[0];
  A[2][1] = edgez[1];
  A[2][2] = edgez[2];

  helper::Decompose< Real >::polarDecomposition(A, r);
}

template < class DataTypes>
void HexaFEMForceField< DataTypes >::accumulateForcePolar(Vector& force, const Vector& p, Volume w, const Index& a,
                                                          const Index& b, const Index& c, const Index& d, const Index& e,
                                                          const Index& f, const Index& g, const Index& h)
{
  TransformationMatrix& rot = this->m_rotation[w];
  Coord horizontal;
  horizontal = (p[b] - p[a] + p[c] - p[d] + p[f] - p[e] + p[g] - p[h]) * .25;
  Coord vertical;
  vertical = (p[d] - p[a] + p[c] - p[b] + p[h] - p[e] + p[g] - p[f]) * .25;
  Coord normal2handv;
  normal2handv = (p[e] - p[a] + p[f] - p[b] + p[h] - p[d] + p[g] - p[c]) * .25;
  computeRotationPolar(rot, horizontal, vertical, normal2handv);

  // positions of the deformed and displaced tetrahedron in its frame
  helper::fixed_array< Coord, 8 > deformed;
  deformed[0] = rot * p[a];
  deformed[1] = rot * p[b];
  deformed[2] = rot * p[c];
  deformed[3] = rot * p[d];
  deformed[4] = rot * p[e];
  deformed[5] = rot * p[f];
  deformed[6] = rot * p[g];
  deformed[7] = rot * p[h];

  // displacement
  Displacement U;
  for (int k = 0; k < 8; ++k)
  {
    int indice = k * 3;
    for (int l = 0; l < 3; ++l)
      U[indice + l] = this->m_rotatedInitialElements[w][k][l] - deformed[k][l];
  }

  Displacement F;
  if (this->f_updateStiffnessMatrix.getValue())
  {
    computeElementStiffnessMatrix(w, this->m_materielStiffness[w], deformed[0], deformed[1], deformed[2], deformed[3],
                                  deformed[4], deformed[5], deformed[6], deformed[7]);
  }

  if (!this->f_assembling.getValue())
  {
    // compute element force
    F = this->m_stiffnessElementMatrix[w] * U;

    force[a] += rot.multTranspose(Deriv(F[0], F[1], F[2]));
    force[b] += rot.multTranspose(Deriv(F[3], F[4], F[5]));
    force[c] += rot.multTranspose(Deriv(F[6], F[7], F[8]));
    force[d] += rot.multTranspose(Deriv(F[9], F[10], F[11]));
    force[e] += rot.multTranspose(Deriv(F[12], F[13], F[14]));
    force[f] += rot.multTranspose(Deriv(F[15], F[16], F[17]));
    force[g] += rot.multTranspose(Deriv(F[18], F[19], F[20]));
    force[h] += rot.multTranspose(Deriv(F[21], F[22], F[23]));
  }
  else
  {
    serr << "TODO(HexaFEMForceField): support for assembling system matrix when using polar method." << sendl;
  }
}

template < class DataTypes> void HexaFEMForceField< DataTypes >::draw(const core::visual::VisualParams* vparams)
{
  using defaulttype::Vec;

  if (!this->getTopology())
    return;
  if (!vparams->displayFlags().getShowForceFields() && l_vonMisesStressColorMap)
  {
      l_vonMisesStressColorMap->f_showLegend.setValue(false);
      return;
  }
  if( d_computeVonMisesStress.getValue() && l_vonMisesStressColorMap)
      l_vonMisesStressColorMap->f_showLegend.setValue(true);

  if (!this->mstate)
    return;
  if (!this->f_drawing.getValue())
    return;

  helper::ReadAccessor< Data< VecCoord > > x = this->mstate->read(core::VecCoordId::position());

  if (vparams->displayFlags().getShowWireFrame())
  {
    vparams->drawTool()->setPolygonMode(0, true);
  }

    /// vonMises stress
    Real minVM = (Real)1e20, maxVM = (Real)-1e20;
    if (d_computeVonMisesStress.getValue() > 0 ) {

		this->m_topology->foreach_cell([&,this](Volume w)
		{
			minVM = (m_sxx[w] < minVM) ? m_sxx[w] : minVM;
			maxVM = (m_sxx[w] > maxVM) ? m_sxx[w] : maxVM;
		});
		if (maxVM < prevMaxStress)
			maxVM = prevMaxStress;
	}

//  std::vector< defaulttype::Vector3 > points[6];
	std::vector< defaulttype::Vector3 > hexaPoints[6];

	this->m_topology->foreach_cell([&,this](Volume w)
	{
    assert(this->getTopology()->is_hexa(w));
    {
      hexaPoints[0].clear();
      hexaPoints[1].clear();
      hexaPoints[2].clear();
      hexaPoints[3].clear();
      hexaPoints[4].clear();
      hexaPoints[5].clear();

      const auto& t = this->getTopology()->get_dofs(w);
      Coord center =
          (x[t[HexaFF_Order::ind[0]]] + x[t[HexaFF_Order::ind[1]]] + x[t[HexaFF_Order::ind[2]]] + x[t[HexaFF_Order::ind[3]]] +
           x[t[HexaFF_Order::ind[4]]] + x[t[HexaFF_Order::ind[5]]] + x[t[HexaFF_Order::ind[6]]] + x[t[HexaFF_Order::ind[7]]]) *
          0.125;
      Coord pa = x[t[HexaFF_Order::ind[0]]] - (x[t[HexaFF_Order::ind[0]]] - center) * 0.15;
      Coord pb = x[t[HexaFF_Order::ind[1]]] - (x[t[HexaFF_Order::ind[1]]] - center) * 0.15;
      Coord pc = x[t[HexaFF_Order::ind[3]]] - (x[t[HexaFF_Order::ind[3]]] - center) * 0.15;
      Coord pd = x[t[HexaFF_Order::ind[2]]] - (x[t[HexaFF_Order::ind[2]]] - center) * 0.15;
      Coord pe = x[t[HexaFF_Order::ind[4]]] - (x[t[HexaFF_Order::ind[4]]] - center) * 0.15;
      Coord pf = x[t[HexaFF_Order::ind[5]]] - (x[t[HexaFF_Order::ind[5]]] - center) * 0.15;
      Coord pg = x[t[HexaFF_Order::ind[7]]] - (x[t[HexaFF_Order::ind[7]]] - center) * 0.15;
      Coord ph = x[t[HexaFF_Order::ind[6]]] - (x[t[HexaFF_Order::ind[6]]] - center) * 0.15;

//      points[0].push_back(pa);
//      points[0].push_back(pb);
//      points[0].push_back(pf);
//      points[0].push_back(pa);
//      points[0].push_back(pf);
//      points[0].push_back(pe);
      hexaPoints[0].push_back(pa);
      hexaPoints[0].push_back(pb);
      hexaPoints[0].push_back(pf);
      hexaPoints[0].push_back(pe);

//      points[1].push_back(pa);
//      points[1].push_back(pc);
//      points[1].push_back(pg);
//      points[1].push_back(pa);
//      points[1].push_back(pg);
//      points[1].push_back(pe);
      hexaPoints[1].push_back(pa);
      hexaPoints[1].push_back(pc);
      hexaPoints[1].push_back(pg);
      hexaPoints[1].push_back(pe);

//      points[2].push_back(pc);
//      points[2].push_back(pd);
//      points[2].push_back(ph);
//      points[2].push_back(pc);
//      points[2].push_back(ph);
//      points[2].push_back(pg);
      hexaPoints[2].push_back(pc);
      hexaPoints[2].push_back(pd);
      hexaPoints[2].push_back(ph);
      hexaPoints[2].push_back(pg);

//      points[3].push_back(pd);
//      points[3].push_back(pb);
//      points[3].push_back(pf);
//      points[3].push_back(pd);
//      points[3].push_back(pf);
//      points[3].push_back(ph);
      hexaPoints[3].push_back(pd);
      hexaPoints[3].push_back(pb);
      hexaPoints[3].push_back(pf);
      hexaPoints[3].push_back(ph);


//      points[4].push_back(pa);
//      points[4].push_back(pb);
//      points[4].push_back(pd);
//      points[4].push_back(pa);
//      points[4].push_back(pd);
//      points[4].push_back(pc);
      hexaPoints[4].push_back(pa);
      hexaPoints[4].push_back(pb);
      hexaPoints[4].push_back(pd);
      hexaPoints[4].push_back(pc);

//      points[5].push_back(pe);
//      points[5].push_back(pg);
//      points[5].push_back(ph);
//      points[5].push_back(pe);
//      points[5].push_back(ph);
//      points[5].push_back(pf);
      hexaPoints[5].push_back(pe);
      hexaPoints[5].push_back(pg);
      hexaPoints[5].push_back(ph);
      hexaPoints[5].push_back(pf);
    }
//  }
//  vparams->drawTool()->drawQuads(hexaPoints[0], Vec< 4, float >(0.1, 0.0, 1.0, 1.0));
//  vparams->drawTool()->drawQuads(hexaPoints[1], Vec< 4, float >(0.1, 0.1, 1.0, 1.0));
//  vparams->drawTool()->drawQuads(hexaPoints[2], Vec< 4, float >(0.1, 0.2, 1.0, 1.0));
//  vparams->drawTool()->drawQuads(hexaPoints[3], Vec< 4, float >(0.1, 0.3, 1.0, 1.0));
//  vparams->drawTool()->drawQuads(hexaPoints[4], Vec< 4, float >(0.1, 0.4, 1.0, 1.0));
//  vparams->drawTool()->drawQuads(hexaPoints[5], Vec< 4, float >(0.1, 0.5, 1.0, 1.0));

    defaulttype::Vec4f col;
    if (d_computeVonMisesStress.getValue() > 0 ) {
        sofa::helper::ColorMap::evaluator<Real> evalColor = l_vonMisesStressColorMap->getEvaluator(minVM, maxVM);
        col = evalColor(m_sxx[w]);
        col[3] = d_vonMisesStressAlpha.getValue();//1.0f;
        vparams->drawTool()->drawQuads(hexaPoints[0], col); col[1]+=0.025;
        vparams->drawTool()->drawQuads(hexaPoints[1], col); col[1]+=0.025;
        vparams->drawTool()->drawQuads(hexaPoints[2], col); col[1]+=0.025;
        vparams->drawTool()->drawQuads(hexaPoints[3], col); col[1]+=0.025;
        vparams->drawTool()->drawQuads(hexaPoints[4], col); col[1]+=0.025;
        vparams->drawTool()->drawQuads(hexaPoints[5], col);
    }
    else
    {
      vparams->drawTool()->drawQuads(hexaPoints[0], Vec< 4, float >(0.1, 0.0, 1.0, 1.0));
      vparams->drawTool()->drawQuads(hexaPoints[1], Vec< 4, float >(0.1, 0.1, 1.0, 1.0));
      vparams->drawTool()->drawQuads(hexaPoints[2], Vec< 4, float >(0.1, 0.2, 1.0, 1.0));
      vparams->drawTool()->drawQuads(hexaPoints[3], Vec< 4, float >(0.1, 0.3, 1.0, 1.0));
      vparams->drawTool()->drawQuads(hexaPoints[4], Vec< 4, float >(0.1, 0.4, 1.0, 1.0));
      vparams->drawTool()->drawQuads(hexaPoints[5], Vec< 4, float >(0.1, 0.5, 1.0, 1.0));
    }
    });
//  std::cout << "minVM : " << minVM << std::endl;
//  std::cout << "maxVM : " << maxVM << std::endl;
//  std::cout << "m_vonMises[0] : " << m_vonMisesPerElement[0] << std::endl;

//  vparams->drawTool()->drawTriangles(points[0], Vec< 4, float >(0.1, 0.0, 1.0, 1.0));
//  vparams->drawTool()->drawTriangles(points[1], Vec< 4, float >(0.1, 0.1, 1.0, 1.0));
//  vparams->drawTool()->drawTriangles(points[2], Vec< 4, float >(0.1, 0.2, 1.0, 1.0));
//  vparams->drawTool()->drawTriangles(points[3], Vec< 4, float >(0.1, 0.3, 1.0, 1.0));
//  vparams->drawTool()->drawTriangles(points[4], Vec< 4, float >(0.1, 0.4, 1.0, 1.0));
//  vparams->drawTool()->drawTriangles(points[5], Vec< 4, float >(0.1, 0.5, 1.0, 1.0));

  if (vparams->displayFlags().getShowWireFrame())
  {
    vparams->drawTool()->setPolygonMode(0, false);
  }
}

} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_HEXA_FEM_FF_INL
