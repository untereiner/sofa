/*******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 1       *
*                (c) 2006-2007 MGH, INRIA, USTL, UJF, CNRS                     *
*                                                                              *
* This library is free software; you can redistribute it and/or modify it      *
* under the terms of the GNU Lesser General Public License as published by the *
* Free Software Foundation; either version 2.1 of the License, or (at your     *
* option) any later version.                                                   *
*                                                                              *
* This library is distributed in the hope that it will be useful, but WITHOUT  *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
* for more details.                                                            *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this library; if not, write to the Free Software Foundation,      *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
*                                                                              *
* Contact information: contact@sofa-framework.org                              *
*                                                                              *
* Authors: J. Allard, P-J. Bensoussan, S. Cotin, C. Duriez, H. Delingette,     *
* F. Faure, S. Fonteneau, L. Heigeas, C. Mendoza, M. Nesme, P. Neumann,        *
* and F. Poyer                                                                 *
*******************************************************************************/
#include <sofa/component/visualmodel/VisualModelImpl.h>
#include <sofa/helper/gl/RAII.h>
#include <sofa/helper/vector.h>
#include <sofa/defaulttype/Quat.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/io/Mesh.h>
#include <sofa/helper/rmath.h>
#include <sstream>

namespace sofa
{

namespace component
{

namespace visualmodel
{

using namespace sofa::defaulttype;
using namespace sofa::core::componentmodel::topology;

void VisualModelImpl::parse(core::objectmodel::BaseObjectDescription* arg)
{
    this->core::VisualModel::parse(arg);
    VisualModelImpl* obj = this;

    if (arg->getAttribute("normals")!=NULL)
        obj->setUseNormals(atoi(arg->getAttribute("normals"))!=0);

    if (arg->getAttribute("castshadow")!=NULL)
        obj->setCastShadow(atoi(arg->getAttribute("castshadow"))!=0);

    std::string filename = arg->getAttribute("filename","");
    std::string loader = arg->getAttribute("loader","");
    std::string texturename = arg->getAttribute("texturename","");
    obj->load(filename, loader, texturename);
    if (arg->getAttribute("flip")!=NULL)
    {
        obj->flipFaces();
    }
    if (arg->getAttribute("color"))
    {
        obj->setColor(arg->getAttribute("color"));
    }
    if (arg->getAttribute("scaleTex")!=NULL)
    {
        obj->applyUVScale(atof(arg->getAttribute("scaleTex","1.0")), atof(arg->getAttribute("scaleTex","1.0")));
    }
    if (arg->getAttribute("du")!=NULL || arg->getAttribute("dv")!=NULL)
    {
        obj->applyUVTranslation(atof(arg->getAttribute("du","0.0")), atof(arg->getAttribute("dv","0.0")));
    }
    if (arg->getAttribute("scale")!=NULL)
    {
        obj->applyScale(atof(arg->getAttribute("scale","1.0")));
    }
    if (arg->getAttribute("dx")!=NULL || arg->getAttribute("dy")!=NULL || arg->getAttribute("dz")!=NULL)
    {
        obj->applyTranslation(atof(arg->getAttribute("dx","0.0")),atof(arg->getAttribute("dy","0.0")),atof(arg->getAttribute("dz","0.0")));
    }
    if (arg->getAttribute("rx")!=NULL)
    {
        obj->applyRotation(Quat(Vec3d(1,0,0), atof(arg->getAttribute("rx","0.0"))*R_PI/180));
    }
    if (arg->getAttribute("ry")!=NULL)
    {
        obj->applyRotation(Quat(Vec3d(0,1,0), atof(arg->getAttribute("ry","0.0"))*R_PI/180));
    }
    if (arg->getAttribute("rz")!=NULL)
    {
        obj->applyRotation(Quat(Vec3d(0,0,1), atof(arg->getAttribute("rz","0.0"))*R_PI/180));
    }

    sofa::helper::io::Mesh::Material M=material.getValue();
    float r,g,b,a;
    if (arg->getAttribute("ambient")!=NULL)
    {
        sscanf(arg->getAttribute("ambient"),"%f %f %f %f", &r, &g, &b, &a);
        M.ambient = Vec4f(r,g,b,a);
    }
    if (arg->getAttribute("useambient")!=NULL)
    {
        M.useAmbient = atoi(arg->getAttribute("useambient"))!=0;
    }
    if (arg->getAttribute("diffuse")!=NULL)
    {
        sscanf(arg->getAttribute("diffuse"),"%f %f %f %f", &r, &g, &b, &a);
        M.diffuse = Vec4f(r,g,b,a);
    }
    if (arg->getAttribute("usediffuse")!=NULL)
    {
        M.useDiffuse = atoi(arg->getAttribute("usediffuse"))!=0;
    }
    if (arg->getAttribute("emissive")!=NULL)
    {
        sscanf(arg->getAttribute("emissive"),"%f %f %f %f", &r, &g, &b, &a);
        M.emissive = Vec4f(r,g,b,a);
    }
    if (arg->getAttribute("useemissive")!=NULL)
    {
        M.useEmissive = atoi(arg->getAttribute("useemissive"))!=0;
    }
    if (arg->getAttribute("specular")!=NULL)
    {
        sscanf(arg->getAttribute("specular"),"%f %f %f %f", &r, &g, &b, &a);
        M.specular = Vec4f(r,g,b,a);
    }
    if (arg->getAttribute("usespecular")!=NULL)
    {
        M.useSpecular = atoi(arg->getAttribute("usespecular"))!=0;
    }

    if (arg->getAttribute("shininess")!=NULL)
    {
        sscanf(arg->getAttribute("shininess"),"%f", &a);
        M.shininess = a;
    }
    if (arg->getAttribute("useshininess")!=NULL)
    {
        M.useShininess = atoi(arg->getAttribute("useshininess"))!=0;
    }
    material.setValue(M);
}

SOFA_DECL_CLASS(VisualModelImpl)

int VisualModelImplClass = core::RegisterObject("Generic visual model. If a viewer is active it will replace the VisualModel alias, otherwise nothing will be displayed.")
        .add< VisualModelImpl >()
        .addAlias("VisualModel")
        ;

VisualModelImpl::VisualModelImpl() //const std::string &name, std::string filename, std::string loader, std::string textureName)
    :  useTopology(false), lastMeshRev(-1), useNormals(true), castShadow(true),material(dataField(&material,"material","Material")) //, tex(NULL)
{
    inputVertices = &vertices;
}

VisualModelImpl::~VisualModelImpl()
{
    if (inputVertices != &vertices) delete inputVertices;
}

bool VisualModelImpl::isTransparent()
{
    return (material.getValue().useDiffuse && material.getValue().diffuse[3] < 1.0);
}

void VisualModelImpl::draw()
{
    if (!isTransparent())
        internalDraw();
}

void VisualModelImpl::drawTransparent()
{
    if (isTransparent())
        internalDraw();
}

void VisualModelImpl::drawShadow()
{
    if (!isTransparent() && getCastShadow())
    {
        //std::cout << "drawShadow for "<<getName()<<std::endl;
        internalDraw();
    }
}

bool VisualModelImpl::load(const std::string& filename, const std::string& loader, const std::string& textureName)
{
    bool tex = false;
    if (!textureName.empty())
    {
        tex = loadTexture(textureName);
    }

    if (!filename.empty())
    {
        //name = filename;
        helper::io::Mesh *objLoader;
        if (loader.empty())
            objLoader = helper::io::Mesh::Create(filename);
        else
            objLoader = helper::io::Mesh::Factory::CreateObject(loader, filename);

        if (!objLoader)
        {
            return false;
        }
        else
        {
            const vector< vector< vector<int> > > &facetsImport = objLoader->getFacets();
            const vector<Vector3> &verticesImport = objLoader->getVertices();
            const vector<Vector3> &normalsImport = objLoader->getNormals();
            const vector<Vector3> &texCoordsImport = objLoader->getTexCoords();

            const helper::io::Mesh::Material &materialImport = objLoader->getMaterial();

            if (materialImport.activated)
            {
                helper::io::Mesh::Material M;
                M = materialImport;
                material.setValue(M);
            }

            std::cout << "Vertices Import size : " << verticesImport.size() << " (" << normalsImport.size() << " normals)." << std::endl;

            int nbVIn = verticesImport.size();
            // First we compute for each point how many pair of normal/texcoord indices are used
            // The map store the final index of each combinaison
            vector< std::map< std::pair<int,int>, int > > vertTexNormMap;
            vertTexNormMap.resize(nbVIn);
            for (unsigned int i = 0; i < facetsImport.size(); i++)
            {
                const vector<vector <int> >& vertNormTexIndex = facetsImport[i];
                if (vertNormTexIndex[0].size() < 3) continue; // ignore lines
                const vector<int>& verts = vertNormTexIndex[0];
                const vector<int>& texs = vertNormTexIndex[1];
                const vector<int>& norms = vertNormTexIndex[2];
                for (unsigned int j = 0; j < verts.size(); j++)
                {
                    vertTexNormMap[verts[j]][std::make_pair((tex?texs[j]:-1), (useNormals?norms[j]:0))] = 0;
                }
            }

            // Then we can compute how many vertices are created
            int nbVOut = 0;
            bool vsplit = false;
            for (int i = 0; i < nbVIn; i++)
            {
                int s = vertTexNormMap[i].size();
                nbVOut += s;
                if (s!=1)
                    vsplit = true;
            }

            // Then we can create the final arrays

            vertices.resize(nbVOut);
            vnormals.resize(nbVOut);

            if (tex)
                vtexcoords.resize(nbVOut);

            if (vsplit)
            {
                inputVertices = new ResizableExtVector<Coord>;
                inputVertices->resize(nbVIn);
                vertPosIdx.resize(nbVOut);
                vertNormIdx.resize(nbVOut);
            }
            else
                inputVertices = &vertices;

            int nbNOut = 0; /// Number of different normals
            for (int i = 0, j = 0; i < nbVIn; i++)
            {
                if (vsplit)
                    (*inputVertices)[i] = verticesImport[i];
                std::map<int, int> normMap;
                for (std::map<std::pair<int, int>, int>::iterator it = vertTexNormMap[i].begin();
                        it != vertTexNormMap[i].end(); ++it)
                {
                    vertices[j] = verticesImport[i];
                    int t = it->first.first;
                    int n = it->first.second;
                    if ((unsigned)n < normalsImport.size())
                        vnormals[j] = normalsImport[n];
                    if ((unsigned)t < texCoordsImport.size())
                        vtexcoords[j] = texCoordsImport[t];
                    if (vsplit)
                    {
                        vertPosIdx[j] = i;
                        if (normMap.count(n))
                            vertNormIdx[j] = normMap[n];
                        else
                            vertNormIdx[j] = normMap[n] = nbNOut++;
                    }
                    it->second = j++;
                }
            }
            if (!vsplit) nbNOut = nbVOut;
            else if (nbNOut == nbVOut) vertNormIdx.resize(0);

            std::cout << "Vertices Export size : " << nbVOut << " (" << nbNOut << " normals)." << std::endl;

            std::cout << "Facets Import size : " << facetsImport.size() << std::endl;

            // Then we create the triangles and quads

            for (unsigned int i = 0; i < facetsImport.size(); i++)
            {
                const vector<vector <int> >& vertNormTexIndex = facetsImport[i];
                if (vertNormTexIndex[0].size() < 3) continue; // ignore lines
                const vector<int>& verts = vertNormTexIndex[0];
                const vector<int>& texs = vertNormTexIndex[1];
                const vector<int>& norms = vertNormTexIndex[2];
                vector<int> idxs;
                idxs.resize(verts.size());
                for (unsigned int j = 0; j < verts.size(); j++)
                {
                    idxs[j] = vertTexNormMap[verts[j]][std::make_pair((tex?texs[j]:-1), (useNormals?norms[j]:0))];
                    if ((unsigned)idxs[j] >= (unsigned)nbVOut)
                    {
                        std::cerr << "ERROR(VisualModelImpl): index "<<idxs[j]<<" out of range\n";
                        idxs[j] = 0;
                    }
                }

                if (verts.size() == 4)
                {
                    // quad
                    quads.push_back(helper::make_array(idxs[0],idxs[1],idxs[2],idxs[3]));
                }
                else
                {
                    // triangle(s)
                    for (unsigned int j = 2; j < verts.size(); j++)
                    {
                        triangles.push_back(helper::make_array(idxs[0],idxs[j-1],idxs[j]));
                    }
                }
            }

            std::cout << "Facets Export size : ";
            if (!triangles.empty())
                std::cout << triangles.size() << " triangles";
            if (!quads.empty())
                std::cout << quads.size() << " quads";
            std::cout << "." << std::endl;

            //for (unsigned int i = 0; i < triangles.size() ; i++)
            //    std::cout << "T"<<i<<": "<<triangles[i][0]<<" "<<triangles[i][1]<<" "<<triangles[i][2]<<std::endl;

            computeNormals();
            computeBBox();
        }
    }
    else
    {
        useTopology = true;
        modified = true;
    }

    if (!xformsModified)
    {
        // add one identity matrix
        xforms.resize(1);
    }
    return true;
}

void VisualModelImpl::applyTranslation(double dx, double dy, double dz)
{
    Vector3 d((GLfloat)dx,(GLfloat)dy,(GLfloat)dz);
    VecCoord& x = *getVecX();
    for (unsigned int i = 0; i < x.size(); i++)
    {
        x[i] += d;
    }
    update();
}

void VisualModelImpl::applyRotation(Quat q)
{
    VecCoord& x = *getVecX();
    for (unsigned int i = 0; i < x.size(); i++)
    {
        x[i] = q.rotate(x[i]);
    }
    update();
}

void VisualModelImpl::applyScale(double scale)
{
    VecCoord& x = *getVecX();
    for (unsigned int i = 0; i < x.size(); i++)
    {
        x[i] *= (GLfloat) scale;
    }
    update();
}

void VisualModelImpl::applyUVTranslation(double dU, double dV)
{
    for (unsigned int i = 0; i < vtexcoords.size(); i++)
    {
        vtexcoords[i][0] += (GLfloat) dU;
        vtexcoords[i][1] += (GLfloat) dV;
    }
}

void VisualModelImpl::applyUVScale(double scaleU, double scaleV)
{
    for (unsigned int i = 0; i < vtexcoords.size(); i++)
    {
        vtexcoords[i][0] *= (GLfloat) scaleU;
        vtexcoords[i][1] *= (GLfloat) scaleV;
    }
}

void VisualModelImpl::init()
{
    VisualModel::init();
    update();
}

void VisualModelImpl::computeNormals()
{
    if (vertNormIdx.empty())
    {
        int nbn = vertices.size();
        /*		std::cerr << "nb of visual vertices"<<nbn<<std::endl;
        		std::cerr << "nb of visual triangles"<<triangles.size()<<std::endl; */

        ResizableExtVector<Coord>& normals = vnormals;

        normals.resize(nbn);
        for (int i = 0; i < nbn; i++)
            normals[i].clear();

        for (unsigned int i = 0; i < triangles.size() ; i++)
        {
            const Coord  v1 = vertices[triangles[i][0]];
            const Coord  v2 = vertices[triangles[i][1]];
            const Coord  v3 = vertices[triangles[i][2]];
            Coord n = cross(v2-v1, v3-v1);
            n.normalize();
            normals[triangles[i][0]] += n;
            normals[triangles[i][1]] += n;
            normals[triangles[i][2]] += n;
        }

        for (unsigned int i = 0; i < quads.size() ; i++)
        {
            const Coord & v1 = vertices[quads[i][0]];
            const Coord & v2 = vertices[quads[i][1]];
            const Coord & v3 = vertices[quads[i][3]];
            Coord n = cross(v2-v1, v3-v1);
            n.normalize();
            normals[quads[i][0]] += n;
            normals[quads[i][1]] += n;
            normals[quads[i][2]] += n;
            normals[quads[i][3]] += n;
        }

        for (unsigned int i = 0; i < normals.size(); i++)
        {
            normals[i].normalize();
        }
    }
    else
    {
        vector<Coord> normals;
        int nbn = 0;
        for (unsigned int i = 0; i < vertNormIdx.size(); i++)
            if (vertNormIdx[i] >= nbn)
                nbn = vertNormIdx[i]+1;

        normals.resize(nbn);
        for (int i = 0; i < nbn; i++)
            normals[i].clear();

        for (unsigned int i = 0; i < triangles.size() ; i++)
        {
            const Coord & v1 = vertices[triangles[i][0]];
            const Coord & v2 = vertices[triangles[i][1]];
            const Coord & v3 = vertices[triangles[i][2]];
            Coord n = cross(v2-v1, v3-v1);
            n.normalize();
            normals[vertNormIdx[triangles[i][0]]] += n;
            normals[vertNormIdx[triangles[i][1]]] += n;
            normals[vertNormIdx[triangles[i][2]]] += n;
        }

        for (unsigned int i = 0; i < quads.size() ; i++)
        {
            const Coord & v1 = vertices[quads[i][0]];
            const Coord & v2 = vertices[quads[i][1]];
            const Coord & v3 = vertices[quads[i][3]];
            Coord n = cross(v2-v1, v3-v1);
            n.normalize();
            normals[vertNormIdx[quads[i][0]]] += n;
            normals[vertNormIdx[quads[i][1]]] += n;
            normals[vertNormIdx[quads[i][2]]] += n;
            normals[vertNormIdx[quads[i][3]]] += n;
        }

        for (unsigned int i = 0; i < normals.size(); i++)
        {
            normals[i].normalize();
        }
        vnormals.resize(vertices.size());
        for (unsigned int i = 0; i < vertices.size(); i++)
        {
            vnormals[i] = normals[vertNormIdx[i]];
        }
    }
}

void VisualModelImpl::computeBBox()
{
    const VecCoord& x = vertices;
    Vec3f minBBox(1e10,1e10,1e10);
    Vec3f maxBBox(-1e10,-1e10,-1e10);
    for (unsigned int i = 0; i < x.size(); i++)
    {
        const Coord& p = x[i];
        for (int c=0; c<3; c++)
        {
            if (p[c] > maxBBox[c]) maxBBox[c] = p[c];
            if (p[c] < minBBox[c]) minBBox[c] = p[c];
        }
    }
    bbox[0] = minBBox;
    bbox[1] = maxBBox;
}

bool VisualModelImpl::addBBox(double* minBBox, double* maxBBox)
{
    if (bbox[0][0] > bbox[1][0]) return false;
    for (unsigned int i=0; i<xforms.size(); i++)
    {
        const Vec3f& center = xforms[i].getCenter();
        const Quat& orientation = xforms[i].getOrientation();
        for (int j=0; j<8; j++)
        {
            Coord p ( bbox[(j>>0)&1][0], bbox[(j>>1)&1][1], bbox[(j>>2)&1][2]);
            p = orientation.rotate(p) + center;
            for (int c=0; c<3; c++)
            {
                if (p[c] > maxBBox[c]) maxBBox[c] = p[c];
                if (p[c] < minBBox[c]) minBBox[c] = p[c];
            }
        }
    }
    return true;
}

void VisualModelImpl::flipFaces()
{
    for (unsigned int i = 0; i < triangles.size() ; i++)
    {
        int temp = triangles[i][1];
        triangles[i][1] = triangles[i][2];
        triangles[i][2] = temp;
    }
    for (unsigned int i = 0; i < vnormals.size(); i++)
    {
        vnormals[i] = -vnormals[i];
    }
}

void VisualModelImpl::setColor(float r, float g, float b, float a)
{
    helper::io::Mesh::Material M = material.getValue();
    M.setColor(r,g,b,a);
    material.setValue(M);
}

static int hexval(char c)
{
    if (c>='0' && c<='9') return c-'0';
    else if (c>='a' && c<='f') return (c-'a')+10;
    else if (c>='A' && c<='F') return (c-'A')+10;
    else return 0;
}

void VisualModelImpl::setColor(std::string color)
{
    if (color.empty()) return;
    float r = 1.0f;
    float g = 1.0f;
    float b = 1.0f;
    float a = 1.0f;
    if (color[0]>='0' && color[0]<='9')
    {
        sscanf(color.c_str(),"%f %f %f %f", &r, &g, &b, &a);
    }
    else if (color[0]=='#' && color.length()>=7)
    {
        r = (hexval(color[1])*16+hexval(color[2]))/255.0f;
        g = (hexval(color[3])*16+hexval(color[4]))/255.0f;
        b = (hexval(color[5])*16+hexval(color[6]))/255.0f;
        if (color.length()>=9)
            a = (hexval(color[7])*16+hexval(color[8]))/255.0f;
    }
    else if (color[0]=='#' && color.length()>=4)
    {
        r = (hexval(color[1])*17)/255.0f;
        g = (hexval(color[2])*17)/255.0f;
        b = (hexval(color[3])*17)/255.0f;
        if (color.length()>=5)
            a = (hexval(color[4])*17)/255.0f;
    }
    else if (color == "white")    { r = 1.0f; g = 1.0f; b = 1.0f; }
    else if (color == "black")    { r = 0.0f; g = 0.0f; b = 0.0f; }
    else if (color == "red")      { r = 1.0f; g = 0.0f; b = 0.0f; }
    else if (color == "green")    { r = 0.0f; g = 1.0f; b = 0.0f; }
    else if (color == "blue")     { r = 0.0f; g = 0.0f; b = 1.0f; }
    else if (color == "cyan")     { r = 0.0f; g = 1.0f; b = 1.0f; }
    else if (color == "magenta")  { r = 1.0f; g = 0.0f; b = 1.0f; }
    else if (color == "yellow")   { r = 1.0f; g = 1.0f; b = 0.0f; }
    else if (color == "gray")     { r = 0.5f; g = 0.5f; b = 0.5f; }
    else
    {
        std::cerr << "Unknown color "<<color<<std::endl;
        return;
    }
    setColor(r,g,b,a);
}

void VisualModelImpl::update()
{
    if (modified && !vertices.empty() || useTopology)
    {
        if (useTopology)
        {
            /** HD : build also a Ogl description from main Topology. But it needs to be build only since the topology update
            is taken care of by the handleTopologyChange() routine */
            sofa::core::componentmodel::topology::BaseTopology* pst = dynamic_cast<sofa::core::componentmodel::topology::BaseTopology *>(getContext()->getMainTopology());
            if (pst)
            {
                useTopology=false;
                computeMeshFromTopology(pst);
            }
            else
            {
                topology::MeshTopology* topology = dynamic_cast<topology::MeshTopology*>(getContext()->getTopology());
                if (topology != NULL && topology->getRevision() != lastMeshRev)
                {
                    computeMesh(topology);
                }
            }
        }
        computePositions();
        computeNormals();
        computeBBox();
        modified = false;
    }
}


void VisualModelImpl::computePositions()
{
    if (!vertPosIdx.empty())
    {
        // Need to transfer positions
        for (unsigned int i=0 ; i < vertices.size(); ++i)
            vertices[i] = (*inputVertices)[vertPosIdx[i]];
    }
}

void VisualModelImpl::computeMesh(topology::MeshTopology* topology)
{
    if (vertices.empty())
    {
        if (!topology->hasPos()) return;
        vertices.resize(topology->getNbPoints());
        for (unsigned int i=0; i<vertices.size(); i++)
        {
            vertices[i][0] = (Real)topology->getPX(i);
            vertices[i][1] = (Real)topology->getPY(i);
            vertices[i][2] = (Real)topology->getPZ(i);
        }
    }

    lastMeshRev = topology->getRevision();
    const vector<topology::MeshTopology::Triangle>& inputTriangles = topology->getTriangles();
    triangles.resize(inputTriangles.size());
    for (unsigned int i=0; i<triangles.size(); ++i)
        triangles[i] = inputTriangles[i];
    const vector<topology::MeshTopology::Quad>& inputQuads = topology->getQuads();
    quads.resize(inputQuads.size());
    for (unsigned int i=0; i<quads.size(); ++i)
        quads[i] = inputQuads[i];
}

void VisualModelImpl::computeMeshFromTopology(sofa::core::componentmodel::topology::BaseTopology* bt)
{
    if (vertices.empty())
    {
        BaseMechanicalState *bs= dynamic_cast<BaseMechanicalState *>(getContext()->getMechanicalState());
        assert(bs);
        assert(bs->getSize());
        vertices.resize(bs->getSize());
        for (unsigned int i=0; i<vertices.size(); i++)
        {
            vertices[i][0] = (Real)bs->getPX(i);
            vertices[i][1] = (Real)bs->getPY(i);
            vertices[i][2] = (Real)bs->getPZ(i);
        }
    }
    TopologyContainer *container=bt->getTopologyContainer();
    sofa::component::topology::TriangleSetTopologyContainer *tstc= dynamic_cast<sofa::component::topology::TriangleSetTopologyContainer *>(container);
    if (tstc)
    {
        const std::vector<sofa::component::topology::Triangle> &triangleArray=tstc->getTriangleArray();
        triangles.resize(triangleArray.size());
        for (unsigned int i=0; i<tstc->getNumberOfTriangles(); ++i)
            triangles[i] = triangleArray[i];
    }
}
void VisualModelImpl::handleTopologyChange()
{
    sofa::core::componentmodel::topology::BaseTopology *topology = static_cast<sofa::core::componentmodel::topology::BaseTopology *>(getContext()->getMainTopology());

    std::list<const TopologyChange *>::const_iterator itBegin=topology->firstChange();
    std::list<const TopologyChange *>::const_iterator itEnd=topology->lastChange();

    while( itBegin != itEnd )
    {
        core::componentmodel::topology::TopologyChangeType changeType = (*itBegin)->getChangeType();
        // Since we are using identifier, we can safely use C type casts.
        switch( changeType )
        {

        case core::componentmodel::topology::TRIANGLESADDED:
        {
            const sofa::component::topology::TrianglesAdded *ta=dynamic_cast< const sofa::component::topology::TrianglesAdded * >( *itBegin );
            Triangle t;
            for (unsigned int i=0; i<ta->getNbAddedTriangles(); ++i)
            {
                t[0]=(int)(ta->triangleArray[i])[0];
                t[1]=(int)(ta->triangleArray[i])[1];
                t[2]=(int)(ta->triangleArray[i])[2];
                triangles.push_back(t);
            }
            break;
        }

        case core::componentmodel::topology::TRIANGLESREMOVED:
        {
            const std::vector<unsigned int> &tab = ( dynamic_cast< const sofa::component::topology::TrianglesRemoved *>( *itBegin ) )->getArray();
            unsigned int last = triangles.size() -1;
            Triangle tmp;
            for (unsigned int i = 0; i < tab.size(); ++i)
            {
                tmp = triangles[tab[i]];
                triangles[tab[i]] = triangles[last];
                triangles[last] = tmp;

                --last;
            }
            triangles.resize( triangles.size() - tab.size() );
            break;
        }

        default:
            // Ignore events that are not Triangle  related.
            break;
        }; // switch( changeType )

        ++itBegin;
    } // while( changeIt != last; )
}
void VisualModelImpl::initTextures()
{
    //if (tex)
    //{
    //    tex->init();
    //}
}

void VisualModelImpl::exportOBJ(std::string name, std::ostream* out, std::ostream* mtl, int& vindex, int& nindex, int& tindex)
{
    *out << "g "<<name<<"\n";

    if (mtl != NULL) // && !material.name.empty())
    {
        std::string name; // = material.name;
        if (name.empty())
        {
            static int count = 0;
            std::ostringstream o; o << "mat" << ++count;
            name = o.str();
        }
        *mtl << "newmtl "<<name<<"\n";
        *mtl << "illum 4\n";
        if (material.getValue().useAmbient)
            *mtl << "Ka "<<material.getValue().ambient[0]<<' '<<material.getValue().ambient[1]<<' '<<material.getValue().ambient[2]<<"\n";
        if (material.getValue().useDiffuse)
            *mtl << "Kd "<<material.getValue().diffuse[0]<<' '<<material.getValue().diffuse[1]<<' '<<material.getValue().diffuse[2]<<"\n";
        *mtl << "Tf 1.00 1.00 1.00\n";
        *mtl << "Ni 1.00\n";
        if (material.getValue().useSpecular)
            *mtl << "Ks "<<material.getValue().specular[0]<<' '<<material.getValue().specular[1]<<' '<<material.getValue().specular[2]<<"\n";
        if (material.getValue().useShininess)
            *mtl << "Ns "<<material.getValue().shininess<<"\n";
        if (material.getValue().useDiffuse && material.getValue().diffuse[3]<1.0)
            *mtl << "Tf "<<material.getValue().diffuse[3]<<' '<<material.getValue().diffuse[3]<<' '<<material.getValue().diffuse[3]<<"\n";

        *out << "usemtl "<<name<<'\n';
    }
    const ResizableExtVector<Coord>& x = *inputVertices;

    int nbv = x.size();

    for (unsigned int i=0; i<x.size(); i++)
    {
        *out << "v "<< std::fixed << x[i][0]<<' '<< std::fixed <<x[i][1]<<' '<< std::fixed <<x[i][2]<<'\n';
    }

    int nbn = 0;

    if (vertNormIdx.empty())
    {
        nbn = vnormals.size();
        for (unsigned int i=0; i<vnormals.size(); i++)
        {
            *out << "vn "<< std::fixed << vnormals[i][0]<<' '<< std::fixed <<vnormals[i][1]<<' '<< std::fixed <<vnormals[i][2]<<'\n';
        }
    }
    else
    {
        for (unsigned int i = 0; i < vertNormIdx.size(); i++)
        {
            if (vertNormIdx[i] >= nbn)
                nbn = vertNormIdx[i]+1;
        }
        vector<int> normVertIdx(nbn);
        for (unsigned int i = 0; i < vertNormIdx.size(); i++)
        {
            normVertIdx[vertNormIdx[i]]=i;
        }
        for (int i = 0; i < nbn; i++)
        {
            int j = normVertIdx[i];
            *out << "vn "<< std::fixed << vnormals[j][0]<<' '<< std::fixed <<vnormals[j][1]<<' '<< std::fixed <<vnormals[j][2]<<'\n';
        }
    }

    int nbt = 0;
    if (!vtexcoords.empty())
    {
        nbt = vtexcoords.size();
        for (unsigned int i=0; i<vtexcoords.size(); i++)
        {
            *out << "vt "<< std::fixed << vtexcoords[i][0]<<' '<< std::fixed <<vtexcoords[i][1]<<'\n';
        }
    }

    for (unsigned int i = 0; i < triangles.size() ; i++)
    {
        *out << "f";
        for (int j=0; j<3; j++)
        {
            int i0 = triangles[i][j];
            int i_p = vertPosIdx.empty() ? i0 : vertPosIdx[i0];
            int i_n = vertNormIdx.empty() ? i0 : vertNormIdx[i0];
            if (vtexcoords.empty())
                *out << ' ' << i_p+vindex+1 << "//" << i_n+nindex+1;
            else
                *out << ' ' << i_p+vindex+1 << '/' << i0+tindex+1 << '/' << i_n+nindex+1;
        }
        *out << '\n';
    }
    for (unsigned int i = 0; i < quads.size() ; i++)
    {
        *out << "f";
        for (int j=0; j<4; j++)
        {
            int i0 = quads[i][j];
            int i_p = vertPosIdx.empty() ? i0 : vertPosIdx[i0];
            int i_n = vertNormIdx.empty() ? i0 : vertNormIdx[i0];
            if (vtexcoords.empty())
                *out << ' ' << i_p+vindex+1 << "//" << i_n+nindex+1;
            else
                *out << ' ' << i_p+vindex+1 << '/' << i0+tindex+1 << '/' << i_n+nindex+1;
        }
        *out << '\n';
    }
    *out << std::endl;
    vindex+=nbv;
    nindex+=nbn;
    tindex+=nbt;
}

} // namespace visualmodel

} // namespace component

} // namespace sofa

