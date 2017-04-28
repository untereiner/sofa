/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, development version     *
*                (c) 2006-2017 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_PARALLELTYPES_H
#define SOFA_PARALLELTYPES_H

#include <sofa/defaulttype/Vec.h>

#include <sofa/helper/accessor.h>
#include <sofa/helper/vector.h>
#include <sofa/helper/random.h>
#include <sofa/defaulttype/MapMapSparseMatrix.h>
#include <iostream>
#include <algorithm>
#include <memory>
#include <sofa/helper/logging/Messaging.h>
#include "ParallelMemoryManager.h"
#include "vector_parallel.h"

namespace sofa
{

namespace defaulttype
{

template<class T>
class ParallelVector : public helper::vector<T,ParallelMemoryManager<T> >
{
public :
    typedef size_t size_type;

    ParallelVector() : helper::vector<T,ParallelMemoryManager<T> >() {}

    ParallelVector(size_type n) : helper::vector<T,ParallelMemoryManager<T> >(n) {}

    ParallelVector(const helper::vector<T,ParallelMemoryManager< T > >& v) : helper::vector<T,ParallelMemoryManager<T> >(v) {}

};

template<class TCoord, class TDeriv, class TReal = typename TCoord::value_type>
class ParallelVectorTypes
{
public:
    typedef TCoord Coord;
    typedef TDeriv Deriv;
    typedef TReal Real;
    typedef helper::vector<Coord> VecCoord;
    typedef helper::vector<Deriv> VecDeriv;
    typedef helper::vector<Real> VecReal;

    enum { spatial_dimensions = Coord::spatial_dimensions };
    enum { coord_total_size = Coord::total_size };
    enum { deriv_total_size = Deriv::total_size };

    typedef Coord CPos;
    static const CPos& getCPos(const Coord& c) { return c; }
    static void setCPos(Coord& c, const CPos& v) { c = v; }
    typedef Deriv DPos;
    static const DPos& getDPos(const Deriv& d) { return d; }
    static void setDPos(Deriv& d, const DPos& v) { d = v; }

    typedef MapMapSparseMatrix<Deriv> MatrixDeriv;


protected:

    /// @internal size dependant specializations
    /// @{

    /// default implementation for size >= 3
    template<int N, class T>
    struct Impl
    {
        static void set( Coord& c, T x, T y, T z )
        {
            c[0] = (Real)x;
            c[1] = (Real)y;
            c[2] = (Real)z;
        }

        static void get( T& x, T& y, T& z, const Coord& c )
        {
            x = (T) c[0];
            y = (T) c[1];
            z = (T) c[2];
        }

        static void add( Coord& c, T x, T y, T z )
        {
            c[0] += (Real)x;
            c[1] += (Real)y;
            c[2] += (Real)z;
        }
    };

    /// specialization for size == 2
    template<class T>
    struct Impl<2,T>
    {
        static void set( Coord& c, T x, T y, T )
        {
            c[0] = (Real)x;
            c[1] = (Real)y;
        }

        static void get( T& x, T& y, T& z, const Coord& c )
        {
            x = (T) c[0];
            y = (T) c[1];
            z = (T) 0;
        }

        static void add( Coord& c, T x, T y, T )
        {
            c[0] += (Real)x;
            c[1] += (Real)y;
        }
    };

    /// specialization for size == 1
    template<class T>
    struct Impl<1,T>
    {
        static void set( Coord& c, T x, T, T )
        {
            c[0] = (Real)x;
        }

        static void get( T& x, T& y, T& z, const Coord& c )
        {
            x = (T) c[0];
            y = (T) 0;
            z = (T) 0;
        }

        static void add( Coord& c, T x, T, T )
        {
            c[0] += (Real)x;
        }
    };

    ///@}



public:

    template<typename T>
    static void set(Coord& c, T x, T y, T z)
    {
        Impl<spatial_dimensions,T>::set(c,x,y,z);
    }

    template<typename T>
    static void get(T& x, T& y, T& z, const Coord& c)
    {
        Impl<spatial_dimensions,T>::get(x,y,z,c);
    }

    /// Return a Deriv with random value. Each entry with magnitude smaller than the given value.
    static Deriv randomDeriv( Real minMagnitude, Real maxMagnitude )
    {
        Deriv result;
        set( result, Real(helper::drand(minMagnitude,maxMagnitude)), Real(helper::drand(minMagnitude,maxMagnitude)), Real(helper::drand(minMagnitude,maxMagnitude)) );
        return result;
    }

    static Deriv coordDifference(const Coord& c1, const Coord& c2)
    {
        return (Deriv)(c1-c2);
    }

    template<typename T>
    static void add(Coord& c, T x, T y, T z)
    {
        Impl<spatial_dimensions,T>::add(c,x,y,z);
    }

    static const char* Name();

    static Coord interpolate(const helper::vector< Coord > &ancestors, const helper::vector< Real > &coefs)
    {
        assert(ancestors.size() == coefs.size());

        Coord c;

        for (std::size_t i = 0; i < ancestors.size(); i++)
        {
            c += ancestors[i] * coefs[i];
        }

        return c;
    }
};


/// Custom vector allocator class allowing data to be allocated at a specific location (such as for transmission through DMA, PCI-Express, Shared Memory, Network)
template<class T>
class ParallelExtVectorAllocator
{
public:
    typedef T              value_type;
    typedef size_t   size_type;
    virtual ~ParallelExtVectorAllocator() {}
    virtual void resize(value_type*& data, size_type size, size_type& maxsize, size_type& cursize)=0;
    virtual void close(value_type*& data)=0;
    virtual void cloneTo( std::unique_ptr<ParallelExtVectorAllocator>& clone ) = 0; ///< clone "this" into given "clone"
};

/// Custom vector class.
///
/// This class allows custom buffer allocation while not having any virtual methods using a bridge pattern with ParallelExtVectorAllocator
template<class T>
class ParallelExtVector
{
public:
    typedef T               value_type;
    typedef size_t  size_type;
    typedef T&              reference;
    typedef const T&        const_reference;
    typedef T*              iterator;
    typedef const T*        const_iterator;

protected:
    value_type* data;
    size_type   maxsize;
    size_type   cursize;
    std::unique_ptr<ParallelExtVectorAllocator<T> > allocator;

public:
    explicit ParallelExtVector(ParallelExtVectorAllocator<T>* alloc = NULL) : data(NULL),  maxsize(0), cursize(0), allocator(alloc) {}
    ParallelExtVector(size_type size, ParallelExtVectorAllocator<T>* alloc) : data(NULL), maxsize(0), cursize(0), allocator(alloc) { resize(size); }
    ~ParallelExtVector() { if (allocator.get()) allocator->close(data); }

    void init() {}

    void setAllocator(ParallelExtVectorAllocator<T>* alloc)
    {
        if (alloc != allocator.get())
        {
            if (cursize)
            {
                value_type* oldData = data;
                size_type size = cursize;

                data = NULL;
                maxsize = 0;
                cursize = 0;
                if(alloc)
                    alloc->resize(data, size, maxsize, cursize);
                if(data != 0 && oldData != 0)
                {
                    std::copy(oldData, oldData + size, data);
                }
                if(allocator.get())
                    allocator->close(oldData);
            }
            allocator.reset(alloc);
        }
    }
    void setData(value_type* d, size_type s) { data=d; maxsize=s; cursize=s; }
    T* getData() { return this->data; }
    const T* getData() const { return this->data; }

    value_type& operator[](size_type i) { return data[i]; }
    const value_type& operator[](size_type i) const { return data[i]; }
    size_type size() const { return cursize; }
    bool empty() const { return cursize==0; }
    void reserve(size_type size)
    {
        if (size <= maxsize)
            return;
        size_type temp = cursize;
        if (allocator.get())
            allocator->resize(data, size, maxsize, temp);
        else
        {
            msg_error("VecTypes") << "reserve: invalid reserve request ("<<size<<">"<<maxsize<<") on external vector without allocator.";
        }
    }
    void resize(size_type size)
    {
        if (size <= maxsize)
            cursize = size;
        else if (allocator.get())
            allocator->resize(data, size, maxsize, cursize);
        else
        {
            cursize = maxsize;
            msg_error("VecTypes") << "resize: invalid resize request ("<<size<<">"<<maxsize<<") on external vector without allocator.";
        }
    }
    void clear()
    {
        resize(0);
    }
    void push_back(const T& v)
    {
        size_type i = this->size();
        resize(i+1);
        (*this)[i] = v;
    }
    T* begin() { return getData(); }
    const T* begin() const { return getData(); }
    T* end() { return getData()+size(); }
    const T* end() const { return getData()+size(); }

    ParallelExtVector& operator=(const ParallelExtVector& ev)
    {
        if(allocator.get())
        {
            allocator->close(data);
        }
        ev.allocator->cloneTo( allocator );
        if(allocator.get())
        {
            allocator->resize(data, ev.cursize, maxsize, cursize);
            if(data != 0)
            {
                std::copy(ev.begin(), ev.end(), data);
            }
        }
        return *this;
    }

    ParallelExtVector(const ParallelExtVector& ev)
        : data(0), maxsize(0), cursize(0)
    {
        ev.allocator->cloneTo( allocator );
        allocator->resize(data, ev.cursize, maxsize, cursize);
        if(data != 0)
        {
            std::copy(ev.begin(), ev.end(), data);
        }
    }


/// Output stream
    inline friend std::ostream& operator<< ( std::ostream& os, const ParallelExtVector<T>& vec )
    {
        if( vec.size()>0 )
        {
            for( std::size_t i=0; i<vec.size()-1; ++i ) os<<vec[i]<<" ";
            os<<vec[vec.size()-1];
        }
        return os;
    }

/// Input stream
    inline friend std::istream& operator>> ( std::istream& in, ParallelExtVector<T>& vec )
    {
        T t;
        vec.clear();
        while(in>>t)
        {
            vec.push_back(t);
        }
        if( in.rdstate() & std::ios_base::eofbit ) { in.clear(); }
        return in;
    }

};

template<class T>
class ParallelDefaultAllocator : public ParallelExtVectorAllocator<T>
{
public:
    typedef typename ParallelExtVectorAllocator<T>::value_type value_type;
    typedef typename ParallelExtVectorAllocator<T>::size_type size_type;
    virtual void close(value_type*& data)
    {
        delete[] data;
        data = 0;
    }
    virtual void resize(value_type*& data, size_type size, size_type& maxsize, size_type& cursize)
    {
        if (size > maxsize)
        {
            T* oldData = data;
            maxsize = (size > 2*maxsize ? size : 2*maxsize);
            data = new T[maxsize];
            if(oldData)
            {
                std::copy(oldData, oldData+cursize, data);
                delete[] oldData;
            }
        }
        cursize = size;
    }
   virtual void cloneTo( std::unique_ptr< ParallelExtVectorAllocator<T> >& clone )
    {
        clone.reset( new ParallelDefaultAllocator<T> );
    }
};

/// Resizable custom vector class using ParallelDefaultAllocator
template<class T>
class ResizableParallelExtVector : public ParallelExtVector<T>
{
public:
    typedef typename ParallelExtVector<T>::value_type value_type;
    typedef typename ParallelExtVector<T>::size_type size_type;
    ResizableParallelExtVector()
        : ParallelExtVector<T>(new ParallelDefaultAllocator<T>)
    {
    }

    ResizableParallelExtVector(const ResizableParallelExtVector& ev)
        :ParallelExtVector<T>(ev)
    {
        this->setAllocator(new ParallelDefaultAllocator<T>);
    }
};

template<class TCoord, class TDeriv, class TReal = typename TCoord::value_type>
class ParallelExtVectorTypes
{
public:
    typedef TCoord Coord;
    typedef TDeriv Deriv;
    typedef TReal Real;
    typedef ResizableParallelExtVector<Coord> VecCoord;
    typedef ResizableParallelExtVector<Deriv> VecDeriv;
    typedef ResizableParallelExtVector<Real> VecReal;

    enum { spatial_dimensions = Coord::spatial_dimensions };
    enum { coord_total_size = Coord::total_size };
    enum { deriv_total_size = Deriv::total_size };

    typedef Coord CPos;
    static const CPos& getCPos(const Coord& c) { return c; }
    static void setCPos(Coord& c, const CPos& v) { c = v; }
    typedef Deriv DPos;
    static const DPos& getDPos(const Deriv& d) { return d; }
    static void setDPos(Deriv& d, const DPos& v) { d = v; }

    typedef MapMapSparseMatrix<Deriv> MatrixDeriv;


protected:

    /// @internal size dependant specializations
    /// @{

    /// default implementation for size >= 3
    template<int N, class T>
    struct Impl
    {
        static void set( Coord& c, T x, T y, T z )
        {
            c[0] = (Real)x;
            c[1] = (Real)y;
            c[2] = (Real)z;
        }

        static void get( T& x, T& y, T& z, const Coord& c )
        {
            x = (T) c[0];
            y = (T) c[1];
            z = (T) c[2];
        }

        static void add( Coord& c, T x, T y, T z )
        {
            c[0] += (Real)x;
            c[1] += (Real)y;
            c[2] += (Real)z;
        }
    };

    /// specialization for size == 2
    template<class T>
    struct Impl<2,T>
    {
        static void set( Coord& c, T x, T y, T )
        {
            c[0] = (Real)x;
            c[1] = (Real)y;
        }

        static void get( T& x, T& y, T& z, const Coord& c )
        {
            x = (T) c[0];
            y = (T) c[1];
            z = (T) 0;
        }

        static void add( Coord& c, T x, T y, T )
        {
            c[0] += (Real)x;
            c[1] += (Real)y;
        }
    };

    /// specialization for size == 1
    template<class T>
    struct Impl<1,T>
    {
        static void set( Coord& c, T x, T, T )
        {
            c[0] = (Real)x;
        }

        static void get( T& x, T& y, T& z, const Coord& c )
        {
            x = (T) c[0];
            y = (T) 0;
            z = (T) 0;
        }

        static void add( Coord& c, T x, T, T )
        {
            c[0] += (Real)x;
        }
    };

    ///@}

public:


    template<typename T>
    static void set(Coord& c, T x, T y, T z)
    {
        Impl<spatial_dimensions,T>::set(c,x,y,z);
    }

    template<typename T>
    static void get(T& x, T& y, T& z, const Coord& c)
    {
        Impl<spatial_dimensions,T>::get(x,y,z,c);
    }

    template<typename T>
    static void add(Coord& c, T x, T y, T z)
    {
        Impl<spatial_dimensions,T>::add(c,x,y,z);
    }

    static const char* Name();

    static Coord interpolate(const helper::vector< Coord > & ancestors, const helper::vector< Real > & coefs)
    {
        assert(ancestors.size() == coefs.size());

        Coord c;

        for (std::size_t i = 0; i < ancestors.size(); i++)
        {
            c += ancestors[i] * coefs[i];
        }

        return c;
    }
};



#ifndef SOFA_FLOAT

/// 3D DOFs, double precision
typedef ParallelVectorTypes<Vec3d,Vec3d,double> ParallelVec3dTypes;
template<> inline const char* ParallelVec3dTypes::Name() { return "ParallelVec3d"; }
/// 3D external DOFs, double precision
typedef ParallelExtVectorTypes<Vec3d,Vec3d,double> ParallelExtVec3dTypes;
template<> inline const char* ParallelExtVec3dTypes::Name() { return "ParallelExtVec3d"; }

/// 2D DOFs, double precision
typedef ParallelVectorTypes<Vec2d,Vec2d,double> ParallelVec2dTypes;
template<> inline const char* ParallelVec2dTypes::Name() { return "ParallelVec2d"; }
/// 2D external DOFs, double precision
typedef ParallelExtVectorTypes<Vec2d,Vec2d,double> ParallelExtVec2dTypes;
template<> inline const char* ParallelExtVec2dTypes::Name() { return "ParallelExtVec2d"; }

/// 1D DOFs, double precision
typedef ParallelVectorTypes<Vec1d,Vec1d,double> ParallelVec1dTypes;
template<> inline const char* ParallelVec1dTypes::Name() { return "ParallelVec1d"; }
/// 1D external DOFs, double precision
typedef ParallelExtVectorTypes<Vec1d,Vec1d,double> ParallelExtVec1dTypes;
template<> inline const char* ParallelExtVec1dTypes::Name() { return "ParallelExtVec1d"; }

/// 6D DOFs, double precision
typedef ParallelVectorTypes<Vec6d,Vec6d,double> ParallelVec6dTypes;
template<> inline const char* ParallelVec6dTypes::Name() { return "ParallelVec6d"; }
/// 6D external DOFs, double precision
typedef ParallelExtVectorTypes<Vec6d,Vec6d,double> ParallelExtVec6dTypes;
template<> inline const char* ParallelExtVec6dTypes::Name() { return "ParallelExtVec6d"; }


#endif

/*#ifndef SOFA_DOUBLE*/

/// 3f DOFs, single precision
typedef ParallelVectorTypes<Vec3f,Vec3f,float> ParallelVec3fTypes;
template<> inline const char* ParallelVec3fTypes::Name() { return "ParallelVec3f"; }
/// 3f external DOFs, single precision
typedef ParallelExtVectorTypes<Vec3f,Vec3f,float> ParallelExtVec3fTypes;
template<> inline const char* ParallelExtVec3fTypes::Name() { return "ParallelExtVec3f"; }

/// 2f DOFs, single precision
typedef ParallelVectorTypes<Vec2f,Vec2f,float> ParallelVec2fTypes;
template<> inline const char* ParallelVec2fTypes::Name() { return "ParallelVec2f"; }
/// 2f external DOFs, single precision
typedef ParallelExtVectorTypes<Vec2f,Vec2f,float> ParallelExtVec2fTypes;
template<> inline const char* ParallelExtVec2fTypes::Name() { return "ParallelExtVec2f"; }

/// 1f DOFs, single precision
typedef ParallelVectorTypes<Vec1f,Vec1f,float> ParallelVec1fTypes;
template<> inline const char* ParallelVec1fTypes::Name() { return "ParallelVec1f"; }
/// 1f external DOFs, single precision
typedef ParallelExtVectorTypes<Vec1f,Vec1f,float> ParallelExtVec1fTypes;
template<> inline const char* ParallelExtVec1fTypes::Name() { return "ParallelExtVec1f"; }

/// 6f DOFs, single precision
typedef ParallelVectorTypes<Vec6f,Vec6f,float> ParallelVec6fTypes;
template<> inline const char* ParallelVec6fTypes::Name() { return "ParallelVec6f"; }
/// 6f external DOFs, single precision
typedef ParallelExtVectorTypes<Vec6f,Vec6f,float> ParallelExtVec6fTypes;
template<> inline const char* ParallelExtVec6fTypes::Name() { return "ParallelExtVec6f"; }

//#endif



#ifdef SOFA_FLOAT
/// 6D DOFs, single precision (default)
typedef ParallelVec6fTypes ParallelVec6Types;
/// 3D DOFs, single precision (default)
typedef ParallelVec3fTypes ParallelVec3Types;
/// 2D DOFs, single precision (default)
typedef ParallelVec2fTypes ParallelVec2Types;
/// 1D DOFs, single precision (default)
typedef ParallelVec1fTypes ParallelVec1Types;
/// 6D external DOFs, single precision (default)
typedef ExtVec6fTypes ParallelExtVec6Types;
/// 3D external DOFs, single precision (default)
typedef ParallelExtVec3fTypes ParallelExtVec3Types;
/// 2D external DOFs, single precision (default)
typedef ParallelExtVec2fTypes ParallelExtVec2Types;
/// 1D external DOFs, single precision (default)
typedef ParallelExtVec1fTypes ParallelExtVec1Types;
#else
/// 6D DOFs, double precision (default)
typedef ParallelVec6dTypes ParallelVec6Types;
/// 3D DOFs, double precision (default)
typedef ParallelVec3dTypes ParallelVec3Types;
/// 2D DOFs, double precision (default)
typedef ParallelVec2dTypes ParallelVec2Types;
/// 1D DOFs, double precision (default)
typedef ParallelVec1dTypes ParallelVec1Types;
/// 6D external DOFs, double precision (default)
typedef ParallelExtVec6dTypes ParallelExtVec6Types;
/// 3D external DOFs, double precision (default)
typedef ParallelExtVec3dTypes ParallelExtVec3Types;
/// 2D external DOFs, double precision (default)
typedef ParallelExtVec2dTypes ParallelExtVec2Types;
/// 1D external DOFs, double precision (default)
typedef ParallelExtVec1dTypes ParallelExtVec1Types;
#endif


// Specialization of the defaulttype::DataTypeInfo type traits template

template<class T>
struct DataTypeInfo< sofa::defaulttype::ParallelExtVector<T> > : public VectorTypeInfo<sofa::defaulttype::ParallelExtVector<T> >
{
    // Remove copy-on-write behavior which is normally activated for vectors
//    enum { CopyOnWrite     = 0 };

    static std::string name() { std::ostringstream o; o << "ParallelExtVector<" << DataTypeName<T>::name() << ">"; return o.str(); }
};

template<class T>
struct DataTypeInfo< sofa::defaulttype::ResizableParallelExtVector<T> > : public VectorTypeInfo<sofa::defaulttype::ResizableParallelExtVector<T> >
{
    // Remove copy-on-write behavior which is normally activated for vectors
//    enum { CopyOnWrite     = 0 };

    static std::string name() { std::ostringstream o; o << "ResizableParallelExtVector<" << DataTypeName<T>::name() << ">"; return o.str(); }
};

} // namespace defaulttype


namespace helper
{

template<class T>
class ReadAccessor< defaulttype::ParallelExtVector<T> > : public ReadAccessorVector< defaulttype::ParallelExtVector<T> >
{
public:
    typedef ReadAccessorVector< defaulttype::ParallelExtVector<T> > Inherit;
    typedef typename Inherit::container_type container_type;
    ReadAccessor(const container_type& c) : Inherit(c) {}
};

template<class T>
class WriteAccessor< defaulttype::ParallelExtVector<T> > : public WriteAccessorVector< defaulttype::ParallelExtVector<T> >
{
public:
    typedef WriteAccessorVector< defaulttype::ParallelExtVector<T> > Inherit;
    typedef typename Inherit::container_type container_type;
    WriteAccessor(container_type& c) : Inherit(c) {}
};

template<class T>
class ReadAccessor< defaulttype::ResizableParallelExtVector<T> > : public ReadAccessorVector< defaulttype::ResizableParallelExtVector<T> >
{
public:
    typedef ReadAccessorVector< defaulttype::ResizableParallelExtVector<T> > Inherit;
    typedef typename Inherit::container_type container_type;
    ReadAccessor(const container_type& c) : Inherit(c) {}
};

template<class T>
class WriteAccessor< defaulttype::ResizableParallelExtVector<T> > : public WriteAccessorVector< defaulttype::ResizableParallelExtVector<T> >
{
public:
    typedef WriteAccessorVector< defaulttype::ResizableParallelExtVector<T> > Inherit;
    typedef typename Inherit::container_type container_type;
    WriteAccessor(container_type& c) : Inherit(c) {}
};

} // namespace helper

} // namespace sofa

#endif
