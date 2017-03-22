#ifndef FIFO_H
#define FIFO_H

#ifdef LOG_THREADS
//#include <omp.h>
#endif //LOG_THREADS
#include <sofa/core/core.h>

#include <thread>
#include <mutex>
#include <condition_variable>

//#include <sofa/core/objectmodel/DDGNode.h>
//#include <sofa/core/objectmodel/BaseData.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/core/objectmodel/Handle.h>
 
namespace sofa
{

namespace core
{

namespace objectmodel
{

//class DDGNode;
//class BaseData;
template < class T >
class Data;
template < class T >
class Handle;

template < class T = void* >
class ThreadSafeQueue
{
private:
    Data<T>* handledData;
//    DDGNode* handledData;
    std::mutex externalMut;
    std::mutex internalMut;
    std::condition_variable cond;
//    Handle* tail;
#ifdef LOG_THREADS
    std::string name;
#endif //LOG_THREADS
public:
    Handle<T>* head;
    Handle<T>* tail;

#ifndef LOG_THREADS
    ThreadSafeQueue( Data<T>* d );
//    ThreadSafeQueue( DDGNode* d );
#else //LOG_THREADS
    ThreadSafeQueue( Data<T>* d, std::string s );
//    ThreadSafeQueue( DDGNode* d, std::string s );
#endif //LOG_THREADS
    ~ThreadSafeQueue( );

    void push(Handle<T>* elt);
    Handle<T>* pop();
    
    Data<T>* getData();
//    BaseData* getData();
    std::mutex* getExternalMutex();

    //the following methods are not thread-safe
    //queue->mut should be locked by the user before calling getTail, getHead, isHead or isempty
    Handle<T>* getTail();
    Handle<T>* getHead();
    bool isHead(Handle<T>* elt);
    bool isempty();
};

} //namespace sofa
} //namespace core
} //namespace objectmodel

#endif //FIFO_H
