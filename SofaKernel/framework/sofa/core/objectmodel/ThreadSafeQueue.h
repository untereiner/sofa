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
//#include <sofa/core/objectmodel/Data.h>
#include <sofa/core/objectmodel/Handle.h>
 
namespace sofa
{

namespace core
{

namespace objectmodel
{

//class DDGNode;
//class BaseData;
//template < class T >
//class Data;
class Handle;

class ThreadSafeQueue
{
private:
//    Data* handledData;
//    DDGNode* handledData;
    std::mutex externalMut;
    std::mutex internalMut;
    std::condition_variable cond;
//    Handle* tail;
#ifdef LOG_THREADS
    std::string name;
#endif //LOG_THREADS
public:
    Handle* head;
    Handle* tail;

#ifndef LOG_THREADS
    ThreadSafeQueue( );
//    ThreadSafeQueue( DDGNode* d );
#else //LOG_THREADS
    ThreadSafeQueue( std::string s );
#endif //LOG_THREADS
    ~ThreadSafeQueue( );

    void push(Handle* elt);
    Handle* pop();
    
//    Data* getData()
//    {
//        return handledData;
//    }
//    BaseData* getData();
    std::mutex* getExternalMutex();
    //the following methods are not thread-safe
    //queue->mut should be locked by the user before calling getTail, getHead, isHead or isempty
    Handle* getTail();
    Handle* getHead();
    bool isHead(Handle* elt);
    bool isempty();
};

} //namespace sofa
} //namespace core
} //namespace objectmodel

#endif //FIFO_H
