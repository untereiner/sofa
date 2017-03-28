#ifndef HANDLE_H
#define HANDLE_H

#include <iostream>
#include <thread>
#include <mutex>
#include <condition_variable>

 
#include <sofa/core/objectmodel/BaseData.h>
 
namespace sofa
{

namespace core
{

namespace objectmodel
{
class BaseData;

class Handle
{
private:
    BaseData* handledData;
    const bool accessType;    //true if Read access, false if Write Access
    std::mutex mut;     //mutex used to wait on cond
    std::condition_variable cond;     //condition variable used to wait in acquire for the handle to be the head of accessQueue
    Handle* next;
    Handle* previousRead;
    bool isAcquired;
    bool isHead;
    unsigned int nbReaders; //numbers of following read accesses attached/fusionned with this (0 if no following reads)
//    std::list< Handle* > consecutiveReads; //if this follows other read handle in accessQueue, fusionnedHandle point to the first of the consecutive read handles
//    //if fusionnedHandle!=NULL, then this is attached to a previous read handle


public:
    Handle( BaseData* d, bool aType );
    ~Handle();

    //if write access : pushs the handle in queue
    //if read access : queue->tail->nbReaders++ if queue->tail is a read, else pushs the handle in queue
    void request();
    //if write access : pops the handle from the queue
    //if read access : pops the handle from the queue if nbReaders==1, else nbReaders--
    //it is the user's responsability to do acquire before release
    void release();
    //waits for the data access to be granted
    //should be called after request and before release, might cause deadlocks otherwise
    BaseData* acquire();
    //if isAcquired, returns a pointer to the handledData
    //should be called after request and before release, might cause deadlocks otherwise
    BaseData* map();


    bool getType();
    //we need to lock the mutex before calling any of the following methods on a Handle contained in the queue
    void setNext( Handle* elt );
    void setPrevious( Handle* elt );
    Handle* getNext();
    void setHead( bool b );
    unsigned int getNbReaders();
    void incrementNbReaders();
    void decrementNbReaders();
    bool acquired();


};

} // namespace objectmodel

} // namespace core

} // namespace sofa

#endif //HANDLE_H
