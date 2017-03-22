#ifndef HANDLE_H
#define HANDLE_H

#include <thread>
#include <mutex>
#include <condition_variable>

#include <sofa/core/objectmodel/ThreadSafeQueue.h>
#include <sofa/core/objectmodel/Data.h>
 
namespace sofa
{

namespace core
{

namespace objectmodel
{

template < class T >
class ThreadSafeQueue;
template < class T >
class Data;

/*
 * @brief Handling of the Data access requests via a thread-safe FIFO
 */
template < class T = void* >
class Handle
{
private:
    ThreadSafeQueue<T>* accessQueue;  ///<FIFO
    bool accessType;    ///<true if Read access, false if Write Access
    std::mutex mut;     ///<mutex used to wait on cond
    std::condition_variable cond;     ///<condition variable used to wait in acquire for the handle to be the head of accessQueue
    Handle<T>* next;    ///<pointer to the next Handle in the FIFO
    unsigned int nbReaders; ///<numbers of following read accesses attached/fusionned with this (0 if no following reads)
    Handle<T>* fusionnedHandle; ///<if this follows other read handle in accessQueue, fusionnedHandle point to the first of the consecutive read handles
                             ///<if fusionnedHandle!=NULL, then this is attached to a previous read handle
public:
    /// Constructor for thread-safe accesses via a FIFO
    Handle( ThreadSafeQueue<T>* fifo, bool aType );
    /// Destructor.  -> Does it need to be explicitely defined ??????
    ~Handle();

    /**
     * @brief if write access : pushs the handle in queue
     *        if read access : queue->tail->nbReaders++ if queue->tail is a read, else pushs the handle in queue
     */
    void request();
    /**
     * @brief if write access : pops the handle from the queue
     *        if read access : pops the handle from the queue if nbReaders==1, else nbReaders--
     *        it is the user's responsability to do acquire before release
     * */
    void release();
    /**
     * @brief waits for the data access to be granted
     *        should be called after request and before release, might cause deadlocks otherwise
     */
    Data<T>* acquire();


    void setNext( Handle<T>* elt );
    Handle<T>* getNext();
    bool getType();
    //we need to lock the mutex before calling any of the following methods on a Handle contained in the queue
    //which methods should I put here ?


};

} // namespace objectmodel

} // namespace core

} // namespace sofa


#endif //HANDLE_H
