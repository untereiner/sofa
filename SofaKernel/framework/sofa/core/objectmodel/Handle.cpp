/*
 * Handle.cpp
 *
 *  Created on: 22 mars 2017
 *      Author: maxime
 */

#include <sofa/core/objectmodel/Handle.h>

namespace sofa
{

namespace core
{

namespace objectmodel
{
//Handle::Handle( BaseData* d, bool aType ) : handledData(d), handledData(d->m_handledData), accessType(aType), mut(), cond(), next(NULL), nbReaders(0), fusionnedHandle(NULL)
Handle::Handle( BaseData* d, bool aType ) : handledData(d), accessType(aType), mut(), cond(), next(NULL), nbReaders(0), fusionnedHandle(NULL)
{
}

Handle::~Handle( )
{
}




/**
 * @brief if write access : pushs the handle in queue
 *        if read access : queue->tail->nbReaders++ if queue->tail is a read, else pushs the handle in queue
 */
void Handle::request()
{
    PRINT_LOG_MSG("......calling Handle::request with accessType = " << accessType);
    //we need to lock the fifo mutex to access its components
    std::unique_lock<std::mutex> lock_q(*handledData->getExternalMutex());
    // if read access and  the last request in the fifo is a read, we "fusion" this with the last request in the fifo
    if ( !handledData->isempty() && (handledData->getTail())->accessType ) {
      fusionnedHandle = handledData->tail;
      //we need to lock the mutex of handledData->tail to avoid race condition on handledData->tail->nbReaders
      std::unique_lock<std::mutex> lock_h(fusionnedHandle->mut);
      fusionnedHandle->nbReaders++;
      lock_h.unlock();
//        nbReaders++;  //this has no attached readers, so nbReaders is still 0
      PRINT_LOG_MSG(" from Handle::request, fusionnedHandle->nbReaders = " << fusionnedHandle->nbReaders);
    }
    else {
      handledData->push( this );
    }
    lock_q.unlock();
  }
/**
 * @brief if write access : pops the handle from the queue
 *        if read access : pops the handle from the queue if nbReaders==1, else nbReaders--
 *        it is the user's responsability to do acquire before release
 * */
void Handle::release()
{
  PRINT_LOG_MSG("......calling Handle::release with accessType = " << accessType);
  //we need to lock the fifo mutex to access its components
  std::unique_lock<std::mutex> lock_q(*handledData->getExternalMutex());

//      bool doPop = false;  //should we pop the queue ?

//      //we need to lock the mutex of the handle to avoid race condition on nbReaders and multiple handledData->pop() calls
//      std::unique_lock<std::mutex> lock_h(mut);

  if ( accessType && fusionnedHandle!=NULL ) {  //this is a read access attached to a previous read handle
//          //we need to lock the mutexes of both handles to avoid race condition on nbReaders and multiple handledData->pop() calls
//          std::unique_lock<std::mutex> lock_h(mut, std::defer_lock);
//          std::unique_lock<std::mutex> lock_fh(fusionnedHandle->mut, std::defer_lock);
//          std::lock(lock_h, lock_fh);
      // not needed since we use the fifo external mutex
      if ( fusionnedHandle->nbReaders > 0 ) {  //this is not the last remaining read attached to fusionnedHandle (including fusionnedHandle itself)
                                               //if fusionnedHandle->nbReaders == 0, then it means that fusionnedHandle has released its access, and this is the last handle attached to fusionnedHandle
        fusionnedHandle->nbReaders--;  //this is not attached to fusionnedHandle anymore
        fusionnedHandle = NULL;
      }
      else {  //this is the last remaining read attached to fusionnedHandle, so we can pop the handle
            handledData->pop();
            if( (this->next) != NULL ) {
              next->cond.notify_all();    //we notify_all since multiple threads can be waiting on the same condition var if next is a read access
              PRINT_LOG_MSG(" from Handle::release, next handle(s) have been notified !");
            }
      }
//          std::unlock(lock_h, lock_fh);
  }
  else if ( accessType && fusionnedHandle==NULL )  {  //this is a read access not attached to a previous read handle
//        std::unique_lock<std::mutex> lock_h(mut);
    if ( nbReaders > 0 ) {  //this still has attached read handles, so we can not pop the handle
      nbReaders--;
    }
    else {  //this is a solitary read access, so we can pop the handle
          handledData->pop();
          if( (this->next) != NULL ) {
            next->cond.notify_all();    //we notify_all since multiple threads can be waiting on the same condition var if next is a read access
            PRINT_LOG_MSG(" from Handle::release, next handle(s) have been notified !");
          }
    }
//        lock_h.unlock();
  }
  else if ( !accessType ) {  //this a write access, so we can pop the handle
        handledData->pop();
        if( (this->next) != NULL ) {
          next->cond.notify_all();    //we notify_all since multiple threads can be waiting on the same condition var if next is a read access
          PRINT_LOG_MSG("from Handle::release, next handle(s) have been notified !");
        }
  }
  else {
    std::cerr << "ERROR !!! this case should not exist" << std::endl;
  }
  //unlock the fifo mutex
  lock_q.unlock();
}
/**
 * @brief waits for the data access to be granted
 *        should be called after request and before release, might cause deadlocks otherwise
 */
BaseData* Handle::acquire()
{
  PRINT_LOG_MSG("......calling Handle::acquire with accessType = " << accessType);
  //we need to lock the fifo mutex to access its components
  std::unique_lock<std::mutex> lock_q(*handledData->getExternalMutex());
  //wait for this (or the read handle we "fusionned" this with) to be head of the queue
  while( !handledData->isHead(this) && !handledData->isHead(fusionnedHandle) ) {
    PRINT_LOG_MSG(" from Handle::acquire with accessType = " << accessType << " calling wait....");

    cond.wait(lock_q);

    if( handledData->isHead(this) ) {
        PRINT_LOG_MSG(" from Handle::acquire with accessType = " << accessType << " wait done !");
    }
  }
  //return the data
  return handledData->getData();
//      return handledData;
//      return handledData->getData();
}


void Handle::setNext( Handle* elt )
{
  next = elt;
}
Handle* Handle::getNext()
{
  return next;
}
bool Handle::getType()
{
  return accessType;
}



} // namespace objectmodel

} // namespace core

} // namespace sofa


