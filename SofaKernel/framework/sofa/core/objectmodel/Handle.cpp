#include <sofa/core/objectmodel/Handle.h>

namespace sofa
{

namespace core
{

namespace objectmodel
{
    Handle::Handle( ThreadSafeQueue<T>* fifo, bool aType ) : accessQueue(fifo), accessType(aType), mut(), cond(), next(NULL), nbReaders(0), fusionnedHandle(NULL)
    {
    }

    Handle::~Handle( )
    {
    }

    /*
     *
     */
    void Handle::request() {
      PRINT_LOG_MSG("......calling Handle::request with accessType = " << accessType);
      //we need to lock the fifo mutex to access its components
      std::unique_lock<std::mutex> lock_q(*accessQueue->getExternalMutex());
      // if read access and  the last request in the fifo is a read, we "fusion" this with the last request in the fifo
      if ( !accessQueue->isempty() && (accessQueue->getTail())->accessType ) {
        fusionnedHandle = accessQueue->tail;
        //we need to lock the mutex of accessQueue->tail to avoid race condition on accessQueue->tail->nbReaders
        std::unique_lock<std::mutex> lock_h(fusionnedHandle->mut);
        fusionnedHandle->nbReaders++;
        lock_h.unlock();
//        nbReaders++;  //this has no attached readers, so nbReaders is still 0
        PRINT_LOG_MSG(" from Handle::request, fusionnedHandle->nbReaders = " << fusionnedHandle->nbReaders);
      }
      else {
        accessQueue->push( this );
      }
      lock_q.unlock();
    }

    /*
     *
     */
    void Handle::release()
    {
      PRINT_LOG_MSG("......calling Handle::release with accessType = " << accessType);
      //we need to lock the fifo mutex to access its components
      std::unique_lock<std::mutex> lock_q(*accessQueue->getExternalMutex());

//      bool doPop = false;  //should we pop the queue ?

//      //we need to lock the mutex of the handle to avoid race condition on nbReaders and multiple accessQueue->pop() calls
//      std::unique_lock<std::mutex> lock_h(mut);

      if ( accessType && fusionnedHandle!=NULL ) {  //this is a read access attached to a previous read handle
//          //we need to lock the mutexes of both handles to avoid race condition on nbReaders and multiple accessQueue->pop() calls
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
                accessQueue->pop();
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
              accessQueue->pop();
              if( (this->next) != NULL ) {
                next->cond.notify_all();    //we notify_all since multiple threads can be waiting on the same condition var if next is a read access
                PRINT_LOG_MSG(" from Handle::release, next handle(s) have been notified !");
              }
        }
//        lock_h.unlock();
      }
      else if ( !accessType ) {  //this a write access, so we can pop the handle
            accessQueue->pop();
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

    /*
     *
     */
    Data<T>* Handle::acquire()
    {
      PRINT_LOG_MSG("......calling Handle::acquire with accessType = " << accessType);
      //we need to lock the fifo mutex to access its components
      std::unique_lock<std::mutex> lock_q(*accessQueue->getExternalMutex());
      //wait for this (or the read handle we "fusionned" this with) to be head of the queue
      while( !accessQueue->isHead(this) && !accessQueue->isHead(fusionnedHandle) ) {
        PRINT_LOG_MSG(" from Handle::acquire with accessType = " << accessType << " calling wait....");

        cond.wait(lock_q);

        if( accessQueue->isHead(this) ) {
            PRINT_LOG_MSG(" from Handle::acquire with accessType = " << accessType << " wait done !");
        }
      }
      //return the data
      return accessQueue->getData();
    }


    void Handle::setNext( Handle<T>* elt )
    {
      next = elt;
    }

    Handle<T>* Handle::getNext()
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




