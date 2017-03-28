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

    Handle::Handle( BaseData* d, bool aType )
    : handledData(d),
      accessType(aType),
      mut(),
      cond(),
      next(NULL),
      previousRead(NULL),
      isAcquired(false),
      isHead(false),
      nbReaders(0)
    {
    }

    Handle::~Handle( )
    {
    }

    /**
     * @brief if write access : pushs the handle in queue
     *        if read access : queue->tail->nbReaders++ if queue->tail is a read, else pushs the handle in queue
     */
    void Handle::request() {
      PRINT_LOG_MSG("......calling Handle::request with accessType = " << accessType);
      //we need to lock the fifo mutex to access its components
      std::unique_lock<std::mutex> lock_q(*handledData->getExternalMutex());
      // if read access and the last request in the fifo is a read, we "fusion" this with the last request in the fifo
      if ( !handledData->isEmpty() && (handledData->getTail())->accessType )
      {
          (handledData->getTail())->setNext(this);
          previousRead = handledData->getTail();
          //if the last handle in the fifo is a read and is acquired, we can perform the read immediately
          if ( previousRead->isAcquired )
          {
              isAcquired = true;
              isHead = true;
          }
          //increment the number of consecutive readers
          nbReaders = 0;
          previousRead->incrementNbReaders(); //recursive call : increments the nbReaders of all consecutive previous reads
//          //we need to lock the mutex of previousRead to avoid race condition on previousRead->nbReaders
//          std::unique_lock<std::mutex> lock_h(fusionnedHandle->mut);
//          fusionnedHandle->nbReaders++;
//          lock_h.unlock();
////          nbReaders++;  //this has no attached readers, so nbReaders is still 0
          PRINT_LOG_MSG(" from Handle::request, previousRead->nbReaders = " << previousRead->nbReaders);
      }
      else if ( !handledData->isEmpty() )
      {
          (handledData->getTail())->setNext(this);
      }
      else //fifo is empty
      {
          isHead = true;
          isAcquired = true;
      }
//      else {  // write access or  the last request in the fifo is a write
//        handledData->push( this );
//      }
      handledData->push( this );
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
        //pop this from the fifo
        handledData->pop();
        //
        if ( accessType && previousRead!=NULL ) //this is a read access attached to a previous read handle
        {
            next->setPrevious(previousRead);
            previousRead->decrementNbReaders();  //recursive call : decrements the nbReaders of all consecutive previous reads
            if ( nbReaders > 0 ) //there are following reads
            {
              previousRead->setNext(next);
            }
        }
        else if ( accessType && nbReaders > 0 ) //this is the first Handle of a set of consecutive read accesses
        {
            //sets the next Handle as Head of the fifo
            next->setHead( true );
        }
        else if( (this->next) != NULL ) //this is a solitary read or a write access
        {
            Handle* handleIterator = next;
            PRINT_LOG_MSG(" from Handle::release, next has " << next->getNbReaders() << " readers");
            for ( unsigned int i = 0; i <= next->getNbReaders(); ++i ) // we notify all consecutive reads, or just the next Handle
            {
              PRINT_LOG_MSG(" from Handle::release, notifying 1 handle...");
              handleIterator->setHead( true );
              handleIterator->cond.notify_one();
              handleIterator = handleIterator->next;
            }
            //sets the next Handle as Head of the fifo
            PRINT_LOG_MSG(" from Handle::release, next handle(s) have been notified !");
        }
        else //this was the last Handle in the fifo
        {
        }
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
        if ( isAcquired )
        {
            return this->map();
        }
        else //wait for this to be head of the queue (consecutive read accesses can be head simultaneously)
        {
            while( !isHead ) {
                PRINT_LOG_MSG(" from Handle::acquire with accessType = " << accessType << " calling wait....");
                cond.wait(lock_q);
                PRINT_LOG_MSG(" from Handle::acquire with accessType = " << accessType << " thread woke up with isHead = " << isHead);
                if( isHead ) {
                    PRINT_LOG_MSG(" from Handle::acquire with accessType = " << accessType << " wait done !");
                }
            }
            isAcquired = true;
        }
        //return the data
        return this->map();
    }

    BaseData* Handle::map()
    {
        return handledData;
    }




    void Handle::setNext( Handle* elt )
    {
      next = elt;
    }
    void Handle::setPrevious( Handle* elt )
    {
      previousRead = elt;
    }

    Handle* Handle::getNext()
    {
      return next;
    }

    bool Handle::getType()
    {
      return accessType;
    }

    void Handle::setHead( bool b )
    {
        isHead = b;
    }

    unsigned int Handle::getNbReaders()
    {
        return nbReaders;
    }

    void Handle::incrementNbReaders()
    {
        nbReaders++;
        if ( previousRead != NULL ) //recursive call
        {
            previousRead->incrementNbReaders();
        }
    }

    void Handle::decrementNbReaders()
    {
        nbReaders--;
        if ( previousRead != NULL ) //recursive call
        {
            previousRead->decrementNbReaders();
        }
    }

    bool Handle::acquired()
    {
        return isAcquired;
    }


} // namespace objectmodel

} // namespace core

} // namespace sofa


