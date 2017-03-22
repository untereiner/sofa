#include <sofa/core/objectmodel/ThreadSafeQueue.h>
//#include <sofa/core/objectmodel/Handle.h>

namespace sofa
{

namespace core
{

namespace objectmodel
{


#ifndef LOG_THREADS
    ThreadSafeQueue::ThreadSafeQueue( Data<T>* d ) : handledData(d), externalMut(), internalMut(), cond(), head(NULL), tail(NULL)
    {
    }
#else
    ThreadSafeQueue::ThreadSafeQueue( Data<T>* d, std::string s ) : handledData(d), externalMut(), internalMut(), cond(), name(s), head(NULL), tail(NULL)
    {
    }
#endif //LOG_THREADS
    ThreadSafeQueue::~ThreadSafeQueue( )
    {
    }

    void ThreadSafeQueue::push(Handle* elt)
    {
      std::unique_lock<std::mutex> lock(internalMut);
      if ( !isempty() ) {
        tail->setNext(elt);
      }
      else {
        head = elt;
      }
      tail = elt;
      lock.unlock();
      cond.notify_one();  //notifies a waiting pop (might happen is the queue was empty)

      PRINT_LOG_MSG(" pushed one " << elt->getType() << " Handle to queue " << name);
    }

    Handle<T>* ThreadSafeQueue::pop()
    {
        std::unique_lock<std::mutex> lock(internalMut);
        while( isempty() )
        {
            cond.wait(lock);
        }

        Handle* elt = head;
        head = head->getNext();
        if( head == NULL ) {
          tail = NULL;
          PRINT_LOG_MSG(" popped the last elt of queue " << name);
        }
        else {
          PRINT_LOG_MSG(" popped one elt of queue " << name);
        }
        return elt;
    }

    Data<T>* ThreadSafeQueue::getData()
    {
        return handledData;
    }

    std::mutex* ThreadSafeQueue::getExternalMutex()
    {
      return &externalMut;
    }

    Handle<T>* ThreadSafeQueue::getTail()
    {
      return tail;
    }
    Handle<T>* ThreadSafeQueue::getHead()
    {
      return head;
    }
    bool ThreadSafeQueue::isHead(Handle* elt)
    {
      return (elt ==  head);
    }

    bool ThreadSafeQueue::isempty()
    {
      return tail == NULL;
    }

} //namespace sofa
} //namespace core
} //namespace objectmodel


