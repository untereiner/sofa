/*
 * ThreadSafeQueue.cpp
 *
 *  Created on: 22 mars 2017
 *      Author: maxime
 */

#include <sofa/core/objectmodel/ThreadSafeQueue.h>

namespace sofa
{

namespace core
{

namespace objectmodel
{
#ifndef LOG_THREADS
ThreadSafeQueue::ThreadSafeQueue( ) : externalMut(), internalMut(), cond(), head(NULL), tail(NULL)
//    ThreadSafeQueue( Data* d ) : handledData(d), externalMut(), internalMut(), cond(), head(NULL), tail(NULL)
{
}
//    ThreadSafeQueue( DDGNode* d );
#else //LOG_THREADS
ThreadSafeQueue::ThreadSafeQueue( std::string s ) : externalMut(), internalMut(), cond(), name(s), head(NULL), tail(NULL)
{
//    ThreadSafeQueue( DDGNode* d, std::string s );
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
Handle* ThreadSafeQueue::pop()
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

//    Data* getData()
//    {
//        return handledData;
//    }
//    BaseData* getData();
std::mutex* ThreadSafeQueue::getExternalMutex()
{
  return &externalMut;
}
//the following methods are not thread-safe
//queue->mut should be locked by the user before calling getTail, getHead, isHead or isempty
Handle* ThreadSafeQueue::getTail()
{
  return tail;
}
Handle* ThreadSafeQueue::getHead()
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



