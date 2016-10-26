#include "elementcache.hpp" 


ElementCacheClass::ElementCacheClass(long ele): HashClassData(), parameter(0)
{
  diffusion=-1.0;
}
 

ElementCacheClass::~ElementCacheClass()
{
  delete [] parameter;
}
