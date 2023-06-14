#include <per/prPhysicalElement.h>

prPhysicalElement& prPhysicalElement::operator=(const prPhysicalElement& pe)
{
    p = pe.p;
    
    return *this;
}

unsigned int prPhysicalElement::size()
{
  return p.getRows();
}

