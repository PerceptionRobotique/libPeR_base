#include <per/prPhysicalElement.h>

prPhysicalElement& prPhysicalElement::operator=(const prPhysicalElement& pe)
{
    p = pe.p;
    
    return *this;
}
