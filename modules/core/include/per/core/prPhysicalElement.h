/*!
 \file prPhyscialElement.h
 \brief Header file for the base prPhysicalElement class
 \author Guillaume CARON
 \version 0.2
 \date september 2020
 */

#if !defined(_PRPHYSICALELEMENT_H)
#define _PRPHYSICALELEMENT_H

#include <per/prcommon.h>
#include <visp/vpColVector.h>

/*!
 \class PER_EXPORT prPhysicalElement prPhysicalElement.h <per/prPhysicalElement.h>
 \brief Base class for physical elements (geometric, kinematic, dynamic, spectral, and so on)
 */
class PER_EXPORT prPhysicalElement
{
public:
    prPhysicalElement() {}
    virtual ~prPhysicalElement() = default;
    
    prPhysicalElement& operator=(const prPhysicalElement& pe);

    unsigned int size();
    
protected:
    vpColVector p; /*!< parameters vector of the physical element */

};

#endif  //_PRPHYSICALELEMENT_H
