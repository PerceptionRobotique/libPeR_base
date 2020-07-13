/*!
 \file prPhyscialElement.h
 \brief Header file for the base prPhysicalElement class
 \author Guillaume CARON
 \version 0.1
 \date february 2017
 */

#if !defined(_PRPHYSICALELEMENT_H)
#define _PRPHYSICALELEMENT_H

#include <visp/vpColVector.h>

/*!
 \class prPhysicalElement prPhysicalElement.h <per/prPhysicalElement.h>
 \brief Base class for physical elements (geometric, kinematic, dynamic, spectral, and so on)
 */
class prPhysicalElement
{
public:
    prPhysicalElement() {}
    ~prPhysicalElement() {}
    
    prPhysicalElement& operator=(const prPhysicalElement& pe);
    
protected:
    vpColVector p; /*!< parameters vector of the physical element */

};

#endif  //_PRPHYSICALELEMENT_H
