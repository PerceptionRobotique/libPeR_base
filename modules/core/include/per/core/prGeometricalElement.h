/*!
 \file prGeometricElement.h
 \brief Header file for the base prGeometricElement class
 \author Guillaume CARON
 \version 0.1
 \date february 2017
 */

#if !defined(_PRGEOMETRICELEMENT_H)
#define _PRGEOMETRICELEMENT_H

#include <per/prPhysicalElement.h>

/*!
 \class PER_EXPORT prGeometricElement prGeometricElement.h <per/prGeometricElement.h>
 \brief Base class defining a geometric element (point, line, plane, ...)
 */
class PER_EXPORT prGeometricElement : public prPhysicalElement
{
public:
	virtual ~prGeometricElement() override = default;
};

#endif  //_PRGEOMETRICELEMENT_H
