/*!
 \file prPointVec.h
 \brief Header file for the prPoint class
 \author Guillaume CARON
 \version 0.1
 \date february 2017
 */

#if !defined(_PRPOINTVEC_H)
#define _PRPOINTVEC_H

#include <per/prGeometricalElement.h>

/*!
 \class PER_EXPORT prPointVec prPointVec.h <per/prPointVec.h>
 \brief Base class defining a generic point 
 */
class PER_EXPORT prPointVec : public prGeometricElement
{
public:
	virtual ~prPointVec() override = default;
};

#endif  //_PRPOINTVEC_H
