/*!
 \file prSensorModel.h
 \brief Header file for the base class of sensor models
 \author Guillaume CARON
 \version 0.1
 \date january 2017
 */

#if !defined(_PRSENSORMODEL_H)
#define _PRSENSORMODEL_H

#include <per/prcommon.h>

/*!
 \class PER_EXPORT prSensorModel prSensorModel.h <per/prSensorModel.h>
 \brief Base class for dealing with geometric functions of a sensor.
 */
class PER_EXPORT prSensorModel
{
public:
	virtual ~prSensorModel() = default;
};

#endif  //_PRSENSORMODEL_H
