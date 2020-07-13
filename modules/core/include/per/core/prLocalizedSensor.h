/*!
 \file prLocalizedSensor.h
 \brief Definition of a localized sensor = sensor model + sensor pose
 \author Guillaume CARON
 \version 0.1
 \date january 2017
 */

#if !defined(_PRLOCALIZEDSENSOR_H)
#define _PRLOCALIZEDSENSOR_H

#include <per/prSensorModel.h>
#include <visp/vpHomogeneousMatrix.h>

/*!
 \class prLocalizedSensor prLocalizedSensor.h <per/prLocalizedSensor.h>
 \brief Container class for dealing with a sensor having a pose in the world.
 */
class prLocalizedSensor
{
public:
    prSensorModel *sensor; /*!< Pointer to a sensor model */
    vpHomogeneousMatrix sMw; /*!< The world to sensor frame change */

};

#endif  //_PRLOCALIZEDSENSOR_H
