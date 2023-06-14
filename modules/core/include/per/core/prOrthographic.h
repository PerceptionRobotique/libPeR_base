/*!
 \file prOrthographic.h
 \brief Definition of the orthographic camera model
 \author Guillaume CARON
 \version 0.1
 \date january 2017
 */

#if !defined(_PRORTHOGRAPHIC_H)
#define _PRORTHOGRAPHIC_H

#include <per/prCameraModel.h>

/*!
 \class PER_EXPORT prOrthographic prOrthographic.h <per/prOrthographic.h>
 \brief Class for dealing with the projection functions of an orthographic camera.
 */
class PER_EXPORT prOrthographic : public prCameraModel
{
public:
	virtual ~prOrthographic() override = default;
};

#endif  //_PRORTHOGRAPHIC_H
