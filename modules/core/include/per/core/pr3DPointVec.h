/*!
 \file pr3DPointVec.h
 \brief Header file for the pr3DPointVec class
 \author Guillaume CARON
 \version 0.1
 \date february 2017
 */

#if !defined(_PR3DPOINTVEC_H)
#define _PR3DPOINTVEC_H

#include <per/prPointVec.h>

/*!
 \class PER_EXPORT pr3DPointVec pr3DPointVec.h <per/pr3DPointVec.h>
 \brief Class defining a feature point existing in 2D pixels, normalized 2D, camera frame 3D, world frame 3D
 */
class PER_EXPORT pr3DPointVec : public prPointVec
{
public:
	virtual ~pr3DPointVec() override = default;
};

#endif  //_PR3DPOINTVEC_H
