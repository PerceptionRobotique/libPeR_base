/*!
 \file pr2DPointVec.h
 \brief Header file for the pr2DPointVec class
 \author Guillaume CARON
 \version 0.1
 \date february 2017
 */

#if !defined(_PR2DPOINTVEC_H)
#define _PR2DPOINTVEC_H

#include <per/prPointVec.h>

/*!
 \class PER_EXPORT pr2DPointVec pr2DPointVec.h <per/pr2DPointVec.h>
 \brief Class defining a 2D point or vector
 */
class PER_EXPORT pr2DPointVec : public prPointVec
{
public:
	virtual ~pr2DPointVec() override = default;
};

#endif  //_PR2DPOINTVEC_H
