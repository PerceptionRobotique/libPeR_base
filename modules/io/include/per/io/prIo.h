/*!
 \file prIo.h
 \brief Header file for the prIo class
 \author Guillaume CARON
 \version 0.1
 \date january 2017
 */

#if !defined(_PRIO_H)
#define _PRIO_H

#include <per/prcommon.h>

/*!
 \class PER_EXPORT prIo prIo.h <per/prIo.h>
 \brief Base class for writing to files and reading from files.
 */
class PER_EXPORT prIo {
public:
	virtual ~prIo() = default;
};

#endif  //_PRIO_H
