/*!
 \file prFeature.h
 \brief Definition of the base class of features
 \author Guillaume CARON
 \version 0.1
 \date january 2017
 */

#if !defined(_PRFEATURE_H)
#define _PRFEATURE_H

#include <per/prcommon.h>

/*!
 \class PER_EXPORT prFeature prFeature.h <per/prFeature.h>
 \brief Base class for a computer sensing feature 
*/
class PER_EXPORT prFeature
{
public:
    /*!
     * \fn prFeature()
     * \brief Default constructor of a prFeature object
     */
    prFeature();

    /*!
     * \fn ~prFeature()
     * \brief Default destructor of a prFeature object
     */
    virtual ~prFeature();
    
    /*!
     * \fn virtual double *toDouble(unsigned int place = 0) const = 0
     * \brief conversion of the feature to an array of double
     *
     * \param place the place of defintion of the feature
     * \return the array (pointer kept for freeing memory) of double values
     */
    virtual double *toDouble(unsigned int place = 0) = 0;
    
protected:
    double *doubleFeature;
};

#endif  //_PRFEATURE_H
