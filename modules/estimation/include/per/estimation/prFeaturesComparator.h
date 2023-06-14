/*!
 \file prFeaturesComparator.h
 \brief Definition of the base class of feature sets comparator
 \author Guillaume CARON
 \version 0.1
 \date august 2017
 */

#if !defined(_PRFEATURESCOMPARATOR_H)
#define _PRFEATURESCOMPARATOR_H

#include <per/prFeaturesSet.h>

/*!
 \class PER_EXPORT prFeaturesComparator<T> prFeaturesComparator.h <per/prFeaturesComparator.h>
 \brief Base class for comparing sets of features
 */
template<typename T>
class PER_EXPORT prFeaturesComparator
{
public:
    /*!
     * \fn prFeaturesComparator()
     * \brief Default constructor of a prFeaturesComparator object
     */
    prFeaturesComparator()
    {
        
    }
    
    /*!
     * \fn ~prFeaturesComparator()
     * \brief Default Destructor of a prFeaturesComparator object
     */
    virtual ~prFeaturesComparator()
    {
        
    }
};

#endif  //_PRFEATURESCOMPARATOR_H
