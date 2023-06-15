/*!
 \file prGaussianMixtureSample.h
 \brief Header file for the prGaussianMixtureSample class
 \author Guillaume CARON
 \version 0.1
 \date august 2017
 */
#if !defined(_PRGAUSSIANMIXTURESAMPLE_H)
#define _PRGAUSSIANMIXTURESAMPLE_H

#include <per/prIntensity.h>

/*!
 \class PER_EXPORT prGaussianMixtureSample prGaussianMixtureSample.h <per/prGaussianMixtureSample.h>
 \brief Class defining a gaussian mixture sample from a general intensity
 */

template<typename Ts, class Tcam = prCameraModel>
class PER_EXPORT prGaussianMixtureSample : public prIntensity<Ts, Tcam>
{
public:
    /*!
     * \fn prGaussianMixtureSample(float lambda_g = 1.f)
     * \brief Contructor of the prGaussianMixtureSample generic Gaussian mixture sample feature type
     * \param lambda_g the Gaussian expansion parameter (\lambda_g)
     */
    prGaussianMixtureSample(float _lambda_g = 1.f) : lambda_g(_lambda_g)
    {
        
    }
    
    /*!
     * \fn ~prPhotometricGMS
     * \brief Destructor of a prPhotometricGMS photometric Gaussian mixture sample feature object
     */
    virtual ~prGaussianMixtureSample() override
    {
        
    }

//private:
    float lambda_g; /*!< Gaussian expansion parameter (\lambda_g)*/
};

#endif  //_PRGAUSSIANMIXTURESAMPLE_H
