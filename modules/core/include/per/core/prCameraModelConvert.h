/*!
 \file prCameraModelConvert.h
 \brief Conversion class between camera models
 \author Guillaume CARON
 \version 0.1
 \date march 2022, june 2022
 */

#if !defined(_PRCAMERAMODELCONVERT_H)
#define _PRCAMERAMODELCONVERT_H

#include <per/prOmni.h>
#include <per/prFisheyeEquidistant.h>
#include <per/prPolyCart.h>

/*!
 \class PER_EXPORT prCameraModelConvert prCameraModelConvert.h <per/prCameraModelConvert.h>
 \brief Class for dealing with converting a camera model to another
 */
class PER_EXPORT prCameraModelConvert {
public:
    /*!
     * \fn static double convert(prFisheyeEquidistant &fecam, prOmni &ocam, double theoreticalFOV = 180)
     * \brief Converter of a prFisheyeEquidistant camera model to a prOmni camera model
     * TODO: generalize to au different to av (currently, only the horizontal parameters are considered to output au=av)
     * \param fecam the input fisheye equidistant camera model
     * \param ocam the output omni camera model
     * \param theoreticalFOV (optional) theoretical FOV [degrees] (as close to the reality for a better approximation). Default 180 degrees.
     * \return the average residual (digital image plane) of the conversion that is an approximation
     */
    static double convert(prFisheyeEquidistant &fecam, prOmni &ocam, double theoreticalFOV = 180);

    virtual ~prCameraModelConvert() = default;

		static double convert(prOmni &ocam, prPolyCart &pccam, double theoreticalFOV = 180);

    /*!
     * \fn static double convert_distortions(prCameraModel *incam, prCameraModel *outcam)
     * \brief Converter of distorsion parameters to a smaller number of distorsion parameters: Plumb-Bob model to Rational-polynomial model
     * Assumes outcam is already allocated and only the distorsion parameters will be changed
     * TODO: add an enum for other cases of conversion as method parameter
     * \param incam the input camera model
     * \param outcam the output camera model
     * \param theoreticalFOV (optional) theoretical FOV [degrees] (as close to the reality for a better approximation). Default 180 degrees.
     * \return the average residual (digital image plane) of the conversion that is an approximation
     */
    static double convert_distortions(prCameraModel *incam, prCameraModel *outcam, double theoreticalFOV = 90);
    
    /*!
     * \fn static double convertOptBearing(prOmni &ocam, prPolyCart &pccam, unsigned imWidth, unsigned imHeight)
     * \brief Converter of a prOmni camera model to a prPolyCart camera model
     * TODO: generalize to prCameraModel pointers
     * \param ocam the input omni camera model
     * \param fecam the output polynomial cartesian camera model
     * \param imWidth width in pixels of the digital image plane
     * \param imHeight height in pixels of the digital image plane
     * \return the average residual (digital image plane) of the conversion that is an approximation
     */
    //static double convertOptBearing(prOmni &ocam, prPolyCart &pccam, unsigned imWidth, unsigned imHeight);
    
};

#endif  //_PRCAMERAMODELCONVERT_H
