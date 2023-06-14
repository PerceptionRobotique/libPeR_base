/*!
 \file prAcquisitionModel.h
 \brief Header file for the base prAcquisitionModel class
 \author Guillaume CARON
 \version 0.2
 \date September 2020
 */

#if !defined(_PRACQUISITIONMODEL_H)
#define _PRACQUISITIONMODEL_H

#include <per/prcommon.h>

#include <per/prGeometricalElement.h>

/*!
 \class PER_EXPORT prAcquisitionModel prAcquisitionModel.h <per/prAcquisitionModel.h>
 \brief Class defining the base class of a data acquisition model
 */
template<typename T> class PER_EXPORT prAcquisitionModel
{
public:

    /*!
     * \fn prAcquisitionModel()
     * \brief Constructor of a base prAcquisitionModel object.
     */
    prAcquisitionModel() : nbSamples(0), ge(NULL), bitmap(NULL), bitmapf(NULL), isBitmapfSet(false)
    {
        
    }

    /*!
     * \fn ~prAcquisitionModel()
     * \brief Destructor of the prAcquisitionModel object : frees memory of grayscale intensities (raw and interpolated). Geometric samples of the signal are not freed. Such free must be handled by daughter classes.
     */
    virtual ~prAcquisitionModel()
    {
        if(bitmap != NULL)
        {
            delete [] bitmap;
            bitmap = NULL;
        }

        if(ge != NULL)
        {
            delete [] ge;
            ge = NULL;
        }
        
        if(bitmapf != NULL)
        {
            delete [] bitmapf;
            bitmapf = NULL;
        }
    }
    
    /*!
     * \fn void setInterpType(prInterpType inttyp_)
     * \brief Initialize the type of interpolation to consider for pixel intensities
     * \param inttyp_ the interpolation type
     * \return Nothing
     */
    void setInterpType(prInterpType inttyp_)
    {
        inttyp = inttyp_;
    }

    /*!
     * \fn unsigned long getNbSamples()
     * \brief Accessor to the number of spherical samples in the spherical image
     * \return the number of samples
     */
    unsigned long getNbSamples()
    {
        return nbSamples;
    }

    /*!
     * \fn int toAbsZN()
     * \brief Transform the pixel intensities so that they are zero-mean, normalized, absoluted and normalized again so that their sum is equal to one
     * \return  0 if the spherical image is well built
     *         -1 if there is no spherical image points loaded
     */
    int toAbsZN()
    {
        if(nbSamples == 0)
            return -1;

        unsigned long ns;
        double mean;
        float *pt_bitmapf;
        if(isBitmapfSet == false)
        {
            if(bitmapf != NULL)
            {
                delete [] bitmapf;
                bitmapf = NULL;
            }
            bitmapf = new float[nbSamples];
            isBitmapfSet = true;
            
            T *pt_bitmap = bitmap;
            mean = 0.0;
            //mean intensity computation
            for(ns = 0 ; ns < nbSamples ; ns++, pt_bitmap++)
                mean += *pt_bitmap;
            mean /= (double)nbSamples;
            
            //intensities centering
            pt_bitmap = bitmap;
            pt_bitmapf = bitmapf;
            for(ns = 0 ; ns < nbSamples ; ns++, pt_bitmap++, pt_bitmapf++)
                *pt_bitmapf = (double)(*pt_bitmap) - mean;
        }
        else
        {
            mean = 0.0;
            pt_bitmapf = bitmapf;
            //mean intensity computation
            for(ns = 0 ; ns < nbSamples ; ns++, pt_bitmapf++)
                mean += *pt_bitmapf;
            mean /= (double)nbSamples;
            
            //centering intensities
            pt_bitmapf = bitmapf;
            for(ns = 0 ; ns < nbSamples ; ns++, pt_bitmapf++)
                *pt_bitmapf = (double)(*pt_bitmapf) - mean;
        }
        
        //computation of the mean of absolute values of centered intensities
        mean = 0.0;
        pt_bitmapf = bitmapf;
        for(ns = 0 ; ns < nbSamples ; ns++, pt_bitmapf++)
        {
            if((*pt_bitmapf) < 0)
                *pt_bitmapf = -(*pt_bitmapf);

            mean += *pt_bitmapf;
        }
        mean /= (double)nbSamples;
        
        //normalization of intensities by the mean
        pt_bitmapf = bitmapf;
        double sum = 0.0;
        for(ns = 0 ; ns < nbSamples ; ns++, pt_bitmapf++)
        {
            *pt_bitmapf /= mean;
            sum += *pt_bitmapf;
        }
        
        //normalization of intensities by their sum
        pt_bitmapf = bitmapf;
        for(ns = 0 ; ns < nbSamples ; ns++, pt_bitmapf++)
            *pt_bitmapf /= sum;
        
        return 0;
    }

//protected:
    prInterpType inttyp;  /*!< interpolation type */
    
    prGeometricElement *ge; /*!< array of signal sampling elements (usually points) */

    unsigned long nbSamples;  /*!< number of sample points on the sphere */
    T *bitmap;  /*!< linear array of pixels of the spherical image */
    float *bitmapf; /*!< linear array of interpolated pixels of the sperical image */
    bool isBitmapfSet; /*!< boolean to know if the bitmap of interpolated pixels is allocated and filled */
};

#endif  //_PRACQUISITIONMODEL_H
