/*!
 \file prRegularlySampledCPImage.h
 \brief Header file for the prRegularlySampledCPImage class
 \author Guillaume CARON
 \version 0.1
 \date september 2020
 */

#if !defined(_PRREGULARLYSAMPLEDCPIMAGE_H)
#define _PRREGULARLYSAMPLEDCPIMAGE_H

#include <per/prAcquisitionModel.h>
#include <per/prCartesian2DPointVec.h>

#include <per/prCameraModel.h>
#include <per/prPointFeature.h>

//#include <per/prcommon.h>

#include <visp/vpImage.h>

/*!
 \class PER_EXPORT prRegularlySampledCPImage prRegularlySampledCPImage.h <per/prRegularlySampledCPImage.h>
 \brief Class defining the image place as a grid (row major) of planar points in Cartesian coordinates lying 
 *      and regularly sampled
 */
template<typename T> class PER_EXPORT prRegularlySampledCPImage : public prAcquisitionModel<T>
{
public:
    
    /*!
     * \fn prRegularlySampledCPImage(unsigned int _height, unsigned int _width)
     * \brief Constructor of a prRegularlySampledCPImage object (CP stands for Cartesian Plane) : initializes the size of the grid, i.e. the resolution of the planar image, and initializes planar points coordinates.
     */
    prRegularlySampledCPImage(unsigned int _height = 0, unsigned int _width = 0) : height(_height), width(_width) //, up(NULL) 
    {
        initializeGrid();
    }

    /*!
     * \fn ~prRegularlySampledCPImage()
     * \brief Destructor of the prRegularlySampledCPImage object : frees memory of planar points coordinates array.
     */
    virtual ~prRegularlySampledCPImage() override
    {
        if(ge != NULL)
        {
            delete [] ge;
            ge = NULL;
        }
    }
    
    /*!
     * \fn int initializeGrid()
     * \brief Allocates and fill the array of points coordinates as well as the array of "pixels" (not set in that function since needs data)
     * \return  0 if the planar points coordinates are well initialized
     */
    int initializeGrid()
    {
        nbSamples = width*height;
        bitmap = new T[nbSamples];
        ge = new prCartesian2DPointVec[nbSamples];
        
        //compute samples coordinates
        prCartesian2DPointVec *pt_up = (prCartesian2DPointVec *)ge;
        for(unsigned int i = 0 ; i < height ; i++)
          for(unsigned int j = 0 ; j < width ; j++, pt_up++)
          {
            pt_up->setPoint(j, i);
          }
    
        return 0;
    }
    
    
    /*!
     * \fn int buildFrom(vpImage<T> & I, prSensorModel *camera = NULL, vpImage<unsigned char> *Mask = NULL)
     * \brief Builds the pixels map of the planar image from the acquired image data
     * \param I the input image of which pixels are of type T
     * \param camera the pointer to intrinsic camera parameters (not used)
     * \param Mask the pointer to a Mask image in order to ignore some pixels (if NULL, every pixel of the input image is considered)
     * \return  0 if the pixel grid is well built
     *         -1 if the pixel grid was not initialized
     */
    int buildFrom(vpImage<T> & I, prSensorModel *camera = NULL, vpImage<unsigned char> *Mask = NULL)
    {
        if(nbSamples == 0)
            return -1;
        
        if(inttyp == IMAGEPLANE_BILINEAR)
        {
            if(bitmapf != NULL)
            {
                delete [] bitmapf;
                bitmapf = NULL;
            }
            bitmapf = new float[nbSamples];
            isBitmapfSet = true;
        }
        
        T *pt_bitmap = bitmap;
        float *pt_bitmapf = bitmapf;

        prCartesian2DPointVec *pt_up = (prCartesian2DPointVec *)ge;
        unsigned int imWidth = I.getWidth(), imHeight = I.getHeight(), i, j;
        float imHeightScale = imHeight/(float)height;
        float imWidthScale = imWidth/(float)width;

        float u, v, du, dv, unmdu, unmdv;
        for(unsigned long ns = 0 ; ns < nbSamples ; ns++, pt_up++, pt_bitmap++)
        {
            //bool consider = false;
            *pt_bitmap = 0;
            if(inttyp == IMAGEPLANE_BILINEAR)
                *pt_bitmapf = 0;
            
            u = pt_up->get_x()*imWidthScale;
            v = pt_up->get_y()*imHeightScale;
                
            if( (u >= 0) && (v >= 0) && (u < (imWidth-1)) && (v < (imHeight-1))  )
            {
                switch(inttyp)
                {
                    case IMAGEPLANE_BILINEAR:
                        i = (int)v; dv = v-i; unmdv = 1.0-dv;
                        j = (int)u; du = u-j; unmdu = 1.0-du;
                        
                        if((Mask != NULL))
                        {
                            if ((*Mask)[i][j] != 0)
                            {
                                *pt_bitmapf += I[i][j]*unmdv*unmdu;
                                //consider = true;
                            }

                            if ((*Mask)[i+1][j] != 0)
                            {
                                *pt_bitmapf += I[i+1][j]*dv*unmdu;
                                //consider = true;
                            }

                            if ((*Mask)[i][j+1] != 0)
                            {
                                *pt_bitmapf += I[i][j+1]*unmdv*du;
                                //consider = true;
                            }

                            if ((*Mask)[i+1][j+1] != 0)
                            {
                                *pt_bitmapf += I[i+1][j+1]*dv*du;
                                //consider = true;
                            }
                        }
                        else
                        {
                            *pt_bitmapf = I[i][j]*unmdv*unmdu + I[i+1][j]*dv*unmdu + I[i][j+1]*unmdv*du + I[i+1][j+1]*dv*du;
                            //consider = true;
                        }
                        
                        *pt_bitmap = *pt_bitmapf;
                        
                        break;
                    case IMAGEPLANE_NEARESTNEIGH:
                    default:
                        
                        i = vpMath::round(v);
                        j = vpMath::round(u);
                        
                        if(Mask != NULL)
                        {
                            if((*Mask)[i][j] != 0)
                            {
                                *pt_bitmap = I[i][j];
                                //consider = true;
                            }
                        }
                        else
                        {
                            *pt_bitmap = I[i][j];
                            //consider = true;
                        }
                        
                        break;
                }
            }
            if(inttyp == IMAGEPLANE_BILINEAR)
                pt_bitmapf++;
        }
        
        return 0;
    }
    
    /*!
     * \fn int toImage(vpImage<T> & I_r, vpPoseVector & r, prSensorModel *camera, vpImage<unsigned char> *Mask = NULL)
     * \brief Maps the planar grid to an image (mainly for visualization purpose)
     * \param I_r (output) the output image of which pixels are of type T
     * \param r the pose vector of the frame change to apply to the grid before mapping to the image
     * \param camera the pointer to the camera intrinsic parameters (not used)
     * \param Mask the pointer to a Mask image in order to ignore some pixels (if NULL, every pixel of the input image is considered)
     * \return  0 if the output image image is well built
     *         -1 if there is no initialized planar grid
     */
    int toImage(vpImage<T> & I_r, vpPoseVector & r, prSensorModel *camera, vpImage<unsigned char> *Mask = NULL)
    {
        if(nbSamples == 0)
            return -1;
        
        T *pt_bitmap = bitmap;
        float *pt_bitmapf = bitmapf;

        prCartesian2DPointVec *pt_up = (prCartesian2DPointVec *)ge;
        unsigned int imWidth = I_r.getWidth(), imHeight = I_r.getHeight(), i, j;
        float imHeightScale = imHeight/(float)height;
        float imWidthScale = imWidth/(float)width;

        double u, v;
        
        for(unsigned long ns = 0 ; ns < nbSamples ; ns++, pt_up++, pt_bitmap++)
        {
            u = pt_up->get_x()*imWidthScale;
            v = pt_up->get_y()*imHeightScale;
                
            if( (u >= 0) && (v >= 0) && (u < (imWidth)) && (v < (imHeight))  )
            {
                i = (unsigned int)v;
                j = (unsigned int)u;
                
                /*if((Mask != NULL))
                {
                    if ((*Mask)[i][j] != 0)
                    {
                        I_r[i][j] = *pt_bitmap;
                    }
                    else
                    {
                        I_r[i][j] = *pt_bitmap;
                    }
                    
                }
                else*/
                    I_r[i][j] = *pt_bitmap;
            }
        }
        
        return 0;
    }
    

    /*!
     * \fn int getRawSample(unsigned long i, prCartesian2DPointVec & pXS, T & val)
     * \brief Extract the sample position and raw value at given index
     * \param i the index at which to extract the raw
     * \param pXS (output) the output coordinates of the sample
     * \param val (output) the output raw value at the sample coordinates
     * \return  0 if the raw sample could be got
     *         -1 if the index is out of range
     */
    int getRawSample(unsigned long i, prCartesian2DPointVec & pup, T & val)
    {
        if(i >= nbSamples)
            return -1;
        
        pup = ((prCartesian2DPointVec *)ge)[i];
        val = bitmap[i];
        
        return 0;
    }
    
    /*!
     * \fn int getSample(unsigned long i, prCartesian2DPointVec & pXS, T & val)
     * \brief Extract the sample position and value (that might have been interpolated) at given index
     * \param i the index at which to extract the raw
     * \param pXS (output) the output coordinates of the sample
     * \param val (output) the output value at the sample coordinates
     * \return  0 if the sample could be got
     *         -1 if the index is out of range
     */
    int getSample(long i, prCartesian2DPointVec & pup, float & val)
    {
        if( (i >= nbSamples) || (i < 0) )
        {
            pup = ((prCartesian2DPointVec *)ge)[0];
            val = 0.f;
            return -1;
        }
        
        if(!isBitmapfSet)
            return -2;
        
        pup = ((prCartesian2DPointVec *)ge)[i];
        val = bitmapf[i];
        
        return 0;
    }

    /*!
     * \fn int getSample(unsigned long i, T & val)
     * \brief Extract the sample position and value (that might have been interpolated) at given index
     * \param i the index at which to extract the raw
     * \param pXS (output) the output coordinates of the sample
     * \param val (output) the output value at the sample coordinates
     * \return  0 if the sample could be got
     *         -1 if the index is out of range
     */
    int getSample(long i, float & val)
    {
        if( (i >= nbSamples) || (i < 0) )
        {
            val = 0.f;
            return -1;
        }
        
        if(!isBitmapfSet)
            return -2;
        
        val = bitmapf[i];
        
        return 0;
    }
    
    /*!
     * \fn prRegularlySampledCPImage<T> &operator=(const prRegularlySampledCPImage<T> &RSCPI)
     * \brief = operator override to copy prRegularlySampledCPImage object
     * TODO: better deal with the memcpy when pixels are not of a basic type
     * TODO: check if ge[j] = RSCPI.ge[j];  needs pointer casts
     * \param RSCPI the prRegularlySampledCPImage to copy to the current
     * \return the prRegularlySampledCPImage<T> self reference object
     */
    prRegularlySampledCPImage<T> &operator=(const prRegularlySampledCPImage<T> &RSCPI)
    {
        height = RSCPI.height;
        width = RSCPI.width;

        inttyp = RSCPI.inttyp;
        nbSamples = RSCPI.nbSamples;
        isBitmapfSet = RSCPI.isBitmapfSet;
        
        bitmap = new T[nbSamples];
        std::memcpy(bitmap, RSCPI.bitmap, nbSamples*sizeof(T)); // OK si T est un type de base
        ge = new prCartesian2DPointVec[nbSamples];
        for(int j = 0 ; j < nbSamples ; j++)
            ge[j] = RSCPI.ge[j];
        if(isBitmapfSet)
        {
            bitmapf = new float[nbSamples];
            std::memcpy(bitmapf, RSCPI.bitmapf, nbSamples*sizeof(float));
        }
        
        return *this;
    }

    /*!
     * \fn unsigned int getSampleDim()
     * \brief accessor to the number of dimensions of a point sample
     * TODO: make it more generic
     * \return the number of dimensions of a point sample
     */
    unsigned int getSampleDim()
    {
      return 2;
    }
    
    int changeSampleFrame(unsigned long & ns, vpHomogeneousMatrix & dMc, prCartesian2DPointVec & x, double & rho)
		{
			x = ((prCartesian2DPointVec *)ge)[ns].changeFrame(dMc);

			rho = 1.;

			return 0;
		}
    
    /*!
     * \fn int getSampleRho(unsigned long & ns, double & rho)
     * \brief accessor to the depth of a pixel (here only for compatibility) 
     * \param ns 
     * \param rho 
     * \return 0
     */
		int getSampleRho(unsigned long & ns, double & rho)
		{
			//TODO: test ns?
			
				rho = 1.;
				
			
				return 0;
		}

        /*!
     * \fn int getSampleDepth(unsigned long & ns, double & Z)
     * \brief accessor to the depth of a pixel (here only for compatibility) 
     * \param ns 
     * \param Z 
     * \return 0
     */
		int getSampleDepth(unsigned long & ns, double & Z)
		{
			//TODO: test ns?
			/*if(isGe3DSpaceSet)
			{*/
				Z = 1.;
			/*}
			else
			{
				rho = 1.;
			}*/

			return 0;
		}
    
//private:
    unsigned int height; /*!< height of the planar image grid */
    unsigned int width; /*!< width of the planar image grid */

    //prCartesian2DPointVec *up; /*!< array of Cartesian planar points coordinates */

    using prAcquisitionModel<T>::ge;
    using prAcquisitionModel<T>::nbSamples;
    using prAcquisitionModel<T>::bitmap;
    using prAcquisitionModel<T>::bitmapf;
    using prAcquisitionModel<T>::isBitmapfSet;
    using prAcquisitionModel<T>::inttyp;
};

#endif  //_PRREGULARLYSAMPLEDCPIMAGE_H
