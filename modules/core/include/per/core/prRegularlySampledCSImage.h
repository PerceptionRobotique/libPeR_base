/*!
 \file prRegularlySampledCSImage.h
 \brief Header file for the prRegularlySampledCSImage class
 \author Guillaume CARON
 \version 0.1
 \date april 2017
 */

#if !defined(_PRREGULARLYSAMPLEDCSIMAGE_H)
#define _PRREGULARLYSAMPLEDCSIMAGE_H

#include <per/prAcquisitionModel.h>
#include <per/prCartesian3DPointVec.h>

#include <per/prStereoModel.h>
#include <per/prOmni.h>
#include <per/prEquirectangular.h>
#include <per/prPointFeature.h>

#include <per/prRegularSphereData.h>

//#include <per/prcommon.h>

#include <visp/vpImage.h>

/*!
 \class PER_EXPORT prRegularlySampledCSImage prRegularlySampledCSImage.h <per/prRegularlySampledCSImage.h>
 \brief Class defining the image place as a set of spherical points in Cartesian coordinates lying 
 *      on a unit sphere and regularly sampled, in the geodesic meaning
 */
template<typename T> class PER_EXPORT prRegularlySampledCSImage : public prAcquisitionModel<T>
{
public:
    
    /*!
     * \fn prRegularlySampledCSImage(unsigned int _subdivLevels = 3)
     * \brief Constructor of a prRegularlySampledCSImage object (CS stands for Cartesian Spherical) : initializes the number of subdivision level of the sphere, i.e. the resolution of the spherical image, and loads the spherical points.
     */
    prRegularlySampledCSImage(unsigned int _subdivLevels = 3) : subdivLevels(_subdivLevels)//, XS(NULL) //, nbSamples(0), bitmap(NULL), bitmapf(NULL), isBitmapfSet(false)
    {
        loadSphere();
    }

    /*!
     * \fn ~prRegularlySampledCSImage()
     * \brief Destructor of the prRegularlySampledCSImage object : frees memory of grayscale intensities (raw and interpolated) and spherical points coordinates array.
     */
    virtual ~prRegularlySampledCSImage() override
    {
        /*
        //already done by the parent
        if(bitmap != NULL)
            delete [] bitmap;
        
        if(bitmapf != NULL)
            delete [] bitmapf;
        */
        /*
        //not to be done as kept in the array of prRegularSphereData
        if(XS != NULL)
            delete [] XS;*/
    }
    
    /*!
     * \fn int loadSphere()
     * \brief Allocates and fill the array of points coordinates as well as the array of "pixels" (not set in that function since needs data)
     * \return  0 if the spherical points coordinates are well loaded
     *         -1 if the number of subdivision level is too high
     */
    int loadSphere()
    {
        if(subdivLevels > prData::maxSubdivLevel)
            return -1;
        
        nbSamples = prData::nbSphericalPts[subdivLevels];
        bitmap = new T[nbSamples];
        ge = new prCartesian3DPointVec[nbSamples];
        
        //copy samples coordinates
        prCartesian3DPointVec *pt_XS = (prCartesian3DPointVec *)ge;
        const prCartesian3DPointVec *pt_XS_src = prData::samplesCoordinates[subdivLevels];
        float inorme = 1.0f;
        for(unsigned long i = 0 ; i < nbSamples ; i++, pt_XS++, pt_XS_src++)
        {
            //renormalize coordinates (might be necessary due to the text exportation
            inorme = 1.0f/sqrt(pt_XS_src->get_X()*pt_XS_src->get_X() + pt_XS_src->get_Y()*pt_XS_src->get_Y() + pt_XS_src->get_Z()*pt_XS_src->get_Z());
            pt_XS->set(pt_XS_src->get_X()*inorme, pt_XS_src->get_Y()*inorme, pt_XS_src->get_Z()*inorme, pt_XS_src->get_W());
        }
    
        return 0;
    }
    
    
    /*!
     * \fn int buildFromTwinOmni(vpImage<T> & I, prStereoModel & stereoCam, vpImage<unsigned char> *Mask = NULL)
     * \brief Builds the pixels map of the spherical image from the acquired dual omni image data
     * TODO: to generalize for other sensor types than the dual omni
     * \param I the input image of which pixels are of type T
     * \param stereoCam the dual omni camera intrinsic and extrinsic parameters
     * \param Mask the pointer to a Mask image in order to ignore some pixels (if NULL, every pixel of the input image is considered)
     * \return  0 if the spherical image is well built
     *         -1 if there is no spherical image points loaded
     */
    int buildFromTwinOmni(vpImage<T> & I, prStereoModel & stereoCam, vpImage<unsigned char> *Mask = NULL)
    {
        I_req = I;
        MaskS.resize(Mask->getHeight(), Mask->getWidth());
        std::memcpy(MaskS.bitmap, Mask->bitmap, Mask->getHeight() * Mask->getWidth() * sizeof(unsigned char));

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
        vpHomogeneousMatrix c2Rc1 = stereoCam.sjMr[1];
        prCartesian3DPointVec XSs;
        prCartesian3DPointVec *pt_XS = (prCartesian3DPointVec *)ge;
        unsigned int icam = 0, imWidth = I.getWidth(), imHeight = I.getHeight(), i, j;
        prPointFeature P;
        double u, v, du, dv, unmdu, unmdv;

        for(unsigned long ns = 0 ; ns < nbSamples ; ns++, pt_XS++, pt_bitmap++)
        {
            //bool consider = false;
            *pt_bitmap = 0;
            if(inttyp == IMAGEPLANE_BILINEAR)
                *pt_bitmapf = 0;
            
            P.sX = *pt_XS;
            //All that follows should be put in a function
            if(P.get_Z() > 0)
            {
                icam = 0;
            }
            else
            {
                icam = 1;
                P.sX = P.sX.changeFrame(c2Rc1);
            }
            
            //if(P.get_Z() > 0.0)
            {
                ((prOmni *)(stereoCam.sen[icam]))->project3DImage(P);
                
                ((prOmni *)(stereoCam.sen[icam]))->meterPixelConversion(P);
            
                u = P.get_u();
                v = P.get_v();
                
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
            }
            if(inttyp == IMAGEPLANE_BILINEAR)
                pt_bitmapf++;
        }
        
        return 0;
    }
    
    /*!
     * \fn int toTwinOmni(vpImage<T> & I_r, vpPoseVector & r, prStereoModel & stereoCam, vpImage<unsigned char> *Mask = NULL)
     * \brief Projects the spherical image to a dual omni image (mainly for visualization purpose)
     * TODO: to generalize for other sensor types than the dual omni
     * \param I_r (output) the output image of which pixels are of type T
     * \param r the pose vector of the frame change to apply to the spherical image before projection in the dual omni image plane
     * \param stereoCam the dual omni camera intrinsic and extrinsic parameters
     * \param Mask the pointer to a Mask image in order to ignore some pixels (if NULL, every pixel of the input image is considered ; if not NULL, masked pixels are not changed)
     * \return  0 if the spherical image is well built
     *         -1 if there is no spherical image points loaded
     */
    int toTwinOmni(vpImage<T> & I_r, vpPoseVector & r, prStereoModel & stereoCam, vpImage<unsigned char> *Mask = NULL)
    {
        if(nbSamples == 0)
            return -1;
        
        T *pt_bitmap = bitmap;
        float *pt_bitmapf = bitmapf;
        vpHomogeneousMatrix c2Rc1 = stereoCam.sjMr[1];
        prCartesian3DPointVec XSs;
        prCartesian3DPointVec *pt_XS = (prCartesian3DPointVec *)ge;
        unsigned int icam = 0, imWidth = I_r.getWidth(), imHeight = I_r.getHeight(), i, j;
        prPointFeature P;
        double u, v, du, dv, unmdu, unmdv;
        
        vpHomogeneousMatrix dMc;
        dMc.buildFrom(r);
        dMc = dMc.inverse(); //a justifier clairement dans la doc
        
        for(unsigned long ns = 0 ; ns < nbSamples ; ns++, pt_XS++, pt_bitmap++)
        {
            P.sX = *pt_XS;
            P.sX = P.sX.changeFrame(dMc);
            if(P.get_Z() > 0)
            {
                icam = 0;
            }
            else
            {
                icam = 1;
                P.sX = P.sX.changeFrame(c2Rc1);
            }
            
            if(P.get_Z() > 0.0)
            {
                ((prOmni *)(stereoCam.sen[icam]))->project3DImage(P);
                
                ((prOmni *)(stereoCam.sen[icam]))->meterPixelConversion(P);
                
                u = P.get_u();
                v = P.get_v();
                
                if( (u >= 0) && (v >= 0) && (u < (imWidth-1)) && (v < (imHeight-1))  )
                {
                    i = (unsigned int)v;
                    j = (unsigned int)u;
                    
				            if((Mask != NULL))
				            {
				                if ((*Mask)[i][j] != 0)
				                   I_r[i][j] = *pt_bitmap;
				            }
				            else
				            {
				                I_r[i][j] = *pt_bitmap;
				            }
                }
            }
        }
        
        return 0;
    }

    /*!
     * \fn int buildFromEquiRect(vpImage<T> & I, prEquirectangular & equiCam, vpImage<unsigned char> *Mask = NULL)
     * \brief Builds the pixels map of the spherical image from the acquired equirectangular image data
     * TODO: to generalize for other sensor types than the equirectangular by merging with buildFromTwinOmni
     * \param I the input image of which pixels are of type T
     * \param equiCam the dual omni camera intrinsic and extrinsic parameters
     * \param Mask the pointer to a Mask image in order to ignore some pixels (if NULL, every pixel of the input image is considered)
     * \return  0 if the spherical image is well built
     *         -1 if there is no spherical image points loaded
     */
    int buildFromEquiRect(vpImage<T> & I, prEquirectangular & equiCam, vpImage<unsigned char> *Mask = NULL)
    {
        I_req = I;
        MaskS.resize(Mask->getHeight(), Mask->getWidth());
        std::memcpy(MaskS.bitmap, Mask->bitmap, Mask->getHeight() * Mask->getWidth() * sizeof(unsigned char));

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

        prCartesian3DPointVec XSs;
        prCartesian3DPointVec *pt_XS = (prCartesian3DPointVec *)ge;
        unsigned int icam = 0, imWidth = I.getWidth(), imHeight = I.getHeight(), i, j;
        prPointFeature P;
        double u, v, du, dv, unmdu, unmdv;

        for(unsigned long ns = 0 ; ns < nbSamples ; ns++, pt_XS++, pt_bitmap++)
        {
            //bool consider = false;
            *pt_bitmap = 0;
            if(inttyp == IMAGEPLANE_BILINEAR)
                *pt_bitmapf = 0;
            
            P.sX = *pt_XS;
            
            equiCam.project3DImage(P);
            
            equiCam.meterPixelConversion(P);
            
            u = P.get_u();
            v = P.get_v();
                
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
     * \fn int toEquiRect(vpImage<T> & I_r, vpPoseVector & r, prEquirectangular & equiCam, vpImage<unsigned char> *Mask = NULL)
     * \brief Projects the spherical image to an equirectangular image (mainly for visualization purpose)
     * TODO: to generalize for other sensor types than the equirectangular by merging with toTwinOmni
     * \param I_r (output) the output image of which pixels are of type T
     * \param r the pose vector of the frame change to apply to the spherical image before projection in the dual omni image plane
     * \param equiCam the equirectangular camera intrinsic and extrinsic parameters
     * \param Mask the pointer to a Mask image in order to ignore some pixels (if NULL, every pixel of the input image is considered ; if not NULL, masked pixels are not changed)
     * \return  0 if the spherical image is well built
     *         -1 if there is no spherical image points loaded
     */
    int toEquiRect(vpImage<T> & I_r, vpPoseVector & r, prEquirectangular & equiCam, vpImage<unsigned char> *Mask = NULL)
    {
        if(nbSamples == 0)
            return -1;
        
        T *pt_bitmap = bitmap;
        float *pt_bitmapf = bitmapf;

        prCartesian3DPointVec XSs;
        prCartesian3DPointVec *pt_XS = (prCartesian3DPointVec *)ge;
        unsigned int icam = 0, imWidth = I_r.getWidth(), imHeight = I_r.getHeight(), i, j;
        prPointFeature P;
        double u, v, du, dv, unmdu, unmdv;
        
        vpHomogeneousMatrix dMc;
        dMc.buildFrom(r);
        dMc = dMc.inverse(); //a justifier clairement dans la doc
        
        for(unsigned long ns = 0 ; ns < nbSamples ; ns++, pt_XS++, pt_bitmap++)
        {
            P.sX = *pt_XS;
            P.sX = P.sX.changeFrame(dMc);

            equiCam.project3DImage(P);
            
            equiCam.meterPixelConversion(P);
                
            u = P.get_u();
            v = P.get_v();
            
            if( (u >= 0) && (v >= 0) && (u < (imWidth-1)) && (v < (imHeight-1))  )
            {
                i = (unsigned int)v;
                j = (unsigned int)u;
                
               if((Mask != NULL))
                {
                    if ((*Mask)[i][j] != 0)
                       I_r[i][j] = *pt_bitmap;
                }
                else
                {
                    I_r[i][j] = *pt_bitmap;
                }
            }
        }
        
        return 0;
    }

    

    /*!
     * \fn int getRawSample(unsigned long i, prCartesian3DPointVec & pXS, T & val)
     * \brief Extract the sample position and raw value at given index
     * \param i the index at which to extract the raw
     * \param pXS (output) the output coordinates of the sample
     * \param val (output) the output raw value at the sample coordinates
     * \return  0 if the raw sample could be got
     *         -1 if the index is out of range
     */
    int getRawSample(unsigned long i, prCartesian3DPointVec & pXS, T & val)
    {
        if(i >= nbSamples)
            return -1;
        
        pXS = ((prCartesian3DPointVec *)ge)[i];
        val = bitmap[i];
        
        return 0;
    }
    
    /*!
     * \fn int getSample(unsigned long i, prCartesian3DPointVec & pXS, T & val)
     * \brief Extract the sample position and value (that might have been interpolated) at given index
     * \param i the index at which to extract the raw
     * \param pXS (output) the output coordinates of the sample
     * \param val (output) the output value at the sample coordinates
     * \return  0 if the sample could be got
     *         -1 if the index is out of range
     */
    int getSample(unsigned long i, prCartesian3DPointVec & pXS, float & val)
    {
        if(i >= nbSamples)
            return -1;
        
        if(!isBitmapfSet)
            return -2;
        
        pXS = ((prCartesian3DPointVec *)ge)[i];
        val = bitmapf[i];
        
        return 0;
    }
    
    /*!
     * \fn prRegularlySampledCSImage<T> &operator=(const prRegularlySampledCSImage<T> &RSCSI)
     * \brief = operator override to copy prRegularlySampledCSImage object
     * TODO: better deal with the memcpy when pixels are not of a basic type
     * TODO: check if ge[j] = RSCSI.ge[j];  needs pointer casts
     * \param RSCSI the prRegularlySampledCSImage to copy to the current
     * \return the prRegularlySampledCSImage<T> self reference object
     */
    prRegularlySampledCSImage<T> &operator=(const prRegularlySampledCSImage<T> &RSCSI)
    {
        subdivLevels = RSCSI.subdivLevels;
        inttyp = RSCSI.inttyp;
        nbSamples = RSCSI.nbSamples;
        isBitmapfSet = RSCSI.isBitmapfSet;
        
        bitmap = new T[nbSamples];
        std::memcpy(bitmap, RSCSI.bitmap, nbSamples*sizeof(T)); // OK si T est un type de base
        ge = new prCartesian3DPointVec[nbSamples];
        for(int j = 0 ; j < nbSamples ; j++)
            ge[j] = RSCSI.ge[j]; 
        if(isBitmapfSet)
        {
            bitmapf = new float[nbSamples];
            std::memcpy(bitmapf, RSCSI.bitmapf, nbSamples*sizeof(float));
        }
        I_req = RSCSI.I_req;

        MaskS.resize(RSCSI.MaskS.getHeight(), RSCSI.MaskS.getWidth());
        std::memcpy(MaskS.bitmap, RSCSI.MaskS.bitmap, RSCSI.MaskS.getHeight() * RSCSI.MaskS.getWidth() * sizeof(unsigned char));
        
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
      return 3;
    }

		int changeSampleFrame(unsigned long & ns, vpHomogeneousMatrix & dMc, prCartesian3DPointVec & XS, double & rho)
		{
			XS = ((prCartesian3DPointVec *)ge)[ns].changeFrame(dMc);

			rho = 1.;

			return 0;
		}
		
		 /*!
     * \fn int getSampleRho(unsigned long & ns, double & rho)
     * \brief accessor to the depth of a spherical point 
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

    vpImage<T> getSampledImage()
    {
        return I_req;
    }
    
//private:
    unsigned int subdivLevels; /*!< number of subdivision levels of the sphere */
    vpImage<T> I_req;
    vpImage<unsigned char> MaskS;

    //prCartesian3DPointVec *XS; /*!< array of Cartesian spherical points coordinates */

    using prAcquisitionModel<T>::ge;
    using prAcquisitionModel<T>::nbSamples;
    using prAcquisitionModel<T>::bitmap;
    using prAcquisitionModel<T>::bitmapf;
    using prAcquisitionModel<T>::isBitmapfSet;
    using prAcquisitionModel<T>::inttyp;
};

#endif  //_PRREGULARLYSAMPLEDCSIMAGE_H
