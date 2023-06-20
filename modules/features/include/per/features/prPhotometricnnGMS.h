/*!
 \file prPhotometricnnGMS.h
 \brief Header file for the prPhotometricnnGMS class
 \author Guillaume CARON
 \version 0.1
 \date september 2020
 */

#if !defined(_PRPHOTOMETRICnnGMS_H)
#define _PRPHOTOMETRICnnGMS_H

#include <per/prGaussianMixtureSample.h>
//#include <per/prRegularlySampledCPImage.h>

#define PIPIP3S2 15.7496099457

/*!
 \class PER_EXPORT prPhotometricnnGMS prPhotometricnnGMS.h <per/prPhotometricnnGMS.h>
 \brief Class defining a photometric based non-normalized gaussian mixture sample
 */
template<typename Ts>
class PER_EXPORT prPhotometricnnGMS : public prGaussianMixtureSample<Ts> {
public:
    /*!
     * \fn prPhotometricnnGMS<Ts>(float lambda_g = 1.f)
     * \brief Contructor of the prPhotometricnnGMS<Ts> photometric non-noramlized Gaussian mixture sample feature type
     * \param lambda_g the non-noramlized Gaussian expansion parameter (\lambda_g)
     */
    prPhotometricnnGMS<Ts>(float lambda_g = 1.f): /*prGaussianMixtureSample<Ts>(lambda_g),*/ gms(0.f), L(6), featureJacobianComputed(false), E_precomp(NULL), Eug_precomp_interlaced(NULL), nbSamples(0), nbDim(0), iSampleStart(0), iSampleJump(0)
    {
        setLambda(lambda_g);
        halfSide = (unsigned int)(3.f*lambda_g+0.5f);
        iSampleBlock = 2*halfSide+1;
        iSampleEnd = 2*halfSide+1;
    }
    
    /*!
     * \fn prPhotometricnnGMS<Ts>(prPhotometricnnGMS<Ts> const &Src)
     * \brief Copy contructor of a prPhotometricnnGMS<Ts> photometric non-normalized Gaussian mixture sample feature type
     * \param Src source object to copy
     */
    prPhotometricnnGMS<Ts>(prPhotometricnnGMS<Ts> const &Src): gms(Src.getGMS()), featureJacobianComputed(Src.featureJacobianComputed), halfSide(Src.halfSide), iSampleStart(Src.iSampleStart), iSampleBlock(Src.iSampleBlock), iSampleJump(Src.iSampleJump), iSampleEnd(Src.iSampleEnd)
    {
        setLambda(Src.getLambda());
        s_g = Src.getS_g();
        nbDim = Src.nbDim;
        L = Src.L;

        if(Src.E_precomp == NULL)
        {
          E_precomp = NULL;
          Eug_precomp_interlaced = NULL;
          nbSamples = 0;
        }
        else
        {
          nbSamples = Src.nbSamples;

          E_precomp = new double[nbSamples];
          Eug_precomp_interlaced = new double[nbDim*nbSamples];

          memcpy(E_precomp, Src.E_precomp, nbSamples*sizeof(double));

          memcpy(Eug_precomp_interlaced, Src.Eug_precomp_interlaced, nbDim*nbSamples*sizeof(double));
        }
    }
    
    /*!
     * \fn ~prPhotometricnnGMS
     * \brief Destructor of a prPhotometricnnGMS photometric non-normalized Gaussian mixture sample feature object
     */
    ~prPhotometricnnGMS()
    {
        if(E_precomp != NULL)
        {
          delete [] E_precomp;
          E_precomp = NULL;
          delete [] Eug_precomp_interlaced;
          Eug_precomp_interlaced = NULL;
        }
    }
    
    /*!
     * \fn void setLambda(float _lambda_g)
     * \brief Sets the Gaussian spread parameter lambda_g and computes "constants" used from it
     * \param _lambda_g the Gaussian spread parameter
     * \return Nothing
     */
    void setLambda(float _lambda_g)
    {
        if(_lambda_g > 0.f)
        {
            prGaussianMixtureSample<Ts>::lambda_g = _lambda_g;
            oneOl2 = 1.f/(prGaussianMixtureSample<Ts>::lambda_g*prGaussianMixtureSample<Ts>::lambda_g);
            mOneO2l2 = -oneOl2*0.5;
            normFact = 1.0;
        }
        else
        {
            prGaussianMixtureSample<Ts>::lambda_g = 1.f;
            mOneO2l2 = 0.f;
            oneOl2 = 0.;
            normFact = 0.f;
        }
        /*
        halfSide = (unsigned int)(3.f*prGaussianMixtureSample<Ts>::lambda_g+0.5f);
        iSampleBlock = 2*halfSide+1;
        iSampleEnd = 2*halfSide+1;*/
    }
    
    /*!
     * \fn float getLambda() const
     * \brief Accessor to the Gaussian spread parameter lambda_g
     * \return the Gaussian spread parameter (\lambda_g)
     */
    float getLambda() const
    {
        return prGaussianMixtureSample<Ts>::lambda_g;
    }
    
    /*!
     * \fn int cartesianBuildFrom(prRegularlySampledCPImage<unsigned char> & IS_req, Ts & _s_g, bool computeFeatureDerivatives = false, double rho = 1., unsigned long sourceIndex = 0)
     * \brief Computes the Gaussian mixture sample value from the set of all Gaussians (one for each "pixel") , considering a Cartesian straight distance between samples, and, optionally (iff computeFeatureDerivatives == true), computes the photometric Gms drivatives with respect to the 6 DOFs pose
     * TODO: take degrees of freedom into account to save processing time
     * TODO: take the "surface constraint" into account (the fact that, for spherical images, the point set is on a sphere, leads to distances computed on the sphere surface, not straight is space) ; a parameter of the mother class of prRegularlySampledCPImage ?..
     * \param IS_req the input data from which to compute the Gaussian mixture sample value
     * \param _s_g the Gaussian mixture smaple parameters (essentially, the location)
     * \param computeFeatureDerivatives boolean to let, if true, the method to compute the derivates with respect to the pose or not, if false
     * \return 0 if the computation succeeds
     */
		template<typename T, template <typename TT> class Tacquisition>
    int buildFrom(Tacquisition<T> & IS_req, Ts & _s_g, bool computeFeatureDerivatives = false, double rho = 1., unsigned long sourceIndex = 0)
    {
        s_g = _s_g;
        gms = 0.;
        double gs = 0.;

        //vpRowVector dgmsXg(s_g.size());
        //dgmsXg = 0.;
        nbDim = s_g.size();
        L.resize(nbDim);
               
        Ts XS;
        float val;
        
        //double XmXg, YmYg, ZmZg;
        vpColVector XmXg(nbDim);
//        vpColVector Lt(nbDim);
        double dist2;

        iSampleStart = sourceIndex-halfSide-IS_req.width*halfSide;
        iSampleBlock = 2*halfSide+1;
        iSampleJump = IS_req.width-iSampleBlock;
        iSampleEnd = iSampleStart+iSampleBlock*IS_req.width;//+iSampleBlock;

//      std::cout << sourceIndex << " " << iSampleStart << " " << iSampleEnd << " " << iSampleBlock << " " << iSampleJump << std::endl;

        unsigned long nbSamples_ = iSampleBlock*iSampleBlock;
        unsigned long totalSamples = IS_req.getNbSamples();
        unsigned int iDim;

        if( (E_precomp != NULL) && (nbSamples_ != nbSamples) )
        {
          delete [] E_precomp;
          E_precomp = NULL;
          delete [] Eug_precomp_interlaced;
          Eug_precomp_interlaced = NULL;
        }
        nbSamples = nbSamples_;

        if(E_precomp == NULL)
        {
          E_precomp = new double[nbSamples];
          Eug_precomp_interlaced = new double[nbDim*nbSamples];
        }
        
        pt_E_precomp = E_precomp;
        pt_Eug_precomp_interlaced = Eug_precomp_interlaced;

        if(computeFeatureDerivatives)
        {
//          unsigned int count = 0;
          for(long iSample = iSampleStart ; iSample < iSampleEnd ; iSample+=iSampleBlock+iSampleJump)
          {
            /*if(iSample < 0)
              continue;*/
    
            for(long i = iSample; (i < iSample+iSampleBlock) && (i < totalSamples) ; i++, pt_E_precomp++)//, count++)
//          for(unsigned long i = 0 ; i < nbSamples ; i++, pt_E_precomp++)
            {
             /* if(i < 0)
                continue;*/

//              if(i == sourceIndex)
//                std::cout << count << std::endl;

//            std::cout << i << std::endl;

                IS_req.getSample(i, XS, val);

                XmXg = XS - s_g;
                //XmXg = XS.get_X()-s_g.get_X();
                //YmYg = XS.get_Y()-s_g.get_Y();
                //ZmZg = XS.get_Z()-s_g.get_Z();

                dist2 = XmXg.sumSquare();
                //dist2 = XmXg*XmXg + YmYg*YmYg + ZmZg*ZmZg;

                *pt_E_precomp = exp(dist2*mOneO2l2);
                gs = val*(*pt_E_precomp); //*normFact

    //            gs = val*exp(dist2*mOneO2l2);
                gms += gs;


                    //Lt += (*pt_Eug_precomp)*gs;

                    for(iDim = 0 ; iDim < nbDim ; iDim++, pt_Eug_precomp_interlaced++)
                    {
                      *pt_Eug_precomp_interlaced = XmXg[iDim]*oneOl2;
                      L.data[iDim] += (*pt_Eug_precomp_interlaced)*gs;//*XmXg*oneOl2;
                    }

    //                Lt = Lt + gs*XmXg*oneOl2;
                    //dgmsXg[0] = dgmsXg[0] + XmXg*gs*oneOl2;
                    //dgmsXg[1] = dgmsXg[1] + YmYg*gs*oneOl2;
                    //dgmsXg[2] = dgmsXg[2] + ZmZg*gs*oneOl2;

            }
          }
//          std::cout << count << " / " << nbSamples << std::endl;
        
        //la primitive ne calcul sa derivee que par rapport a son sampler
        //le reste calcule par l'estimateur selon son espace

            featureJacobianComputed = true;
        }
        else        
        {
          for(long iSample = iSampleStart ; iSample < iSampleEnd ; iSample+=iSampleBlock+iSampleJump)
          {
            /*if(iSample < 0)
              continue;*/
    
            for(long i = iSample; (i < iSample+iSampleBlock) && (i < totalSamples) ; i++, pt_E_precomp++)
//          for(unsigned long i = 0 ; i < nbSamples ; i++, pt_E_precomp++)
            {
              IS_req.getSample(i, XS, val);

              XmXg = XS - s_g;
              //XmXg = XS.get_X()-s_g.get_X();
              //YmYg = XS.get_Y()-s_g.get_Y();
              //ZmZg = XS.get_Z()-s_g.get_Z();

              dist2 = XmXg.sumSquare();
              //dist2 = XmXg*XmXg + YmYg*YmYg + ZmZg*ZmZg;

              *pt_E_precomp = exp(dist2*mOneO2l2);
              gs = val*(*pt_E_precomp); //*normFact

  //            gs = val*exp(dist2*mOneO2l2);
              gms += gs;
            }
          }
        }

        return 0;
    }
    
		//template<typename T, template <typename TT> class Tacquisition>
    //int updateFrom(Tacquisition<T> & IS_req, bool computeFeatureDerivatives = false, unsigned long sourceIndex = 0)
		template<typename T, template <typename TT> class Tacquisition>
    int updateFrom(Tacquisition<T> & IS_req, bool computeFeatureDerivatives = false, unsigned long sourceIndex = 0)
    {
        //gs = 0.;

        unsigned long totalSamples = IS_req.getNbSamples();
        
        //useful later
        //IS_req.getSampleRho(sourceIndex, rho);


        /*
        if(IS_req.getNbSamples() != nbSamples)
        {
          return -1;
        }
*/
        //vpColVector XmXg(nbDim);
        //vpColVector Lt(nbDim);
        //double dist2;

        pt_E_precomp = E_precomp;
        pt_Eug_precomp_interlaced = Eug_precomp_interlaced;

        //s_g = _s_g;
        gms = 0.;
        L.resize(nbDim);
        //L=0;

        pt_L0 = L.data;
//      std::cout << iSampleStart << " " << iSampleEnd << " " << iSampleBlock << " " << iSampleJump << std::endl;

        if(computeFeatureDerivatives)
        {
          //pt_bitmapf = IS_req.bitmapf+iSampleStart;
          //for(iSample = 0 ; iSample < nbSamples ; iSample++, pt_E_precomp++, pt_bitmapf++)
          for(long iSample = iSampleStart ; iSample < iSampleEnd ; iSample+=iSampleBlock+iSampleJump)//, pt_bitmapf+=iSampleJump)
          {
          /*  if(iSample < 0)
              continue;*/

            for(long i = iSample; (i < iSample+iSampleBlock) && (i < totalSamples) ; i++, pt_E_precomp++)//, pt_bitmapf++)
            {
/*              if(i < 0)
                continue;*/

//              std::cout << i << std::endl;

              IS_req.getSample(i, val);

              //XmXg = XS - s_g;

              //dist2 = XmXg.sumSquare();
              gs = val*(*pt_E_precomp);//exp(-dist2*oneO2l2); //*normFact
              //gs = (*pt_bitmapf)*(*pt_E_precomp);//exp(-dist2*oneO2l2); //*normFact

              gms += gs;


              //Lt += (*pt_Eug_precomp)*gs;//*XmXg*oneOl2;
              pt_L = pt_L0;
              for(iDim = 0 ; iDim < nbDim ; iDim++, pt_L++, pt_Eug_precomp_interlaced++)
              {
                *pt_L += (*pt_Eug_precomp_interlaced)*gs;//*XmXg*oneOl2;
              }

            }
          }
       
          //la primitive ne calcul sa derivee que par rapport a son sampler
          //le reste calcule par l'estimateur selon son espace
          featureJacobianComputed = true;
        }
        else
        {
          //for(iSample = 0 ; iSample < nbSamples ; iSample++, pt_E_precomp++)
          for(long iSample = iSampleStart ; iSample < iSampleEnd ; iSample+=iSampleBlock+iSampleJump)//, pt_bitmapf+=iSampleJump)
          {
          /*  if(iSample < 0)
              continue;*/

            for(long i = iSample; (i < iSample+iSampleBlock) && (i < totalSamples) ; i++, pt_E_precomp++)//, pt_bitmapf++)
            {
            IS_req.getSample(iSample, XS, val);

            //XmXg = XS - s_g;

            //dist2 = XmXg.sumSquare();

            gs = val*(*pt_E_precomp);//exp(-dist2*oneO2l2); //*normFact

            gms += gs;
            }
          }
        }
        
        return 0;
    }
    
    /*!
     * \fn int getFeatureJacobian(vpRowVector & Ls)
     * \brief Accessor to the feature Jacobian of the non-normalized Gaussian mixture sample
     * \param Ls (output) the output row vector that receives the photometric nnGMS interaction matrix row
     * \return  0 if the Jacobian could be get (if it was computed, actually)
     *         -1 if the Jacobian was not computed (no cannot be got)
     */
    int getFeatureJacobian(vpRowVector & Ls)
    {
        if(!featureJacobianComputed)
            return -1;
        
        Ls = L;
        
        return 0;
    }
    
    /*!
     * \fn prPhotometricnnGMS<Ts> operator-(const prPhotometricnnGMS<Ts> &f)
     * \brief - operator override to substract two photometric nnGms values
     * TODO : check line fOut.s_g = s_g - f.s_g;
     * \param f the nnGms to substract from the current
     * \return the prPhotometricnnGMS<Ts> object resulting from the difference computation
     */
    prPhotometricnnGMS<Ts> operator-(const prPhotometricnnGMS<Ts> &f)
    {
        prPhotometricnnGMS<Ts> fOut(getLambda());
        
        fOut.gms = gms - f.getGMS();
        //fOut.s_g = s_g - f.s_g;
        
        return fOut;
    }

    /*!
     * \fn prPhotometricnnGMS<Ts> operator*(const prPhotometricnnGMS<Ts> &f)
     * \brief * operator override to multiply two photometric nnGms values
     * TODO: check line fOut.s_g = s_g * f.s_g;
     * TODO: move in the mother class?
     * \param f the nnGms to multiply with the current
     * \return the prPhotometricnnGMS<Ts> object resulting from the product computation
     */
    prPhotometricnnGMS<Ts> operator*(const prPhotometricnnGMS<Ts> &f)
    {
        prPhotometricnnGMS<Ts> fOut(getLambda());

        fOut.gms = gms * f.getGMS();
        //fOut.s_g = s_g * f.s_g;
        
        return fOut;
    }
    
    /*!
     * \fn prPhotometricnnGMS<Ts> operator*(const double w)
     * \brief * operator override to multiply a nnGms value with a scalar
     * TODO : check line fOut.s_g = s_g * w;
     * \param w the scalar to multiply with the nnGms
     * \return the prPhotometricnnGMS<Ts> object resulting from the product with a scalar computation
     */
    prPhotometricnnGMS<Ts> operator*(const double w)
    {
        prPhotometricnnGMS<Ts> fOut(getLambda());

        fOut.gms = gms * w;
        //fOut.s_g = s_g * w;
        
        return fOut;
    }
    
    /*!
     * \fn prPhotometricnnGMS<Ts> &operator+=(const prPhotometricnnGMS<Ts> &f)
     * \brief += operator override to add two nnGms values, updating the current one
     * TODO : check line s_g += f.s_g;
     * \param f the nnGms to add to the current
     * \return the prPhotometricnnGMS<Ts> self reference object
     */
    prPhotometricnnGMS<Ts> &operator+=(const prPhotometricnnGMS<Ts> &f)
    {
        gms += f.getGMS();
        //s_g += f.s_g;
        
        return *this;
    }
    
    /*!
     * \fn prPhotometricnnGMS<Ts> &operator=(const prPhotometricnnGMS<Ts> &f)
     * \brief = operator override to copy nnGms object
     * \param f the nnGms to copy to the current
     * \return the prPhotometricnnGMS<Ts> self reference object
     */
    prPhotometricnnGMS<Ts> &operator=(const prPhotometricnnGMS<Ts> &f)
    {
        gms = f.getGMS();
        s_g = f.getS_g();
        nbDim = f.nbDim;
        mOneO2l2 = f.mOneO2l2;
        normFact = f.normFact;
        oneOl2 = f.oneOl2;
        prGaussianMixtureSample<Ts>::lambda_g = f.getLambda();
        if(f.doubleFeature != NULL)
        {
            if(prFeature::doubleFeature != NULL)
                delete [] prFeature::doubleFeature;
            prFeature::doubleFeature = new double[1];
            
            prFeature::doubleFeature[0] = f.doubleFeature[0];
        }
        featureJacobianComputed = f.featureJacobianComputed;

        halfSide = f.halfSide;
        iSampleStart = f.iSampleStart;
        iSampleBlock = f.iSampleBlock;
        iSampleJump = f.iSampleJump;
        iSampleEnd = f.iSampleEnd;
        L = f.L;

        if(f.E_precomp != NULL)
        {
          if( (E_precomp != NULL) && (f.nbSamples != nbSamples) )
          {
            delete [] E_precomp;
            E_precomp = NULL;
            delete [] Eug_precomp_interlaced;
            Eug_precomp_interlaced = NULL;
          }
          nbSamples = f.nbSamples;

          if(E_precomp == NULL)
          {
            E_precomp = new double[nbSamples];
            Eug_precomp_interlaced = new double[nbDim*nbSamples];
          }

          memcpy(E_precomp, f.E_precomp, nbSamples*sizeof(double));

          memcpy(Eug_precomp_interlaced, f.Eug_precomp_interlaced, nbDim*nbSamples*sizeof(double));
        }

        return *this;
    }
    
    /*!
     * \fn double getGMS() const
     * \brief Accessor to the Gms value
     * \return the Gms value
     */
    double getGMS() const { return gms;}
    
    /*!
     * \fn double *toDouble(unsigned int place = 0)
     * \brief implementation of the toDouble method of the mother class converting the photometric Gms rich feature type to a common double value
     * TODO: move in the mother class?
     * \param place ?
     * \return the pointer to the double type vlaue of the Gms
     */
    double *toDouble(unsigned int place = 0)
    {
        if(prFeature::doubleFeature != NULL)
            delete [] prFeature::doubleFeature;
        
        prFeature::doubleFeature = new double[1];
        prFeature::doubleFeature[0] = gms;
        
        return prFeature::doubleFeature;
    }
    
    /*!
     * \fn Ts getS_g() const
     * \brief Accessor to the Gaussian mixture sample location type parameter
     * TODO: move in the mother class?
     * \return the Gms location
     */
    Ts getS_g() const
    {
        return s_g;
    }
    
private:
    Ts s_g; /*!< !!!!!! feature ? Check the code using the class to remind */
    double gms; /*!< gaussian mixture sample response value */
    vpRowVector L; /*!< interaction matrix of the gaussian mixture sample over the pose*/
    
    double mOneO2l2; /*!< -1.0 / (2.0 * lambda_g^2) */
    double oneOl2; /*!< 1.0 / (lambda_g^2) */
    double normFact; /*!< 1.0 / (lambda_g^3 * (2.0 * PI)^3.0/2.0) */
    
    bool featureJacobianComputed; /*!< boolean to test if the pose Jacobian has already been computed to avoid useless computation */

    double *E_precomp;
    double *Eug_precomp_interlaced;
    unsigned long nbSamples;
    unsigned int nbDim;

    unsigned int halfSide;
    long iSampleStart;
    unsigned long iSampleBlock;
    unsigned long iSampleJump;
    long iSampleEnd;

    double gs;
    double *pt_E_precomp;
    double *pt_Eug_precomp_interlaced;
    unsigned int iDim;
    double *pt_L, *pt_L0;
    float *pt_bitmapf;
    unsigned long iSample;
    float val;
    Ts XS;
};

#endif  //_PRPHOTOMETRICnnGMS_H
