/*!
 \file prPhotometricGMS.h
 \brief Header file for the prPhotometricGMS class
 \author Guillaume CARON
 \version 0.1
 \date august 2017
 */

#if !defined(_PRPHOTOMETRICGMS_H)
#define _PRPHOTOMETRICGMS_H

#include <per/prGaussianMixtureSample.h>
//#include <per/prRegularlySampledCSImage.h>

#define PIPIP3S2 15.7496099457
#define PIS3 1.0471975512

/*!
 \class PER_EXPORT prPhotometricGMS prPhotometricGMS.h <per/prPhotometricGMS.h>
 \brief Class defining a photometric based gaussian mixture sample
 */
template<typename Ts>
class PER_EXPORT prPhotometricGMS : public prGaussianMixtureSample<Ts> {
public:
    /*!
     * \fn prPhotometricGMS<Ts>(float lambda_g = 1.f)
     * \brief Contructor of the prPhotometricGMS<Ts> photometric Gaussian mixture sample feature type
     * \param lambda_g the Gaussian expansion parameter (\lambda_g)
     */
    prPhotometricGMS<Ts>(float lambda_g = 1.f, bool boundedGaussian_ = false): /*prGaussianMixtureSample<Ts>(lambda_g),*/ gms(0.f), L(6), poseJacobianComputed(false), boundedGaussian(boundedGaussian_), kernelBound(-1)
    {
        setLambda(lambda_g);
    }
    
    /*!
     * \fn prPhotometricGMS<Ts>(prPhotometricGMS<Ts> const &Src)
     * \brief Copy contructor of a prPhotometricGMS<Ts> photometric Gaussian mixture sample feature type
     * \param Src source object to copy
     */
    prPhotometricGMS<Ts>(prPhotometricGMS<Ts> const &Src): gms(Src.getGMS()), poseJacobianComputed(Src.poseJacobianComputed), boundedGaussian(Src.boundedGaussian)
    {
        setLambda(Src.getLambda());
        s_g = Src.getS_g();
        L = Src.L;
    }
    
    /*!
     * \fn ~prPhotometricGMS
     * \brief Destructor of a prPhotometricGMS photometric Gaussian mixture sample feature object
     */
    virtual ~prPhotometricGMS() override
    {
        
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
            oneO2l2 = oneOl2*0.5;
            normFact = oneOl2 / (prGaussianMixtureSample<Ts>::lambda_g * PIPIP3S2);
						if(boundedGaussian && (prGaussianMixtureSample<Ts>::lambda_g < PIS3))
							kernelBound = cos(3.0*prGaussianMixtureSample<Ts>::lambda_g);
        }
        else
        {
            prGaussianMixtureSample<Ts>::lambda_g = 1.f;
            oneO2l2 = 0.f;
            oneOl2 = 0.;
            normFact = 0.f;
        }
    }
    
    /*!
     * \fn float getLambda() const
     * \brief Accessor tp the Gaussian spread parameter lambda_g
     * \return the Gaussian spread parameter (\lambda_g)
     */
    float getLambda() const
    {
        return prGaussianMixtureSample<Ts>::lambda_g;
    }
    
    /*!
     * \fn int cartesianBuildFrom(prRegularlySampledCSImage<unsigned char> & IS_req, Ts & _s_g, bool computePoseDerivatives = false)
     * \brief Computes the Gaussian mixture sample value from the set of all Gaussians (one for each "pixel") , considering a Cartesian straight distance between samples, and, optionally (iff computePoseDerivatives == true), computes the photometric Gms drivatives with respect to the 6 DOFs pose
     * TODO: take degrees of freedom into account to save processing time
     * TODO: take the "surface constraint" into account (the fact that, for spherical images, the point set is on a sphere, leads to distances computed on the sphere surface, not straight is space) ; a parameter of the mother class of prRegularlySampledCSImage ?..
     * \param IS_req the input data from which to compute the Gaussian mixture sample value
     * \param _s_g the Gaussian mixture smaple parameters (essentially, the location)
     * \param computePoseDerivatives boolean to let, if true, the method to compute the derivates with respect to the pose or not, if false
     * \return 0 if the computation succeeds
     */
     template<typename T, template <typename TT> class Tacquisition>
    int cartesianBuildFrom(Tacquisition<T> & IS_req, Ts & _s_g, bool computePoseDerivatives = false)
    {
        s_g = _s_g;
        gms = 0.;
        double gs = 0.;
        vpRowVector dgmsXg(3);
        dgmsXg = 0.;
        
        Ts XS;
        float val;
        
        unsigned long nbSamples = IS_req.getNbSamples();
        double XmXg, YmYg, ZmZg;
        double dist2;
        
        for(unsigned long i = 0 ; i < nbSamples ; i++)
        {
            IS_req.getSample(i, XS, val);
            
            XmXg = XS.get_X()-s_g.get_X();
            YmYg = XS.get_Y()-s_g.get_Y();
            ZmZg = XS.get_Z()-s_g.get_Z();
            
            dist2 = XmXg*XmXg + YmYg*YmYg + ZmZg*ZmZg;

            gs = val*normFact*exp(-dist2*oneO2l2);

            gms += gs;

            if(computePoseDerivatives)
            {
                dgmsXg[0] = dgmsXg[0] + XmXg*gs*oneOl2;
                dgmsXg[1] = dgmsXg[1] + YmYg*gs*oneOl2;
                dgmsXg[2] = dgmsXg[2] + ZmZg*gs*oneOl2;
            }
        }
        
        if(computePoseDerivatives)
        {
            vpMatrix dXgdr(3,6);
            dXgdr[0][0] = -1.; dXgdr[0][1] =  0.; dXgdr[0][2] =  0.;
            dXgdr[1][0] =  0.; dXgdr[1][1] = -1.; dXgdr[1][2] =  0.;
            dXgdr[2][0] =  0.; dXgdr[2][1] =  0.; dXgdr[2][2] = -1.;
            
            dXgdr[0][3] =           0.; dXgdr[0][4] = -s_g.get_Z(); dXgdr[0][5] =  s_g.get_Y();
            dXgdr[1][3] =  s_g.get_Z(); dXgdr[1][4] =           0.; dXgdr[1][5] = -s_g.get_X();
            dXgdr[2][3] = -s_g.get_Y(); dXgdr[2][4] =  s_g.get_X(); dXgdr[2][5] =           0.;
            
            L = dgmsXg*dXgdr;

            poseJacobianComputed = true;
        }
        
        return 0;
    }
    
    /*!
     * \fn int buildFrom(prRegularlySampledCSImage<unsigned char> & IS_req, Ts & _s_g, bool computePoseDerivatives = false, double rho = 1., unsigned long sourceIndex = 0)
     * \brief Computes the Gaussian mixture sample value from the set of all Gaussians (one for each "pixel") , considering a circle arch distance between samples, since they are lying on a spherical surface, and, optionally (iff computePoseDerivatives == true), computes the photometric Gms drivatives with respect to the 6 DOFs pose
     * TODO: take degrees of freedom into account to save processing time
     * TODO: take the "surface constraint" into account (the fact that, for spherical images, the point set is on a sphere, leads to distances computed on the sphere surface, not straight is space) ; a parameter of the mother class of prRegularlySampledCSImage ?..
     * \param IS_req the input data from which to compute the Gaussian mixture sample value
     * \param _s_g the Gaussian mixture sample parameters (essentially, the location)
     * \param computePoseDerivatives boolean to let, if true, the method to compute the derivates with respect to the pose or not, if false
		 * \param rho the distance of 3D point to the camera (1 by default if no 3D data is used)
     * \return 0 if the computation succeeds
     */
    template<typename T, template <typename TT> class Tacquisition>
    int buildFrom(Tacquisition<T> & IS_req, Ts & _s_g, bool computePoseDerivatives = false, double rho = 1., unsigned long sourceIndex = 0)
    {
	      s_g = _s_g;
			
				//ensure s_g coordinates are of unit norm (OK for SImage only)
				double iNorm = 1.0/sqrt(s_g.get_X()*s_g.get_X()+s_g.get_Y()*s_g.get_Y()+s_g.get_Z()*s_g.get_Z());
				s_g.setPoint(s_g.get_X()*iNorm, s_g.get_Y()*iNorm, s_g.get_Z()*iNorm);

	      gms = 0.;
	      double gs = 0.;
	      vpRowVector dgmsXg(3);
	      dgmsXg = 0.;
	      
	      Ts XS;
	      float val;
	      
	      unsigned long nbSamples = IS_req.getNbSamples();
	      double dot;
	      double dist;
	      double fact_dHdXs;
	      double sqrt1mdotdot;
	      
				if(boundedGaussian)
				{
			    for(unsigned long i = 0 ; i < nbSamples ; i++)
			    {
			        IS_req.getSample(i, XS, val);

			        dot = XS.get_X()*s_g.get_X() + XS.get_Y()*s_g.get_Y() + XS.get_Z()*s_g.get_Z();
							if(dot > kernelBound)
							{
			        	dot = (dot < -1.)?-1.:(dot > 1.)?1.:dot;
			        	dist = acos(dot);
			        
			        	gs = val*normFact*exp(-dist*dist*oneO2l2);
			        
			        	gms += gs;
			        	if(computePoseDerivatives)
			        	{
			            sqrt1mdotdot = sqrt(1.-dot*dot);
			            if(sqrt1mdotdot != 0.)
			                fact_dHdXs = gs*oneOl2*dist/sqrt1mdotdot;
			            else
			                fact_dHdXs = 0.;
			            
			            dgmsXg[0] = dgmsXg[0] + fact_dHdXs*XS.get_X();
			            dgmsXg[1] = dgmsXg[1] + fact_dHdXs*XS.get_Y();
			            dgmsXg[2] = dgmsXg[2] + fact_dHdXs*XS.get_Z();
			        	}
			    		}
					}
				}
				else
				{
			    for(unsigned long i = 0 ; i < nbSamples ; i++)
			    {
			        IS_req.getSample(i, XS, val);

			        dot = XS.get_X()*s_g.get_X() + XS.get_Y()*s_g.get_Y() + XS.get_Z()*s_g.get_Z();
			        dot = (dot < -1.)?-1.:(dot > 1.)?1.:dot;
			        dist = acos(dot);
			        
			        gs = val*normFact*exp(-dist*dist*oneO2l2);
			        
			        /*std::cout << val << " " << normFact << " " << exp(-dist*dist*oneO2l2) << " " << dist << " " << dot << " " << s_g.get_X() << " " << s_g.get_Y() << " " << s_g.get_Z() << " " << oneO2l2 << std::endl;
			        
			        if(isnan(gs))
			        	exit(-13);
			        	*/		        
			        gms += gs;
			        if(computePoseDerivatives)
			        {
			            sqrt1mdotdot = sqrt(1.-dot*dot);
			            if(sqrt1mdotdot != 0.)
			                fact_dHdXs = gs*oneOl2*dist/sqrt1mdotdot;
			            else
			                fact_dHdXs = 0.;
			            
			            dgmsXg[0] = dgmsXg[0] + fact_dHdXs*XS.get_X();
			            dgmsXg[1] = dgmsXg[1] + fact_dHdXs*XS.get_Y();
			            dgmsXg[2] = dgmsXg[2] + fact_dHdXs*XS.get_Z();
			        }
			    }
				}
					
					//std::cout << gs << std::endl;
		    if(rho > 0.)
        {
		      if(computePoseDerivatives)
		      {
		          vpMatrix dXgdr(3,6);
							double iRho = 1./rho;
							
							//std::cout << "iRho: " << iRho  << " " << s_g.get_X() << " " << s_g.get_Y() << " " << s_g.get_Z() << std::endl;
							
							dXgdr[0][0] = (s_g.get_X()*s_g.get_X()-1.)*iRho; dXgdr[0][1] =      s_g.get_X()*s_g.get_Y()*iRho; dXgdr[0][2] =      s_g.get_X()*s_g.get_Z()*iRho;
		          dXgdr[1][0] =                       dXgdr[0][1]; dXgdr[1][1] = (s_g.get_Y()*s_g.get_Y()-1.)*iRho; dXgdr[1][2] =      s_g.get_Y()*s_g.get_Z()*iRho;
		          dXgdr[2][0] =                       dXgdr[0][2]; dXgdr[2][1] =                       dXgdr[1][2]; dXgdr[2][2] = (s_g.get_Z()*s_g.get_Z()-1.)*iRho;
							//why was it like below before correcting with the actual columns of translational dofs?
		          //dXgdr[0][0] = -1.; dXgdr[0][1] =  0.; dXgdr[0][2] =  0.;
		          //dXgdr[1][0] =  0.; dXgdr[1][1] = -1.; dXgdr[1][2] =  0.;
		          //dXgdr[2][0] =  0.; dXgdr[2][1] =  0.; dXgdr[2][2] = -1.;
		          
		          dXgdr[0][3] =           0.; dXgdr[0][4] = -s_g.get_Z(); dXgdr[0][5] =  s_g.get_Y();
		          dXgdr[1][3] =  s_g.get_Z(); dXgdr[1][4] =           0.; dXgdr[1][5] = -s_g.get_X();
		          dXgdr[2][3] = -s_g.get_Y(); dXgdr[2][4] =  s_g.get_X(); dXgdr[2][5] =           0.;
		          
		          L = dgmsXg*dXgdr;
		          
		          poseJacobianComputed = true;
		      }
        }
        else
        {
        	L.resize(1,6);
        	poseJacobianComputed = true;
        }
        
        return 0;
    }

   // template<typename T>
   template<typename T, template <typename TT> class Tacquisition>
    int updateFrom(Tacquisition<T> & IS_req, bool computeFeatureDerivatives = false, unsigned long sourceIndex = 0)
    {
    	/*
      std::cout << "prPhotometricGMS::updateFrom not implemented yet" << std::endl;
      return -1;
      */
      double rho;
      
      IS_req.getSampleRho(sourceIndex, rho);
      return buildFrom(IS_req, s_g, computeFeatureDerivatives, rho, sourceIndex);
    }
    
    /*!
     * \fn int getPoseJacobian(vpRowVector & Ls)
     * \brief Accessor to the pose Jacobian of the Gaussian mixture sample
     * \param Ls (output) the output row vector that receives the photometric Gms interaction matrix row
     * \return  0 if the Jacobian could be get (if it was computed, actually)
     *         -1 if the Jacobian was not computed (no cannot be got)
     */
    int getPoseJacobian(vpRowVector & Ls)
    {
        if(!poseJacobianComputed)
            return -1;
        
        Ls = L;
        
        return 0;
    }
    
    /*!
     * \fn prPhotometricGMS<Ts> operator-(const prPhotometricGMS<Ts> &f)
     * \brief - operator override to substract two photometric Gms values
     * TODO : check line fOut.s_g = s_g - f.s_g;
     * \param f the Gms to substract from the current
     * \return the prPhotometricGMS<Ts> object resulting from the difference computation
     */
    prPhotometricGMS<Ts> operator-(const prPhotometricGMS<Ts> &f)
    {
        prPhotometricGMS<Ts> fOut(getLambda());
        
        fOut.gms = gms - f.getGMS();
        //fOut.s_g = s_g - f.s_g;
        
        return fOut;
    }

    /*!
     * \fn prPhotometricGMS<Ts> operator*(const prPhotometricGMS<Ts> &f)
     * \brief * operator override to multiply two photometric Gms values
     * TODO: check line fOut.s_g = s_g * f.s_g;
     * TODO: move in the mother class?
     * \param f the Gms to multiply with the current
     * \return the prPhotometricGMS<Ts> object resulting from the product computation
     */
    prPhotometricGMS<Ts> operator*(const prPhotometricGMS<Ts> &f)
    {
        prPhotometricGMS<Ts> fOut(getLambda());

        fOut.gms = gms * f.getGMS();
        //fOut.s_g = s_g * f.s_g;
        
        return fOut;
    }
    
    /*!
     * \fn prPhotometricGMS<Ts> operator*(const double w)
     * \brief * operator override to multiply a Gms value with a scalar
     * TODO : check line fOut.s_g = s_g * w;
     * \param w the scalar to multiply with the Gms
     * \return the prPhotometricGMS<Ts> object resulting from the product with a scalar computation
     */
    prPhotometricGMS<Ts> operator*(const double w)
    {
        prPhotometricGMS<Ts> fOut(getLambda());

        fOut.gms = gms * w;
        //fOut.s_g = s_g * w;
        
        return fOut;
    }
    
    /*!
     * \fn prPhotometricGMS<Ts> &operator+=(const prPhotometricGMS<Ts> &f)
     * \brief += operator override to add two Gms values, updating the current one
     * TODO : check line s_g += f.s_g;
     * \param f the Gms to add to the current
     * \return the prPhotometricGMS<Ts> self reference object
     */
    prPhotometricGMS<Ts> &operator+=(const prPhotometricGMS<Ts> &f)
    {
        gms += f.getGMS();
        //s_g += f.s_g;
        
        return *this;
    }
    
    /*!
     * \fn prPhotometricGMS<Ts> &operator=(const prPhotometricGMS<Ts> &f)
     * \brief = operator override to copy Gms object
     * \param f the Gms to copy to the current
     * \return the prPhotometricGMS<Ts> self reference object
     */
    prPhotometricGMS<Ts> &operator=(const prPhotometricGMS<Ts> &f)
    {
        gms = f.getGMS();
        s_g = f.getS_g();
        oneO2l2 = f.oneO2l2;
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
        poseJacobianComputed = f.poseJacobianComputed;
				boundedGaussian = f.boundedGaussian;
				kernelBound = f.kernelBound;

        L = f.L;
        
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
    
    double oneO2l2; /*!< 1.0 / (2.0 * lambda_g^2) */
    double oneOl2; /*!< 1.0 / (lambda_g^2) */
    double normFact; /*!< 1.0 / (lambda_g^3 * (2.0 * PI)^3.0/2.0) */
    
    bool poseJacobianComputed; /*!< boolean to test if the pose Jacobian has already been computed to avoid useless computation */

		bool boundedGaussian; /*!< boolean to test if GMS is computed on the whole domain or a subpart of it */
		double kernelBound; /*!< cos(3*lambda_g) */
};

#endif  //_PRPHOTOMETRICGMS_H
