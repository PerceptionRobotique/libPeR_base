/*!
 \file prPoseSphericalEstim.h
 \brief Definition of the class to estimate the pose of a spherical camera
 \author Guillaume CARON
 \version 0.2
 \date June 2017, January 2022
 */

#if !defined(_PRPOSESPHERICALESTIM_H)
#define _PRPOSESPHERICALESTIM_H

#include <per/prCameraPoseEstim.h>

#include <per/prRegularlySampledCSImage.h>
#include <per/prSensorModel.h>
#include <per/prFeaturesComparator.h>

#include <per/prFeaturesSet.h>

#include <visp/vpExponentialMap.h>

//#define FICOUT

// #define LM


/*!
 \class PER_EXPORT prPoseSphericalEstim prPoseSphericalEstim.h <per/prPoseSphericalEstim.h>
 \brief Class implementing pose estimation algorithms for spherical cameras
 */
template<typename FSetT, typename FSetCmpT>
class PER_EXPORT prPoseSphericalEstim : public prCameraPoseEstim {
    
public:
    
    /*!
     * \fn prPoseSphericalEstim()
     * \brief Contructor of the prPoseSphericalEstim
     */
    prPoseSphericalEstim(double seuilResidu_ = 1e-9) : nb_s(0), iterMax(10), nbIncrMax(5), seuilResidu(seuilResidu_) //1e-6 en gyro pur ; 1e-9 en odometrie
    {
        setdof();
    }
    
    /*!
     * \fn ~prPoseSphericalEstim()
     * \brief Destructor of a prPoseSphericalEstim object
     */
    virtual ~prPoseSphericalEstim() override
    {
        if(ficIter.is_open())
            ficIter.close();
    }
    
    /*!
     * \fn void setddl(bool tX, bool tY, bool tZ, bool rX, bool rY, bool rZ)
     * \brief Sets which degrees of freedom must be considered in the pose estimation
     * \param tX activates or deactivate the tX degree of freedom
     * \param tY activates or deactivate the tY degree of freedom
     * \param tZ activates or deactivate the tZ degree of freedom
     * \param rX activates or deactivate the rX degree of freedom
     * \param rY activates or deactivate the rY degree of freedom
     * \param rZ activates or deactivate the rZ degree of freedom
     * \return Nothing
     */
    void setdof(bool tX = true, bool tY = true, bool tZ = true, bool rX = true, bool rY = true, bool rZ = true)
    {
        dof[0] = tX;
        dof[1] = tY;
        dof[2] = tZ;
        dof[3] = rX;
        dof[4] = rY;
        dof[5] = rZ;
    }
    
    /*!
     * \fn void buildFrom(T & I)
     * \brief Build the features set of reference for the pose estimator
     * \param I_ the source image
     * \return Nothing
     */
    void buildFrom(FSetT & FeatureSet_)
    {
        FeatureSet = FeatureSet_;
    }

    /*!
     * \fn void startSaveIterations(char *fileName)
     * \brief Opens the file in which estimation iterations will be saved
     * \param fileName the filename
     * \return Nothing
     */
    void startSaveIterations(char *fileName)
    {
        if(ficIter.is_open())
            ficIter.close();
        ficIter.open(fileName);
    }
    
    /*!
     * \fn void stopSaveIterations(char *fileName)
     * \brief Closes the file in which estimation iterations are saved
     * \return Nothing
     */
    void stopSaveIterations()
    {
        if(ficIter.is_open())
            ficIter.close();
    }

    /**
     * \fn void setIterMax(unsigned int maxIter)
     * \brief Set the Iter Max 
     * \param maxIter set max iter 
     * \return Nothing
     */
    inline void setIterMax(unsigned int maxIter){iterMax = maxIter;}

    /**
     * \fn unsigned int getIterMax()
     * \brief Get the Iter Max number
     * \return unsigned int 
     */
    unsigned int getIterMax(){return iterMax;}
    
    /*double detect(vpImage<unsigned char> & I, double & tux_optim, double & tuy_optim, double & tuz_optim, double step = M_PI/180.0, double low_range = 0.0, double high_range = 2.0*M_PI, vpImage<unsigned char> *Id_rot = NULL);*/

    /*!
     * \fn double track(T & Id, double & tux_optim, double & tuy_optim, double & tuz_optim, vpImage<unsigned char> *Id_rot = NULL)
     * \brief Update the pose of Id so that it matches I
     * \param Id the input image on which degrees of freedom will have an impact
     * \param r the 6-vector that contains pose elements values
     * \return the features comparator at convergence
     */
    //FSetCmpT & track(FSetT & FeatureSetd, vpPoseVector & r_optim, bool robust = false)
    double track(FSetT & FeatureSetd, vpPoseVector & r_optim, double lambda = 1.0, bool robust = false)
    {
        //Optimisation loop parameters
        double residual = 1e20, residual_1, residual_1_back;
        residual_1 = 2*residual;
#ifdef LM
        double mu = 0.01;
#endif
        vpColVector v, v6(6), e;
        vpPoseVector r, r_back;
        
        //degrees of freedom parameters
        vpHomogeneousMatrix dMc;
        unsigned int nbDOF = 0, numDOF, indDOF;
        for (numDOF = 0 ; numDOF < 6 ; numDOF++)
            if (dof[numDOF])
                nbDOF++;
        
        // Interaction matrix, Hessien, error,...
        vpMatrix Ls, LsRaw; // desired interaction matrix
#ifdef LM
        vpMatrix Hs, H;
        vpMatrix diagHs(nbDOF,nbDOF) ;
#endif
        vpColVector error, error_back, weighted_error, weighted_error_back ; // Error FSet-FSetd
        vpColVector weights, weights_back;
        
        // Compute the interaction matrix
        // link the variation of image intensity to camera motion
        // here it is computed at the desired position
        FeatureSetd.computeFeaturePoseJacobian(LsRaw, nbDOF, dof); //L_feat_r
        
        /*vpMatrix LsRaw_des = LsRaw, LsRaw_cur;
        LsRaw_cur.resize(LsRaw.getRows(), LsRaw.getCols(), false);*/
        
        r = r_optim; // init si "odometrie"
        dMc.buildFrom(r);
        
        unsigned int iter = 0;
        unsigned int nbIncr = 0;
#ifdef FICOUT
        std::ofstream ficOut("opt_r.txt");
#endif //FICOUT
        while( (iter < iterMax) && (nbIncr < nbIncrMax) && (fabs(residual - residual_1) > seuilResidu) )
        {
            residual_1 = residual;
            residual = 0;
            
            //if L_desired
            //FeatureSet.update(dMc, false);
            
            //if L_current (best for 3dofs)
            FeatureSet.update(dMc, true); 
            FeatureSet.computeFeaturePoseJacobian(LsRaw, nbDOF, dof);
						
						Ls.resize(LsRaw.getRows(), LsRaw.getCols(), false); // does only a test if Ls is already of the correct size
						//If Ls' size is 0, it means none of the current or desired feature sets was initialized with computing the pose Jacobian --> TODO test

            //if (L_desired + L_current) / 2
            /*FeatureSet.update(dMc, true);
            FeatureSet.computeFeaturePoseJacobian(LsRaw_cur, nbDOF, dof);
            LsRaw = LsRaw_cur + LsRaw_des;
            LsRaw *= 0.5;*/
            
            FSetCmpT errorComputer(FeatureSet, FeatureSetd, robust);
            error = errorComputer.getFeaturesCmpVpColVector();
            if(robust)
            {
                weights = errorComputer.getRobustWeigths();
                residual = errorComputer.getRobustCost().toDouble()[0]; // a revoir
                if(residual == 0)
                {
                    weights = 1;
                    residual = errorComputer.getCost().toDouble()[0]; // a revoir
                }
            }
            else
                residual = errorComputer.getCost().toDouble()[0]; // a revoir
            
            std::cout << "residual : " << residual << std::endl;
            
            if((iter < 1) || (residual < residual_1)) // le precedent TC a permis de reduire l'erreur
            {
                std::cout << "dim : " << r.t() << std::endl;
                
                nbIncr = 0;
#ifdef LM
                //mu *= 0.5;
                mu *= 0.1;
#endif
                residual_1_back = residual_1;
                error_back = error;
                r_back = r;
                
                /*                Ls = LsFull;
                 //Mettre a zero les lignes pour les pixels notToBeProcessed...
                 sId.checkInteractionMatrixLines(Ls);
                 //A considerer pour la calcul des poids de robustesse
                 */
                if(robust)
                {
                    weighted_error.resize(error.getRows(), error.getCols(), false);
                    
                    for(int i = 0 ; i < error.size() ; i++)
                    {
                        double wi = weights[i];
                        for(int j = 0 ; j < nbDOF ; j++)
                            Ls[i][j] = wi*LsRaw[i][j];

                        weighted_error[i] = wi*error[i];
                    }
#ifdef LM
                    // Compute the Hessian H = L^TL
                    Hs = Ls.AtA() ;
#endif
                    weights_back = weights;
                    weighted_error_back = weighted_error;
                }
                else
                {
#ifdef LM
                    // Compute the Hessian H = L^TL
                    Hs = LsRaw.AtA() ;
#endif
                }

#ifdef LM
                // Compute the Hessian diagonal for the Levenberg-Marquartd
                // optimization process
                for(int j = 0 ; j < nbDOF ; j++) diagHs[j][j] = Hs[j][j];
#endif
#ifdef FICOUT
                ficOut << r.t() << std::endl << residual << std::endl;
#endif //FICOUT
            }
            else // le precedent TC a fait augmenter l'erreur
            {
                std::cout << "aug : " << r.t() << std::endl;
#ifdef LM
                nbIncr++;
                
                //on augmente mu et on revient dans la configuration de l'iteration precedente
                mu *= 10.0;
                //mu *= 2.0;

                residual = residual_1;
                residual_1 = residual_1_back;
                
                error = error_back;
                if(robust)
                {
                    weights = weights_back;
                    weighted_error = weighted_error_back;
                }
                
                r = r_back;
#endif
#ifndef LM
                //break;
#endif
            }
            
#ifdef LM
            //std::cout << "LM term" << std::endl;
            // Compute the levenberg Marquartd term
            {
                H = ((mu * diagHs) + Hs).inverseByLU();
            }
#endif
            if(robust)
            {
#ifdef LM
                            //std::cout << "control law" << std::endl;
                //  compute the control law
                e = H * Ls.t() * weighted_error ;
#else //GN
                e = Ls.pseudoInverseEigen3() * weighted_error;
#endif
            }
            else
            {
#ifdef LM
                //  compute the control law
                e = H * LsRaw.t() * error ;
#else //GN
                e = LsRaw.pseudoInverseEigen3() * error;
#endif
            }

            //
            v = -lambda*e;
            
            //update the DOFs
            indDOF = 0;
            for (numDOF = 0 ; numDOF < 6 ; numDOF++)
                if (dof[numDOF])
                {
                    v6[numDOF] = v[indDOF];
                    indDOF++;
                }
                else
                    v6[numDOF] = 0;

            dMc = vpExponentialMap::direct(v6).inverse()*dMc;
            r.buildFrom(dMc);
            
            iter++;
        }
#ifdef FICOUT
        ficOut.close();
#endif //FICOUT
        
        if(ficIter.is_open())
            ficIter << iter << std::endl;
        
        std::cout << iter << " residu : " << residual << " (new) | " << error.sumSquare() << " (Feat) ; nbIncr : " << nbIncr << std::endl;
        
        r_optim = r_back;
        
        return residual;
    }
    
    /*!
     * \fn double initControl(float lambda_ = 1.f, bool current_ = true)
     * \brief Initialize the matrices and the control law parameters (gain, ...)
     * \param lambda_ gain
     * \param current_ boolean to set current (true) or desired (false) interaction matrix
     */
    void initControl(float lambda_ = 1.f, bool current_ = true)
    {
        current = current_;

        //Optimisation loop parameters
        residual = 1e20;

#ifdef LM
        mu = 0.01;
#endif
        lambda = lambda_;
        
        //degrees of freedom parameters
        unsigned int numDOF, indDOF;
        nbDOF = 0;
        for (numDOF = 0 ; numDOF < 6 ; numDOF++)
            if (dof[numDOF])
                nbDOF++;
        
#ifdef LM
        diagHs.resize(nbDOF,nbDOF) ;
#endif
      
      	
      	if(!current)
      	{  
		      // Compute the interaction matrix
		      // link the variation of image feature to camera motion
		      // here it is computed at the desired position
		      FeatureSet.computeFeaturePoseJacobian(LsRaw, nbDOF, dof); 
        }
        else
        {
        	
        }
       	

        iter = 0;
        nbIncr = 0;

        controlInit = true;
    }

/*!
     * \fn double control(FSetT & FeatureSetc, vpColVector & v, bool robust = false)
     * \brief Compute the pose increment to move FeatureSetc toward FeatureSet
     * \param FeatureSetc the current feature set
     * \param v the 6-vector that contains pose increment values
     * \return the residual at convergence or -1 if control not initialized
     */
//    template<class Tcmp>
    double control(FSetT & FeatureSetc, vpColVector & v, bool robust = false)
    {
      if(!controlInit)
        return -1;

      double residual_back = residual;

      if(current)
      {
				FeatureSetc.computeFeaturePoseJacobian(LsRaw, nbDOF, dof);
						
				Ls.resize(LsRaw.getRows(), LsRaw.getCols(), false); // does only a test if Ls is already of the correct size
				//If Ls' size is 0, it means none of the current or desired feature sets was initialized with computing the pose Jacobian --> TODO test
      }      
            
      FSetCmpT errorComputer(FeatureSetc, FeatureSet, robust);
      error = errorComputer.getFeaturesCmpVpColVector();
      if(robust)
      {
        weights = errorComputer.getRobustWeigths();
        residual = errorComputer.getRobustCost().toDouble()[0]; // a revoir
        if(residual == 0)
        {
          weights = 1;
          residual = errorComputer.getCost().toDouble()[0]; // a revoir
        }
      }
      else
        residual = errorComputer.getCost().toDouble()[0]; // a revoir
            
//      std::cout << "residual : " << residual << std::endl;
            
      if((iter < 1) || (residual < residual_back)) // le precedent TC a permis de reduire l'erreur
      {
#ifdef VERBOSE
          std::cout << "dim : " << std::endl;
#endif
          nbIncr = 0;
#ifdef LM
          //mu *= 0.5;
          mu *= 0.1;
#endif
          residual_back = residual;
      }
      else
      {
#ifdef VERBOSE
          std::cout << "aug : " << std::endl;
#endif

#ifdef LM
          nbIncr++;
                
          //on augmente mu et on revient dans la configuration de l'iteration precedente
          mu *= 10.0;
#endif
      }                

      if(robust)
      {
          weighted_error.resize(error.getRows(), error.getCols(), false);
                    
          for(int i = 0 ; i < error.size() ; i++)
          {
            double wi = weights[i];
            for(int j = 0 ; j < nbDOF ; j++)
              Ls[i][j] = wi*LsRaw[i][j];

            weighted_error[i] = wi*error[i];
          }
#ifdef LM
          // Compute the Hessian H = L^TL
          Hs = Ls.AtA() ;
#endif

      }
      else
      {
#ifdef LM
        // Compute the Hessian H = L^TL
        Hs = LsRaw.AtA() ;
#endif
      }

#ifdef LM
      // Compute the Hessian diagonal for the Levenberg-Marquartd
      // optimization process
      for(int j = 0 ; j < nbDOF ; j++) diagHs[j][j] = Hs[j][j];

      //std::cout << "LM term" << std::endl;
      // Compute the levenberg Marquartd term
      {
          H = ((mu * diagHs) + Hs).inverseByLU();
      }
#endif
      if(robust)
      {
#ifdef LM
        //std::cout << "control law" << std::endl;
        //  compute the control law
        e = H * Ls.t() * weighted_error ;
#else //GN
        e = Ls.pseudoInverse() * weighted_error;
#endif
      }
      else
      {
#ifdef LM
        //  compute the control law
        e = H * LsRaw.t() * error ;
#else //GN
        e = LsRaw.pseudoInverseEigen3() * error;
#endif
      }

      //
      v = -lambda*e;
            
      iter++;
      
      return residual;
    }
    
private:
    bool dof[6]; /*!< commentaire */
    
    FSetT FeatureSet;

    unsigned int nb_s;
    
    unsigned int iterMax, nbIncrMax;
    double seuilResidu;
    
    std::ofstream ficIter;
    
    // TODO: delete below variable declarations in function track
		double residual;
#ifdef LM
    double mu;
#endif
    double lambda;
    vpColVector v, e;
        
    unsigned int nbDOF;

    // Interaction matrix, Hessian, error,...
    vpMatrix Ls, LsRaw; // interaction matrix
#ifdef LM
    vpMatrix Hs, H;
    vpMatrix diagHs;
#endif

    vpColVector error, error_back, weighted_error, weighted_error_back ; // Error FSet-FSetd
    vpColVector weights, weights_back;

    bool current;
    unsigned int iter, nbIncr;

    bool controlInit;
        
};

#endif  //_PRPOSESPHERICALESTIM_H
