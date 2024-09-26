/*!
 \file prCameraPoseEstim.h
 \brief Definition of the class to estimate the pose of a camera
 \author Guillaume CARON
 \version 0.2
 \date December 2023
 */


#if !defined(_PRCAMERAPOSEESTIM_H)
#define _PRCAMERAPOSEESTIM_H

#include <per/prcommon.h>
#include <per/prPoseEstim.h>

#include <per/prCameraModel.h>

#include <per/prPointFeature.h>

//#define FICOUT

//#define LM

//#define VERBOSE

template<class FSetT, class FSetCmpT>
class PER_EXPORT prCameraPoseEstim : public prPoseEstim {
public:

/*!
     * \fn prCameraPoseEstim()
     * \brief Contructor of the prCameraPoseEstim
     */
    prCameraPoseEstim() : nb_s(0), iterMax(40), nbIncrMax(5), seuilResidu(1e-6), current(true), controlInit(false), pointFeatures(NULL), Depth(1.0), cam(nullptr)
    {
        setdof();
    }
    
    /*!
     * \fn ~prCameraPoseEstim()
     * \brief Destructor of a prCameraPoseEstim object
     */
    virtual ~prCameraPoseEstim() override
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

	/*!
     * \fn double track(FSetT & FeatureSetd, vpPoseVector & r_optim, bool robust = false)
     * \brief Update the pose of FeatureSet so that it matches FeatureSetd
     * \param FeatureSetd the input data toward which degrees of freedom will have an impact on FeatureSet
     * \param r the 6-vector that contains pose elements values
     * \return the residual error at convergence
     */
    double track(FSetT & FeatureSetd, vpPoseVector & r_optim, bool robust = false)
    {
        //Optimisation loop parameters
        double residual = -1;

		//to do
     
	    return residual;
    }

	/*!
     * \fn double initControl(FSetT & FeatureSetc, vpColVector & v, bool robust = false)
     * \brief Compute the pose increment to move FeatureSetc toward FeatureSet
     * \param FeatureSetc the current feature set
     * \param v the 6-vector that contains pose increment values
     * \return the residual at convergence
     */
    void initControl(float lambda_ = 1.f, float cst_Depth = 1.f, bool current_ = true)
    {
        current = current_;
        Depth = cst_Depth;

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
        
        initPointFeatures(); // actually the variable sampling the acquisition model

        iter = 0;
        nbIncr = 0;

        controlInit = true;
    }

	/*!
     * \fn void initPointFeatures()
     * \brief e Initialize the table of feature points (digital image plane and normalized image plane) associated to the sampler of the feature set 
     */
	void initPointFeatures()
  {
  
    if(pointFeatures != NULL)
    {
      delete [] pointFeatures;
      pointFeatures = NULL;
    }

    unsigned int n_samples = FeatureSet.sampler.nbSamples;
    pointFeatures = new prPointFeature[n_samples];

    prPointFeature *pt_pointFeatures = pointFeatures;
    for(unsigned int i = 0 ; i < n_samples ; i++, pt_pointFeatures++)
    {
      pt_pointFeatures->setPixUV(((prCartesian2DPointVec *)FeatureSet.sampler.ge)[i].get_x(), ((prCartesian2DPointVec *)FeatureSet.sampler.ge)[i].get_y());
      cam->pixelMeterConversion(*pt_pointFeatures);
    }

  }

	/*!
     * \fn void computeScenePointJacobian(prPointFeature & P, vpMatrix & LX)
     * \brief Compute the Jacobian of a 3D point with repsect to the sensor pose
     * \param P the point feature
     * \param LX the output interaction matrix of a 3D point (size: 3x6)
     */
	void computeScenePointJacobian(prPointFeature & P, vpMatrix & LX)
    {
    	double X = P.get_X(), Y = P.get_Y(), Z = P.get_Z();
	
    	LX.resize(3,6,false);
	
	    LX[0][0] = -1.;
	    LX[0][1] = 0.;
	    LX[0][2] = 0.;
	    LX[0][3] = 0.;
	    LX[0][4] = -Z;
	    LX[0][5] = Y;
	
	    LX[1][0] = 0.;
	    LX[1][1] = -1.;
	    LX[1][2] = 0.;
	    LX[1][3] = Z;
	    LX[1][4] = 0.;
	    LX[1][5] = -X;

	    LX[2][0] = 0.;
	    LX[2][1] = 0.;
	    LX[2][2] = -1.;
	    LX[2][3] = -Y;
	    LX[2][4] = X;
	    LX[2][5] = 0.;
    }

	/*!
     * \fn void setSensor(prSensorModel *sensor)
     * \brief sets the camera sensor model
     * \param sensor the camera model to consider
     */
	void setSensor(prSensorModel *sensor)
    {
      cam = (prCameraModel *)sensor;
    }

	/*!
     * \fn void computePoseJacobian(FSetT& FeatureSetc)
     * \brief computes the Jacobian of the variable sampling the acquisition model with respect to the sensor pose
     * \param FeatureSetc the current feature set
     */
    void computePoseJacobian(FSetT& FeatureSetc)
    {
      unsigned int nbDim = 2; //!!! specific
      unsigned int numDOF;
      unsigned long nbSamples = FeatureSet.set.size();
      vpMatrix Ls(nbDim, nbDOF), LuX(nbDim,3), LXP(3,6);
      double *pt_Lug, *pt_Ls;

      Lug.resize(nbDim*nbSamples, nbDOF, false);

      unsigned int lig = 0;
      prPointFeature *pt_pointFeatures = pointFeatures;
      double Depth = 1.0;
      for(unsigned long ns = 0 ; ns < nbSamples ; ns++, pt_pointFeatures++) //, it_set++
      {
        //computePointJacobian(it_set->s_g, Ls, au, av);// L_u_g si s_g est un point 2D numerique cartesien observe par une camera perspective
        //R�cup�ration d'un Depth coh�rent sinon utilisation du Depth voisin
        double Depth_c; //!!! specific
        FeatureSetc.input.getSampleDepth(ns, Depth_c); //!!! specific
        if (Depth_c > 0) Depth = Depth_c; //!!! specific

        bool to3DOK = true;
		
        if(cam != nullptr)
        {
          //std::cout << "unproject " << ns << std::endl;
          to3DOK = cam->unProject(*pt_pointFeatures, Depth);
        }

        if(to3DOK)
        {

          //std::cout << "computeScenePointJacobian" << std::endl;
          computeScenePointJacobian(*pt_pointFeatures, LXP); 

          if(cam != nullptr)
          {
            //std::cout << "computeSensorJacobian" << std::endl;
            cam->computeSensorJacobian(*pt_pointFeatures, LuX); //!!! needs to depend on the acquisition model sampler space
            Ls = LuX*LXP; // L_u_g si s_g est un point 2D numerique cartesien observe par une camera perspective
          }
          else // useless for the moment but useful later to align just point clouds together
          {
            Ls = LXP;
          }

        }
        else
          Ls=0;

        for (lig = 0; lig < nbDim; lig++)
        {
            pt_Lug = Lug[nbDim * ns + lig]; // get pointer on line ns of L
            pt_Ls = Ls[lig];
            for (numDOF = 0; numDOF < 6; numDOF++, pt_Ls++)
                if (dof[numDOF])
                {
                    *pt_Lug = *pt_Ls;
                    pt_Lug++;
                }
        }
      }

    }

	/*!
     * \fn void composeFeaturesPoseJacobians(vpMatrix & LfRaw, vpMatrix & Lug, vpMatrix & LsRaw)
     * \brief computes the Jacobian of the variable sampling the acquisition model with respect to the sensor pose
     * \param LfRaw features Jacobian of size nb features \times nb dimensions of the sampler elements
     * \param Lug sampler's Jacobian of size nb dimensions of the sampler elements * nb features \time nb degrees of freedom
     * \param LsRaw output Jacobian matrix resulting of per feature composition of LfRaw and Lug
     */
    void composeFeaturesPoseJacobians(vpMatrix & LfRaw, vpMatrix & Lug, vpMatrix & LsRaw)
    {
      unsigned int nbCol_Lf = LfRaw.getCols();
      unsigned int nbCol_Lug = Lug.getCols();
      vpRowVector LfRaw_l;
      double *pt_LfRaw_l;
      vpMatrix Lfug(1,nbCol_Lug);

      LsRaw.resize(LfRaw.getRows(), nbCol_Lug);//, false);
      double *pt_LfRaw0 = LfRaw.data;
      double *pt_LsRaw = LsRaw.data;
      double *pt_Lug0 = Lug.data;
      double *pt_Lug, *pt_LfRaw;

      unsigned int iCol_Lug;
      unsigned int iCol_Lf;
      unsigned int incr_pt_Lug0 = nbCol_Lug*(nbCol_Lf-1);

      for(unsigned int lig = 0 ; lig < LfRaw.getRows() ; lig++)
      {
        for(iCol_Lug = 0 ; iCol_Lug < nbCol_Lug ; iCol_Lug++, pt_LsRaw++, pt_Lug0++)
        {
          pt_LfRaw = pt_LfRaw0;
          pt_Lug = pt_Lug0;
          for(iCol_Lf = 0 ; iCol_Lf < nbCol_Lf ; iCol_Lf++, pt_LfRaw++, pt_Lug+=nbCol_Lug)
          {
            *pt_LsRaw += (*pt_LfRaw)*(*pt_Lug);
          }
          
        }
        pt_LfRaw0 = pt_LfRaw;
        pt_Lug0 += incr_pt_Lug0;
      }
    }

    /*!
     * \fn double interactionAndError(FSetT & FeatureSetc, vpMatrix & Lc, vpColVector &errc, bool robust = false)
     * \brief Compute the interaction matrix and error vector of FeatureSetc
     * \param FeatureSetc the current feature set
     * \param Lc the matrix that contains the interaction matrix
     * \param errc the column vector that contains the error vector
     * \return the residual or -1 if control not initialized
     */
    double interactionAndError(FSetT & FeatureSetc, vpMatrix & Lc, vpColVector &errc, bool robust = false)
    {
      if(!controlInit)
        return -1;

      double residual_back = residual;

      if(current)
      {
        //FeatureSetc.computeFeaturePoseJacobian(LsRaw, nbDOF, dof);
        //std::cout << "computeFeaturesJacobian" << std::endl;
        FeatureSetc.computeFeaturesJacobian(LfRaw); //, nbDOF, dof);
        //LfRaw: n_samples lig ; sample_dims cols
        /*computePoseJacobians(Lug);
        //Lug: n_samples*sample_dims lig ; nbDOF cols
        LsRaw = composeJacobians(LfRaw, Lug);
        */
        //maybe updatePoseJacobian(), in case Depth are not constant or samples moved/changed
        //LfRaw: n_samples lig ; nbDOF cols
        //std::cout << "computePoseJacobian" << std::endl;
        computePoseJacobian(FeatureSetc);
        //std::cout << "composeFeaturesPoseJacobians" << std::endl;
        composeFeaturesPoseJacobians(LfRaw, Lug, LsRaw);
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
        //e = H * Ls.t() * weighted_error ;
        Lc = H * Ls.t();
        errc = weighted_error;
#else //GN
        //e = Ls.pseudoInverse() * weighted_error;
        Lc = Ls;
        errc = weighted_error;
#endif
      }
      else
      {
#ifdef LM
        //  compute the control law
        //e = H * LsRaw.t() * error ;
        Lc = H * LsRaw.t();
        errc = error;
#else //GN
        //e = LsRaw.pseudoInverseEigen3() * error;
        Lc = LsRaw;
        errc = error;
#endif
      }
            
      iter++;
      
      return residual;
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
            //FeatureSetc.computeFeaturePoseJacobian(LsRaw, nbDOF, dof);
            //std::cout << "computeFeaturesJacobian" << std::endl;
            FeatureSetc.computeFeaturesJacobian(LfRaw); //, nbDOF, dof);
            //LfRaw: n_samples lig ; sample_dims cols
            /*computePoseJacobians(Lug);
            //Lug: n_samples*sample_dims lig ; nbDOF cols
            LsRaw = composeJacobians(LfRaw, Lug);
            */
            //maybe updatePoseJacobian(), in case Depth are not constant or samples moved/changed
            //LfRaw: n_samples lig ; nbDOF cols
            //std::cout << "computePoseJacobian" << std::endl;
            computePoseJacobian(FeatureSetc);
            //std::cout << "composeFeaturesPoseJacobians" << std::endl;
            composeFeaturesPoseJacobians(LfRaw, Lug, LsRaw);
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

    /*!
     * \fn double getCond(bool robust = false)
     * \brief Compute the condition number of the interaction matrix
     * \param robust to consider the raw of the robustified interaction matrix
     * \return the condition number
     */
//    template<class Tcmp>
    double getCond(bool robust = false)
    {
      vpMatrix V;
      vpColVector w;

      if(robust && Ls.getRows() == 0)
        robust = false;

      if(robust)
      {
        Ls.svdEigen3(w, V);
      }
      else
        LsRaw.svdEigen3(w, V);

      return w[0]/w[w.getRows()-1];
    }

private:
    bool dof[6]; /*!< commentaire */
    
    FSetT FeatureSet;
    prPointFeature *pointFeatures;
    double Depth; //!!! either Z or \rho or else depending on the camera model 

    unsigned int nb_s;

    unsigned int iterMax, nbIncrMax;
    double seuilResidu;
    
    std::ofstream ficIter;

    double residual;
#ifdef LM
    double mu;
#endif
    double lambda;
    vpColVector v, e;
        
    unsigned int nbDOF;

    // Interaction matrix, Hessian, error,...
    vpMatrix Lug, Ls, LfRaw, LsRaw; // desired interaction matrix
#ifdef LM
    vpMatrix Hs, H;
    vpMatrix diagHs;
#endif

    vpColVector error, error_back, weighted_error, weighted_error_back ; // Error FSet-FSetd
    vpColVector weights, weights_back;

    bool current;
    unsigned int iter, nbIncr;

    bool controlInit;

    prCameraModel *cam;

};

#undef VERBOSE

#endif  //_PRCAMERAPOSEESTIM_H
