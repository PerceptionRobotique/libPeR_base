/*!
 \file prPosePerspectiveEstim.h
 \brief Definition of the class to estimate the pose of a perspective camera
 \author Guillaume CARON
 \version 0.2
 \date September 2020
 */


#if !defined(_PRPOSEPERSPECTIVEESTIM_H)
#define _PRPOSEPERSPECTIVEESTIM_H

#include <per/prCameraPoseEstim.h>

#include <per/prCameraModel.h>
#include <per/prFeaturesComparator.h>

#include <per/prFeaturesSet.h>
#include <per/prPointFeature.h>

//#include <per/prCartesian2DPointVec.h>

//#define FICOUT

//#define LM

//#define VERBOSE

template<class FSetT, class FSetCmpT>
//template<class Tfeat, class FSetT, class FSetCmpT> //Tfeat useless?
class PER_EXPORT prPosePerspectiveEstim : public prCameraPoseEstim<FSetT, FSetCmpT> {

public:
    
    /*!
     * \fn prPosePerspectiveEstim()
     * \brief Contructor of the prPosePerspectiveEstim
     */
    prPosePerspectiveEstim() : nb_s(0), iterMax(40), nbIncrMax(5), seuilResidu(1e-6), current(true), controlInit(false), pointFeatures(NULL), Z(1.0), iZ(1.0)
    {
        setdof();
    }
    
    /*!
     * \fn ~prPosePerspectiveEstim()
     * \brief Destructor of a prPosePerspectiveEstim object
     */
    virtual ~prPosePerspectiveEstim() override
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
     * \param FeatureSetd the input image toward which degrees of freedom will have an impact on FeatureSet
     * \param r the 6-vector that contains pose elements values
     * \return the residual error at convergence
     */
    double track(FSetT & FeatureSetd, vpPoseVector & r_optim, bool robust = false)
    {
        //Optimisation loop parameters
        double residual = -1;
/* = 1e20, residual_1, residual_1_back;
        residual_1 = 2*residual;
#ifdef LM
        double mu = 0.01;
#endif
        double lambda = 1.0;
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
        Ls.resize(LsRaw.getRows(), LsRaw.getCols(), false);
        
        //vpMatrix LsRaw_des = LsRaw, LsRaw_cur;
        //LsRaw_cur.resize(LsRaw.getRows(), LsRaw.getCols(), false);
        
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
            //if (L_desired + L_current) / 2
            //FeatureSet.update(dMc, true);
            //FeatureSet.computeFeaturePoseJacobian(LsRaw_cur, nbDOF, dof);
            //LsRaw = LsRaw_cur + LsRaw_des;
            //LsRaw *= 0.5;
            
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
                
                //                Ls = LsFull;
                // //Mettre a zero les lignes pour les pixels notToBeProcessed...
                // sId.checkInteractionMatrixLines(Ls);
                // //A considerer pour la calcul des poids de robustesse
                 
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
                e = Ls.pseudoInverse() * weighted_error;
#endif
            }
            else
            {
#ifdef LM
                //  compute the control law
                e = H * LsRaw.t() * error ;
#else //GN
                e = LsRaw.pseudoInverse() * error;
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
        */
        return residual;
    }

    /*!
     * \fn double initControl(FSetT & FeatureSetc, vpColVector & v, bool robust = false)
     * \brief Compute the pose increment to move FeatureSetc toward FeatureSet
     * \param FeatureSetc the current feature set
     * \param v the 6-vector that contains pose increment values
     * \return the residual at convergence
     */
    void initControl(float lambda_ = 1.f, float cst_Z = 1.f, bool current_ = true)
    {
        current = current_;
        Z = cst_Z;

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
        
        // Compute the interaction matrix
        // link the variation of image intensity to camera motion
        // here it is computed at the desired position
        //FeatureSet.computeFeaturePoseJacobian(LsRaw, nbDOF, dof); //L_feat_r
        // non, c'est ici :
        initPointFeatures();
        //computePoseJacobian();
        //avec le set de poisitions de sampling dans FeatureSet --> Lug: n_samples*sample_dims lig ; nbDOF cols
        //Ls.resize(LsRaw.getRows(), LsRaw.getCols(), false);

        iter = 0;
        nbIncr = 0;

        controlInit = true;
    }

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

    //void CPosePerspective::computeJacobianForVVS(CPoint & P, CModel *cam, vpMatrix & Ls)
    void computePointJacobian(prPointFeature & P, vpMatrix & Ls, double au = 1.0, double av = 1.0, double* Z = nullptr)
    {
    	double x = P.get_x(), y = P.get_y();
    	/*double Z = P.get_Z(), iZ;*/

        if (!Z) Z = &(this->Z);
	    iZ = 1.0/(*Z);
	
    	Ls.resize(2,6,false);
	
	    Ls[0][0] = -au * iZ;
	    Ls[0][1] = 0;
	    Ls[0][2] = au * x * iZ;
	    Ls[0][3] = au * x * y;
	    Ls[0][4] = -au * (1.0+x*x);
	    Ls[0][5] = au * y;
	
	    Ls[1][0] = av * 0;
	    Ls[1][1] = -av * iZ;
	    Ls[1][2] = av * y * iZ;
	    Ls[1][3] = av * (1+y*y);
	    Ls[1][4] = -av * x * y;
	    Ls[1][5] = -av * x;
    }


    void setSensor(prSensorModel *sensor)
    {
      cam = (prCameraModel *)sensor;
    }

    void computePoseJacobian(FSetT& FeatureSetc)
    {
//      typename std::vector<Tfeat>::iterator it_set = FeatureSet.set.begin();

      //T u_g;
        unsigned int nbDim = 2;
      unsigned int numDOF;
      unsigned long nbSamples = FeatureSet.set.size();
      vpMatrix Ls(nbDim, nbDOF);
      double *pt_Lug, *pt_Ls;

      Lug.resize(nbDim*nbSamples, nbDOF, false);
        
      double au = cam->getau(), av = cam->getav();

      unsigned int lig = 0;
      prPointFeature *pt_pointFeatures = pointFeatures;
      double Z = 1.0;
      for(unsigned long ns = 0 ; ns < nbSamples ; ns++, pt_pointFeatures++) //, it_set++
      {
        //computePointJacobian(it_set->s_g, Ls, au, av);// L_u_g si s_g est un point 2D numerique cartesien observe par une camera perspective
        //R�cup�ration d'un Z coh�rent sinon utilisation du Z voisin
        double Zc;
        FeatureSetc.input.getSampleDepth(ns, Zc);
        if (Zc > 0) Z = Zc;
        computePointJacobian(*pt_pointFeatures, Ls, au, av, &Z);// L_u_g si s_g est un point 2D numerique cartesien observe par une camera perspective

        for (lig = 0; lig < nbDim; lig++)
        {
            pt_Lug = Lug[nbDim * ns + lig]; // get pointer on line ns of L
            pt_Ls = Ls[lig];//.data;
            for (numDOF = 0; numDOF < 6; numDOF++, pt_Ls++)
                if (dof[numDOF])
                {
                    *pt_Lug = *pt_Ls;
                    //std::cout << *pt_Lug << " ";
                    pt_Lug++;
                }
            //std::cout << std::endl;
        }
      }

    }

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
//        std::cout << "LfRaw " << LfRaw.extract(lig, 0, 1, nbCol_Lf) << std::endl;
//        std::cout << "Lug " << Lug.extract(nbCol_Lf*lig, 0, nbCol_Lf, nbCol_Lug) << std::endl;

/*        Lfug = LfRaw.extract(lig, 0, 1, nbCol_Lf)*Lug.extract(nbCol_Lf*lig, 0, nbCol_Lf, nbCol_Lug);
        memcpy(LsRaw[lig], Lfug.data, nbCol_Lug*sizeof(double));*/

        for(iCol_Lug = 0 ; iCol_Lug < nbCol_Lug ; iCol_Lug++, pt_LsRaw++, pt_Lug0++)
        {
//          std::cout << std::endl << iCol_Lug;
          pt_LfRaw = pt_LfRaw0;
          pt_Lug = pt_Lug0;
          for(iCol_Lf = 0 ; iCol_Lf < nbCol_Lf ; iCol_Lf++, pt_LfRaw++, pt_Lug+=nbCol_Lug)
          {
//            std::cout << " " << iCol_Lf;
            *pt_LsRaw += (*pt_LfRaw)*(*pt_Lug);
          }
          
        }
        pt_LfRaw0 = pt_LfRaw;
        pt_Lug0 += incr_pt_Lug0;
/*
        std::cout << "Prod ";
        for(int col = 0 ; col < 6 ; col++)
          std::cout << LsRaw[lig][col] << " ";
        std::cout << std::endl;*/
      }
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
            FeatureSetc.computeFeaturesJacobian(LfRaw); //, nbDOF, dof);
            //LfRaw: n_samples lig ; sample_dims cols
            /*computePoseJacobians(Lug);
            //Lug: n_samples*sample_dims lig ; nbDOF cols
            LsRaw = composeJacobians(LfRaw, Lug);
            */
            //maybe updatePoseJacobian(), in case Z are not constant or samples moved/changed
            //LfRaw: n_samples lig ; nbDOF cols
            computePoseJacobian(FeatureSetc);
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
     * \fn double control(FSetT & FeatureSetc, vpColVector & v, bool robust = false)
     * \brief Compute the pose increment to move FeatureSetc toward FeatureSet
     * \param FeatureSetc the current feature set
     * \param v the 6-vector that contains pose increment values
     * \return the residual at convergence or -1 if control not initialized
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
    double Z, iZ;

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

#endif  //_PRPOSEPERSPECTIVEESTIM_H
