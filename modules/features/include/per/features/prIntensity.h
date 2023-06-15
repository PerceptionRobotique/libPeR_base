/*!
 \file prIntensity.h
 \brief Header file for the generic prIntensity class
 \author Guillaume CARON
 \version 0.1
 \date january 2017
 */

#if !defined(_PRINTENSITY_H)
#define _PRINTENSITY_H

#include <per/prFeature.h>
#include <per/prRegularlySampledCSImage.h>

#include <per/prNeighborhood.h>

/*!
 \class PER_EXPORT prIntensity prIntensity.h <per/prIntensity.h>
 \brief Class for the generic intensity feature
 */
template <typename Ts, class Tcam>
class PER_EXPORT prIntensity : public prFeature
{
public:
    /*!
     * \fn prIntensity<Ts>()
     * \brief Contructor of the prIntensity<Ts> photometric sample feature type
     */
    prIntensity<Ts, Tcam>() : intensityVal(0.f), L(6), poseJacobianComputed(false)
    {
        isNeighSet = false;
        nbDof = 6;
    }

    /*!
     * \fn prIntensity<Ts>(prIntensity<Ts> const &Src)
     * \brief Copy contructor of a prIntensity<Ts> photometric sample feature type
     * \param Src source object to copy
     */
    prIntensity<Ts, Tcam>(prIntensity<Ts, Tcam> const &Src) : intensityVal(Src.getGMS()), poseJacobianComputed(Src.poseJacobianComputed)
    {
        isNeighSet = false;
        s_g = Src.getS_g();
        L = Src.L;
        camera = Src.camera;
    }

    /*!
     * \fn ~prIntensity
     * \brief Destructor of a prIntensity photometric sample feature object
     */
    ~prIntensity()
    {
    }
    //virtual ~prIntensity() override = default;

    template <typename T, template <typename TT> class Tacquisition>
    int buildFrom(Tacquisition<T> &IS_req, Ts &_s_g, bool computePoseDerivatives = false, double rho = 1., unsigned long sourceIndex = 0)
    {
        return buildFromCam(camera, IS_req, _s_g, computePoseDerivatives, rho, sourceIndex);
    }

    /*!
     * \fn int buildFrom(prEquirectangular *equiCam, prRegularlySampledCSImage<unsigned char> & IS_req, Ts & _s_g, bool computePoseDerivatives = false, double rho = 1., unsigned long sourceIndex = 0)
     * \brief Computes the sample value from the set of a sample (one for each "pixel") from an equirectangular given camera, and, optionally (iff computePoseDerivatives == true), computes the photometric drivatives with respect to the 6 DOFs pose
     * TODO: take degrees of freedom into account to save processing time
     * TODO: find how to take properly account of the image type, if unsigned char the photometric value is put in the [0,1] interval
     * \param equiCam pointer to the equirectangular camera, used to project the feature on the image plane and likewise in the other direction
     * \param IS_req the input data from which to compute the Gaussian mixture sample value
     * \param _s_g the Gaussian mixture sample parameters (essentially, the location)
     * \param computePoseDerivatives boolean to let, if true, the method to compute the derivates with respect to the pose or not, if false
     * \param rho the distance of 3D point to the camera (1 by default if no 3D data is used)
     * \param sourceIndex the index of the sample to consider
     * \return 0 if the computation succeeds
     */
    template <typename T, template <typename TT> class Tacquisition>
    int buildFromCam(prEquirectangular *equiCam, Tacquisition<T> &IS_req, Ts &_s_g, bool computePoseDerivatives = false, double rho = 1., unsigned long sourceIndex = 0)
    {
        s_g = _s_g;

        int imHeight = IS_req.I_req.getHeight();
        int imWidth = IS_req.I_req.getWidth();

        // CALCUL DE L'INTENSITE DE LA FEATURE MISE A JOUR
        float val = 0.;

        vpRowVector derivativeVec(3);
        derivativeVec = 0.;

        if (isNeighSet == false)
        {
            prNeighborhood neighborhoodEquiCam;
            neighborhoodEquiCam.buildNeighborsCartSphere(imWidth, imHeight, equiCam);
            deltaCSSampling = neighborhoodEquiCam.deltaCSSampling;

            Neigh = neighborhoodEquiCam.Neigh;

            isNeighSet = true;

            std::cout << "neigh ok" << std::endl;
        }

        prPointFeature P_temp;
        P_temp.set_X(s_g.get_X());
        P_temp.set_Y(s_g.get_Y());
        P_temp.set_Z(s_g.get_Z());

        equiCam->project3DImage(P_temp);
        equiCam->meterPixelConversion(P_temp);

        int u = P_temp.get_u();
        int v = P_temp.get_v();

        if ((u > 5) && (v > 5) && (u < (imWidth - 6)) && (v < (imHeight - 6)))
        {
            indexU = (unsigned int)v;
            indexV = (unsigned int)u;
            T valTemp = IS_req.I_req[indexU][indexV];
            intensityVal = (float)valTemp;

            if (IS_req.MaskS[indexU][indexV] != 0)
            {
                IXs = derivativeFilter(IS_req.I_req, indexU, indexV, Neigh[0]) / deltaCSSampling; // VoisXs
                IYs = derivativeFilter(IS_req.I_req, indexU, indexV, Neigh[1]) / deltaCSSampling; // VoisYs
                IZs = derivativeFilter(IS_req.I_req, indexU, indexV, Neigh[2]) / deltaCSSampling; // VoisZs
            }
            else
            {
                IXs = IYs = IZs = 0.;
            }

            if (computePoseDerivatives)
            {
                derivativeVec[0] = IXs;
                derivativeVec[1] = IYs;
                derivativeVec[2] = IZs;
            }
        }
        computePoseDerivative(derivativeVec, rho, computePoseDerivatives);

        return 0;
    }

    /*!
     * \fn int buildFrom(prStereoModel *stereoCam, prRegularlySampledCSImage<unsigned char> & IS_req, Ts & _s_g, bool computePoseDerivatives = false, double rho = 1., unsigned long sourceIndex = 0)
     * \brief Computes the sample value from the set of a sample (one for each "pixel") from an equirectangular given camera, and, optionally (iff computePoseDerivatives == true), computes the photometric drivatives with respect to the 6 DOFs pose
     * TODO: take degrees of freedom into account to save processing time
     * TODO: find how to take properly account of the image type, if unsigned char the photometric value is put in the [0,1] interval
     * \param stereoCam pointer to the dual fisheye camera, used to project the feature on the image plane and likewise in the other direction
     * \param IS_req the input data from which to compute the Gaussian mixture sample value
     * \param _s_g the Gaussian mixture sample parameters (essentially, the location)
     * \param computePoseDerivatives boolean to let, if true, the method to compute the derivates with respect to the pose or not, if false
     * \param rho the distance of 3D point to the camera (1 by default if no 3D data is used)
     * \param sourceIndex the index of the sample to consider
     * \return 0 if the computation succeeds
     */
    template <typename T, template <typename TT> class Tacquisition>
    int buildFromCam(prStereoModel *stereoCam, Tacquisition<T> &IS_req, Ts &_s_g, bool computePoseDerivatives = false, double rho = 1., unsigned long sourceIndex = 0)
    {
        s_g = _s_g;


        inttyp = IS_req.inttyp;

        int imHeight = IS_req.I_req.getHeight();
        int imWidth = IS_req.I_req.getWidth();


        // CALCUL DE L'INTENSITE DE LA FEATURE MISE A JOUR
        float val = 0.;

        vpRowVector derivativeVec(3);
        derivativeVec = 0.;

        if (isNeighSet == false)
        {
            for (int numSen = 0; numSen < stereoCam->get_nbsens(); numSen++)
            {
                prNeighborhood neighborhoodStereo;
                neighborhoodStereo.buildNeighborsCartSphere(imWidth, imHeight, (prOmni *)(stereoCam->sen[numSen]));
                neighborVec.push_back(neighborhoodStereo.Neigh);
                deltaCSSamplingVec.push_back(neighborhoodStereo.deltaCSSampling);
            }

            isNeighSet = true;
        }

        prPointFeature P_temp;
        P_temp.set_X(s_g.get_X());
        P_temp.set_Y(s_g.get_Y());
        P_temp.set_Z(s_g.get_Z());

        int icam;

        if (P_temp.get_Z() > 0)
        {
            icam = 0;
        }
        else
        {
            icam = 1;
            P_temp.sX = P_temp.sX.changeFrame(stereoCam->sjMr[1]);
        }

        ((prOmni *)(stereoCam->sen[icam]))->project3DImage(P_temp);
        ((prOmni *)(stereoCam->sen[icam]))->meterPixelConversion(P_temp);

        int u = P_temp.get_u();
        int v = P_temp.get_v();

        if ((u >= 0) && (v >= 0) && (u < (imWidth - 1)) && (v < (imHeight - 1)))
        {
            indexU = (unsigned int)v;
            indexV = (unsigned int)u;

            T valTemp = IS_req.I_req[indexU][indexV];
            intensityVal = (float)valTemp;

            if (IS_req.MaskS[indexU][indexV] != 0)
            {
                IXs = derivativeFilter(IS_req.I_req, indexU, indexV, neighborVec[icam][0]) / deltaCSSamplingVec[icam]; // VoisXs
                IYs = derivativeFilter(IS_req.I_req, indexU, indexV, neighborVec[icam][1]) / deltaCSSamplingVec[icam]; // VoisYs
                IZs = derivativeFilter(IS_req.I_req, indexU, indexV, neighborVec[icam][2]) / deltaCSSamplingVec[icam]; // VoisZs
            }
            else
            {
                IXs = IYs = IZs = 0.;
            }

            // To take account of both sides of the sphere
            if (computePoseDerivatives)
            {
                if (icam == 0)
                {
                    derivativeVec[0] = IXs;
                    derivativeVec[1] = IYs;
                    derivativeVec[2] = IZs;
                }
                else
                {
                    derivativeVec[0] = -IXs;
                    derivativeVec[1] = -IYs;
                    derivativeVec[2] = -IZs;
                }
            }
        }
        computePoseDerivative(derivativeVec, rho, computePoseDerivatives);

        return 0;
    }

    /*!
     * \fn void computePoseDerivative(vpRowVector &derivativeVec, double rho, bool computePoseDerivatives)
     * \brief Computes the photometric derivatives with respect to the 6 DOFs pose
     * \param derivativeVec the computed vector containing the derivatives of the spherical image along each direction
     * \param rho the distance of 3D point to the camera (1 by default if no 3D data is used)
     * \param computePoseDerivatives boolean to let, if true, the method to compute the derivates with respect to the pose or not, if false
     * \return 0 if the computation succeeds
     */
    void computePoseDerivative(vpRowVector &derivativeVec, double rho, bool computePoseDerivatives)
    {

        vpMatrix dXgdr(3, 6);
        if (rho > 0.)
        {
            if (computePoseDerivatives)
            {
                double iRho = 1. / rho;

                dXgdr[0][0] = (s_g.get_X() * s_g.get_X() - 1.) * iRho;
                dXgdr[0][1] = s_g.get_X() * s_g.get_Y() * iRho;
                dXgdr[0][2] = s_g.get_X() * s_g.get_Z() * iRho;
                dXgdr[1][0] = dXgdr[0][1];
                dXgdr[1][1] = (s_g.get_Y() * s_g.get_Y() - 1.) * iRho;
                dXgdr[1][2] = s_g.get_Y() * s_g.get_Z() * iRho;
                dXgdr[2][0] = dXgdr[0][2];
                dXgdr[2][1] = dXgdr[1][2];
                dXgdr[2][2] = (s_g.get_Z() * s_g.get_Z() - 1.) * iRho;

                dXgdr[0][3] = 0.;
                dXgdr[0][4] = -s_g.get_Z();
                dXgdr[0][5] = s_g.get_Y();
                dXgdr[1][3] = s_g.get_Z();
                dXgdr[1][4] = 0.;
                dXgdr[1][5] = -s_g.get_X();
                dXgdr[2][3] = -s_g.get_Y();
                dXgdr[2][4] = s_g.get_X();
                dXgdr[2][5] = 0.;

                L = derivativeVec * dXgdr;

                poseJacobianComputed = true;
            }
        }
        else
        {
            L.resize(1, 6);
            poseJacobianComputed = true;
        }
    }

    template <typename T, template <typename TT> class Tacquisition>
    int updateFrom(Tacquisition<T> &IS_req, bool computeFeatureDerivatives = false, unsigned long sourceIndex = 0)
    {

        std::cout << "prIntensity::updateFrom not implemented yet" << std::endl;
        // return -1;

        double rho;

        IS_req.getSampleRho(sourceIndex, rho);
        return buildFrom(IS_req, s_g, computeFeatureDerivatives, rho, sourceIndex);
    }

    /*!
     * \fn double derivativeFilter(vpImage<T> &fr, int r, int c, int ****N)
     * \brief Computes the derivative of a sample based on its neighbors
     * \TODO: Take account of the image type, if unsigned char, divide by 255, else, already in the [0,1] interval so do nothing
     * \param fr The considered image whose derivative will be computed
     * \param r row index where to compute the derivative
     * \param c col index where to compute the derivative
     * \param N The neighboorhood already computed that gives the indexes of the neighbors used for the derivative computation
     * \return The derivative value of the image at a given index
     */
    template <typename T>
    double derivativeFilter(vpImage<T> &fr, int r, int c, int ****N)
    {
        double derTot;
        int *rN = N[r][c][0], *cN = N[r][c][1];

        inttyp = IMAGEPLANE_NEARESTNEIGH;

        switch (inttyp)
        {
        case IMAGEPLANE_BILINEAR:
            float rowI[6];
            int u, v;
            float du, dv, unmdu, unmdv;
            for (int i = 0; i < 6; i++)
            {
                v = rN[i];
                dv = rN[i] - v;
                unmdv = 1.0f - dv;
                u = cN[i];
                du = cN[i] - u;
                unmdu = 1.0f - du;

                rowI[i] = (float)fr[v][u] * unmdv * unmdu + (float)fr[v + 1][u] * dv * unmdu + (float)fr[v][u + 1] * unmdv * du + (float)fr[v + 1][u + 1] * dv * du;
            }
            derTot = (2047.0 * (rowI[3] - rowI[2]) + 913.0 * (rowI[4] - rowI[1]) + 112.0 * (rowI[5] - rowI[0])) * 0.00011879306;
            break;
        case IMAGEPLANE_NEARESTNEIGH:
        default:

            double der1 = 2047.0 * ((float)fr[rN[3]][cN[3]] - (float)fr[rN[2]][cN[2]]); // / 255.0;
            double der2 = 913.0 * ((float)fr[rN[4]][cN[4]] - (float)fr[rN[1]][cN[1]]);  // / 255.0;
            double der3 = 112.0 * ((float)fr[rN[5]][cN[5]] - (float)fr[rN[0]][cN[0]]);  // / 255.0;
            derTot = (der1 + der2 + der3) * 0.00011879306;

            break;
        }

        return derTot;
    }

    /*!
     * \fn int getPoseJacobian(vpRowVector & Ls)
     * \brief Accessor to the pose Jacobian of the Gaussian mixture sample
     * \param Ls (output) the output row vector that receives the photometric Gms interaction matrix row
     * \return  0 if the Jacobian could be get (if it was computed, actually)
     *         -1 if the Jacobian was not computed (no cannot be got)
     */
    int getPoseJacobian(vpRowVector &Ls)
    {
        if (!poseJacobianComputed)
            return -1;

        Ls = L;

        return 0;
    }

    /*!
     * \fn prIntensity<Ts> operator-(const prIntensity<Ts> &f)
     * \brief - operator override to substract two photometric photometric values
     * TODO : check line fOut.s_g = s_g - f.s_g;
     * \param f the photometric to substract from the current
     * \return the prIntensity<Ts> object resulting from the difference computation
     */
    prIntensity<Ts, Tcam> operator-(const prIntensity<Ts, Tcam> &f)
    {
        prIntensity<Ts, Tcam> fOut;

        fOut.intensityVal = intensityVal - f.getGMS();

        return fOut;
    }

    /*!
     * \fn prIntensity<Ts> operator*(const prIntensity<Ts> &f)
     * \brief * operator override to multiply two photometric photometric values
     * TODO: check line fOut.s_g = s_g * f.s_g;
     * TODO: move in the mother class?
     * \param f the photometric to multiply with the current
     * \return the prIntensity<Ts> object resulting from the product computation
     */
    prIntensity<Ts, Tcam> operator*(const prIntensity<Ts, Tcam> &f)
    {
        prIntensity<Ts, Tcam> fOut;

        fOut.intensityVal = intensityVal * f.getGMS();

        return fOut;
    }

    /*!
     * \fn prIntensity<Ts> operator*(const double w)
     * \brief * operator override to multiply a photometric value with a scalar
     * TODO : check line fOut.s_g = s_g * w;
     * \param w the scalar to multiply with the photometric
     * \return the prIntensity<Ts> object resulting from the product with a scalar computation
     */
    prIntensity<Ts, Tcam> operator*(const double w)
    {
        prIntensity<Ts, Tcam> fOut;

        fOut.intensityVal = intensityVal * w;

        return fOut;
    }

    /*!
     * \fn prIntensity<Ts> &operator+=(const prIntensity<Ts> &f)
     * \brief += operator override to add two photometric values, updating the current one
     * TODO : check line s_g += f.s_g;
     * \param f the photometric to add to the current
     * \return the prIntensity<Ts> self reference object
     */
    prIntensity<Ts, Tcam> &operator+=(const prIntensity<Ts, Tcam> &f)
    {
        // s_g += f.s_g;

        intensityVal += f.getGMS();

        return *this;
    }

    /*!
     * \fn prIntensity<Ts> &operator=(const prIntensity<Ts> &f)
     * \brief = operator override to copy photometric object
     * \TODO: make a clean memory copy of the neighbors to avoid computing them again
     * \param f the photometric to copy to the current
     * \return the prIntensity<Ts> self reference object
     */
    prIntensity<Ts, Tcam> &operator=(const prIntensity<Ts, Tcam> &f)
    {
        intensityVal = f.getGMS();
        s_g = f.getS_g();
        if (f.doubleFeature != NULL)
        {
            if (prFeature::doubleFeature != NULL)
                delete[] prFeature::doubleFeature;
            prFeature::doubleFeature = new double[1];

            prFeature::doubleFeature[0] = f.doubleFeature[0];
        }
        poseJacobianComputed = f.poseJacobianComputed;

        camera = f.camera;
        indexU = f.indexU;
        indexV = f.indexV;
        IXs = f.IXs;
        IYs = f.IYs;
        IZs = f.IZs;

        inttyp = f.inttyp;

        L = f.L;

        return *this;
    }

    /*!
     * \fn double getGMS() const (written as getGMS() for now to avoid some conflicts)
     * \brief Accessor to the photometric value
     * \return the photometric value
     */
    float getGMS() const
    {
        return intensityVal;
    }

    /*!
     * \fn double *toDouble(unsigned int place = 0)
     * \brief implementation of the toDouble method of the mother class converting the photometric feature type to a common double value
     * TODO: move in the mother class?
     * \param place ?
     * \return the pointer to the double type vlaue of the photometric
     */
    double *toDouble(unsigned int place = 0)
    {
        if (prFeature::doubleFeature != NULL)
            delete[] prFeature::doubleFeature;

        prFeature::doubleFeature = new double[1];
        prFeature::doubleFeature[0] = (double)intensityVal;

        return prFeature::doubleFeature;
    }

    /*!
     * \fn Ts getS_g() const
     * \brief Accessor to the Gaussian mixture sample location type parameter
     * TODO: move in the mother class?
     * \return the photometric location
     */
    Ts getS_g() const
    {
        return s_g;
    }

    void setSensor(Tcam *camera_)
    {
        camera = camera_;
    }

    unsigned int indexU, indexV;
    double IXs, IYs, IZs;
    std::vector<int *****> neighborVec;
    Tcam *camera;
    prNeighborhood neighborhood;

private:
    Ts s_g; /*!< !!!!!! feature ? Check the code using the class to remind */
    float intensityVal;
    vpRowVector L; /*!< interaction matrix of the gaussian mixture sample over the pose*/

    bool poseJacobianComputed; /*!< boolean to test if the pose Jacobian has already been computed to avoid useless computation */

    int *****Neigh;  // Neighborhood
    bool isNeighSet; // bool to check if the neighborhood is already set or not
    int dim_s;
    double deltaCSSampling;
    int nbDof;

    std::vector<double> deltaCSSamplingVec;

    prInterpType inttyp;
};

#endif //_PRINTENSITY_H
