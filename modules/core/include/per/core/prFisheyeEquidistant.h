/*!
 \file prFisheyeEquidistant.h
 \brief Definition of an equidistant fisheye camera model (f-theta model)
 \author Guillaume CARON
 \version 0.1
 \date march 2022
 */

#if !defined(_PRFISHEYEEQUIDISTANT_H)
#define _PRFISHEYEEQUIDISTANT_H

#include <per/prCameraModel.h>

/*!
 \class PER_EXPORT prFisheyeEquidistant prFisheyeEquidistant.h <per/prFisheyeEquidistant.h>
 \brief Class for dealing with the projection functions of an equidistant fisheye camera.
 */
class PER_EXPORT prFisheyeEquidistant : public prCameraModel {
public:
    /*!
     * \fn prFisheyeEquidistant(double au=0,double av=0,double u0=0,double v0=0, double k1=0,double k2=0,double k3=0,double k4=0,double k5=0)
     * \brief Contructor of the prFisheyeEquidistant camera model initializing intrinsic parameters
     * \param au the horizontal scale factor from the normalized image plane to the digital image (\alpha_u)
     * \param av the vertical scale factor from the normalized image plane to the digital image (\alpha_v)
     * \param u0 the horizontal coordinate of the principal point (u_0)
     * \param v0 the vertical coordinate of the principal point (v_0)
     * \param k1 the second order radial distorsion parameter (k_1)
     * \param k2 the fourth order radial distorsion parameter (k_2)
     * \param k3 the sixth order radial distorsion parameter (k_3)
     * \param k4 the first tangential distorsion parameter (k_4)
     * \param k5 the second tangential distorsion parameter (k_5)
     */
    prFisheyeEquidistant(double au=0,double av=0,double u0=0,double v0=0, double k1=0,double k2=0,double k3=0,double k4=0,double k5=0);
    
    /*!
     * \fn ~prFisheyeEquidistant()
     * \brief Destructor of a prFisheyeEquidistant object
     */
    virtual ~prFisheyeEquidistant() override = default;

    
    /*!
     * \fn void init(double au,double av,double u0,double v0, double k1=0,double k2=0,double k3=0,double k4=0,double k5=0)
     * \brief Initializing intrinsic parameters
     * \param au the horizontal scale factor from the normalized image plane to the digital image (\alpha_u)
     * \param av the vertical scale factor from the normalized image plane to the digital image (\alpha_v)
     * \param u0 the horizontal coordinate of the principal point (u_0)
     * \param v0 the vertical coordinate of the principal point (v_0)
     * \param k1 the second order radial distorsion parameter (k_1)
     * \param k2 the fourth order radial distorsion parameter (k_2)
     * \param k3 the sixth order radial distorsion parameter (k_3)
     * \param k4 the first tangential distorsion parameter (k_4)
     * \param k5 the second tangential distorsion parameter (k_5)
     * \return Nothing
     */
    void init(double au,double av,double u0,double v0, double k1=0,double k2=0,double k3=0,double k4=0,double k5=0);
    
    /*!
     * \fn void project3DImage(prPointFeature & P)
     * \brief Project a 3D point expressed in the equidistant fisheye camera frame to the normalized image plane
     * \param P the point to project
     * \return Nothing
     */
    void project3DImage(prPointFeature & P);

    /*!
     * \fn virtual void unProject(prPointFeature & P, double & Depth)=0
     * \brief Compute the 3D point expressed in the sensor frame from the normalized image plane and the depth (must be implemented in a specified sensor model)
     * \param P the point to unproject
     * \param Depth the Depth used to unproject the P (/rho here...)
     * \return true if unprojection could be done
     */
    bool unProject(prPointFeature & P, double & Depth);
    
    /*!
     * \fn void project3DSphere(prPointFeature & P, double & Xs, double & Ys, double & Zs)
     * \brief Project a 3D point expressed in the equidistant fisheye camera frame to the unit sphere centered at the camera origin
     * \param P the point to project
     * \param Xs the Cartesian X coordinate of the projected point on the unit sphere
     * \param Ys the Cartesian Y coordinate of the projected point on the unit sphere
     * \param Zs the Cartesian Z coordinate of the projected point on the unit sphere
     * \return Nothing
     */
    void project3DSphere(prPointFeature & P, double & Xs, double & Ys, double & Zs);
    
    /*!
     * \fn virtual void projectImageSphere(prPointFeature & P, double & Xs, double & Ys, double & Zs)
     * \brief Inverse projection from the normalized image plane to the unit sphere centered at the camera origin
     * \param P the point to backproject
     * \param Xs the Cartesian X coordinate of the projected point on the unit sphere
     * \param Ys the Cartesian Y coordinate of the projected point on the unit sphere
     * \param Zs the Cartesian Z coordinate of the projected point on the unit sphere
     * \return Nothing
     */
    virtual bool projectImageSphere(prPointFeature & P, double & Xs, double & Ys, double & Zs);
    
    /*!
     * \fn prFisheyeEquidistant& operator=(const prFisheyeEquidistant& cam)
     * \brief = operator overload for copying a prFisheyeEquidistant camera model
     *
     * \param cam the source prFisheyeEquidistant
     * \return the affected prFisheyeEquidistant.
     */
    prFisheyeEquidistant& operator=(const prFisheyeEquidistant& cam);

    /*!
     * \fn void computeSensorJacobian(prPointFeature & P, vpMatrix & LuX)
     * \brief computes the Jacobian of the sensor digital image point coordinates with respect to the 3D coordinates in the camera frame
     * \param P the point for which computing the Jacobian
     * \param LuX the output Jacobian Matrix
     * \return Nothing
     */
    void computeSensorJacobian(prPointFeature & P, vpMatrix & LuX);
};

#endif  //_PRFISHEYEEQUIDISTANT_H
