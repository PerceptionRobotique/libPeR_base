/*!
 \file prEquirectangular.h
 \brief Definition of the equirectangular imaging camera
 \author Guillaume CARON
 \version 0.1
 \date january 2017
 */

#if !defined(_PREQUIRECTANGULAR_H)
#define _PREQUIRECTANGULAR_H

#include <per/prCameraModel.h>

/*!
 \class PER_EXPORT prEquirectangular prEquirectangular.h <per/prEquirectangular.h>
 \brief Class for dealing with the projection functions of an equirectangular imaging camera.
 */
class PER_EXPORT prEquirectangular : public prCameraModel
{
public:
    
    /*!
     * \fn prEquirectangular(double au=0,double av=0,double u0=0,double v0=0, double k1=0,double k2=0,double k3=0,double k4=0,double k5=0)
     * \brief Contructor of the prEquirectangular camera model initializing intrinsic parameters
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
    prEquirectangular(double au=0,double av=0,double u0=0,double v0=0,double k1=0,double k2=0,double k3=0,double k4=0,double k5=0);
    
    /*!
     * \fn ~prEquirectangular()
     * \brief Destructor of a prEquirectangular object
     */
    virtual ~prEquirectangular() override = default;
    
    
    /*!
     * \fn void init(double au,double av,double u0,double v0,double k1=0,double k2=0,double k3=0,double k4=0,double k5=0)
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
     * \brief Project a 3D point expressed in the spherical camera frame to the normalized image plane
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
     * \brief Project a 3D point expressed in the spherical camera frame to the unit sphere
     * \param P the point to project
     * \param Xs the Cartesian X coordinate of the projected point on the unit sphere
     * \param Ys the Cartesian Y coordinate of the projected point on the unit sphere
     * \param Zs the Cartesian Z coordinate of the projected point on the unit sphere
     * \return Nothing
     */
    void project3DSphere(prPointFeature & P, double & Xs, double & Ys, double & Zs);
    
    /*!
     * \fn virtual void projectImageSphere(prPointFeature & P, double & Xs, double & Ys, double & Zs)
     * \brief Inverse projection from the normalized image plane to the unit sphere
     * \param P the point to backproject
     * \param Xs the Cartesian X coordinate of the projected point on the unit sphere
     * \param Ys the Cartesian Y coordinate of the projected point on the unit sphere
     * \param Zs the Cartesian Z coordinate of the projected point on the unit sphere
     * \return Nothing
     */
    virtual bool projectImageSphere(prPointFeature & P, double & Xs, double & Ys, double & Zs);
    
    /*!
     * \fn prEquirectangular& operator=(const prEquirectangular& cam)
     * \brief = operator overload for copying a prEquirectangular camera model
     *
     * \param cam the source prEquirectangular
     * \return the affected prEquirectangular.
     */
    prEquirectangular& operator=(const prEquirectangular& cam);

    /*!
     * \fn void computeSensorJacobian(prPointFeature & P, vpMatrix & LuX)
     * \brief computes the Jacobian of the sensor digital image point coordinates with respect to the 3D coordinates in the camera frame
     * \param P the point for which computing the Jacobian
     * \param LuX the output Jacobian Matrix
     * \return Nothing
     */
    void computeSensorJacobian(prPointFeature & P, vpMatrix & LuX);
};

#endif  //_PREQUIRECTANGULAR_H
