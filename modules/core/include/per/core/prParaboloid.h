/*!
 \file prParaboloid.h
 \brief Definition of the adhoc model for paracatadioptric omnidirectional cameras
 \author Guillaume CARON
 \version 0.1
 \date january 2017
 */

#if !defined(_PRPARABOLOID_H)
#define _PRPARABOLOID_H

#include <per/prCameraModel.h>

/*!
 \class prParaboloid prParaboloid.h <per/prParaboloid.h>
 \brief Class for dealing with the projection functions of a paracatadioptric omnidirectional camera (adhoc model).
 */
class prParaboloid : public prCameraModel {
public:
    double h; /*!< Paracatadioptric camera model parameter: latus rectum of the paraboloidal shaped mirror */

    /*!
     * \fn prParaboloid(double au=0,double av=0,double u0=0,double v0=0,double h=0.5, double k1=0,double k2=0,double k3=0,double k4=0,double k5=0)
     * \brief Contructor of the prParaboloid camera model initializing intrinsic parameters
     * \param au the horizontal scale factor from the normalized image plane to the digital image (\alpha_u)
     * \param av the vertical scale factor from the normalized image plane to the digital image (\alpha_v)
     * \param u0 the horizontal coordinate of the principal point (u_0)
     * \param v0 the vertical coordinate of the principal point (v_0)
     * \param h  latus rectum of the paraboloidal shaped mirror (h)
     * \param k1 the second order radial distorsion parameter (k_1)
     * \param k2 the fourth order radial distorsion parameter (k_2)
     * \param k3 the sixth order radial distorsion parameter (k_3)
     * \param k4 the first tangential distorsion parameter (k_4)
     * \param k5 the second tangential distorsion parameter (k_5)
     */
    prParaboloid(double au=0,double av=0,double u0=0,double v0=0,double h=0.5, double k1=0,double k2=0,double k3=0,double k4=0,double k5=0);
    
    /*!
     * \fn ~prParaboloid()
     * \brief Destructor of a prParaboloid object
     */
    ~prParaboloid();
    
    
    /*!
     * \fn void init(double au,double av,double u0,double v0,double h, double k1=0,double k2=0,double k3=0,double k4=0,double k5=0)
     * \brief Initializes intrinsic parameters
     * \param au the horizontal scale factor from the normalized image plane to the digital image (\alpha_u)
     * \param av the vertical scale factor from the normalized image plane to the digital image (\alpha_v)
     * \param u0 the horizontal coordinate of the principal point (u_0)
     * \param v0 the vertical coordinate of the principal point (v_0)
     * \param h  latus rectum of the paraboloidal shaped mirror (h)
     * \param k1 the second order radial distorsion parameter (k_1)
     * \param k2 the fourth order radial distorsion parameter (k_2)
     * \param k3 the sixth order radial distorsion parameter (k_3)
     * \param k4 the first tangential distorsion parameter (k_4)
     * \param k5 the second tangential distorsion parameter (k_5)
     * \return Nothing
     */
    void init(double au,double av,double u0,double v0,double h, double k1=0,double k2=0,double k3=0,double k4=0,double k5=0);
    
    /*!
     * \fn double geth()
     * \brief Accessor to the h parameter of the prParaboloid sensor model
     * \return xi
     */
    double geth(){return h;};
    
    /*!
     * \fn void project3DImage(prPointFeature & P)
     * \brief Project a 3D point expressed in the paracatadioptric camera frame to the normalized image plane
     * \param P the point to project
     * \return Nothing
     */
    void project3DImage(prPointFeature & P);
    
    /*!
     * \fn void projectImageMiroir(prPointFeature & P, double & Xm, double & Ym, double & Zm)
     * \brief Inverse projection from the normalized image plane to the paraboloidal mirror
     * \param P the point to backproject
     * \param Xm the Cartesian X coordinate of the projected point on the mirror
     * \param Ym the Cartesian Y coordinate of the projected point on the mirror
     * \param Zm the Cartesian Z coordinate of the projected point on the mirror
     * \return Nothing
     */
    void projectImageMiroir(prPointFeature & P, double & Xm, double & Ym, double & Zm);
    
    /*!
     * \fn void project3DMiroir(prPointFeature & P, double & Xm, double & Ym, double & Zm)
     * \brief Project a 3D point expressed in the paracatadioptric camera frame to its mirror
     * \param P the point to project
     * \param Xm the Cartesian X coordinate of the projected point on the mirror
     * \param Ym the Cartesian Y coordinate of the projected point on the mirror
     * \param Zm the Cartesian Z coordinate of the projected point on the mirror
     * \return Nothing
     */
    void project3DMiroir(prPointFeature & P, double & Xm, double & Ym, double & Zm);
    
    
    /*!
     * \fn prParaboloid& operator=(const prParaboloid& cam)
     * \brief = operator overload for copying a prParaboloid camera model
     *
     * \param cam the source prParaboloid
     * \return the affected prParaboloid.
     */
    prParaboloid& operator=(const prParaboloid & cam);
};

#endif  //_PRPARABOLOID_H
